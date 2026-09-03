"""Planar polygon-frame DNA origami built on Curved Design routing.

The module deliberately treats Curved Design's DNAxiS/AutoCS result as a
frozen topology.  It only relocates the already-budgeted insertion/deletion
count into vertex-local bend windows and replaces the embedded circular
target geometry with a rounded polygon target geometry.
"""

import base64
import gzip
import hashlib
import itertools
import json
import math
import os
import struct
import zlib

from .athena import encode_geometry_payload
from .curved import (BP_RISE_NM, DNA_HELIX_RADIUS_NM, MIN_RING_BP,
                     build_rings, create_curved_project, safe_name)
from .indelanalysis import (optimize_frame_pair_curvature,
                            short_staple_deletion_protection)
from .twistbend import (TwistBendError, equal_partition_indel_sites,
                        estimate_global_twist)


FRAME_METADATA_VERSION = 1
MAX_INDEL_PER_DOMAIN = 3
OXDNA_LENGTH_NM = 0.8518

FRAME_PRESETS = (
    ("triangle", "Equilateral triangle"),
    ("square", "Square"),
    ("pentagon", "Regular pentagon"),
    ("hexagon", "Regular hexagon"),
    ("rectangle", "Rectangle"),
    ("parallelogram", "Parallelogram"),
    ("image", "Polygon from image"))


def _distance(first, second):
    return math.hypot(second[0] - first[0], second[1] - first[1])


def _signed_area(vertices):
    return 0.5 * sum(
        first[0] * second[1] - second[0] * first[1]
        for first, second in zip(vertices, vertices[1:] + vertices[:1]))


def _unit(vector):
    length = math.hypot(vector[0], vector[1])
    if length <= 1e-12:
        raise ValueError("Polygon contains a zero-length edge.")
    return vector[0] / length, vector[1] / length


def regular_polygon(vertex_count, side_nm):
    vertex_count = int(vertex_count)
    side_nm = float(side_nm)
    if vertex_count < 3 or side_nm <= 0:
        raise ValueError("A regular polygon needs at least three sides.")
    radius = side_nm / (2.0 * math.sin(math.pi / vertex_count))
    phase = math.pi / 2.0
    vertices = [(radius * math.cos(phase + 2.0 * math.pi * index /
                                  vertex_count),
                 radius * math.sin(phase + 2.0 * math.pi * index /
                                  vertex_count))
                for index in range(vertex_count)]
    # Use counter-clockwise order throughout the geometry code.
    return list(reversed(vertices)) if _signed_area(vertices) < 0 else vertices


def rectangle_polygon(width_nm, height_nm):
    width_nm, height_nm = float(width_nm), float(height_nm)
    if min(width_nm, height_nm) <= 0:
        raise ValueError("Rectangle side lengths must be positive.")
    return [(0.0, 0.0), (width_nm, 0.0),
            (width_nm, height_nm), (0.0, height_nm)]


def parallelogram_polygon(first_nm, second_nm, angle_degrees):
    first_nm, second_nm = float(first_nm), float(second_nm)
    angle = math.radians(float(angle_degrees))
    if min(first_nm, second_nm) <= 0 or not 15.0 <= angle_degrees <= 165.0:
        raise ValueError("Parallelogram sides/angle are outside safe limits.")
    offset = (second_nm * math.cos(angle), second_nm * math.sin(angle))
    return [(0.0, 0.0), (first_nm, 0.0),
            (first_nm + offset[0], offset[1]), offset]


def polygon_from_spec(spec):
    shape = str(spec.get("frame_shape", "square")).lower()
    if shape == "triangle":
        return regular_polygon(3, spec["side_nm"])
    if shape == "square":
        return regular_polygon(4, spec["side_nm"])
    if shape == "pentagon":
        return regular_polygon(5, spec["side_nm"])
    if shape == "hexagon":
        return regular_polygon(6, spec["side_nm"])
    if shape == "rectangle":
        return rectangle_polygon(spec["first_side_nm"],
                                 spec["second_side_nm"])
    if shape == "parallelogram":
        return parallelogram_polygon(
            spec["first_side_nm"], spec["second_side_nm"],
            spec.get("corner_angle_degrees", 70.0))
    if shape == "image":
        vertices = [(float(row[0]), float(row[1]))
                    for row in spec.get("image_vertices", ())]
        if len(vertices) < 3:
            raise ValueError("No valid polygon was detected in the image.")
        supplied = [float(value) for value in
                    spec.get("image_side_lengths_nm", ())]
        if supplied:
            if len(supplied) != len(vertices) or min(supplied) <= 0:
                raise ValueError(
                    "Image polygon needs one positive length for every edge.")
            rebuilt = [(0.0, 0.0)]
            for index, length in enumerate(supplied[:-1]):
                first, second = vertices[index], vertices[index + 1]
                direction = _unit((second[0] - first[0],
                                   second[1] - first[1]))
                rebuilt.append((rebuilt[-1][0] + direction[0] * length,
                                rebuilt[-1][1] + direction[1] * length))
            final_direction = _unit((vertices[0][0] - vertices[-1][0],
                                     vertices[0][1] - vertices[-1][1]))
            closure = (rebuilt[-1][0] + final_direction[0] * supplied[-1],
                       rebuilt[-1][1] + final_direction[1] * supplied[-1])
            if math.hypot(*closure) > 0.05 * sum(supplied):
                raise ValueError(
                    "The entered edge lengths do not close the detected "
                    "polygon (closure error exceeds 5%).")
            return rebuilt
        reference = float(spec.get("image_reference_side_nm", 50.0))
        source = _distance(vertices[0], vertices[1])
        if source <= 0:
            raise ValueError("The detected reference edge has zero length.")
        scale = reference / source
        origin = vertices[0]
        return [((point[0] - origin[0]) * scale,
                 (point[1] - origin[1]) * scale) for point in vertices]
    raise ValueError("Unknown Frame Design shape: %s" % shape)


def polygon_metrics(vertices):
    vertices = [tuple(map(float, point)) for point in vertices]
    if len(vertices) < 3:
        raise ValueError("A frame needs at least three vertices.")
    if _signed_area(vertices) < 0:
        vertices.reverse()
    side_lengths = [_distance(first, second) for first, second in
                    zip(vertices, vertices[1:] + vertices[:1])]
    if min(side_lengths) <= 1e-9:
        raise ValueError("Polygon contains a zero-length edge.")
    turns = []
    for index, vertex in enumerate(vertices):
        previous = vertices[index - 1]
        following = vertices[(index + 1) % len(vertices)]
        incoming = _unit((vertex[0] - previous[0],
                          vertex[1] - previous[1]))
        outgoing = _unit((following[0] - vertex[0],
                          following[1] - vertex[1]))
        cross = incoming[0] * outgoing[1] - incoming[1] * outgoing[0]
        dot = max(-1.0, min(1.0,
                           incoming[0] * outgoing[0] +
                           incoming[1] * outgoing[1]))
        turn = math.atan2(cross, dot)
        if turn <= 1e-6:
            raise ValueError(
                "Frame Design currently requires a simple convex polygon.")
        turns.append(turn)
    return {"vertices": vertices, "side_lengths_nm": side_lengths,
            "turn_angles_rad": turns,
            "interior_angles_degrees":
                [180.0 - math.degrees(value) for value in turns],
            "perimeter_sharp_nm": sum(side_lengths)}


def _ring_normal_offsets(spec, neutral_perimeter_nm):
    lattice = str(spec.get("lattice", "square")).lower()
    layers = int(spec.get("layers", 1))
    outer_diameter = neutral_perimeter_nm / math.pi + \
        2.0 * DNA_HELIX_RADIUS_NM
    rings = build_rings(
        "cylinder", float(spec["cross_section_height_nm"]),
        outer_diameter, outer_diameter, layers, lattice=lattice)
    radii = [float(row.get("indel_radius_nm", row["radius_nm"]))
             for row in rings]
    radius_middle = 0.5 * (min(radii) + max(radii))
    heights = [float(row["height_nm"]) for row in rings]
    height_middle = 0.5 * (min(heights) + max(heights))
    for ring, radius, height in zip(rings, radii, heights):
        # Preserve Curved Design's physical bend face.  Its radial coordinate
        # is the coordinate that receives differential helix lengths, so it
        # must remain the in-plane Frame normal (inside/outside).  The older
        # Frame prototype swapped these two axes and silently rotated the
        # Square/Honeycomb cross-section by 90 degrees.
        ring["frame_normal_offset_nm"] = radius - radius_middle
        ring["frame_binormal_offset_nm"] = height - height_middle
    return rings


def _largest_remainder(total, weights):
    total = int(total)
    sign = 1 if total >= 0 else -1
    magnitude = abs(total)
    weight_sum = sum(weights)
    raw = [magnitude * value / weight_sum for value in weights]
    result = [int(math.floor(value)) for value in raw]
    remaining = magnitude - sum(result)
    order = sorted(range(len(raw)), key=lambda index:
                   (-(raw[index] - result[index]), index))
    for index in order[:remaining]:
        result[index] += 1
    return [sign * value for value in result]


def plan_frame(spec):
    """Return a complete rounded-polygon and vertex-indel plan."""
    lattice = str(spec.get("lattice", "square")).lower()
    if lattice not in ("square", "honeycomb"):
        raise ValueError("Frame lattice must be square or honeycomb.")
    domain_size = 7 if lattice == "honeycomb" else 8
    maximum_indel_per_domain = max(1, min(
        MAX_INDEL_PER_DOMAIN,
        int(spec.get("maximum_indel_per_domain", MAX_INDEL_PER_DOMAIN))))
    # Honeycomb has a 21-bp lattice/crossover phase, but the closed common
    # parent used by Curved Design is an alternating 42-bp unit.  Keep these
    # concepts separate: density and phase continue to use 21, while the
    # neutral-axis parent is snapped to 42N.  Square uses 32 for both.
    native_period = 21 if lattice == "honeycomb" else 32
    parent_period = 42 if lattice == "honeycomb" else 32
    layers = int(spec.get("layers", 1))
    base_minimum_parent = 84 if lattice == "honeycomb" else 96
    # create_curved_project receives the base reinforcement layer as its
    # cylinder diameter.  Its minimax common parent is one complete parent
    # period farther out for every reinforcement layer.  Ensure that the
    # inward-shifted base layer also remains above Curved's ring minimum.
    minimum_parent = max(
        base_minimum_parent,
        int(math.ceil((MIN_RING_BP+(layers-1)*parent_period) /
                      float(parent_period))) * parent_period)
    metrics = polygon_metrics(polygon_from_spec(spec))
    requested_perimeter = metrics["perimeter_sharp_nm"]
    # Cross-section offsets do not depend on the provisional parent radius.
    # The exact native parent length is solved after the corner windows are
    # known, because replacing a sharp vertex by an arc shortens the path.
    rings = _ring_normal_offsets(spec, requested_perimeter)

    maximum_delta_by_vertex = []
    per_ring_vertex_delta = []
    for ring in rings:
        offset = float(ring["frame_normal_offset_nm"])
        values = [int(round(angle * offset / BP_RISE_NM))
                  for angle in metrics["turn_angles_rad"]]
        per_ring_vertex_delta.append(values)
    for vertex_index in range(len(metrics["vertices"])):
        maximum_delta_by_vertex.append(max(
            abs(values[vertex_index]) for values in per_ring_vertex_delta))

    bend_mode = str(spec.get("bend_length_mode", "auto")).lower()
    custom_bp = int(spec.get("bend_length_bp", domain_size))
    bend_lengths = []
    domain_counts = []
    average_densities = []
    realized_densities = []
    radii = []
    tangencies = []
    for angle, maximum_delta in zip(
            metrics["turn_angles_rad"], maximum_delta_by_vertex):
        minimum_domains = max(1, int(math.ceil(
            maximum_delta / float(maximum_indel_per_domain))))
        if bend_mode == "custom":
            domains = max(1, int(round(custom_bp / float(domain_size))))
        else:
            domains = minimum_domains
        bend_bp = domains * domain_size
        average = maximum_delta / float(domains)
        realized = int(math.ceil(average))
        bend_length_nm = bend_bp * BP_RISE_NM
        radius = bend_length_nm / angle
        tangent = radius * math.tan(angle / 2.0)
        bend_lengths.append(bend_bp)
        domain_counts.append(domains)
        average_densities.append(average)
        realized_densities.append(realized)
        radii.append(radius)
        tangencies.append(tangent)

    bend_lengths_nm = [value * BP_RISE_NM for value in bend_lengths]
    rounded_requested = (requested_perimeter - 2.0 * sum(tangencies) +
                         sum(bend_lengths_nm))
    nominal_bp = max(minimum_parent, int(round(
        rounded_requested / BP_RISE_NM / parent_period)) * parent_period)
    target_perimeter = nominal_bp * BP_RISE_NM
    curved_base_perimeter_bp = nominal_bp-(layers-1)*parent_period
    curved_outer_diameter_nm = (
        curved_base_perimeter_bp*BP_RISE_NM/math.pi+
        2.0*DNA_HELIX_RADIUS_NM)
    # Solve one uniform scale for the sharp polygon so that its rounded path,
    # not its unrounded outline, is exactly the native closed-loop length.
    scale = ((target_perimeter + 2.0 * sum(tangencies) -
              sum(bend_lengths_nm)) / requested_perimeter)
    vertices = [(point[0] * scale, point[1] * scale)
                for point in metrics["vertices"]]
    metrics = polygon_metrics(vertices)
    rings = _ring_normal_offsets(spec, target_perimeter)

    fit_remaining = []
    for index, side in enumerate(metrics["side_lengths_nm"]):
        remaining = side - tangencies[index] - \
            tangencies[(index + 1) % len(tangencies)]
        fit_remaining.append(remaining)
    feasible = (max(realized_densities or [0]) <=
                maximum_indel_per_domain and
                min(fit_remaining or [0]) >= -1e-8)
    reasons = []
    if max(realized_densities or [0]) > maximum_indel_per_domain:
        reasons.append(
            "vertex indel density exceeds +/-%d per domain" %
            maximum_indel_per_domain)
    if min(fit_remaining or [0]) < -1e-8:
        reasons.append("rounded corners overlap on at least one edge")

    # Native coordinate zero is the exit of the final rounded corner.  Use
    # exact domain-sized arc windows and allocate only the remaining bases to
    # straight edges.  Vertex centres are therefore the true arc midpoints.
    straight_total_bp = nominal_bp - sum(bend_lengths)
    straight_total_nm = sum(max(0.0, value) for value in fit_remaining)
    straight_bp_raw = [
        (straight_total_bp * max(0.0, value) / straight_total_nm
         if straight_total_nm else 0.0) for value in fit_remaining]
    straight_bp = [int(math.floor(value)) for value in straight_bp_raw]
    for index in sorted(range(len(straight_bp)), key=lambda item:
                        (-(straight_bp_raw[item] - straight_bp[item]), item))[
                            :straight_total_bp - sum(straight_bp)]:
        straight_bp[index] += 1
    raw_vertex_centres = []
    cursor = 0.0
    for straight, bend in zip(straight_bp, bend_lengths):
        cursor += straight
        raw_vertex_centres.append(cursor + bend / 2.0)
        cursor += bend
    # Place base 0 at the midpoint of the longest straight edge.  This keeps
    # every vertex-local bend window away from the linear caDNAno array ends
    # without changing the closed path or any AutoCS topology.
    origin_edge = max(range(len(straight_bp)),
                      key=lambda index: (straight_bp[index], -index))
    origin_shift = sum(
        straight_bp[index] + bend_lengths[index]
        for index in range(origin_edge)) + straight_bp[origin_edge] / 2.0
    vertex_centres = [
        (value - origin_shift) % nominal_bp
        for value in raw_vertex_centres]

    # The parent cylinder calculates the same closed-loop residual.  Store a
    # predicted target for UI/reporting; the postprocessor preserves the exact
    # residual emitted by the stable Curved pipeline.
    for ring, values in zip(rings, per_ring_vertex_delta):
        ring["frame_vertex_indels"] = values
        ring["frame_total_indel"] = sum(values)
        ring["frame_target_bp"] = nominal_bp + sum(values)

    return {
        "version": FRAME_METADATA_VERSION, "lattice": lattice,
        "domain_size_bp": domain_size, "native_period_bp": native_period,
        "parent_period_bp": parent_period,
        "minimum_parent_bp": minimum_parent,
        "curved_base_perimeter_bp": curved_base_perimeter_bp,
        "curved_parent_outer_diameter_nm": curved_outer_diameter_nm,
        "nominal_perimeter_bp": nominal_bp,
        "nominal_perimeter_nm": target_perimeter,
        "input_scale_factor": scale, "vertices_nm": vertices,
        "side_lengths_nm": metrics["side_lengths_nm"],
        "straight_native_bp": straight_bp,
        "turn_angles_degrees": [math.degrees(value)
                                 for value in metrics["turn_angles_rad"]],
        "interior_angles_degrees": metrics["interior_angles_degrees"],
        "bend_length_bp": bend_lengths,
        "bend_domain_count": domain_counts,
        "bend_radius_nm": radii,
        "tangent_trim_nm": tangencies,
        "straight_edge_nm": fit_remaining,
        "vertex_native_centres": vertex_centres,
        "native_origin_shift_bp": origin_shift,
        "maximum_abs_indel_by_vertex": maximum_delta_by_vertex,
        "average_max_indel_per_domain": average_densities,
        "realized_max_indel_per_domain": realized_densities,
        "maximum_indel_per_domain_allowed": maximum_indel_per_domain,
        "feasible": feasible, "failure_reasons": reasons,
        "rings": rings}


def _external(entry, helix):
    return (entry == [-1, -1, -1, -1] or int(entry[0]) < 0 or
            int(entry[2]) < 0 or int(entry[0]) != helix or
            int(entry[2]) != helix)


def _safe_sites(row, helix, start, end):
    size = min(len(row["scaf"]), len(row["stap"]))
    return [index for index in range(max(1, int(start)),
                                    min(size - 1, int(end)))
            if not _external(row["scaf"][index], helix) and
            not _external(row["stap"][index], helix)]


def _partition_sites(count, start, end, candidates, domain_size,
                     loads, insertion, maximum_per_domain=3):
    """Choose sites using the same equal-partition hierarchy as Curved.

    The vertex-local bend window remains a hard scope.  Inside it, selection
    follows the shared ordering: equal axial partition intersected with the
    preferred absolute 7/8-bp domain, another safe domain in that partition,
    the preferred domain outside the partition, then the remaining nearest
    safe site.  Existing per-domain loads remain authoritative.
    """
    count = int(count)
    if count <= 0:
        return []
    available = sorted(set(int(site) for site in candidates))
    candidates_by_domain = {}
    for site in available:
        candidates_by_domain.setdefault(site // domain_size, []).append(site)
    capacities = {}
    for domain, sites in candidates_by_domain.items():
        remaining = max(0, int(maximum_per_domain)-
                        int(loads.get(domain, 0)))
        capacities[domain] = (remaining if insertion and sites else
                              min(remaining, len(sites)))
    if sum(capacities.values()) < count:
        raise ValueError(
            "No safe indel capacity remains inside a vertex bend window.")
    try:
        records = equal_partition_indel_sites(
            count, int(math.floor(start)), int(math.ceil(end))-1,
            domain_size, candidates_by_domain, capacities,
            allow_repeated_sites=insertion)
    except TwistBendError as error:
        raise ValueError(
            "No safe indel site remains inside a vertex bend window: %s" %
            error)
    for item in records:
        domain = int(item["domain"])
        loads[domain] = int(loads.get(domain, 0))+1
    return [int(item["idx"]) for item in records]


def _decode_geometry(metadata):
    encoded = metadata.get("geometry_data")
    if not encoded:
        return None
    return json.loads(gzip.decompress(base64.b64decode(encoded)))


def topology_fingerprint(design):
    """Fingerprint all scaffold/staple links, excluding indel arrays."""
    topology = [{"num": int(row["num"]), "scaf": row["scaf"],
                 "stap": row["stap"]}
                for row in sorted(design.get("vstrands", ()),
                                  key=lambda value: int(value["num"]))]
    raw = json.dumps(topology, separators=(",", ":"),
                     sort_keys=True).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _crossover_fingerprint(design):
    """Return only inter-helix links; longitudinal nick changes are ignored."""
    result = set()
    for row in design.get("vstrands", ()):
        helix = int(row["num"])
        for strand_type in ("scaf", "stap"):
            for index, entry in enumerate(row.get(strand_type, ())):
                for offset in (0, 2):
                    other, other_index = map(int, entry[offset:offset+2])
                    if other >= 0 and other != helix:
                        result.add((strand_type, helix, index, offset,
                                    other, other_index))
    return frozenset(result)


def _straight_complete_domain_edges(plan):
    """Return physical straight edges and their complete native domains.

    Frame corner windows are circular, while caDNAno's base arrays are
    linear.  Testing each base centre with a circular signed distance keeps a
    corner that crosses base zero intact.  A 7/8-bp domain is admitted only
    when all of its bases are outside every corner window, so the subsequent
    common-mode correction can never share a domain with a curvature indel.

    Crucially, a physical edge that crosses base zero remains *one* edge with
    two linear caDNAno segments.  Treating those two segments independently
    used to turn a regular triangle into four allocation bins and produced
    asymmetric quotas such as 0/3/2/0 instead of one quota per physical edge.
    """
    nominal = int(plan["nominal_perimeter_bp"])
    domain_size = int(plan["domain_size_bp"])
    windows = [(float(centre) % nominal, float(length), index)
               for index, (centre, length) in enumerate(zip(
                   plan.get("vertex_native_centres", ()),
                   plan.get("bend_length_bp", ())))]
    windows.sort(key=lambda item: item[0])

    def in_bend(index):
        base_centre = float(index) + .5
        for centre, length, unused_vertex in windows:
            delta = ((base_centre-float(centre)+nominal/2.0) % nominal -
                     nominal/2.0)
            if abs(delta) < float(length)/2.0-1e-12:
                return True
        return False

    complete = []
    for domain in range(nominal // domain_size):
        start = domain*domain_size
        end = start+domain_size-1
        if not any(in_bend(index) for index in range(start, end+1)):
            complete.append((start, end))
    if not windows:
        return []

    # A physical straight edge runs from the end of one bend window to the
    # start of the next one.  Domains are assigned by their circular distance
    # from that physical start, not by their linear caDNAno coordinate.
    edges = []
    for index, (centre, length, vertex) in enumerate(windows):
        next_centre, next_length, next_vertex = windows[
            (index+1) % len(windows)]
        edge_start = (centre+length/2.0) % nominal
        edge_end = (next_centre-next_length/2.0) % nominal
        native_length = (edge_end-edge_start) % nominal
        edge_domains = []
        for start, end in complete:
            domain_centre = 0.5*(start+end+1)
            distance = (domain_centre-edge_start) % nominal
            if distance < native_length-1e-9:
                edge_domains.append((distance, start, end))
        edge_domains.sort()
        segments = []
        for unused_distance, start, end in edge_domains:
            if segments and start == segments[-1][1]+1:
                segments[-1] = (segments[-1][0], end)
            else:
                segments.append((start, end))
        edges.append({
            "edge_index": int(vertex),
            "from_vertex": int(vertex),
            "to_vertex": int(next_vertex),
            "start_coordinate_bp": float(edge_start),
            "native_length_bp": int(round(native_length)),
            "segments": segments,
            "eligible_length_bp": int(sum(
                end-start+1 for start, end in segments))})
    return edges


def _frame_twist_helix_data(design, lattice, straight_length):
    """Build a synthetic straight-region mechanical cross-section.

    The residual-twist predictor needs the real lattice coordinates and the
    frozen AutoCS connectivity, but it must not include corner differential
    indels.  Therefore the actual Frame design supplies its cross-section and
    neighbor graph, while every helix receives one synthetic, fully paired
    native interval whose length equals the total eligible straight length.
    """
    from cadnano2.model.document import Document
    from cadnano2.model.enum import LatticeType
    from cadnano2.model.io.legacydecoder import import_legacy_dict

    lattice_type = (LatticeType.Honeycomb if lattice == "honeycomb"
                    else LatticeType.Square)
    document = Document()
    part = import_legacy_dict(
        document, design, lattice_type, forceLatticeType=True)
    if part is None:
        raise RuntimeError(
            "Cannot reload Frame design for straight-edge twist prediction.")
    neighbors = {}
    for vh in part.getVirtualHelices():
        linked = set()
        for strand_set in vh.getStrandSets():
            for strand in strand_set:
                for connected in (strand.connection5p(),
                                  strand.connection3p()):
                    if connected is not None and \
                            connected.virtualHelix() is not vh:
                        linked.add(int(connected.virtualHelix().number()))
        neighbors[int(vh.number())] = linked
    midpoint = max(0, int(straight_length)//2)
    result = {}
    for vh in part.getVirtualHelices():
        number = int(vh.number())
        coord = vh.coord()
        x, y = part.latticeCoordToPositionXY(*coord)
        mechanical_scale = 2.8/(2.0*float(part.radius()))
        result[number] = {
            "number": number,
            "coord": (float(x)*mechanical_scale,
                      float(y)*mechanical_scale),
            "lattice_coord": coord,
            "low": 0, "high": max(0, int(straight_length)-1),
            "forbidden": set(), "deletion_protected": set(),
            "scaffold_intervals": [(0, int(straight_length)-1)],
            "staple_intervals": [(0, int(straight_length)-1)],
            "crossovers": sorted((midpoint, partner)
                                 for partner in neighbors.get(number, ())),
            "insertions": {}}
    return result, float(part.helicalPitch())


def _capped_weighted_allocation(total, weights, capacities):
    """Distribute an integer total proportionally without exceeding caps."""
    total = int(total)
    capacities = [max(0, int(value)) for value in capacities]
    if total < 0 or total > sum(capacities):
        raise ValueError("Straight-edge common-mode capacity is insufficient.")
    values = [0]*len(capacities)
    while sum(values) < total:
        remaining = total-sum(values)
        open_indices = [index for index, capacity in enumerate(capacities)
                        if values[index] < capacity]
        if not open_indices:
            raise ValueError(
                "Straight-edge common-mode allocation exhausted its caps.")
        open_weight = sum(max(0.0, float(weights[index]))
                          for index in open_indices)
        raw = dict((index, (remaining*max(0.0, float(weights[index])) /
                            open_weight if open_weight else
                            remaining/float(len(open_indices))))
                   for index in open_indices)
        additions = dict((index, min(capacities[index]-values[index],
                                     int(math.floor(raw[index]))))
                         for index in open_indices)
        if not any(additions.values()):
            order = sorted(open_indices, key=lambda index: (
                -(raw[index]-math.floor(raw[index])), index))
            additions[order[0]] = 1
        for index in open_indices:
            values[index] += additions[index]
    return values


def _fractional_common_capacity(capacities):
    """Return the largest balanced integer total allowed by per-helix caps.

    A physical edge may have a fractional *mean* quota even though each
    caDNAno helix can only receive an integer number of indels.  For a total
    ``T`` over ``N`` helices, every helix therefore receives either
    ``floor(T/N)`` or ``ceil(T/N)`` edits.  This helper finds the largest T
    for which that one-count spread is feasible.
    """
    capacities = [max(0, int(value)) for value in capacities]
    if not capacities:
        return 0
    helix_count = len(capacities)
    for total in range(sum(capacities), -1, -1):
        base, remainder = divmod(total, helix_count)
        if all(capacity >= base for capacity in capacities) and sum(
                capacity >= base+1 for capacity in capacities) >= remainder:
            return total
    return 0


def _balanced_fractional_edge_quotas(helix_numbers, coordinates,
                                     capacities_by_edge, row_total):
    """Allocate one fractional-mean quota across one equal-edge group.

    Each returned edge row sums to ``row_total``; within a row, integer
    quotas differ by at most one.  The helices receiving the ceiling value
    are chosen to minimize the first moment about the cross-section centroid
    (the no-bending condition), then second-moment anisotropy and repeated
    use across equal physical edges.  Thus a fractional mean is implemented
    by a deterministic, symmetry-balanced pattern rather than random edits.
    """
    helix_numbers = list(map(int, helix_numbers))
    edge_count = len(capacities_by_edge)
    helix_count = len(helix_numbers)
    if not edge_count or not helix_count:
        return []
    row_total = int(row_total)
    base, remainder = divmod(row_total, helix_count)
    centred = {}
    centre_x = sum(float(coordinates[number][0])
                   for number in helix_numbers)/float(helix_count)
    centre_y = sum(float(coordinates[number][1])
                   for number in helix_numbers)/float(helix_count)
    for number in helix_numbers:
        centred[number] = (
            float(coordinates[number][0])-centre_x,
            float(coordinates[number][1])-centre_y)
    radius2 = (sum(x*x+y*y for x, y in centred.values()) /
               float(helix_count))
    radius = math.sqrt(max(radius2, 1e-12))
    extra_counts = dict((number, 0) for number in helix_numbers)
    previous_sets = []
    rows = []

    def subset_score(subset):
        subset = set(subset)
        if not subset:
            return (0.0, 0.0, 0.0, 0.0)
        x_sum = sum(centred[number][0] for number in subset)
        y_sum = sum(centred[number][1] for number in subset)
        first_moment = math.hypot(x_sum, y_sum) / (
            radius*max(1, len(subset)))
        xx_minus_yy = sum(
            centred[number][0]**2-centred[number][1]**2
            for number in subset)
        two_xy = sum(2.0*centred[number][0]*centred[number][1]
                     for number in subset)
        anisotropy = math.hypot(xx_minus_yy, two_xy) / (
            max(radius2, 1e-12)*max(1, len(subset)))
        future_counts = [extra_counts[number]+(number in subset)
                         for number in helix_numbers]
        average = sum(future_counts)/float(helix_count)
        variance = sum((value-average)**2 for value in future_counts) / \
            float(helix_count)
        repeat = sum(len(subset & old) for old in previous_sets) / \
            float(max(1, len(subset))*max(1, len(previous_sets)))
        # The first moment is the direct bending residual and is therefore
        # dominant.  The remaining terms rotate the integer rounding pattern
        # across equal edges and avoid a persistent anisotropic bias.
        return (first_moment, variance, anisotropy, repeat)

    def improve(initial, eligible):
        selected = set(initial)
        best_score = subset_score(selected)
        changed = True
        while changed:
            changed = False
            best_swap = None
            for removed in sorted(selected):
                for added in sorted(set(eligible)-selected):
                    candidate = (selected-{removed}) | {added}
                    score = subset_score(candidate)
                    if score < best_score:
                        best_score = score
                        best_swap = (removed, added)
            if best_swap is not None:
                selected.remove(best_swap[0])
                selected.add(best_swap[1])
                changed = True
        return selected, best_score

    for edge_position, capacities in enumerate(capacities_by_edge):
        capacities = list(map(int, capacities))
        if len(capacities) != helix_count or any(
                capacity < base for capacity in capacities):
            raise RuntimeError(
                "Fractional straight-edge allocation exceeds a helix cap.")
        eligible = [number for number, capacity in zip(
            helix_numbers, capacities) if capacity >= base+1]
        if len(eligible) < remainder:
            raise RuntimeError(
                "Fractional straight-edge rounding has insufficient caps.")
        if not remainder:
            selected = set()
        else:
            candidate_sets = []
            combination_count = math.comb(len(eligible), remainder)
            if combination_count <= 50000:
                candidate_sets.extend(
                    set(values) for values in itertools.combinations(
                        eligible, remainder))
            else:
                angular = sorted(eligible, key=lambda number: (
                    math.atan2(centred[number][1], centred[number][0]),
                    number))
                # Phase-shifted, angularly uniform seeds cover the symmetric
                # solutions without an exponential subset search.
                for phase in range(min(len(angular), 64)):
                    chosen = set()
                    for item in range(remainder):
                        target = (phase + (item+.5)*len(angular) /
                                  float(remainder)) % len(angular)
                        for offset in range(len(angular)):
                            candidate = angular[int(round(
                                target+offset)) % len(angular)]
                            if candidate not in chosen:
                                chosen.add(candidate)
                                break
                    candidate_sets.append(chosen)
                candidate_sets.append(set(sorted(
                    eligible, key=lambda number: (
                        extra_counts[number],
                        math.atan2(centred[number][1],
                                   centred[number][0]), number))[
                                       :remainder]))
            best = None
            for candidate in candidate_sets:
                improved, score = improve(candidate, eligible)
                key = (score, tuple(sorted(improved)))
                if best is None or key < best[0]:
                    best = (key, improved)
            selected = best[1]
        row = {}
        for number in helix_numbers:
            row[number] = base+(1 if number in selected else 0)
            extra_counts[number] += (1 if number in selected else 0)
        if sum(row.values()) != row_total or \
                max(row.values())-min(row.values()) > 1:
            raise RuntimeError(
                "Fractional straight-edge allocation lost integer balance.")
        rows.append(row)
        previous_sets.append(set(selected))
    return rows


def _apply_straight_common_mode_remove_twist(json_path, plan):
    """Cancel intrinsic twist on straight edges using common-mode indels.

    Corner curvature indels are immutable.  Only complete native domains
    wholly contained in non-bending straight edges are eligible.  A physical
    edge has one shared *mean* quota, but the mean may be fractional: each
    helix receives one of the two adjacent integers and the integer rounding
    pattern is balanced about the cross-section centroid.  This approaches
    zero twist more closely without introducing an unintended bend.
    """
    with open(json_path, "r", encoding="utf-8") as source:
        design = json.load(source)
    topology_before = topology_fingerprint(design)
    rows = {int(row["num"]): row for row in design.get("vstrands", ())}
    curvature = dict(design.get("curvature_indels") or {})
    metadata = dict(design.get("curved_metadata") or {})
    physical_edges = _straight_complete_domain_edges(plan)
    intervals = [segment for edge in physical_edges
                 for segment in edge["segments"]]
    straight_length = sum(int(edge["native_length_bp"])
                          for edge in physical_edges)
    eligible_length = sum(int(edge["eligible_length_bp"])
                          for edge in physical_edges)
    helix_numbers = sorted(rows)
    domain_size = int(plan["domain_size_bp"])
    maximum = int(plan.get(
        "maximum_indel_per_domain_allowed", MAX_INDEL_PER_DOMAIN))
    audit = {
        "enabled": True,
        "scope": "complete native domains wholly inside straight edges",
        "straight_intervals_bp": [[start, end]
                                  for start, end in intervals],
        "straight_native_bp": int(straight_length),
        "straight_eligible_complete_domain_bp": int(eligible_length),
        "physical_straight_edges": [{
            "edge_index": int(edge["edge_index"]),
            "from_vertex": int(edge["from_vertex"]),
            "to_vertex": int(edge["to_vertex"]),
            "start_coordinate_bp": float(edge["start_coordinate_bp"]),
            "native_length_bp": int(edge["native_length_bp"]),
            "eligible_length_bp": int(edge["eligible_length_bp"]),
            "segments": [[start, end]
                         for start, end in edge["segments"]]}
            for edge in physical_edges],
        "domain_size_bp": domain_size,
        "maximum_indel_per_domain_allowed": maximum,
        "preserved": ["all bend-window curvature indels",
                      "all scaffold links", "all AutoCS crossovers",
                      "all existing nicks"]}
    if not physical_edges or eligible_length <= 0 or \
            straight_length <= 0 or not helix_numbers:
        audit.update({"status": "no eligible straight domain",
                      "common_indels_per_helix": 0,
                      "total_new_indels": 0})
    else:
        helix_data, helical_pitch = _frame_twist_helix_data(
            design, plan["lattice"], straight_length)
        region = {"helices": helix_numbers, "start": 0,
                  "end": straight_length-1}
        baseline = estimate_global_twist(
            region, helix_data, helical_pitch, extra_base_delta=0)
        baseline_rate = float(baseline["twist_per_base_deg"])
        direction = (-1 if baseline_rate > 1e-12 else
                     1 if baseline_rate < -1e-12 else 0)
        deletion_protected = short_staple_deletion_protection(rows)

        edge_helix_capacity = {}
        edge_balanced_capacity = []
        edge_capacities_by_helix = []
        candidate_cache = {}
        load_cache = {}
        for helix in helix_numbers:
            row = rows[helix]
            loads = {}
            nominal = min(int(plan["nominal_perimeter_bp"]),
                          len(row.get("loop", ())),
                          len(row.get("skip", ())))
            for index in range(nominal):
                value = int(row["loop"][index])+int(row["skip"][index])
                if value:
                    domain = index//domain_size
                    loads[domain] = loads.get(domain, 0)+value
            load_cache[helix] = dict(loads)
            capacities = []
            for edge_position, edge in enumerate(physical_edges):
                edge_candidates = []
                edge_domain_caps = {}
                for start, end in edge["segments"]:
                    candidates = [index for index in _safe_sites(
                        row, helix, start, end+1)
                        if int(row["loop"][index]) == 0 and
                        int(row["skip"][index]) == 0]
                    if direction < 0:
                        protected = deletion_protected.get(helix, set())
                        candidates = [index for index in candidates
                                      if index not in protected]
                    for index in candidates:
                        domain = index//domain_size
                        local = ((float(index)+.5-
                                  float(edge["start_coordinate_bp"])) %
                                 float(plan["nominal_perimeter_bp"]))
                        if local < float(edge["native_length_bp"])+1e-9:
                            edge_candidates.append((index, domain, local))
                for domain in set(item[1] for item in edge_candidates):
                    sites = [item for item in edge_candidates
                             if item[1] == domain]
                    load = int(loads.get(domain, 0))
                    room = (maximum-load if direction > 0 else load+maximum)
                    edge_domain_caps[domain] = min(
                        len(sites), max(0, room))
                candidate_cache[(helix, edge_position)] = (
                    edge_candidates, edge_domain_caps)
                capacities.append(sum(edge_domain_caps.values()))
            edge_helix_capacity[helix] = capacities
        for edge_position in range(len(physical_edges)):
            capacities = [edge_helix_capacity[helix][edge_position]
                          for helix in helix_numbers]
            edge_capacities_by_helix.append(capacities)
            edge_balanced_capacity.append(
                _fractional_common_capacity(capacities))

        # Equal physical edges receive the same total over the cross-section,
        # hence the same possibly fractional mean quota per helix.  Search the
        # reachable *global* totals because the calibrated predictor depends
        # on the complete connected cross-section; predicting one isolated
        # edge can disable its connectivity calibration and give a false
        # positive result.  For each reachable total, retain the least-
        # variance length-normalized allocation.
        length_groups = {}
        for edge_position, edge in enumerate(physical_edges):
            length_groups.setdefault(
                int(edge["native_length_bp"]), []).append(edge_position)
        edge_cross_section_totals = [0]*len(physical_edges)
        handedness_guard_applied = False
        group_predictions = []
        if direction:
            groups = []
            for native_length, positions in sorted(length_groups.items()):
                group_capacity = min(edge_balanced_capacity[position]
                                     for position in positions)
                groups.append((native_length, positions, group_capacity))

            # total -> (density-variance cost, quotas by length group)
            states = {0: (0.0, [])}
            for native_length, positions, group_capacity in groups:
                edge_count = len(positions)
                updated = {}
                for previous_total, (previous_cost,
                                     previous_quotas) in states.items():
                    for cross_section_total in range(group_capacity+1):
                        total = (previous_total + edge_count *
                                 cross_section_total)
                        mean_quota = (cross_section_total /
                                      float(len(helix_numbers)))
                        cost = (previous_cost + edge_count*mean_quota**2 /
                                float(max(1, native_length)))
                        prior = updated.get(total)
                        candidate = (cost, previous_quotas+
                                     [cross_section_total])
                        if prior is None or candidate[0] < prior[0]-1e-12:
                            updated[total] = candidate
                states = updated

            predictions = {}
            for total in states:
                predictions[total] = estimate_global_twist(
                    region, helix_data, helical_pitch,
                    extra_base_delta=direction*total)
            unconstrained_total = min(states, key=lambda total: (
                abs(float(predictions[total]["twist_per_base_deg"])),
                states[total][0], total))
            same_side_totals = [
                total for total in states
                if ((baseline_rate > 0 and float(predictions[total][
                    "twist_per_base_deg"]) > 1e-12) or
                    (baseline_rate < 0 and float(predictions[total][
                    "twist_per_base_deg"]) < -1e-12))]
            selected_total = min(
                same_side_totals or [0], key=lambda total: (
                    abs(float(predictions[total]["twist_per_base_deg"])),
                    states[total][0], total))
            selected_group_totals = states[selected_total][1]
            unconstrained_group_totals = states[unconstrained_total][1]
            handedness_guard_applied = (
                selected_total != unconstrained_total)
            for group_index, (native_length, positions,
                              group_capacity) in enumerate(groups):
                cross_section_total = int(
                    selected_group_totals[group_index])
                unconstrained = int(
                    unconstrained_group_totals[group_index])
                for position in positions:
                    edge_cross_section_totals[position] = \
                        cross_section_total
                group_predictions.append({
                    "native_length_bp": int(native_length),
                    "edge_positions": list(map(int, positions)),
                    "mean_quota_per_helix": (
                        cross_section_total/float(len(helix_numbers))),
                    "cross_section_total_per_edge": cross_section_total,
                    "unconstrained_mean_quota_per_helix": (
                        unconstrained/float(len(helix_numbers))),
                    "balanced_capacity_per_edge": int(group_capacity),
                    "global_final_prediction": predictions[selected_total]})

        coordinates = dict((helix, helix_data[helix]["coord"])
                           for helix in helix_numbers)
        edge_helix_quotas = [dict((helix, 0) for helix in helix_numbers)
                             for unused_edge in physical_edges]
        for unused_length, positions in sorted(length_groups.items()):
            cross_section_total = edge_cross_section_totals[positions[0]]
            rows_for_group = _balanced_fractional_edge_quotas(
                helix_numbers, coordinates,
                [edge_capacities_by_helix[position]
                 for position in positions], cross_section_total)
            for position, row_quotas in zip(positions, rows_for_group):
                edge_helix_quotas[position] = row_quotas
        edge_mean_quotas = [
            total/float(len(helix_numbers))
            for total in edge_cross_section_totals]
        per_helix_counts = dict((helix, sum(
            edge_helix_quotas[position][helix]
            for position in range(len(physical_edges))))
            for helix in helix_numbers)
        total_selected = sum(per_helix_counts.values())
        mean_count = total_selected/float(len(helix_numbers))

        def select_edge_sites(helix, edge_position, quota):
            """Select equal-partition sites in one circular physical edge."""
            if not quota:
                return []
            candidates, capacities = candidate_cache[
                (helix, edge_position)]
            edge_length = float(
                physical_edges[edge_position]["native_length_bp"])
            targets = [(index+.5)*edge_length/float(quota)
                       for index in range(quota)]
            used_sites = set()
            used_by_domain = {}
            selected = []
            for target in targets:
                available = [item for item in candidates
                             if item[0] not in used_sites and
                             used_by_domain.get(item[1], 0) <
                             int(capacities.get(item[1], 0))]
                if not available:
                    raise RuntimeError(
                        "Straight-edge physical-edge quota exhausted its "
                        "safe sites.")
                item = min(available, key=lambda value: (
                    abs(float(value[2])-target), value[1], value[0]))
                used_sites.add(item[0])
                used_by_domain[item[1]] = \
                    used_by_domain.get(item[1], 0)+1
                selected.append(item)
            return selected

        selected_by_helix = {}
        for helix in helix_numbers:
            row = rows[helix]
            selected = []
            loads = load_cache[helix]
            for edge_position, edge_quota in enumerate(edge_helix_quotas):
                quota = int(edge_quota[helix])
                if not quota:
                    continue
                records = select_edge_sites(helix, edge_position, quota)
                for index, domain, unused_local in records:
                    index = int(index)
                    domain = int(domain)
                    if direction > 0:
                        row["loop"][index] += 1
                    else:
                        row["skip"][index] = -1
                    loads[domain] = int(loads.get(domain, 0))+direction
                    if abs(loads[domain]) > maximum:
                        raise RuntimeError(
                            "Straight-edge Remove Twist exceeded the hard "
                            "indel/domain limit.")
                    selected.append(index)
            if len(selected) != per_helix_counts[helix]:
                raise RuntimeError(
                    "Straight-edge Remove Twist lost its balanced quota.")
            selected_by_helix[helix] = sorted(selected)

        interval_quotas = []
        for edge_position, edge in enumerate(physical_edges):
            for start, end in edge["segments"]:
                interval_quotas.append(sum(
                    start <= index <= end
                    for sites in selected_by_helix.values()
                    for index in sites)/float(len(helix_numbers)))

        final_prediction = estimate_global_twist(
            region, helix_data, helical_pitch,
            extra_base_delta=direction*total_selected)
        per_edge_ranges = []
        first_moment_residuals = []
        centre_x = sum(coordinates[number][0]
                       for number in helix_numbers)/float(len(helix_numbers))
        centre_y = sum(coordinates[number][1]
                       for number in helix_numbers)/float(len(helix_numbers))
        radius = math.sqrt(max(1e-12, sum(
            (coordinates[number][0]-centre_x)**2+
            (coordinates[number][1]-centre_y)**2
            for number in helix_numbers)/float(len(helix_numbers))))
        for row_quotas in edge_helix_quotas:
            values = list(row_quotas.values())
            mean = sum(values)/float(len(values))
            x_moment = sum((row_quotas[number]-mean) *
                           (coordinates[number][0]-centre_x)
                           for number in helix_numbers)
            y_moment = sum((row_quotas[number]-mean) *
                           (coordinates[number][1]-centre_y)
                           for number in helix_numbers)
            per_edge_ranges.append([min(values), max(values)])
            first_moment_residuals.append(
                math.hypot(x_moment, y_moment) /
                (radius*max(1, sum(values))))
        audit.update({
            "status": "applied" if total_selected else
                      "baseline already optimal or no legal improvement",
            "baseline_prediction": baseline,
            "final_prediction": final_prediction,
            "correction": ("insertion" if direction > 0 else
                           "deletion" if direction < 0 else "none"),
            "signed_mean_indel_per_helix": direction*mean_count,
            "signed_indel_by_helix": dict((str(helix),
                direction*per_helix_counts[helix])
                for helix in helix_numbers),
            "common_indels_per_helix": mean_count,
            "per_helix_indel_count_range": [
                min(per_helix_counts.values()),
                max(per_helix_counts.values())],
            "edge_balanced_capacity": edge_balanced_capacity,
            "edge_cross_section_total": edge_cross_section_totals,
            "edge_common_quota": edge_mean_quotas,
            "edge_helix_quota": [{str(helix): int(row[helix])
                                   for helix in helix_numbers}
                                  for row in edge_helix_quotas],
            "edge_helix_quota_range": per_edge_ranges,
            "edge_first_moment_residual": first_moment_residuals,
            "interval_common_quota": interval_quotas,
            "balanced_capacity_total": int(sum(edge_balanced_capacity)),
            "handedness_guard": (
                "retain original handedness; final twist must remain "
                "strictly on the original side of zero"),
            "handedness_guard_applied": bool(handedness_guard_applied),
            "equal_length_edge_groups": group_predictions,
            "total_new_indels": int(total_selected),
            "per_helix_sites": {str(helix): selected_by_helix[helix]
                                for helix in helix_numbers},
            "distribution_method":
                ("one circular physical-edge identity across base zero; "
                 "equal physical lengths receive equal fractional mean "
                 "quotas; per-edge integer helix quotas differ by at most "
                 "one and minimize the cross-section first moment; rounding "
                 "patterns rotate across equal edges; equal axial "
                 "partitions; nearest safe site inside an eligible absolute "
                 "native domain; no handedness overshoot")})

        record_by_helix = {int(record["helix"]): record
                           for record in curvature.get("rings", ())}
        for helix, sites in selected_by_helix.items():
            record = record_by_helix.get(helix)
            if record is None:
                continue
            key = "insertions" if direction > 0 else "deletions"
            record[key] = sorted(set(map(int, record.get(key, ()))) |
                                 set(sites))
            record["straight_common_mode_remove_twist_sites"] = sites
            record["straight_common_mode_remove_twist_signed_count"] = \
                int(direction*len(sites))
            record["target_bases"] = (int(record["nominal_bases"])+
                                      len(record.get("insertions", ())) -
                                      len(record.get("deletions", ())))
            final_loads = []
            row = rows[helix]
            for domain in range(int(math.ceil(
                    int(record["nominal_bases"])/float(domain_size)))):
                first = domain*domain_size
                last = min(first+domain_size,
                           int(record["nominal_bases"]))
                final_loads.append(sum(
                    int(row["loop"][index])+int(row["skip"][index])
                    for index in range(first, last)))
            record["final_domain_indel_loads"] = final_loads
            record["maximum_indel_in_one_domain"] = max(
                [abs(value) for value in final_loads] or [0])

        all_loads = [abs(value)
                     for record in curvature.get("rings", ())
                     for value in record.get("final_domain_indel_loads", ())]
        curvature["maximum_insertion_per_domain"] = max(
            [value for record in curvature.get("rings", ())
             for value in record.get("final_domain_indel_loads", ())
             if value > 0] or [0])
        curvature["maximum_deletion_per_domain"] = max(
            [-value for record in curvature.get("rings", ())
             for value in record.get("final_domain_indel_loads", ())
             if value < 0] or [0])
        curvature["maximum_abs_indel_per_domain"] = max(all_loads or [0])

    curvature["frame_straight_common_mode_remove_twist"] = audit
    metadata["frame_straight_common_mode_remove_twist"] = audit
    metadata["curvature_encoding"] = \
        "frame-corner-indels+straight-common-mode-remove-twist"
    design["curvature_indels"] = curvature
    design["curved_metadata"] = metadata
    if topology_fingerprint(design) != topology_before:
        raise RuntimeError(
            "Straight-edge Remove Twist attempted to alter AutoCS topology.")
    with open(json_path, "w", encoding="utf-8", newline="\n") as output:
        json.dump(design, output, separators=(",", ":"))
    return curvature, metadata, audit


def _reoptimize_frame_staple_nicks(json_path, lattice):
    """Repartition Frame staples after vertex-local indels are finalized.

    Curved Design initially breaks staples while its indels still occupy the
    full loop.  Frame subsequently moves that fixed indel budget into corner
    windows, so those old nicks no longer optimize the final actual lengths.
    Rejoin only same-helix nick boundaries and rerun Autobreak against the
    final indels.  AutoCS crossovers are immutable throughout this pass.
    """
    from cadnano2.model.document import Document
    from cadnano2.model.enum import LatticeType
    from cadnano2.model.io.legacydecoder import import_legacy_dict
    from cadnano2.model.io.legacyencoder import legacy_dict_from_part
    from cadnano2.model.parts.part import (_existingStapleNickBoundaries,
                                           _stapleOligoBaseRecords)

    with open(json_path, "r", encoding="utf-8") as source:
        original = json.load(source)
    before_crossovers = _crossover_fingerprint(original)
    lattice_type = (LatticeType.Honeycomb if lattice == "honeycomb"
                    else LatticeType.Square)
    document = Document()
    part = import_legacy_dict(
        document, original, lattice_type, forceLatticeType=True)
    if part is None:
        raise RuntimeError("Cannot reload Frame design for staple nicks.")

    def staple_state():
        state = []
        for oligo in part.oligos():
            if not oligo.isStaple() or oligo.isHybrid() or \
                    oligo.strand5p() is None:
                continue
            records = _stapleOligoBaseRecords(oligo)
            state.append((oligo, int(oligo.actualLength()),
                          sum(int(record[3]) == 0 for record in records)))
        return state

    initial_staple_state = staple_state()
    initial_normal_below_preferred = sum(
        deletion_count < 2 and 21 <= length < 30
        for unused_oligo, length, deletion_count in initial_staple_state)
    used_heal_boundaries = set()

    def heal_target_nicks():
        """Join only products that are too short after final Frame indels.

        Autobreak enforces the normal 21--57 nt range with a 30--50 nt
        preference and uses 58--64 nt only as a proven no-solution exception.
        Frame starts from an already broken Curved parent, however, so a
        pre-existing legal 21--29 nt product would otherwise never be
        reconsidered.  Join that product to one same-helix neighbour and let
        the next Autobreak pass repartition only the local component.  A
        deletion-dense product below 40 nt uses the same local repair.
        """
        healed = 0
        while True:
            target_details = {
                oligo: (length, deletion_count)
                for oligo, length, deletion_count in staple_state()
                if (length < 30 or
                    (deletion_count >= 2 and length < 40))}
            targets = set(target_details)
            if not targets:
                break
            candidates = []
            for helix, upper_index in sorted(
                    _existingStapleNickBoundaries(part)):
                if (helix, upper_index) in used_heal_boundaries:
                    continue
                vh = part.virtualHelix(helix)
                strand_set = vh.stapleStrandSet()
                lower = strand_set.getStrand(upper_index-1)
                upper = strand_set.getStrand(upper_index)
                if lower is None or upper is None or lower is upper or \
                        not strand_set.strandsCanBeMerged(lower, upper):
                    continue
                lower_oligo, upper_oligo = lower.oligo(), upper.oligo()
                if lower_oligo is upper_oligo or not (
                        lower_oligo in targets or upper_oligo in targets):
                    continue
                combined = (int(lower_oligo.actualLength()) +
                            int(upper_oligo.actualLength()))
                target_strand = (lower if lower_oligo in targets else upper)
                other = upper if target_strand is lower else lower
                target_oligo = target_strand.oligo()
                unused_length, deletion_count = target_details[target_oligo]
                if deletion_count >= 2:
                    # A deletion-weakened staple is best repaired directly
                    # into 40--60 nt.  A longer component remains acceptable
                    # because the next strict Autobreak pass can repartition
                    # it without changing any crossover.
                    if 40 <= combined <= 60:
                        repair_class, deviation = 0, abs(combined-50)
                    elif 80 <= combined <= 120:
                        repair_class, deviation = 1, abs(combined-100)
                    elif 61 <= combined < 80:
                        repair_class, deviation = 2, abs(combined-60)
                    elif combined < 40:
                        repair_class, deviation = 3, 40-combined
                    else:
                        repair_class, deviation = 4, combined-120
                    score = (repair_class, deviation, combined,
                             helix, upper_index)
                else:
                    # For an ordinary 21--29 nt product, first prefer a
                    # direct 30--50 nt merge.  Next prefer a 60--100 nt local
                    # component, because the following planner can divide it
                    # into two 30--50 nt products.  Retaining 51--57 nt is a
                    # legal but less preferred fallback.
                    if 30 <= combined <= 50:
                        repair_class, deviation = 0, abs(combined-40)
                    elif 60 <= combined <= 100:
                        repair_class, deviation = 1, abs(combined-80)
                    elif 51 <= combined <= 57:
                        repair_class, deviation = 2, abs(combined-50)
                    elif combined < 30:
                        repair_class, deviation = 3, 30-combined
                    else:
                        repair_class, deviation = 4, abs(combined-80)
                    score = (repair_class, deviation, combined,
                             helix, upper_index)
                candidates.append((score, strand_set, target_strand, other,
                                   (helix, upper_index)))
            if not candidates:
                break
            unused_score, strand_set, target_strand, other, boundary = min(
                candidates, key=lambda item: item[0])
            strand_set.mergeStrands(
                target_strand, other, useUndoStack=False)
            used_heal_boundaries.add(boundary)
            healed += 1
        return healed

    # Re-evaluate the final actual lengths without globally joining all old
    # nicks.  Each iteration can (a) split a newly-long component and (b)
    # locally heal one short/deletion-weakened product.  A following pass then
    # repartitions only that local component.  The loop ends when neither
    # operation changes the design.
    merged = total_nicks = total_skipped = 0
    pass_count = 0
    for pass_count in range(1, 21):
        part._autobreakStaplesApplied = False
        result = part.autoBreakStaples(
            preserveCrossovers=True, markUnbreakable=True,
            preferDeletionDense=True)
        created = int(result.get("nicks", 0))
        total_nicks += created
        total_skipped += int(result.get("skipped", 0))
        healed = heal_target_nicks()
        merged += healed
        if not created and not healed:
            break
    lengths = []
    deletion_dense_lengths = []
    for oligo in part.oligos():
        if not oligo.isStaple() or oligo.isHybrid() or \
                oligo.strand5p() is None:
            continue
        records = _stapleOligoBaseRecords(oligo)
        length = int(oligo.actualLength())
        lengths.append(length)
        if sum(int(record[3]) == 0 for record in records) >= 2:
            deletion_dense_lengths.append(length)

    encoded = legacy_dict_from_part(
        part, os.path.basename(json_path), includeSequences=False)
    for key, value in original.items():
        if key not in ("name", "num_bases", "vstrands"):
            encoded[key] = value
    after_crossovers = _crossover_fingerprint(encoded)
    if after_crossovers != before_crossovers:
        raise RuntimeError(
            "Frame staple nick optimization changed AutoCS crossovers.")
    audit = {
        "method": "iterative targeted post-Frame-indel nick optimization",
        "passes": int(pass_count),
        "healed_target_nicks": int(merged),
        "created_nicks": int(total_nicks),
        "unbreakable_staples": sum(
            value < 21 or value > 64 for value in lengths),
        "optimizer_skipped_attempts": int(total_skipped),
        "minimum_staple_nt": min(lengths or [0]),
        "maximum_staple_nt": max(lengths or [0]),
        "normal_staple_range_nt": [21, 57],
        "normal_preferred_range_nt": [30, 50],
        "normal_below_preferred_before_repartition": int(
            initial_normal_below_preferred),
        "normal_below_preferred_after_repartition": sum(
            21 <= value < 30 for value in lengths),
        "deletion_dense_preferred_range_nt": [40, 60],
        "deletion_dense_staples": len(deletion_dense_lengths),
        "deletion_dense_in_40_60": sum(
            40 <= value <= 60 for value in deletion_dense_lengths),
        "deletion_dense_outside_40_60_nt": [
            value for value in deletion_dense_lengths
            if not 40 <= value <= 60],
        "deletion_dense_lengths_nt": deletion_dense_lengths,
        "preserved": ["all scaffold links", "all AutoCS crossovers",
                      "all indels"]}
    encoded.setdefault("curved_metadata", {})[
        "frame_staple_nick_optimization"] = audit
    # The final nick partition is intentionally newer than the Curved-stage
    # topology fingerprint, while every inter-helix AutoCS crossover remains
    # bit-for-bit unchanged.  Publish both facts explicitly.
    encoded["curved_metadata"]["autocs_topology_fingerprint"] = \
        topology_fingerprint(encoded)
    encoded["curved_metadata"]["autocs_crossover_fingerprint"] = \
        hashlib.sha256(repr(sorted(after_crossovers)).encode("utf-8")).hexdigest()
    encoded.setdefault("curvature_indels", {})[
        "frame_staple_nick_optimization"] = audit
    with open(json_path, "w", encoding="utf-8", newline="\n") as output:
        json.dump(encoded, output, separators=(",", ":"))
    return audit, dict(encoded.get("curved_metadata") or {})


def _rounded_segments(plan):
    vertices = plan["vertices_nm"]
    tangencies = plan["tangent_trim_nm"]
    radii = plan["bend_radius_nm"]
    segments = []
    for index, vertex in enumerate(vertices):
        previous = vertices[index - 1]
        following = vertices[(index + 1) % len(vertices)]
        incoming = _unit((vertex[0] - previous[0],
                          vertex[1] - previous[1]))
        outgoing = _unit((following[0] - vertex[0],
                          following[1] - vertex[1]))
        entry = (vertex[0] - incoming[0] * tangencies[index],
                 vertex[1] - incoming[1] * tangencies[index])
        exit_point = (vertex[0] + outgoing[0] * tangencies[index],
                      vertex[1] + outgoing[1] * tangencies[index])
        left_normal = (-incoming[1], incoming[0])
        centre = (entry[0] + left_normal[0] * radii[index],
                  entry[1] + left_normal[1] * radii[index])
        start_angle = math.atan2(entry[1] - centre[1],
                                 entry[0] - centre[0])
        segments.append({"entry": entry, "exit": exit_point,
                         "centre": centre, "radius": radii[index],
                         "start_angle": start_angle,
                         "turn": math.radians(
                             plan["turn_angles_degrees"][index])})
    path = []
    for index, corner in enumerate(segments):
        previous = segments[index - 1]
        straight_start, straight_end = previous["exit"], corner["entry"]
        straight_length = _distance(straight_start, straight_end)
        if straight_length > 1e-9:
            path.append({"kind": "line", "start": straight_start,
                         "end": straight_end, "length": straight_length})
        path.append({"kind": "arc", "centre": corner["centre"],
                     "radius": corner["radius"],
                     "start_angle": corner["start_angle"],
                     "turn": corner["turn"],
                     "length": corner["radius"] * corner["turn"]})
    return path


def _point_on_path(path, fraction):
    total = sum(segment["length"] for segment in path)
    distance = (fraction % 1.0) * total
    for segment in path:
        if distance <= segment["length"] + 1e-12:
            local = 0.0 if not segment["length"] else \
                distance / segment["length"]
            if segment["kind"] == "line":
                start, end = segment["start"], segment["end"]
                point = (start[0] + (end[0] - start[0]) * local,
                         start[1] + (end[1] - start[1]) * local)
                tangent = _unit((end[0] - start[0], end[1] - start[1]))
            else:
                angle = segment["start_angle"] + segment["turn"] * local
                point = (segment["centre"][0] +
                         segment["radius"] * math.cos(angle),
                         segment["centre"][1] +
                         segment["radius"] * math.sin(angle))
                tangent = (-math.sin(angle), math.cos(angle))
            return point, tangent
        distance -= segment["length"]
    return _point_on_path(path, 0.0)


def _sample_frame(plan, sample_count=None):
    """Sample the neutral rounded-polygon centreline in physical units."""
    path = _rounded_segments(plan)
    sample_count = int(sample_count or max(
        128, len(plan["vertices_nm"]) * 64))
    return [_point_on_path(path, index / float(sample_count))
            for index in range(sample_count)]


def _frame_stl_text(plan, rings):
    """Return a closed rectangular tube following the true frame curve."""
    samples = _sample_frame(plan)
    normal_offsets = [float(row.get("frame_normal_offset_nm", 0.0))
                      for row in rings]
    vertical_offsets = [float(row.get("frame_binormal_offset_nm", 0.0))
                        for row in rings]
    half_helix = DNA_HELIX_RADIUS_NM
    low_normal = min(normal_offsets) - half_helix
    high_normal = max(normal_offsets) + half_helix
    low_z = min(vertical_offsets) - half_helix
    high_z = max(vertical_offsets) + half_helix

    def corners(sample):
        point, tangent = sample
        outward = (tangent[1], -tangent[0])
        return [
            (point[0] + outward[0] * low_normal,
             point[1] + outward[1] * low_normal, low_z),
            (point[0] + outward[0] * high_normal,
             point[1] + outward[1] * high_normal, low_z),
            (point[0] + outward[0] * high_normal,
             point[1] + outward[1] * high_normal, high_z),
            (point[0] + outward[0] * low_normal,
             point[1] + outward[1] * low_normal, high_z)]

    sections = [corners(sample) for sample in samples]
    lines = ["solid cadnano_frame_design"]
    for index, current in enumerate(sections):
        following = sections[(index + 1) % len(sections)]
        for edge in range(4):
            next_edge = (edge + 1) % 4
            for triangle in ((current[edge], following[edge],
                              current[next_edge]),
                             (following[edge], following[next_edge],
                              current[next_edge])):
                lines.extend(("  facet normal 0 0 0", "    outer loop"))
                for point in triangle:
                    lines.append("      vertex %.8f %.8f %.8f" % point)
                lines.extend(("    endloop", "  endfacet"))
    lines.extend(("endsolid cadnano_frame_design", ""))
    return "\n".join(lines)


def _png_chunk(kind, data):
    return (struct.pack(">I", len(data)) + kind + data +
            struct.pack(">I", zlib.crc32(kind + data) & 0xffffffff))


def _frame_preview_png(plan, rings, width=900, height=620):
    """Render a dependency-free top-view preview of every frame helix."""
    background = (248, 249, 251)
    pixels = [bytearray(background * width) for unused in range(height)]
    samples = _sample_frame(plan)
    curves = []
    for ring in rings:
        offset = float(ring.get("frame_normal_offset_nm", 0.0))
        points = []
        for point, tangent in samples:
            outward = (tangent[1], -tangent[0])
            points.append((point[0] + outward[0] * offset,
                           point[1] + outward[1] * offset))
        curves.append(points)
    all_points = [point for curve in curves for point in curve]
    minimum_x = min(point[0] for point in all_points)
    maximum_x = max(point[0] for point in all_points)
    minimum_y = min(point[1] for point in all_points)
    maximum_y = max(point[1] for point in all_points)
    span_x = max(1.0, maximum_x - minimum_x)
    span_y = max(1.0, maximum_y - minimum_y)
    scale = min((width - 80.0) / span_x, (height - 80.0) / span_y)
    centre_x = 0.5 * (minimum_x + maximum_x)
    centre_y = 0.5 * (minimum_y + maximum_y)

    def set_pixel(x_value, y_value, color):
        x_value, y_value = int(round(x_value)), int(round(y_value))
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                x, y = x_value + dx, y_value + dy
                if 0 <= x < width and 0 <= y < height:
                    offset = x * 3
                    pixels[y][offset:offset + 3] = bytes(color)

    for index, curve in enumerate(curves):
        color = ((43, 105, 174) if index % 2 == 0 else
                 (104, 151, 197))
        for point in curve:
            x_value = width / 2.0 + (point[0] - centre_x) * scale
            y_value = height / 2.0 - (point[1] - centre_y) * scale
            set_pixel(x_value, y_value, color)
    raw = b"".join(b"\x00" + bytes(row) for row in pixels)
    return (b"\x89PNG\r\n\x1a\n" +
            _png_chunk(b"IHDR", struct.pack(
                ">IIBBBBB", width, height, 8, 2, 0, 0, 0)) +
            _png_chunk(b"IDAT", zlib.compress(raw, 9)) +
            _png_chunk(b"IEND", b""))


def _rewrite_frame_artifacts(result, plan):
    """Replace provisional cylinder artifacts with true frame artifacts."""
    paths = list(result.get("input_paths") or ())
    stl_path = next((path for path in paths if path.endswith("_shape.stl")),
                    None)
    modules_path = next((path for path in paths
                         if path.endswith("_modules.csv")), None)
    preview_path = next((path for path in paths
                         if path.endswith("_preview.png")), None)
    rings = plan["rings"]
    if stl_path:
        with open(stl_path, "w", encoding="ascii", newline="\n") as output:
            output.write(_frame_stl_text(plan, rings))
    if modules_path:
        with open(modules_path, "w", encoding="utf-8",
                  newline="\n") as output:
            output.write(
                "index,layer,slice,frame_normal_offset_nm,"
                "frame_binormal_offset_nm,nominal_bp,target_bp,"
                "vertex_indels\n")
            for ring in rings:
                vertex_indels = ";".join(
                    str(value) for value in ring["frame_vertex_indels"])
                output.write(
                    "{index},{layer},{slice},{normal:.6f},{vertical:.6f},"
                    "{nominal},{target},{indels}\n".format(
                        index=int(ring["index"]), layer=int(ring["layer"]),
                        slice=int(ring["slice"]),
                        normal=float(ring["frame_normal_offset_nm"]),
                        vertical=float(ring["frame_binormal_offset_nm"]),
                        nominal=int(plan["nominal_perimeter_bp"]),
                        target=int(ring["frame_target_bp"]),
                        indels=vertex_indels))
    if preview_path:
        with open(preview_path, "wb") as output:
            output.write(_frame_preview_png(plan, rings))


def _frame_geometry(geometry, plan, ring_rows):
    if not geometry:
        return geometry
    path = _rounded_segments(plan)
    nominal = float(plan["nominal_perimeter_bp"])
    origin_shift = float(plan.get("native_origin_shift_bp", 0.0))
    offsets = dict((int(row["index"]), (
        float(row.get("frame_normal_offset_nm", 0.0)),
        float(row.get("frame_binormal_offset_nm", 0.0))))
                   for row in ring_rows)
    # Keep an explicit duplex centre axis alongside the nucleotide frames.
    # The nucleotide coordinates intentionally retain their helical phase;
    # deriving a BILD cylinder from those coordinates would therefore draw a
    # corkscrew instead of the requested smooth dsDNA rod.  These positions
    # contain only the rounded-polygon centreline plus the helix cross-section
    # offset, in the same oxDNA units as the stored nucleotide frames.
    helix_axes = {}
    for key, frame in geometry.get("frames", {}).items():
        unused_type, helix_text, base_text = key.split(":")
        helix, base = int(helix_text), int(base_text)
        point, tangent2 = _point_on_path(
            path, (base + 0.5 + origin_shift) / nominal)
        outward = (tangent2[1], -tangent2[0], 0.0)
        tangent = (tangent2[0], tangent2[1], 0.0)
        normal_offset, z_offset = offsets.get(helix, (0.0, 0.0))
        centre = ((point[0] + outward[0] * normal_offset) /
                  OXDNA_LENGTH_NM,
                  (point[1] + outward[1] * normal_offset) /
                  OXDNA_LENGTH_NM,
                  z_offset / OXDNA_LENGTH_NM)
        helix_axes["%d:%d" % (helix, base)] = list(centre)

        # Preserve the duplex phase and paired-strand displacement from the
        # circular DNAxiS frame by remapping its local radial/z/tangent basis.
        old_pos = list(map(float, frame["pos"]))
        radial_length = math.hypot(old_pos[0], old_pos[1])
        old_out = ((old_pos[0] / radial_length if radial_length else 1.0),
                   (old_pos[1] / radial_length if radial_length else 0.0),
                   0.0)
        old_tangent = (-old_out[1], old_out[0], 0.0)
        ring = ring_rows[helix]
        old_centre = (float(ring.get("geometry_radius_nm",
                                     ring["radius_nm"])) /
                      OXDNA_LENGTH_NM * old_out[0],
                      float(ring.get("geometry_radius_nm",
                                     ring["radius_nm"])) /
                      OXDNA_LENGTH_NM * old_out[1],
                      float(ring["height_nm"]) / OXDNA_LENGTH_NM)
        displacement = [old_pos[axis] - old_centre[axis]
                        for axis in range(3)]

        def remap(vector):
            radial = sum(vector[axis] * old_out[axis]
                         for axis in range(3))
            longitudinal = sum(vector[axis] * old_tangent[axis]
                               for axis in range(3))
            vertical = vector[2]
            mapped = [radial * outward[axis] +
                      longitudinal * tangent[axis]
                      for axis in range(3)]
            mapped[2] += vertical
            return mapped

        moved = remap(displacement)
        frame["pos"] = [centre[axis] + moved[axis] for axis in range(3)]
        frame["a1"] = remap(list(map(float, frame["a1"])))
        frame["a3"] = remap(list(map(float, frame["a3"])))
    geometry["source"] = "Frame Design rounded-polygon target geometry"
    geometry["frame_geometry"] = {
        "vertices_nm": plan["vertices_nm"],
        "bend_radius_nm": plan["bend_radius_nm"],
        "bend_length_bp": plan["bend_length_bp"],
        "helix_axes": helix_axes}
    return geometry


def _localize_corner_indels(json_path, plan):
    with open(json_path, "r", encoding="utf-8") as source:
        design = json.load(source)
    topology_before = topology_fingerprint(design)
    rows = dict((int(row["num"]), row) for row in design["vstrands"])
    curvature = dict(design.get("curvature_indels") or {})
    records = [dict(row) for row in curvature.get("rings", ())]
    domain_size = int(plan["domain_size_bp"])
    maximum_per_domain = int(plan.get(
        "maximum_indel_per_domain_allowed", MAX_INDEL_PER_DOMAIN))
    centres = plan["vertex_native_centres"]
    bend_lengths = plan["bend_length_bp"]
    ring_plan = dict((int(item["index"]), item)
                     for item in plan["rings"])
    reports = []
    # Remove the full-loop Curved allocation first.  Frame preserves the
    # exact signed budget but relocates it only into vertex bend windows.
    for record in records:
        helix = int(record["helix"])
        row = rows[helix]
        for index in record.get("insertions", ()):
            row["loop"][int(index)] = 0
        for index in record.get("deletions", ()):
            row["skip"][int(index)] = 0
    deletion_protected = short_staple_deletion_protection(rows)
    for record in records:
        helix = int(record["helix"])
        row = rows[helix]
        nominal_bases = int(record["nominal_bases"])
        if nominal_bases != int(plan["nominal_perimeter_bp"]):
            raise RuntimeError(
                "Frame parent length differs from the Curved parent: "
                "%d != %d bp." %
                (nominal_bases, int(plan["nominal_perimeter_bp"])))
        # The target residual comes from the explicit Frame vertex plan.  Its
        # normal coordinate is the same radial bend coordinate used by
        # Curved Design; only the axial positions are localized to vertices.
        residual = int(ring_plan[helix].get("frame_total_indel", 0))
        # plan_frame already rounds every angle*normal-offset contribution.
        # Re-distributing only the total by largest remainder changes
        # unequal-polygon vertex angles and can reverse the intended local
        # bend.  Preserve the exact per-vertex plan instead.
        by_vertex = list(map(int, ring_plan[helix].get(
            "frame_vertex_indels", ())))
        if len(by_vertex) != len(centres) or sum(by_vertex) != residual:
            raise RuntimeError(
                "Frame per-vertex indel plan is incomplete or inconsistent.")
        inserted, deleted = [], []
        loads = {}
        vertex_reports = []
        for vertex, (count, centre, length) in enumerate(
                zip(by_vertex, centres, bend_lengths)):
            start, end = centre - length / 2.0, centre + length / 2.0
            candidates = _safe_sites(row, helix,
                                     int(math.floor(start)),
                                     int(math.ceil(end)))
            if count < 0:
                candidates = [index for index in candidates if index not in
                              deletion_protected.get(helix, set())]
            selected = _partition_sites(
                abs(count), start, end, candidates, domain_size, loads,
                insertion=count > 0,
                maximum_per_domain=maximum_per_domain)
            if count > 0:
                for index in selected:
                    row["loop"][index] += 1
                inserted.extend(selected)
            elif count < 0:
                for index in selected:
                    row["skip"][index] = -1
                deleted.extend(selected)
            vertex_reports.append({
                "vertex": vertex, "signed_indel": count,
                "window_start": start, "window_end": end,
                "sites": selected})
        if len(inserted) - len(deleted) != residual:
            raise RuntimeError("Frame indel relocation changed helix length.")
        record.update({"insertions": sorted(inserted),
                       "deletions": sorted(deleted),
                       "target_bases": int(record["nominal_bases"]) +
                           residual,
                       "frame_vertices": vertex_reports,
                       "domain_indel_quotas": [
                           loads.get(index, 0) for index in range(
                               int(math.ceil(record["nominal_bases"] /
                                             float(domain_size))))],
                       "maximum_indel_in_one_domain":
                           max(loads.values() or [0])})
        reports.append(record)
    curvature.update({"version": 3, "mode": "frame-corner-localized-indels",
                      "rings": reports, "domain_size_bp": domain_size,
                      "maximum_indel_per_domain_allowed":
                          maximum_per_domain,
                      "maximum_insertion_per_domain": max(
                          [record["maximum_indel_in_one_domain"]
                           for record in reports if record["insertions"]] or
                          [0]),
                      "maximum_deletion_per_domain": max(
                          [record["maximum_indel_in_one_domain"]
                           for record in reports if record["deletions"]] or
                          [0])})
    metadata = dict(design.get("curved_metadata") or {})
    geometry = _frame_geometry(
        _decode_geometry(metadata), plan, plan["rings"])
    if geometry:
        metadata["geometry_data"] = encode_geometry_payload(geometry)
    metadata.update({"format": "cadnano-frame-project-v1",
                     "shape": "frame", "frame_plan": dict(
                         (key, value) for key, value in plan.items()
                         if key != "rings"),
                     "curvature_encoding": "frame-corner-localized-indels",
                     "autocs_topology_fingerprint": topology_before,
                     "autocs_topology_frozen": True})
    design["curvature_indels"] = curvature
    design["curved_metadata"] = metadata
    # Curved Design's stable AutoCS topology remains frozen.  Only move the
    # already-planned indels inside their own vertex windows so every
    # physical adjacent-helix pair follows its theoretical floor/ceiling as
    # closely as the legal coordinates permit.
    pair_summary = optimize_frame_pair_curvature(design)
    curvature = design["curvature_indels"]
    curvature["frame_pair_curvature_summary"] = pair_summary
    topology_after = topology_fingerprint(design)
    if topology_after != topology_before:
        raise RuntimeError(
            "Frame postprocessing attempted to alter frozen AutoCS topology.")
    with open(json_path, "w", encoding="utf-8", newline="\n") as output:
        json.dump(design, output, separators=(",", ":"))
    return curvature, metadata


def create_frame_project(spec, progress=None, cancelled=None):
    """Create a frame while preserving Curved Design's AutoCS topology."""
    plan = plan_frame(spec)
    if not plan["feasible"]:
        raise ValueError("Frame geometry is not feasible: %s." %
                         "; ".join(plan["failure_reasons"]))
    frame_name = safe_name(spec.get("name") or spec.get(
        "frame_shape", "frame"))
    equivalent_outer_diameter = float(
        plan["curved_parent_outer_diameter_nm"])
    curved_spec = {
        "shape": "cylinder", "lattice": plan["lattice"],
        "height_nm": float(spec["cross_section_height_nm"]),
        "maximum_diameter_nm": equivalent_outer_diameter,
        "minimum_diameter_nm": equivalent_outer_diameter,
        "layers": int(spec.get("layers", 1)),
        "scaffold_crossover_density_mode": spec.get(
            "scaffold_crossover_density_mode", "periodic"),
        "scaffold_crossover_density_multiple": int(spec.get(
            "scaffold_crossover_density_multiple", 1)),
        "name": frame_name, "project_root": spec["project_root"],
        "output_name_override": frame_name + "_frame_" + plan["lattice"]}
    result = create_curved_project(
        curved_spec, progress=progress, cancelled=cancelled)
    curvature, metadata = _localize_corner_indels(result["json_path"], plan)
    curvature, metadata, twist_audit = \
        _apply_straight_common_mode_remove_twist(
            result["json_path"], plan)
    staple_audit, metadata = _reoptimize_frame_staple_nicks(
        result["json_path"], plan["lattice"])
    _rewrite_frame_artifacts(result, plan)
    result["metadata"] = metadata
    result["frame_plan"] = plan
    result["frame_straight_common_mode_remove_twist"] = twist_audit
    result["frame_staple_nick_optimization"] = staple_audit
    result["indel_summary"] = {
        "domain_size_bp": plan["domain_size_bp"],
        "maximum_insertion_per_domain":
            curvature["maximum_insertion_per_domain"],
        "maximum_deletion_per_domain":
            curvature["maximum_deletion_per_domain"],
        "maximum_indel_per_domain_allowed": int(
            plan["maximum_indel_per_domain_allowed"])}
    # Update the saved design-settings record after the circular routing
    # scaffold has been converted into a frame.
    for path in result.get("input_paths", ()):
        if path.endswith("_design_settings.json") and os.path.isfile(path):
            with open(path, "w", encoding="utf-8", newline="\n") as output:
                json.dump(metadata, output, indent=2, sort_keys=True)
                output.write("\n")
    return result
