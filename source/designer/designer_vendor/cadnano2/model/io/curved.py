"""Parameterised DNAxiS curved-origami project integration."""

import json
import math
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import queue
import threading
import zlib
from collections import defaultdict

from .twistbend import (TwistBendError,
                        equal_partition_indel_sites)
from .indelanalysis import optimize_curved_pair_curvature


CURVED_METADATA_VERSION = 1
RING_SPACING_NM = 2.8
DNA_HELIX_RADIUS_NM = 1.0
BP_RISE_NM = 0.332
MIN_RING_BP = 72
CURVED_SCAFFOLD_MAX_BASES = 25000


PRESETS = (
    ("bowl", "Bowl", "Open concave bowl"),
    ("sphere", "Sphere / Ellipsoid", "Closed rounded shell"),
    ("cone", "Cone / Frustum", "Tapered rotational shell"),
    ("mushroom", "Mushroom", "Stem with a broad rounded cap"),
    ("gourd", "Gourd", "Two joined rounded lobes"),
    ("vase", "Vase", "Bulb, narrow neck and flared rim"),
    ("cylinder", "Cylinder", "Constant-diameter rotational shell"))


def curved_root():
    return os.path.abspath(os.path.join(
        os.path.dirname(__file__), os.pardir, os.pardir,
        "third_party", "dnaxis"))


def safe_name(value):
    cleaned = "".join(char if char.isalnum() or char in "-_" else "_"
                      for char in str(value).strip())
    return cleaned.strip("-_") or "curved-design"


def _dimension_name_token(value):
    """Return a compact filesystem-safe decimal such as 12p2."""
    text = ("%.6f" % float(value)).rstrip("0").rstrip(".")
    return text.replace("-", "m").replace(".", "p")


def curved_output_name(name, shape, lattice=None, layers=None,
                       height_nm=None, maximum_diameter_nm=None,
                       minimum_diameter_nm=None):
    """Build a descriptive Curved project/file stem.

    New projects include shape, lattice, layer count and requested outer
    dimensions.  Optional arguments preserve compatibility when opening an
    older metadata record that contains only a name and shape.
    """
    base = safe_name(name)
    shape_suffix = "_" + safe_name(shape)
    if not base.lower().endswith(shape_suffix.lower()):
        base += shape_suffix
    details = []
    if lattice is not None:
        details.append(safe_name(str(lattice).lower()))
    if layers is not None:
        details.append("%dL" % int(layers))
    if height_nm is not None:
        details.append("H%s" % _dimension_name_token(height_nm))
    if maximum_diameter_nm is not None and \
            minimum_diameter_nm is not None:
        details.append("D%s-%snm" % (
            _dimension_name_token(maximum_diameter_nm),
            _dimension_name_token(minimum_diameter_nm)))
    return base + (("_" + "_".join(details)) if details else "")


def unique_curved_project_target(project_root, output_name):
    """Return a collision-free project root and one shared numeric suffix.

    The root folder, main JSON, and every automatically generated input file
    use the same ``_N`` suffix.  Existing projects are never reused or
    overwritten, including when only a same-named loose file exists beside
    the requested project folder.
    """
    requested_root = os.path.abspath(project_root)
    parent = os.path.dirname(requested_root)
    base_root = os.path.basename(requested_root.rstrip(os.sep))
    base_output = safe_name(output_name)
    extensions = (".json", ".stl", ".csv", ".png", ".pdb", ".cif",
                  ".top", ".dat")

    def collides(root, name):
        if os.path.lexists(root):
            return True
        return any(os.path.lexists(os.path.join(parent, name + extension))
                   for extension in extensions)

    if not collides(requested_root, base_output):
        return requested_root, base_output, "", 0
    version = 1
    while True:
        suffix = "_%d" % version
        candidate_root = os.path.join(parent, base_root + suffix)
        candidate_output = base_output + suffix
        if not collides(candidate_root, candidate_output):
            return candidate_root, candidate_output, suffix, version
        version += 1


def _lerp(points, t):
    if t <= points[0][0]:
        return points[0][1]
    for first, second in zip(points, points[1:]):
        if t <= second[0]:
            fraction = (t - first[0]) / (second[0] - first[0])
            smooth = fraction * fraction * (3.0 - 2.0 * fraction)
            return first[1] + smooth * (second[1] - first[1])
    return points[-1][1]


def _profile_radius(shape, t, minimum, maximum):
    span = max(0.0, maximum - minimum)
    if shape == "cylinder":
        return maximum
    if shape == "cone":
        return maximum - span * t
    if shape == "bowl":
        return minimum + span * (t ** 1.65)
    if shape == "sphere":
        x = 2.0 * t - 1.0
        return minimum + span * math.sqrt(max(0.0, 1.0 - x * x))
    if shape == "mushroom":
        return _lerp(((0.0, minimum), (0.42, minimum),
                      (0.58, minimum * 1.15), (0.78, maximum),
                      (1.0, max(minimum, maximum * 0.62))), t)
    if shape == "gourd":
        return _lerp(((0.0, max(minimum, maximum * 0.48)),
                      (0.23, maximum * 0.86), (0.46, minimum),
                      (0.72, maximum),
                      (1.0, max(minimum, maximum * 0.58))), t)
    if shape == "vase":
        return _lerp(((0.0, maximum * 0.62), (0.18, maximum),
                      (0.44, maximum * 0.76), (0.60, minimum),
                      (0.82, minimum * 1.05), (1.0, maximum * 0.72)), t)
    raise ValueError("Unknown curved preset: %s" % shape)


def build_rings(shape, height_nm, maximum_diameter_nm,
                minimum_diameter_nm, layers=1, lattice="square"):
    shape = str(shape).lower()
    outer_height = float(height_nm)
    outer_maximum = float(maximum_diameter_nm) / 2.0
    outer_minimum = float(minimum_diameter_nm) / 2.0
    height = outer_height - 2.0 * DNA_HELIX_RADIUS_NM
    maximum = outer_maximum - DNA_HELIX_RADIUS_NM
    minimum = outer_minimum - DNA_HELIX_RADIUS_NM
    layers = int(layers)
    lattice = str(lattice).lower()
    if lattice not in ("square", "honeycomb"):
        raise ValueError("Curved lattice must be square or honeycomb.")
    if height < RING_SPACING_NM:
        raise ValueError(
            "Outer height must be at least %.1f nm." %
            (RING_SPACING_NM + 2.0 * DNA_HELIX_RADIUS_NM))
    if maximum <= 0 or minimum <= 0 or minimum > maximum:
        raise ValueError(
            "Outer diameters must exceed the %.1f nm DNA diameter, and the "
            "minimum cannot exceed the maximum." %
            (2.0 * DNA_HELIX_RADIUS_NM))
    if layers < 1 or layers > 3:
        raise ValueError("Reinforcement layers must be between 1 and 3.")

    # Sample the base-layer centreline by meridional arc length.  Uniform-z
    # sampling is correct only for a cylinder; on a sloped profile it makes
    # neighbouring ring helices farther apart than the requested 2.8 nm.
    dense_count = max(1001, int(math.ceil(height * 100.0)) + 1)
    dense = []
    cumulative = [0.0]
    for index in range(dense_count):
        t = index / float(dense_count - 1)
        point = (
            DNA_HELIX_RADIUS_NM + height * t,
            _profile_radius(shape, t, minimum, maximum))
        if dense:
            previous = dense[-1]
            cumulative.append(cumulative[-1] + math.hypot(
                point[0] - previous[0], point[1] - previous[1]))
        dense.append(point)
    meridian_length = cumulative[-1]
    meridian_pitch = (RING_SPACING_NM if lattice == "square" else
                      RING_SPACING_NM * math.sqrt(3.0) / 2.0)
    interval_count = max(1, int(round(meridian_length / meridian_pitch)))
    if lattice == "honeycomb":
        # An even number of rings gives each serpentine layer transition a
        # native brick-wall rung: layer 0->1 at the high end, then 1->2 at
        # the low end.  Pick the closest odd interval count.
        odd_candidates = [value for value in (
            interval_count - 1, interval_count, interval_count + 1)
                          if value >= 1 and value % 2 == 1]
        interval_count = min(
            odd_candidates,
            key=lambda value: abs(
                meridian_length / float(value) - meridian_pitch))
        target_meridian = interval_count * meridian_pitch

        def scaled_meridian(scale):
            return sum(math.hypot(
                (right[0] - left[0]) * scale,
                right[1] - left[1])
                for left, right in zip(dense, dense[1:]))

        radial_only = scaled_meridian(0.0)
        while target_meridian + 1e-9 < radial_only:
            interval_count += 2
            target_meridian = interval_count * meridian_pitch
        low_scale, high_scale = 0.0, 1.0
        while scaled_meridian(high_scale) < target_meridian:
            high_scale *= 2.0
        for unused in range(60):
            middle = (low_scale + high_scale) / 2.0
            if scaled_meridian(middle) < target_meridian:
                low_scale = middle
            else:
                high_scale = middle
        height_scale = (low_scale + high_scale) / 2.0
        dense = [(DNA_HELIX_RADIUS_NM +
                  (point[0] - DNA_HELIX_RADIUS_NM) * height_scale,
                  point[1]) for point in dense]
        cumulative = [0.0]
        for left, right in zip(dense, dense[1:]):
            cumulative.append(cumulative[-1] + math.hypot(
                right[0] - left[0], right[1] - left[1]))
        meridian_length = cumulative[-1]
    else:
        height_scale = 1.0
    slice_count = interval_count + 1
    target_step = meridian_length / float(interval_count)
    profile = []
    cursor = 0
    for slice_index in range(slice_count):
        target_length = min(
            meridian_length, slice_index * target_step)
        while (cursor + 1 < len(cumulative) and
               cumulative[cursor + 1] < target_length):
            cursor += 1
        if cursor + 1 >= len(cumulative):
            z_value, radius = dense[-1]
        else:
            segment = cumulative[cursor + 1] - cumulative[cursor]
            fraction = (0.0 if segment <= 1e-12 else
                        (target_length - cumulative[cursor]) / segment)
            z_value = (dense[cursor][0] + fraction *
                       (dense[cursor + 1][0] - dense[cursor][0]))
            radius = (dense[cursor][1] + fraction *
                      (dense[cursor + 1][1] - dense[cursor][1]))
        profile.append((z_value, radius))

    radius_floor = MIN_RING_BP * BP_RISE_NM / (2.0 * math.pi)
    grid = []
    tangents = []
    for index in range(len(profile)):
        before = profile[max(0, index - 1)]
        after = profile[min(len(profile) - 1, index + 1)]
        dz = after[0] - before[0]
        dr = after[1] - before[1]
        magnitude = math.hypot(dz, dr) or 1.0
        tangents.append((dz / magnitude, dr / magnitude))

    for layer in range(layers):
        layer_rows = []
        for slice_index, (z_value, base_radius) in enumerate(profile):
            tangent_z, tangent_r = tangents[slice_index]
            normal_z, normal_r = -tangent_r, tangent_z
            if lattice == "square":
                normal_offset = layer * RING_SPACING_NM
                actual_z = z_value
                requested_radius = base_radius + normal_offset
                indel_radius = requested_radius
            else:
                # Honeycomb is a brick-wall strip.  Successive rings are
                # separated by sqrt(3)/2*d tangentially and alternate by d/2
                # normally.  Mirrored layers have a 3d/2 centre-line pitch.
                parity_offset = (RING_SPACING_NM / 2.0
                                 if (layer + slice_index) % 2 else 0.0)
                normal_offset = (
                    layer * 1.5 * RING_SPACING_NM + parity_offset -
                    RING_SPACING_NM / 2.0)
                actual_z = z_value + normal_offset * normal_z
                requested_radius = base_radius + normal_offset * normal_r
                # Honeycomb crossover phases define which neighbouring
                # helices can connect, but they do not compensate for the
                # different circumferences of the alternating inner/outer
                # centre-lines.  Use every helix's actual centre-line radius
                # as its indel target.  This includes both the brick-wall
                # parity offset in a single layer and the additional radial
                # offset in reinforced layers.
                indel_radius = requested_radius
            radius = max(radius_floor, requested_radius)
            indel_radius = max(radius_floor, indel_radius)
            bp = max(MIN_RING_BP, int(round(
                2.0 * math.pi * radius / BP_RISE_NM / 2.0)) * 2)
            indel_bp = max(MIN_RING_BP, int(round(
                2.0 * math.pi * indel_radius / BP_RISE_NM / 2.0)) * 2)
            actual_radius = bp * BP_RISE_NM / (2.0 * math.pi)
            layer_rows.append({
                "bp": bp, "height_nm": actual_z,
                "radius_nm": actual_radius,
                "geometry_radius_nm": radius,
                "indel_radius_nm": indel_radius,
                "indel_bp": indel_bp,
                "requested_radius_nm": radius,
                "layer": layer, "slice": slice_index,
                "meridian_spacing_nm": target_step,
                "meridian_pitch_nm": meridian_pitch,
                "normal_offset_nm": normal_offset,
                "lattice_center_spacing_nm": RING_SPACING_NM,
                "lattice": lattice,
                "requested_outer_height_nm": outer_height,
                "actual_outer_height_nm":
                    2.0 * DNA_HELIX_RADIUS_NM + height * height_scale,
                "outer_dimensions": True,
                "profile_shape": shape})
        grid.append(layer_rows)

    # A serpentine order gives DNAxiS one continuous scaffold pathway while
    # keeping every consecutive pair physically adjacent.
    rings = []
    for layer, rows in enumerate(grid):
        ordered = rows if layer % 2 == 0 else list(reversed(rows))
        rings.extend(ordered)
    for index, ring in enumerate(rings):
        ring["direction"] = bool(index % 2 == 0)
        ring["index"] = index
    return rings


def estimated_scaffold_bases(rings):
    return sum(int(ring["bp"]) for ring in rings)


def curved_indel_plan(rings, lattice):
    """Build one equal-length lattice-periodic parent, then encode curvature.

    DNAxiS-style curved shells do not start from independently rounded ring
    lengths.  Every helix first receives the same complete-period nominal
    length.  Insertions/deletions then encode only the *relative* circumference
    change between rings.  The small absolute rounding error is applied as a
    common radius offset.  A Square one-layer cylinder therefore needs no
    indels, while a Honeycomb cylinder can still need them because its native
    brick-wall channels occupy alternating centre-line radii.
    """
    lattice = str(lattice).lower()
    if lattice == "honeycomb":
        period, domain_size = 21, 7
        minimum_nominal = 84
        parent_period = 42
    elif lattice == "square":
        period, domain_size = 32, 8
        minimum_nominal = 96
        parent_period = 32
    else:
        raise ValueError("Curved lattice must be square or honeycomb.")

    if not rings:
        raise ValueError("Curved design requires at least one ring helix.")

    requested = [int(ring.get("indel_bp", ring["bp"])) for ring in rings]
    def snap_parent(value):
        lower = max(parent_period,
                    (value // parent_period) * parent_period)
        upper = lower if lower == value else lower + parent_period
        return max(minimum_nominal, min(
            (lower, upper),
            key=lambda candidate: (abs(candidate - value), -candidate)))

    # Choose the common reference (the target ring that maps to the unmodified
    # parent) by a minimax search.  This is deliberately not the arithmetic
    # mean: an asymmetric cross-section can make the mean leave one edge with
    # a much larger indel load.  The score first minimizes the largest
    # *absolute* helix load per 7/8-bp domain, before considering integer
    # domain peaks or total load.  Thus -4/0/+4 is preferred over 0/+4/+8
    # even when both plans happen to round to the same <=3 indels/domain
    # limit.  This keeps the physical domain-length change as small and as
    # symmetric as possible instead of merely satisfying the hard bound.
    # The final target circumferences retain every pairwise difference; only a
    # common whole-design circumference shift is introduced by period snap.
    mean_requested = sum(requested) / float(len(requested))
    parent_options = []
    for candidate_reference in range(min(requested), max(requested) + 1):
        candidate_nominal = snap_parent(candidate_reference)
        candidate_domains = max(1, candidate_nominal // domain_size)
        residuals = [value - candidate_reference for value in requested]
        insertion_total = sum(max(0, value) for value in residuals)
        deletion_total = sum(max(0, -value) for value in residuals)
        insertion_peak = int(math.ceil(
            max([0] + [max(0, value) for value in residuals]) /
            float(candidate_domains)))
        deletion_peak = int(math.ceil(
            max([0] + [max(0, -value) for value in residuals]) /
            float(candidate_domains)))
        peak = max(insertion_peak, deletion_peak)
        maximum_absolute_load = max([0] + [abs(value)
                                            for value in residuals])
        parent_options.append((
            maximum_absolute_load / float(candidate_domains),
            maximum_absolute_load,
            peak,
            sum(abs(value) for value in residuals),
            abs(insertion_total - deletion_total),
            insertion_peak + deletion_peak,
            abs(candidate_reference - mean_requested),
            -candidate_nominal,
            candidate_reference,
            candidate_nominal))
    unused_normalized_peak, unused_absolute_peak, unused_peak, unused_total, \
        unused_imbalance, unused_peak_sum, unused_mean_distance, \
        unused_negative_nominal, reference, nominal = min(parent_options)
    common_shift = nominal - reference
    crossover_periods = nominal // period
    domain_count = nominal // domain_size

    planned = [dict(ring) for ring in rings]
    groups = {}
    for ring, requested_bp in zip(planned, requested):
        adjusted_bp = requested_bp + common_shift
        if adjusted_bp < period:
            raise ValueError(
                "The common-period radius adjustment makes a ring too short.")
        ring["requested_bp"] = requested_bp
        ring["bp"] = adjusted_bp
        ring["nominal_bp"] = nominal
        # Keep the native lattice geometry while applying the same small
        # whole-design rounding shift to every physical ring.  Honeycomb's
        # alternating radial channels have already contributed their actual
        # centre-line circumference to requested_bp above.
        geometry_radius = float(ring.get(
            "geometry_radius_nm", ring["radius_nm"]))
        geometry_radius += common_shift * BP_RISE_NM / (2.0 * math.pi)
        ring["geometry_radius_nm"] = geometry_radius
        ring["radius_nm"] = geometry_radius
        groups.setdefault(adjusted_bp, []).append(int(ring["index"]))

    maximum_insertion = 0
    maximum_deletion = 0
    maximum_insertion_per_domain = 0
    maximum_deletion_per_domain = 0
    ring_plans = []
    forced = []
    for adjusted_bp, indices in sorted(groups.items()):
        residual = adjusted_bp - nominal
        average_spacing = adjusted_bp / float(crossover_periods)
        peak = int(math.ceil(
            abs(residual) / float(crossover_periods)))
        if residual > 0:
            maximum_insertion = max(maximum_insertion, peak)
            maximum_insertion_per_domain = max(
                maximum_insertion_per_domain,
                int(math.ceil(residual / float(domain_count))))
        elif residual < 0:
            maximum_deletion = max(maximum_deletion, peak)
            maximum_deletion_per_domain = max(
                maximum_deletion_per_domain,
                int(math.ceil(-residual / float(domain_count))))
        domain_peak = int(math.ceil(
            abs(residual) / float(domain_count)))
        option = {
            "target": adjusted_bp,
            "requested_target": adjusted_bp - common_shift,
            "indices": indices,
            "periods": crossover_periods,
            "nominal": nominal,
            "residual": residual,
            "average_spacing": average_spacing,
            "domain_size": domain_size,
            "domain_count": domain_count,
            "maximum_per_domain": domain_peak,
            "domain_limit_ok": domain_peak <= 3,
            "maximum_per_crossover": peak}
        ring_plans.append(option)
        if not option["domain_limit_ok"]:
            forced.append(option)

    actual_minimum = min(
        int(math.floor(option["average_spacing"]))
        for option in ring_plans)
    actual_maximum = max(
        int(math.ceil(option["average_spacing"]))
        for option in ring_plans)
    is_cylinder = all(str(ring.get("profile_shape", "")).lower() ==
                      "cylinder" for ring in rings)
    diameter_adjustment = 2.0 * common_shift * BP_RISE_NM / (2.0 * math.pi)

    summary = {
        "lattice": lattice,
        "period": period,
        "domain_size_bp": domain_size,
        "domain_count": domain_count,
        "maximum_indel_per_domain_allowed": 3,
        "maximum_insertion_per_domain": maximum_insertion_per_domain,
        "maximum_deletion_per_domain": maximum_deletion_per_domain,
        "domain_limit_feasible": not forced,
        "maximum_insertion_per_crossover": maximum_insertion,
        "maximum_deletion_per_crossover": maximum_deletion,
        "effective_crossover_spacing_minimum": actual_minimum,
        "effective_crossover_spacing_maximum": actual_maximum,
        "requires_confirmation": False,
        "infeasible_ring_groups": [dict(option) for option in forced],
        "ring_groups": ring_plans,
        "common_nominal_bp": nominal,
        "optimized_parent_reference_bp": reference,
        "minimum_nominal_bp": minimum_nominal,
        "equal_length_parent_period_bp": parent_period,
        "baseline_requested_bp": reference,
        "circumference_adjustment_bp": common_shift,
        "diameter_adjustment_nm": diameter_adjustment,
        "equal_length_parent": True,
        "parent_selection": "minimax-domain-load",
        "cylinder_period_snap": is_cylinder and len(groups) == 1}
    return planned, summary


def _stl_text(rings, segments=48):
    profile = sorted(
        (ring for ring in rings if int(ring["layer"]) == 0),
        key=lambda ring: int(ring["slice"]))
    lines = ["solid cadnano_curved_design"]
    for lower, upper in zip(profile, profile[1:]):
        for segment in range(segments):
            angles = (2.0 * math.pi * segment / segments,
                      2.0 * math.pi * (segment + 1) / segments)
            points = []
            for ring, angle in ((lower, angles[0]),
                                (lower, angles[1]),
                                (upper, angles[0]),
                                (upper, angles[1])):
                points.append((ring["radius_nm"] * math.cos(angle),
                               ring["radius_nm"] * math.sin(angle),
                               ring["height_nm"]))
            for triangle in ((points[0], points[1], points[2]),
                             (points[1], points[3], points[2])):
                lines.append("  facet normal 0 0 0")
                lines.append("    outer loop")
                for point in triangle:
                    lines.append("      vertex %.8f %.8f %.8f" % point)
                lines.append("    endloop")
                lines.append("  endfacet")
    lines.extend(("endsolid cadnano_curved_design", ""))
    return "\n".join(lines)


def _png_chunk(kind, data):
    return (struct.pack(">I", len(data)) + kind + data +
            struct.pack(">I", zlib.crc32(kind + data) & 0xffffffff))


def _preview_png(rings, width=900, height=620):
    pixels = [bytearray([248, 249, 251] * width) for unused in range(height)]
    max_radius = max(float(ring["radius_nm"]) for ring in rings)
    max_height = max(float(ring["height_nm"]) for ring in rings) or 1.0
    scale = min((width - 90.0) / (2.0 * max_radius),
                (height - 80.0) / max_height)

    def set_pixel(x_value, y_value, color):
        x_value, y_value = int(x_value), int(y_value)
        if 0 <= x_value < width and 0 <= y_value < height:
            offset = x_value * 3
            pixels[y_value][offset:offset + 3] = bytes(color)

    for ring in rings:
        cx = width / 2.0
        cy = height - 40.0 - float(ring["height_nm"]) * scale
        rx = float(ring["radius_nm"]) * scale
        ry = max(1.5, rx * 0.18)
        color = ((54, 102, 166) if int(ring["layer"]) == 0 else
                 (112, 145, 187))
        steps = max(80, int(rx * 5))
        for step in range(steps):
            angle = 2.0 * math.pi * step / steps
            set_pixel(cx + rx * math.cos(angle),
                      cy + ry * math.sin(angle), color)
    raw = b"".join(b"\x00" + bytes(row) for row in pixels)
    return (b"\x89PNG\r\n\x1a\n" +
            _png_chunk(b"IHDR", struct.pack(
                ">IIBBBBB", width, height, 8, 2, 0, 0, 0)) +
            _png_chunk(b"IDAT", zlib.compress(raw, 9)) +
            _png_chunk(b"IEND", b""))


def _dnaxis_python():
    candidates = [
        os.path.expanduser("~/virtualenvs/dnaxis-cadnano/bin/python"),
        sys.executable, "/opt/homebrew/bin/python3",
        shutil.which("python3")]
    checked = set()
    for candidate in candidates:
        if not candidate or candidate in checked or not os.path.isfile(candidate):
            continue
        checked.add(candidate)
        probe = subprocess.run(
            [candidate, "-c", "import numpy, scipy"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if probe.returncode == 0:
            return candidate
    raise RuntimeError(
        "The DNAxiS scientific runtime is not installed. Expected "
        "~/virtualenvs/dnaxis-cadnano with NumPy and SciPy.")


def _run(command, progress=None, cancelled=None):
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, bufsize=1)
    lines = []
    output_queue = queue.Queue()

    def read_output():
        for line in process.stdout:
            output_queue.put(line)
        output_queue.put(None)

    reader = threading.Thread(target=read_output, daemon=True)
    reader.start()
    stream_finished = False
    while not stream_finished or process.poll() is None:
        try:
            line = output_queue.get(timeout=0.10)
        except queue.Empty:
            line = ""
        if line is None:
            stream_finished = True
        elif line:
            lines.append(line)
            if progress is not None:
                progress(line.rstrip())
        elif progress is not None:
            # Let Qt process the Cancel click even while DNAxiS is silent.
            progress("")
        if cancelled is not None and cancelled():
            process.terminate()
            try:
                process.wait(timeout=3)
            except subprocess.TimeoutExpired:
                process.kill()
            raise RuntimeError("Curved Design was cancelled.")
    code = process.wait()
    reader.join(timeout=1)
    if code:
        full_output = "".join(lines)
        tail = "\n".join(full_output.splitlines()[-80:])
        traceback_at = full_output.rfind("Traceback (most recent call last):")
        if traceback_at >= 0:
            traceback_lines = full_output[traceback_at:].splitlines()[:80]
            detail = "\n".join(traceback_lines) + "\n\nLast output:\n" + tail
        else:
            detail = tail
        raise RuntimeError(
            "DNAxiS failed with exit code %d.\n\n%s" %
            (code, detail))


def _final_crossover_indel_statistics(encoded, lattice, curvature_indels):
    """Measure actual circular spacing between final native crossovers."""
    period = 21 if lattice == "honeycomb" else 32
    rows = {int(row["num"]): row for row in encoded.get("vstrands", [])}
    ring_records = {
        int(record["helix"]): record
        for record in curvature_indels.get("rings", [])}
    intervals = []
    unresolved_pairs = []
    for helix_id, row in rows.items():
        record = ring_records.get(helix_id)
        if record is None:
            continue
        nominal_size = int(record["nominal_bases"])
        loops = row.get("loop", [])
        skips = row.get("skip", [])
        for strand_type in ("scaf", "stap"):
            by_neighbour = {}
            for index, connection in enumerate(row.get(strand_type, [])):
                if index >= nominal_size:
                    break
                # Measure one directed phase family at a time.  Combining
                # the incoming and outgoing halves of a reciprocal block
                # creates false 1/31-bp intervals; outgoing sites repeat at
                # the real lattice period, including across max -> 0.
                neighbour = int(connection[2])
                if neighbour >= 0 and neighbour != helix_id:
                    by_neighbour.setdefault(
                        ("outgoing", neighbour), set()).add(index)
                neighbour = int(connection[0])
                if neighbour >= 0 and neighbour != helix_id:
                    by_neighbour.setdefault(
                        ("incoming", neighbour), set()).add(index)
            for (direction, neighbour), position_set in sorted(
                    by_neighbour.items()):
                positions = sorted(position_set)
                if not positions:
                    unresolved_pairs.append(
                        [strand_type, direction, helix_id, neighbour])
                    continue
                # One directed crossover on a circular helix represents one
                # full-circumference interval, not an unmeasurable pair.
                if len(positions) == 1:
                    position_pairs = [(positions[0],
                                       positions[0] + nominal_size)]
                else:
                    position_pairs = zip(
                        positions,
                        positions[1:] + [positions[0] + nominal_size])
                for first, second in position_pairs:
                    nominal_distance = second - first
                    insertions = 0
                    deletions = 0
                    for offset in range(1, nominal_distance + 1):
                        index = (first + offset) % nominal_size
                        loop = (int(loops[index])
                                if index < len(loops) else 0)
                        skip = (int(skips[index])
                                if index < len(skips) else 0)
                        insertions += max(0, loop)
                        deletions += max(0, -skip)
                    crossover_periods = max(
                        1, int(round(
                            nominal_distance / float(period))))
                    insertion_peak = int(math.ceil(
                        insertions / float(crossover_periods)))
                    deletion_peak = int(math.ceil(
                        deletions / float(crossover_periods)))
                    actual_span = (
                        nominal_distance + insertions - deletions)
                    actual_average = (
                        actual_span / float(crossover_periods))
                    intervals.append({
                        "strand_type": strand_type,
                        "direction": direction,
                        "helix": helix_id, "neighbour": neighbour,
                        "start": first, "end": second % nominal_size,
                        "wraps_boundary": second >= nominal_size,
                        "nominal_distance": nominal_distance,
                        "insertions": insertions,
                        "deletions": deletions,
                        "crossover_periods": crossover_periods,
                        "maximum_insertion_per_crossover": insertion_peak,
                        "maximum_deletion_per_crossover": deletion_peak,
                        "actual_span": actual_span,
                        "actual_spacing_average": actual_average,
                        "actual_spacing_minimum": int(math.floor(
                            actual_average)),
                        "actual_spacing_maximum": int(math.ceil(
                            actual_average))})

    if intervals:
        maximum_insertion = max(
            item["maximum_insertion_per_crossover"] for item in intervals)
        maximum_deletion = max(
            item["maximum_deletion_per_crossover"] for item in intervals)
        spacing_minimum = min(
            item["actual_spacing_minimum"] for item in intervals)
        spacing_maximum = max(
            item["actual_spacing_maximum"] for item in intervals)
    else:
        maximum_insertion = int(curvature_indels.get(
            "maximum_insertion_per_crossover", 0))
        maximum_deletion = int(curvature_indels.get(
            "maximum_deletion_per_crossover", 0))
        spacing_minimum = period - maximum_deletion
        spacing_maximum = period + maximum_insertion
    staple_intervals = [item for item in intervals
                        if item["strand_type"] == "stap"]
    staple_interval_minimum = (min(
        item["actual_span"] for item in staple_intervals)
        if staple_intervals else 0)
    staple_interval_maximum = (max(
        item["actual_span"] for item in staple_intervals)
        if staple_intervals else 0)
    domain_size = 7 if lattice == "honeycomb" else 8
    maximum_insertion_per_domain = 0
    maximum_deletion_per_domain = 0
    domain_violations = []
    for helix, record in sorted(ring_records.items()):
        row = rows.get(helix)
        if row is None:
            continue
        nominal_size = int(record["nominal_bases"])
        domain_count = max(1, int(math.ceil(
            nominal_size / float(domain_size))))
        for domain in range(domain_count):
            first = domain * domain_size
            last = min(nominal_size, first + domain_size)
            insertions = sum(max(0, int(row.get("loop", [])[index]))
                             for index in range(first, last))
            deletions = sum(max(0, -int(row.get("skip", [])[index]))
                            for index in range(first, last))
            maximum_insertion_per_domain = max(
                maximum_insertion_per_domain, insertions)
            maximum_deletion_per_domain = max(
                maximum_deletion_per_domain, deletions)
            if max(insertions, deletions) > 3:
                domain_violations.append({
                    "helix": helix, "domain": domain,
                    "start": first, "end": last - 1,
                    "insertions": insertions,
                    "deletions": deletions})
    return {
        "source": "final native crossovers and domain-distributed indels",
        "period": period,
        "domain_size_bp": domain_size,
        "maximum_indel_per_domain_allowed": 3,
        "maximum_insertion_per_domain": maximum_insertion_per_domain,
        "maximum_deletion_per_domain": maximum_deletion_per_domain,
        "domain_limit_feasible": not domain_violations,
        "domain_limit_violations": domain_violations,
        "maximum_insertion_per_crossover": maximum_insertion,
        "maximum_deletion_per_crossover": maximum_deletion,
        "effective_crossover_spacing_minimum": spacing_minimum,
        "effective_crossover_spacing_maximum": spacing_maximum,
        "staple_crossover_interval_minimum": staple_interval_minimum,
        "staple_crossover_interval_maximum": staple_interval_maximum,
        "staple_crossover_density_maximum_denominator":
            staple_interval_minimum,
        "staple_crossover_density_minimum_denominator":
            staple_interval_maximum,
        "measured_interval_count": len(intervals),
        "unresolved_single_crossover_pairs": unresolved_pairs,
        "intervals": intervals}


def _balance_all_staple_crossover_density_unused(encoded, minimum=21):
    """Redistribute indels so every adjacent staple-xover arc is >= minimum.

    Insertion/deletion totals on each helix are conserved.  Indels are moved
    only between crossover-delimited arcs, and every target remains clear of
    scaffold/staple crossovers and nicks.
    """
    blank = [-1, -1, -1, -1]
    moved_insertions = 0
    moved_deletions = 0
    unresolved = []
    final_spacings = []

    for row in encoded.get("vstrands", []):
        helix = int(row["num"])
        size = len(row.get("stap", []))
        if size < 2:
            continue
        positions = sorted(set(
            index for index, entry in enumerate(row["stap"])
            if ((int(entry[0]) >= 0 and int(entry[0]) != helix) or
                (int(entry[2]) >= 0 and int(entry[2]) != helix))))
        if len(positions) < 2:
            continue

        def arc_indices(start, end):
            result = []
            index = (start + 1) % size
            while True:
                result.append(index)
                if index == end:
                    return result
                index = (index + 1) % size

        arcs = [arc_indices(start, end) for start, end in zip(
            positions, positions[1:] + positions[:1])]

        def weight(index):
            return max(0, 1 + int(row["loop"][index]) +
                       int(row["skip"][index]))

        def spans():
            return [sum(weight(index) for index in arc) for arc in arcs]

        safe = set()
        for index in range(size):
            valid = True
            for strand_type in ("scaf", "stap"):
                entry = row[strand_type][index]
                if (entry == blank or int(entry[0]) != helix or
                        int(entry[2]) != helix):
                    valid = False
                    break
            if valid:
                safe.add(index)

        for unused_iteration in range(
                max(1, sum(abs(minimum - value) for value in spans()) * 2)):
            current = spans()
            deficient = next((index for index, value in enumerate(current)
                              if value < minimum), None)
            if deficient is None:
                break
            changed = False

            # First move a deletion out of the deficient arc.  The recipient
            # arc must retain at least the minimum after receiving it.
            deleted = next((index for index in arcs[deficient]
                            if int(row["skip"][index]) < 0), None)
            if deleted is not None:
                donor_arcs = [index for index, value in enumerate(current)
                              if index != deficient and value > minimum]
                for donor in sorted(
                        donor_arcs, key=lambda index: -current[index]):
                    targets = [index for index in arcs[donor]
                               if index in safe and
                               not row["loop"][index] and
                               not row["skip"][index]]
                    if not targets:
                        continue
                    row["skip"][deleted] = 0
                    row["skip"][targets[len(targets) // 2]] = -1
                    moved_deletions += 1
                    changed = True
                    break
            if changed:
                continue

            # Otherwise transfer one insertion from the most generous arc.
            donor_arcs = [index for index, value in enumerate(current)
                          if index != deficient and value > minimum]
            for donor in sorted(
                    donor_arcs, key=lambda index: -current[index]):
                sources = [index for index in arcs[donor]
                           if int(row["loop"][index]) > 0]
                targets = [index for index in arcs[deficient]
                           if index in safe]
                if not sources or not targets:
                    continue
                source = max(sources, key=lambda index: (
                    int(row["loop"][index]), -index))
                # Prefer fusing at an existing insertion, otherwise use the
                # middle safe position of this crossover-free interval.
                existing = [index for index in targets
                            if int(row["loop"][index]) > 0]
                target = (min(existing, key=lambda index: (
                    int(row["loop"][index]), index)) if existing else
                    targets[len(targets) // 2])
                row["loop"][source] = int(row["loop"][source]) - 1
                row["loop"][target] = int(row["loop"][target]) + 1
                moved_insertions += 1
                changed = True
                break
            if not changed:
                break

        for arc_index, value in enumerate(spans()):
            if value < minimum:
                unresolved.append({
                    "helix": helix,
                    "start": positions[arc_index],
                    "end": positions[(arc_index + 1) % len(positions)],
                    "actual_spacing": value})

    return {
        "minimum_spacing": minimum,
        "moved_insertions": moved_insertions,
        "moved_deletions": moved_deletions,
        "unresolved": unresolved}


def _balance_total_staple_crossover_density(
        encoded, minimum=21, lattice=None):
    """Balance actual spacing for each helix-pair/direction crossover family.

    Scaffold and staple crossovers are combined when they connect the same
    neighbouring helices in the same direction.  Crossovers to a different
    neighbour, or in the opposite direction, belong to a different density
    family and do not participate in this 1/21 check.
    """
    blank = [-1, -1, -1, -1]
    moved_insertions = 0
    moved_deletions = 0
    unresolved = []
    final_spacings = []

    for row in encoded.get("vstrands", []):
        helix = int(row["num"])
        size = len(row.get("stap", []))
        if size < 2:
            continue
        groups = defaultdict(set)
        for strand_type in ("scaf", "stap"):
            for index, entry in enumerate(row.get(strand_type, [])):
                incoming = int(entry[0])
                outgoing = int(entry[2])
                if incoming >= 0 and incoming != helix:
                    groups[("incoming", incoming)].add(index)
                if outgoing >= 0 and outgoing != helix:
                    groups[("outgoing", outgoing)].add(index)

        def arc_indices(start, end):
            result = []
            index = (start + 1) % size
            while True:
                result.append(index)
                if index == end:
                    return result
                index = (index + 1) % size

        constraints = []
        for (direction, neighbour), position_set in sorted(groups.items()):
            positions = sorted(position_set)
            if len(positions) < 2:
                continue
            family = (direction, neighbour)
            for start, end in zip(positions, positions[1:] + positions[:1]):
                constraints.append({
                    "direction": direction, "neighbour": neighbour,
                    "family": family, "family_count": len(positions),
                    "start": start, "end": end,
                    "indices": arc_indices(start, end)})
        if not constraints:
            continue

        safe = set()
        for index in range(size):
            valid = True
            for strand_type in ("scaf", "stap"):
                entry = row[strand_type][index]
                if (entry == blank or int(entry[0]) != helix or
                        int(entry[2]) != helix):
                    valid = False
                    break
            if valid:
                safe.add(index)

        def redistribute_two_family_deletions():
            """Solve overlapping incoming/outgoing deletion quotas at once.

            With two phase families, every safe base is an edge between one
            interval in each family.  Selecting deletion sites is therefore
            a small capacitated bipartite-flow problem.  Solving it in one
            pass avoids the plateau where a necessary first move temporarily
            transfers a one-base deficit to the opposite direction.
            """
            family_members = defaultdict(list)
            for constraint_index, item in enumerate(constraints):
                family_members[item["family"]].append(constraint_index)
            if len(family_members) != 2:
                return 0
            deletion_count = sum(
                1 for value in row["skip"] if int(value) < 0)
            if deletion_count <= 0:
                return 0
            deletion_floor = (19 if lattice == "honeycomb" else 24)
            current_invalid = False
            for item in constraints:
                actual = sum(max(
                    0, 1 + int(row["loop"][index]) +
                    int(row["skip"][index]))
                    for index in item["indices"])
                required = (deletion_floor if any(
                    int(row["skip"][index]) < 0
                    for index in item["indices"]) else minimum)
                if actual < required:
                    current_invalid = True
                    break
            if not current_invalid:
                return 0
            families = sorted(family_members)
            left_items = family_members[families[0]]
            right_items = family_members[families[1]]
            def capacity(constraint_index):
                potential = sum(
                    max(0, 1 + int(row["loop"][index]))
                    for index in constraints[constraint_index]["indices"])
                return max(0, potential - deletion_floor)

            left_capacity = dict(
                (index, capacity(index)) for index in left_items)
            right_capacity = dict(
                (index, capacity(index)) for index in right_items)
            if sum(left_capacity.values()) < deletion_count or \
                    sum(right_capacity.values()) < deletion_count:
                return 0

            left_for_position = {}
            right_for_position = {}
            for constraint_index in left_items:
                for index in constraints[constraint_index]["indices"]:
                    left_for_position[index] = constraint_index
            for constraint_index in right_items:
                for index in constraints[constraint_index]["indices"]:
                    right_for_position[index] = constraint_index
            edge_positions = defaultdict(list)
            old_deletions = set(
                index for index, value in enumerate(row["skip"])
                if int(value) < 0)
            for index in sorted(safe):
                if int(row["loop"][index]) > 0 or \
                        index not in left_for_position or \
                        index not in right_for_position:
                    continue
                key = (left_for_position[index], right_for_position[index])
                edge_positions[key].append(index)

            source = ("source",)
            sink = ("sink",)
            residual = defaultdict(dict)

            def add_edge(first, second, amount):
                residual[first][second] = int(amount)
                residual[second].setdefault(first, 0)

            for index, amount in left_capacity.items():
                add_edge(source, ("left", index), amount)
            for (left, right), positions_for_edge in edge_positions.items():
                add_edge(("left", left), ("right", right),
                         len(positions_for_edge))
            for index, amount in right_capacity.items():
                add_edge(("right", index), sink, amount)

            flow = 0
            while flow < deletion_count:
                parent = {source: None}
                queue_nodes = [source]
                for node in queue_nodes:
                    for neighbour, amount in residual[node].items():
                        if amount > 0 and neighbour not in parent:
                            parent[neighbour] = node
                            queue_nodes.append(neighbour)
                            if neighbour == sink:
                                break
                    if sink in parent:
                        break
                if sink not in parent:
                    return 0
                amount = deletion_count - flow
                node = sink
                while parent[node] is not None:
                    amount = min(amount, residual[parent[node]][node])
                    node = parent[node]
                node = sink
                while parent[node] is not None:
                    previous = parent[node]
                    residual[previous][node] -= amount
                    residual[node][previous] += amount
                    node = previous
                flow += amount

            selected = set()
            for (left, right), positions_for_edge in edge_positions.items():
                original_capacity = len(positions_for_edge)
                remaining = residual[("left", left)].get(
                    ("right", right), original_capacity)
                used = original_capacity - remaining
                preferred = sorted(
                    positions_for_edge,
                    key=lambda index: (index not in old_deletions, index))
                selected.update(preferred[:used])
            if len(selected) != deletion_count:
                return 0
            for index in old_deletions:
                row["skip"][index] = 0
            for index in selected:
                row["skip"][index] = -1
            return len(old_deletions - selected)

        moved_deletions += redistribute_two_family_deletions()

        def base_weight(index):
            return max(0, 1 + int(row["loop"][index]) +
                       int(row["skip"][index]))

        def values():
            return [sum(base_weight(index) for index in item["indices"])
                    for item in constraints]

        def requirements():
            # Curvature deletions have a lattice-specific, physically
            # accepted lower spacing: 19 bp Honeycomb and 24 bp Square.
            # Newly split/non-deletion intervals retain the common 21-bp
            # total-crossover density floor.
            deletion_floor = (19 if lattice == "honeycomb" else 24)
            return [
                deletion_floor if any(
                    int(row["skip"][index]) < 0
                    for index in item["indices"]) else minimum
                for item in constraints]

        def score():
            deficits = [max(0, required - value)
                        for required, value in zip(
                            requirements(), values())]
            return sum(deficits), sum(value > 0 for value in deficits)

        # Indel relocation is a local phase adjustment, not a way to stretch
        # a nominal 4--7-base arc into a legal crossover interval.  Repair at
        # most three missing actual bases per interval; larger deficits must
        # be solved by rolling back the newly added crossover.
        local_repair_limit = 3
        seen_indel_states = {(
            tuple(int(value) for value in row["loop"]),
            tuple(int(value) for value in row["skip"]))}
        for unused_iteration in range(max(1, len(constraints) * 3)):
            current_values = values()
            current_requirements = requirements()
            deficient = next((index for index, value in
                              enumerate(current_values)
                              if value < current_requirements[index] and
                              current_requirements[index] - value <=
                              local_repair_limit), None)
            if deficient is None:
                break
            # A crossover family partitions the entire physical helix.  If
            # its actual circumference is already shorter than
            # ``minimum * crossover_count``, merely moving the same indels
            # around can never satisfy every interval.  Detect this before
            # optimization instead of oscillating insertions hundreds of
            # times between mutually impossible arcs.
            deficient_item = constraints[deficient]
            family_indices = [
                index for index, item in enumerate(constraints)
                if item["family"] == deficient_item["family"]]
            family_total = sum(current_values[index]
                               for index in family_indices)
            family_required = sum(current_requirements[index]
                                  for index in family_indices)
            if family_total < family_required:
                # Try another deficient family if one is feasible.  When
                # none is feasible, no indel relocation can improve the
                # hard-rule score on this helix.
                feasible_deficient = next((
                    index for index, value in enumerate(current_values)
                    if value < current_requirements[index] and
                    current_requirements[index] - value <=
                    local_repair_limit and
                    sum(current_values[other] for other, item in
                        enumerate(constraints)
                        if item["family"] == constraints[index]["family"])
                    >= sum(current_requirements[other]
                           for other, item in enumerate(constraints)
                           if item["family"] ==
                           constraints[index]["family"])),
                    None)
                if feasible_deficient is None:
                    break
                deficient = feasible_deficient
            current_score = score()
            bad_indices = set(constraints[deficient]["indices"])
            target_pool = [index for index in
                           constraints[deficient]["indices"]
                           if index in safe]
            if not target_pool:
                break
            existing_targets = [index for index in target_pool
                                if int(row["loop"][index]) > 0]
            insertion_targets = (existing_targets[:3] or
                                 [target_pool[len(target_pool) // 2]])
            best = None

            # Move one deletion out of the deficient family interval.
            deleted_sources = [index for index in bad_indices
                               if int(row["skip"][index]) < 0]
            deletion_targets = [index for index in sorted(safe - bad_indices)
                                if not row["loop"][index] and
                                not row["skip"][index]]
            for source in deleted_sources[:3]:
                # A circular boundary interval can have only a handful of
                # globally safe recipients.  Sampling eight positions may
                # miss all of them (observed as a legal Square 24-bp plan
                # being left at 22 bp).  Helices contain only a few hundred
                # bases, so exhaustive recipient evaluation is both bounded
                # and substantially cheaper than regenerating the design.
                for target in deletion_targets:
                    row["skip"][source] = 0
                    row["skip"][target] = -1
                    trial_score = score()
                    trial_state = (
                        tuple(int(value) for value in row["loop"]),
                        tuple(int(value) for value in row["skip"]))
                    row["skip"][source] = -1
                    row["skip"][target] = 0
                    # Two overlapping direction families can require a
                    # coordinated two-deletion exchange.  The first move may
                    # keep total deficit unchanged while temporarily splitting
                    # it across two arcs; permit that bounded plateau step.
                    if trial_state not in seen_indel_states and \
                            trial_score[0] <= current_score[0] and (
                            best is None or trial_score < best[0]):
                        best = (trial_score, "deletion", source, target,
                                trial_state)

            # Or move one insertion from outside into the deficient arc.
            insertion_sources = [
                index for index in range(size)
                if index not in bad_indices and int(row["loop"][index]) > 0]
            for source in insertion_sources:
                for target in insertion_targets:
                    native_period = 21 if lattice == "honeycomb" else 32
                    if source // native_period != target // native_period:
                        continue
                    row["loop"][source] = int(row["loop"][source]) - 1
                    row["loop"][target] = int(row["loop"][target]) + 1
                    trial_score = score()
                    trial_state = (
                        tuple(int(value) for value in row["loop"]),
                        tuple(int(value) for value in row["skip"]))
                    row["loop"][source] = int(row["loop"][source]) + 1
                    row["loop"][target] = int(row["loop"][target]) - 1
                    if trial_state not in seen_indel_states and \
                            trial_score[0] <= current_score[0] and (
                            best is None or trial_score < best[0]):
                        best = (trial_score, "insertion", source, target,
                                trial_state)

            if best is None:
                break
            unused_score, kind, source, target, accepted_state = best
            if kind == "deletion":
                row["skip"][source] = 0
                row["skip"][target] = -1
                moved_deletions += 1
            else:
                row["loop"][source] = int(row["loop"][source]) - 1
                row["loop"][target] = int(row["loop"][target]) + 1
                moved_insertions += 1
            seen_indel_states.add(accepted_state)

        row_values = values()
        row_requirements = requirements()
        final_spacings.extend(row_values)
        for item, value, required in zip(
                constraints, row_values, row_requirements):
            if value < required:
                unresolved.append({
                    "helix": helix,
                    "neighbour": int(item["neighbour"]),
                    "direction": item["direction"],
                    "start": int(item["start"]),
                    "end": int(item["end"]),
                    "actual_spacing": int(value),
                    "required_spacing": int(required)})

    return {
        "definition": "same helix-pair and crossover direction; "
                      "scaffold plus staple",
        "minimum_spacing": minimum,
        "deletion_interval_minimum_spacing": (
            19 if lattice == "honeycomb" else 24),
        "final_spacing_minimum": (
            min(final_spacings) if final_spacings else 0),
        "final_spacing_maximum": (
            max(final_spacings) if final_spacings else 0),
        "moved_insertions": moved_insertions,
        "moved_deletions": moved_deletions,
        "unresolved": unresolved}


def _remove_added_staple_xover(
        encoded, lattice, chosen, leave_nicks=False,
        merge_previous_records=()):
    """Remove one serialized staple xover.

    ``leave_nicks`` implements half-removal of a local double-crossover site.
    A later full removal passes that first record through
    ``merge_previous_records`` so both removed sites are restored to ordinary
    same-helix staple runs.
    """
    from cadnano2.model.document import Document
    from cadnano2.model.enum import LatticeType
    from cadnano2.model.io.legacydecoder import import_legacy_dict
    from cadnano2.model.io.legacyencoder import legacy_dict_from_part

    lattice_type = (LatticeType.Honeycomb if lattice == "honeycomb"
                    else LatticeType.Square)
    document = Document()
    part = import_legacy_dict(
        document, encoded, lattice_type, forceLatticeType=True)
    if part is None:
        return None
    first, second, first_index, second_index = map(int, chosen)
    strand5p = None
    strand3p = None
    first_vh = part.virtualHelix(first)
    if first_vh is not None:
        for strand in first_vh.stapleStrandSet():
            connection = strand.connection3p()
            if connection is not None and \
                    strand.idx3Prime() == first_index and \
                    connection.virtualHelix().number() == second and \
                    connection.idx5Prime() == second_index:
                strand5p, strand3p = strand, connection
                break
    if strand5p is None:
        return None

    def side_specs(record):
        record_first, record_second, record_first_index, \
            record_second_index = map(int, record)
        result = []
        for helix, index, is_three_prime in (
                (record_first, record_first_index, True),
                (record_second, record_second_index, False)):
            vh = part.virtualHelix(helix)
            if vh is None:
                continue
            strand_set = vh.stapleStrandSet()
            if is_three_prime:
                delta = 1 if strand_set.isDrawn5to3() else -1
            else:
                delta = -1 if strand_set.isDrawn5to3() else 1
            result.append((strand_set, index, delta))
        return result

    sides = [] if leave_nicks else side_specs(chosen)
    for previous in merge_previous_records:
        sides.extend(side_specs(previous))
    part.removeXover(strand5p, strand3p, useUndoStack=False)
    for strand_set, index, delta in sides:
        # After removeXover the two local pieces are represented by bases on
        # either side of the former crossover coordinate.
        left = strand_set.getStrand(index)
        right = strand_set.getStrand(index + delta)
        if left is not None and right is not None and left is not right and \
                strand_set.strandsCanBeMerged(left, right):
            strand_set.mergeStrands(left, right, useUndoStack=False)
    return legacy_dict_from_part(
        part, encoded.get("name", "curved"), includeSequences=False)


def _rollback_unresolved_added_staple_xovers(
        encoded, lattice, curvature_indels, added_records):
    """Retract only newly added xovers whose <21-bp arc cannot be repaired."""
    retained = [tuple(map(int, record)) for record in added_records]
    rolled_back = []
    balance = _balance_total_staple_crossover_density(
        encoded, minimum=21, lattice=lattice)
    while balance["unresolved"]:
        bad = balance["unresolved"][0]
        helix = int(bad["helix"])
        endpoints = {int(bad["start"]), int(bad["end"])}
        chosen = None
        for record in reversed(retained):
            first, second, first_index, second_index = map(int, record)
            if ((first == helix and first_index in endpoints) or
                    (second == helix and second_index in endpoints)):
                chosen = (first, second, first_index, second_index)
                break
        if chosen is None:
            break
        updated = _remove_added_staple_xover(encoded, lattice, chosen)
        if updated is None:
            break
        retained.remove(chosen)
        rolled_back.append(chosen)
        encoded = updated
        _rebalance_indels_against_final_crossovers(
            encoded, lattice, curvature_indels, only_conflicts=False)
        balance = _balance_total_staple_crossover_density(
            encoded, minimum=21, lattice=lattice)
    return encoded, retained, rolled_back, balance


def _rollback_added_xovers_for_autobreak(
        encoded, lattice, curvature_indels, added_records):
    """Make added double-xover sites compatible with staple breaking.

    A nearby pair of newly inserted staple crossovers is treated as one
    double-crossover site.  First remove only one member, which leaves a nick
    beside the surviving crossover.  If an unbreakable staple still touches
    that site, remove the surviving member as well.  Native AutoCS crossover
    records are never candidates for either stage.
    """
    from cadnano2.model.document import Document
    from cadnano2.model.enum import LatticeType, StrandType
    from cadnano2.model.io.legacydecoder import import_legacy_dict
    from cadnano2.model.parts.part import (
        _bestStapleBreakPlan, _crossoverPositionsByHelix,
        _legalStapleNickBoundaries, _stapleOligoBaseRecords)

    lattice_type = (LatticeType.Honeycomb if lattice == "honeycomb"
                    else LatticeType.Square)
    if lattice == "square":
        settings = (7, 16, 8, 32)
        pair_window = 32
    else:
        settings = (6, 14, 7, 28)
        pair_window = 21
    retained = [tuple(map(int, record)) for record in added_records]
    half_removed = []
    fully_removed = []
    pending_pairs = []

    def circular_index_distance(first, second, size):
        return min((first - second) % size, (second - first) % size)

    def companion_for(record, records, size):
        first, second, first_index, unused_second_index = record
        pair = tuple(sorted((first, second)))
        choices = []
        for other in records:
            if other == record or tuple(sorted(other[:2])) != pair:
                continue
            distance = circular_index_distance(
                first_index, int(other[2]), size)
            if 0 < distance <= pair_window:
                choices.append((distance, other))
        return min(choices)[1] if choices else None

    def break_analysis(current):
        document = Document()
        part = import_legacy_dict(
            document, current, lattice_type, forceLatticeType=True)
        if part is None:
            return None, []
        staple_xovers = _crossoverPositionsByHelix(
            part, StrandType.Staple)
        scaffold_xovers = _crossoverPositionsByHelix(
            part, StrandType.Scaffold)
        bad = []
        for oligo in part.oligos():
            if not oligo.isStaple() or oligo.isHybrid():
                continue
            records = _stapleOligoBaseRecords(oligo)
            deletion_dense = sum(
                int(record[3]) == 0 for record in records) >= 2
            plan_arguments = dict(
                preferredMinimum=30,
                preferredMaximum=50,
                targetLength=40,
                terminalMaximum=49,
                preferDeletionDense=deletion_dense)
            if _bestStapleBreakPlan(
                    oligo, staple_xovers, scaffold_xovers,
                    settings[0], settings[1], settings[2],
                    settings[3],
                    hardMaximum=(60 if deletion_dense else 57),
                    **plan_arguments) is not None:
                continue
            # The 64-nt ceiling is a true second-pass exception: it is
            # considered only after the complete normal 21--57 nt problem
            # (or deletion-dense 40--60 nt preference) is unsatisfiable.
            if _bestStapleBreakPlan(
                    oligo, staple_xovers, scaffold_xovers,
                    settings[0], settings[1], settings[2],
                    settings[3],
                    hardMaximum=64, **plan_arguments) is not None:
                continue
            if 58 <= oligo.actualLength() <= 64:
                if not oligo.isLoop() or _legalStapleNickBoundaries(
                        records, staple_xovers, scaffold_xovers,
                        settings[0], settings[2]):
                    continue
            relaxed_plan = _bestStapleBreakPlan(
                oligo, staple_xovers, scaffold_xovers,
                settings[0], settings[1], settings[2],
                settings[3],
                ignoreIndels=True,
                hardMaximum=(60 if deletion_dense else 57),
                **plan_arguments)
            bad.append({
                "bases": set((int(helix), int(index))
                             for helix, index, unused_strand,
                             unused_length in records),
                "relaxed_plan": relaxed_plan})
        return part, bad

    def unbreakable_bases(current):
        part, details = break_analysis(current)
        return part, [item["bases"] for item in details]

    break_optimization = {
        "attempts": 0, "accepted_passes": 0,
        "moved_insertions": 0, "moved_deletions": 0,
        "resolved_staples": 0}

    def optimize_indels_for_break(current):
        """Clear relaxed-plan nick endpoints within native-period budgets.

        The strict plan and the otherwise-identical ``ignoreIndels`` plan
        provide a direct break-feasibility score.  If only indels block a
        legal plan, move those indels inside the same original crossover
        period, retain the exact period quota, and accept the placement only
        when it reduces the number of unbreakable staples without creating a
        crossover-density violation.
        """
        unused_part, details = break_analysis(current)
        actionable = [item for item in details if item["relaxed_plan"]]
        if not actionable:
            return current
        break_optimization["attempts"] += 1
        trial = json.loads(json.dumps(current))
        rows = {int(row["num"]): row
                for row in trial.get("vstrands", [])}
        reserved = defaultdict(set)
        for item in actionable:
            for helix, upper_index in item["relaxed_plan"]:
                reserved[int(helix)].update(
                    (int(upper_index) - 1, int(upper_index)))
        native_period = 21 if lattice == "honeycomb" else 32
        preferred_phase = 7 if lattice == "honeycomb" else 8
        moved_insertions = moved_deletions = 0

        for helix, blocked in reserved.items():
            row = rows.get(helix)
            if row is None:
                continue
            size = len(row.get("stap", []))

            def safe_target(source):
                period_index = source // native_period
                lower = period_index * native_period
                upper = min(size, lower + native_period)
                candidates = []
                for index in range(max(1, lower), min(size - 1, upper)):
                    if index in blocked or row["loop"][index] or \
                            row["skip"][index]:
                        continue
                    valid = True
                    for strand_type in ("scaf", "stap"):
                        entry = row[strand_type][index]
                        if entry == [-1, -1, -1, -1] or \
                                int(entry[0]) != helix or \
                                int(entry[2]) != helix:
                            valid = False
                            break
                    if valid:
                        candidates.append(index)
                if not candidates:
                    return None
                return min(candidates, key=lambda index: (
                    index % preferred_phase in (0, preferred_phase - 1),
                    min(abs(index - source), size - abs(index - source)),
                    index))

            for source in sorted(blocked):
                if not 0 <= source < size:
                    continue
                loop_value = max(0, int(row["loop"][source]))
                skip_value = int(row["skip"][source])
                if not loop_value and skip_value >= 0:
                    continue
                target = safe_target(source)
                if target is None:
                    continue
                if loop_value:
                    row["loop"][source] = 0
                    row["loop"][target] = \
                        int(row["loop"][target]) + loop_value
                    moved_insertions += loop_value
                if skip_value < 0:
                    row["skip"][source] = 0
                    row["skip"][target] = -1
                    moved_deletions += 1

        def synchronize_records(trial_rows):
            """Keep the curvature plan aligned with accepted encoded sites."""
            for record in curvature_indels.get("rings", []):
                helix = int(record["helix"])
                row = trial_rows.get(helix)
                if row is None:
                    continue
                nominal_size = min(
                    int(record.get("nominal_bases", len(row["loop"]))),
                    len(row["loop"]))
                record["insertions"] = [
                    index for index in range(nominal_size)
                    for unused_copy in range(
                        max(0, int(row["loop"][index])))]
                record["deletions"] = [
                    index for index in range(nominal_size)
                    if int(row["skip"][index]) < 0]

        if moved_insertions or moved_deletions:
            density = _balance_total_staple_crossover_density(
                trial, minimum=21, lattice=lattice)
            unused_part, trial_details = break_analysis(trial)
            if not density["unresolved"] and \
                    len(trial_details) < len(details):
                synchronize_records(rows)
                break_optimization["accepted_passes"] += 1
                break_optimization["moved_insertions"] += moved_insertions
                break_optimization["moved_deletions"] += moved_deletions
                break_optimization["resolved_staples"] += (
                    len(details) - len(trial_details))
                return trial

        # If actual segment lengths—not endpoint occupancy—prevent the
        # relaxed plan, score bounded one-base insertion transfers directly.
        # Moving a unit out of an unbreakable staple (or into a <21-nt one)
        # changes break length while remaining inside the same native period.
        base_rows = {int(row["num"]): row
                     for row in current.get("vstrands", [])}
        candidates = []
        for detail in details:
            bad_bases = detail["bases"]
            by_helix = defaultdict(set)
            for helix, index in bad_bases:
                by_helix[int(helix)].add(int(index))
            for helix, indices in by_helix.items():
                row = base_rows.get(helix)
                if row is None:
                    continue
                size = len(row["loop"])
                sources = [index for index in indices
                           if int(row["loop"][index]) > 0]
                for source in sources:
                    lower = (source // native_period) * native_period
                    upper = min(size, lower + native_period)
                    targets = []
                    for target in range(max(1, lower), min(size - 1, upper)):
                        if target in indices or row["loop"][target] or \
                                row["skip"][target]:
                            continue
                        valid = True
                        for strand_type in ("scaf", "stap"):
                            entry = row[strand_type][target]
                            if entry == [-1, -1, -1, -1] or \
                                    int(entry[0]) != helix or \
                                    int(entry[2]) != helix:
                                valid = False
                                break
                        if valid:
                            targets.append(target)
                    for target in sorted(targets, key=lambda index: (
                            index % preferred_phase in
                            (0, preferred_phase - 1),
                            abs(index - source), index))[:4]:
                        candidates.append((helix, source, target))
        best = None
        for helix, source, target in candidates[:32]:
            break_optimization["attempts"] += 1
            candidate = json.loads(json.dumps(current))
            candidate_rows = {int(row["num"]): row
                              for row in candidate.get("vstrands", [])}
            candidate_row = candidate_rows[helix]
            candidate_row["loop"][source] = \
                int(candidate_row["loop"][source]) - 1
            candidate_row["loop"][target] = \
                int(candidate_row["loop"][target]) + 1
            density = _balance_total_staple_crossover_density(
                candidate, minimum=21, lattice=lattice)
            if density["unresolved"]:
                continue
            unused_part, candidate_details = break_analysis(candidate)
            score = len(candidate_details)
            if score < len(details) and (
                    best is None or score < best[0]):
                best = (score, candidate, candidate_rows)
        if best is None:
            return current
        score, candidate, candidate_rows = best
        synchronize_records(candidate_rows)
        break_optimization["accepted_passes"] += 1
        break_optimization["moved_insertions"] += 1
        break_optimization["resolved_staples"] += len(details) - score
        return candidate

    # Complete-crossover topology gets the first break-aware insertion pass.
    encoded = optimize_indels_for_break(encoded)

    maximum_attempts = max(1, len(retained) * 2)
    for unused_attempt in range(maximum_attempts):
        part, bad_sets = unbreakable_bases(encoded)
        if part is None or not bad_sets:
            break
        size = part.maxBaseIdx() + 1

        # Second stage: a previously half-removed double-xover site still
        # participates in an invalid staple.  Remove its surviving half.
        chosen = None
        completed_pair = None
        for removed, companion in pending_pairs:
            if companion not in retained:
                continue
            endpoints = {(companion[0], companion[2]),
                         (companion[1], companion[3])}
            if any(endpoints & bases for bases in bad_sets):
                chosen = companion
                completed_pair = (removed, companion)
                break

        # First stage: break only one half of a newly added local pair.
        if chosen is None:
            for record in reversed(retained):
                endpoints = {(record[0], record[2]),
                             (record[1], record[3])}
                if not any(endpoints & bases for bases in bad_sets):
                    continue
                companion = companion_for(record, retained, size)
                if companion is not None:
                    chosen = record
                    pending_pairs.append((record, companion))
                    break
            # A standalone added crossover has no distinct half to retain;
            # removing it is already the complete fallback.
            if chosen is None:
                for record in reversed(retained):
                    endpoints = {(record[0], record[2]),
                                 (record[1], record[3])}
                    if any(endpoints & bases for bases in bad_sets):
                        chosen = record
                        completed_pair = (record,)
                        break
        if chosen is None:
            break

        updated = _remove_added_staple_xover(
            encoded, lattice, chosen,
            leave_nicks=completed_pair is None,
            merge_previous_records=(
                (completed_pair[0],) if completed_pair is not None and
                len(completed_pair) == 2 else ()))
        if updated is None:
            break
        encoded = updated
        retained.remove(chosen)
        if completed_pair is None:
            half_removed.append(chosen)
        else:
            fully_removed.append(completed_pair)
            if len(completed_pair) == 2:
                pending_pairs.remove(completed_pair)
        _rebalance_indels_against_final_crossovers(
            encoded, lattice, curvature_indels, only_conflicts=False)
        _balance_total_staple_crossover_density(
            encoded, minimum=21, lattice=lattice)
        # The half-site nick is an explicit exception to ordinary crossover
        # clearance.  Re-score break feasibility and redistribute only within
        # original-period budgets before deciding whether its companion must
        # also be removed.
        encoded = optimize_indels_for_break(encoded)

    unused_part, remaining_bad = unbreakable_bases(encoded)
    return (encoded, retained, half_removed, fully_removed,
            len(remaining_bad), break_optimization)


def _rebalance_indels_against_final_crossovers(
        encoded, lattice, curvature_indels, only_conflicts=False):
    """Distribute the fixed indel budget uniformly over 7/8-bp domains.

    Crossovers and nicks are immutable inputs to this pass.  Every helix keeps
    its exact target circumference, while its insertions or deletions are
    spread over the common parent axis.  The Dietz-style local strain guard is
    a hard maximum of three inserted/deleted bases in one native domain.
    """
    domain_size = 7 if lattice == "honeycomb" else 8
    maximum_per_domain = 3
    rows = {int(row["num"]): row for row in encoded.get("vstrands", [])}
    blank = [-1, -1, -1, -1]

    # Deleting even one base from a nominal 21-position staple component
    # makes it violate the hard 21-nt lower bound.  Map those components once
    # and keep deletion sites on longer products; no crossover must move.
    deletion_protected = defaultdict(set)
    staple_nodes = set()
    for helix_id, row in rows.items():
        for index, entry in enumerate(row.get("stap", [])):
            if entry != blank:
                staple_nodes.add((helix_id, index))
    unseen = set(staple_nodes)
    while unseen:
        seed = unseen.pop()
        component = {seed}
        stack = [seed]
        while stack:
            helix_id, index = stack.pop()
            entry = rows[helix_id]["stap"][index]
            for neighbour in ((int(entry[0]), int(entry[1])),
                              (int(entry[2]), int(entry[3]))):
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    component.add(neighbour)
                    stack.append(neighbour)
        if len(component) <= 21:
            for helix_id, index in component:
                deletion_protected[helix_id].add(index)

    relocated_conflicts = 0
    maximum_observed = 0
    records = curvature_indels.get("rings", [])
    for record_order, record in enumerate(records):
        helix_id = int(record["helix"])
        row = rows.get(helix_id)
        if row is None:
            continue
        nominal_size = int(record["nominal_bases"])
        insertion_count = len(record.get("insertions", []))
        deletion_count = len(record.get("deletions", []))
        if not insertion_count and not deletion_count:
            continue
        if insertion_count and deletion_count:
            raise RuntimeError(
                "Helix %d contains both curvature insertions and deletions."
                % helix_id)
        count = insertion_count or deletion_count
        is_insertion = bool(insertion_count)
        domain_count = max(1, int(math.ceil(
            nominal_size / float(domain_size))))
        if count > domain_count * maximum_per_domain:
            raise RuntimeError(
                "Helix %d requires %d indels across %d %d-bp domains; "
                "the hard limit is +/-3 per domain. Reduce curvature, "
                "increase the selected span, or change the cross-section."
                % (helix_id, count, domain_count, domain_size))

        safe_by_domain = defaultdict(list)
        for index in range(1, nominal_size - 1):
            valid = True
            for strand_type in ("scaf", "stap"):
                entry = row[strand_type][index]
                if (entry == blank or int(entry[0]) != helix_id or
                        int(entry[2]) != helix_id):
                    valid = False
                    break
            if valid:
                domain = min(domain_count - 1, index // domain_size)
                if (is_insertion or
                        index not in deletion_protected[helix_id]):
                    safe_by_domain[domain].append(index)

        capacities = []
        for domain in range(domain_count):
            candidates = safe_by_domain.get(domain, [])
            capacities.append(
                maximum_per_domain if is_insertion and candidates else
                min(maximum_per_domain, len(candidates)))
        if sum(capacities) < count:
            raise RuntimeError(
                "Helix %d has only %d safe domain indel positions for %d "
                "required indels after protecting crossovers, nicks and "
                "short staples." % (helix_id, sum(capacities), count))

        previous_insertions = list(record.get("insertions", []))
        previous_deletions = list(record.get("deletions", []))
        for index in range(nominal_size):
            row["loop"][index] = 0
            row["skip"][index] = 0
        try:
            selected_records = equal_partition_indel_sites(
                count, 0, nominal_size - 1, domain_size,
                dict((domain, safe_by_domain.get(domain, []))
                     for domain in range(domain_count)),
                dict(enumerate(capacities)),
                allow_repeated_sites=is_insertion)
        except TwistBendError as error:
            raise RuntimeError(
                "Helix %d cannot satisfy equal-partition indel placement: %s"
                % (helix_id, error))
        selected = [int(item['idx']) for item in selected_records]
        quotas = [0] * domain_count
        for item in selected_records:
            quotas[int(item['domain'])] += 1

        if is_insertion:
            for index in selected:
                row["loop"][index] = int(row["loop"][index]) + 1
            record["insertions"] = sorted(selected)
            record["deletions"] = []
        else:
            for index in selected:
                row["skip"][index] = -1
            record["insertions"] = []
            record["deletions"] = sorted(selected)
        record["domain_size_bp"] = domain_size
        record["domain_indel_quotas"] = quotas
        record["indel_placements"] = [{
            "base": int(item["idx"]),
            "value": 1 if is_insertion else -1,
            "partition": int(item["assignment"]["partition"]),
            "partition_count": int(count),
            "ideal_base": float(item["assignment"]["target"]),
            "partition_start": int(item["assignment"][
                "partition_start"]),
            "partition_end": int(item["assignment"]["partition_end"]),
            "target_domain": int(item["assignment"]["domain"]),
            "domain": int(item["domain"]),
        } for item in selected_records]
        record["indel_distribution_method"] = (
            "equal-partition-and-native-domain-intersection-first; "
            "nearest-safe-site repair; no forced stagger")
        record["maximum_indel_in_one_domain"] = max(quotas or [0])
        maximum_observed = max(maximum_observed, max(quotas or [0]))
        relocated_conflicts += sum(
            1 for old, new in zip(
                sorted(previous_insertions or previous_deletions),
                sorted(selected)) if old != new)

    return {
        "relocated_conflicts": relocated_conflicts,
        "fused_insertions": 0,
        "domain_size_bp": domain_size,
        "maximum_indel_per_domain_allowed": maximum_per_domain,
        "maximum_indel_per_domain_observed": maximum_observed,
        "only_conflicts_requested": bool(only_conflicts)}


def _close_circular_staple_boundaries(encoded, curvature_indels):
    """Join AutoCS' artificial index-0 staple nick on every DNA ring."""
    rows = {int(row["num"]): row for row in encoded.get("vstrands", [])}
    closed = 0
    for record in curvature_indels.get("rings", []):
        helix_id = int(record["helix"])
        nominal_size = int(record["nominal_bases"])
        row = rows.get(helix_id)
        if row is None or nominal_size < 2:
            continue
        strand = row.get("stap", [])
        if len(strand) < nominal_size:
            continue
        first = strand[0]
        last = strand[nominal_size - 1]
        if first == [-1, -1, -1, -1] or last == [-1, -1, -1, -1]:
            continue
        # Forward-drawn staple: N-1 is the 3' end and 0 is the 5' end.
        if (int(last[2]) < 0 and int(first[0]) < 0):
            last[2:4] = [helix_id, 0]
            first[0:2] = [helix_id, nominal_size - 1]
            closed += 1
        # Reverse-drawn staple: 0 is the 3' end and N-1 is the 5' end.
        elif (int(first[2]) < 0 and int(last[0]) < 0):
            first[2:4] = [helix_id, nominal_size - 1]
            last[0:2] = [helix_id, 0]
            closed += 1
    return closed


def _add_curved_long_interval_staple_crossovers(
        part, lattice, allowed_pairs=None, curvature_indels=None):
    """Split long actual-bp staple intervals before Autobreak.

    Candidates are native lattice staple sites.  Same-direction intervals
    retain the 21-bp floor, with the established curvature-deletion exception
    of 19 bp Honeycomb / 24 bp Square.  Opposite-direction sites on one helix
    retain the inclusive 8-bp Square / 7-bp Honeycomb clearance.  Indels are
    intentionally not a rejection criterion: the final rebalance pass
    relocates them after the crossover set is complete.
    """
    from cadnano2.model.enum import StrandType
    from cadnano2.model.parts.part import (
        _bestStapleBreakPlan, _crossoverPositionsByHelix,
        _stapleOligoBaseRecords)

    allowed_pairs = (set(tuple(sorted(pair)) for pair in allowed_pairs)
                     if allowed_pairs is not None else None)
    nominal_size = part.maxBaseIdx() + 1
    opposite_minimum = 8 if lattice == "square" else 7
    threshold = 42

    indel_records = dict(
        (int(record["helix"]), record)
        for record in (curvature_indels or {}).get("rings", []))
    weights = {}
    same_direction_minimum = {}
    deletion_floor = 19 if lattice == "honeycomb" else 24
    for vh in part.getVirtualHelices():
        record = indel_records.get(vh.number(), {})
        insertion_count = len(record.get("insertions", []))
        deletion_count = len(record.get("deletions", []))
        projected_scale = (
            (nominal_size + insertion_count - deletion_count) /
            float(nominal_size))
        # Candidate selection uses the pre-design average insertion budget,
        # not the provisional physical positions.  A native period with e.g.
        # 21 planned insertions is therefore evaluated as 42 actual bases
        # before any insertion site has been chosen.
        weights[vh.number()] = [projected_scale] * nominal_size
        same_direction_minimum[vh.number()] = (
            deletion_floor if deletion_count else 21)

    def arc_distance(helix, first, second):
        values = weights[helix]
        total = 0
        index = (first + 1) % nominal_size
        while True:
            total += values[index]
            if index == second:
                return total
            index = (index + 1) % nominal_size

    def circular_distance(helix, first, second):
        return min(arc_distance(helix, first, second),
                   arc_distance(helix, second, first))

    def role_positions():
        positions = defaultdict(set)
        for vh in part.getVirtualHelices():
            for strand in vh.stapleStrandSet():
                connection = strand.connection3p()
                if connection is None or \
                        connection.virtualHelix() is vh:
                    continue
                first = vh.number()
                second = connection.virtualHelix().number()
                positions[(first, "out")].add(strand.idx3Prime())
                positions[(first, "out", second)].add(
                    strand.idx3Prime())
                positions[(second, "in")].add(connection.idx5Prime())
                positions[(second, "in", first)].add(
                    connection.idx5Prime())
        return positions

    scaffold_positions = defaultdict(set)
    scaffold_role_positions = defaultdict(set)
    for first, second, first_index, second_index in \
            part._existingScaffoldCrossoverRecords():
        scaffold_positions[first].add(first_index)
        scaffold_positions[second].add(second_index)
        scaffold_role_positions[(first, "out", second)].add(first_index)
        scaffold_role_positions[(second, "in", first)].add(second_index)

    candidates = set()
    for vh in part.getVirtualHelices():
        staple_set = vh.stapleStrandSet()
        drawn_5_to_3 = staple_set.isDrawn5to3()
        for neighbour, index, strand_type, is_low in \
                part.potentialCrossoverList(vh):
            if strand_type != StrandType.Staple:
                continue
            pair = tuple(sorted((vh.number(), neighbour.number())))
            if allowed_pairs is not None and pair not in allowed_pairs:
                continue
            from_is_5p = ((is_low and drawn_5_to_3) or
                          (not is_low and not drawn_5_to_3))
            if from_is_5p and \
                    index not in scaffold_positions[vh.number()] and \
                    index not in scaffold_positions[neighbour.number()]:
                candidates.add((vh.number(), neighbour.number(), index))

    added = []
    while candidates:
        positions = role_positions()
        options = []
        invalid = []
        for first, second, index in sorted(candidates):
            role_specs = ((first, "out", second),
                          (second, "in", first))
            split_data = []
            legal = True
            for helix, role, neighbour in role_specs:
                opposite = positions.get(
                    (helix, "in" if role == "out" else "out"), ())
                pair_role_positions = set(
                    scaffold_role_positions.get(
                        (helix, role, neighbour), set()))
                pair_role_positions.update(
                    positions.get((helix, role, neighbour), set()))
                # Long-interval reinforcement is defined by the combined
                # scaffold+staple crossover family for this same helix pair
                # and direction.  Measuring only the existing staple sites
                # can select a midpoint just 4--7 bases from a scaffold
                # crossover and then force a needless rollback.
                same = sorted(pair_role_positions)
                resulting_count = len(pair_role_positions | {index})
                if sum(weights[helix]) < \
                        same_direction_minimum[helix] * resulting_count:
                    legal = False
                    break
                if any(circular_distance(helix, index, other) <
                       opposite_minimum for other in opposite):
                    legal = False
                    break
                if not same:
                    total = sum(weights[helix])
                    if total < threshold:
                        legal = False
                        break
                    split_data.append((total, None, None))
                    continue
                previous = min(
                    same,
                    key=lambda value: arc_distance(helix, value, index))
                following = min(
                    same,
                    key=lambda value: arc_distance(helix, index, value))
                left = arc_distance(helix, previous, index)
                right = arc_distance(helix, index, following)
                total = left + right
                # Do not proportionally assign the period's insertion budget
                # to the two candidate sides here.  A legal native site can
                # split a Honeycomb coordinate interval 7/14; the fixed
                # insertion quota is allowed to be distributed unevenly so
                # both final actual sides reach 21 bases.  Only the complete
                # original interval must have at least 42 actual bases now.
                if total < threshold:
                    legal = False
                    break
                split_data.append((total, left, right))
            if not legal:
                invalid.append((first, second, index))
                continue

            concrete = [(left, right) for unused_total, left, right
                        in split_data if left is not None]
            if lattice == "square" and concrete:
                minimum_side = min(min(left, right)
                                   for left, right in concrete)
                maximum_side = min(max(left, right)
                                   for left, right in concrete)
                tier = (3 if minimum_side >= 32 else
                        2 if minimum_side >= 21 and
                        maximum_side >= 32 else 1)
            else:
                tier = 1
            balance = sum(abs(left - right) for left, right in concrete)
            options.append((
                -tier,
                -min(total for total, unused_left, unused_right
                     in split_data),
                balance, index, first, second, split_data))
        for candidate in invalid:
            candidates.discard(candidate)
        if not options:
            break
        unused_tier, unused_total, unused_balance, index, first, second, \
            unused_split_data = min(options)
        candidates.discard((first, second, index))
        first_vh = part.virtualHelix(first)
        second_vh = part.virtualHelix(second)
        first_strand = first_vh.stapleStrandSet().getStrand(index)
        second_strand = second_vh.stapleStrandSet().getStrand(index)
        if first_strand is None or second_strand is None or \
                first_strand.hasXoverAt(index) or \
                second_strand.hasXoverAt(index):
            continue
        part.createXover(
            first_strand, index, second_strand, index,
            useUndoStack=False)
        created = any(
            strand.idx3Prime() == index and
            strand.connection3p() is not None and
            strand.connection3p().virtualHelix().number() == second
            for strand in first_vh.stapleStrandSet())
        if not created:
            continue
        added.append((first, second, index, index))
    return added


def _relocate_scaffold_nick_away_from_added_xovers(encoded, added_records):
    """Move the unique scaffold nick off a newly added staple crossover."""
    if not added_records:
        return None
    rows = {int(row["num"]): row for row in encoded.get("vstrands", [])}
    conflict_positions = defaultdict(set)
    for first, second, first_index, second_index in added_records:
        conflict_positions[int(first)].add(int(first_index))
        conflict_positions[int(second)].add(int(second_index))
    five_prime = []
    three_prime = []
    for helix, row in rows.items():
        for index, entry in enumerate(row.get("scaf", [])):
            if entry == [-1, -1, -1, -1]:
                continue
            if int(entry[0]) < 0:
                five_prime.append((helix, index))
            if int(entry[2]) < 0:
                three_prime.append((helix, index))
    if len(five_prime) != 1 or len(three_prime) != 1:
        raise RuntimeError(
            "Curved scaffold must have one nick before relocation.")
    old_five, old_three = five_prime[0], three_prime[0]
    if old_five[1] not in conflict_positions[old_five[0]] and \
            old_three[1] not in conflict_positions[old_three[0]]:
        return None

    candidates = []
    for helix, row in rows.items():
        size = len(row.get("scaf", []))
        for index, entry in enumerate(row.get("scaf", [])):
            if int(entry[2]) != helix:
                continue
            target = int(entry[3])
            if target < 0 or target >= size or \
                    min((target - index) % size,
                        (index - target) % size) != 1:
                continue
            endpoints = (index, target)
            unsafe = False
            for endpoint in endpoints:
                if row["loop"][endpoint] or row["skip"][endpoint]:
                    unsafe = True
                    break
                for strand_type in ("scaf", "stap"):
                    connection = row[strand_type][endpoint]
                    if ((int(connection[0]) >= 0 and
                         int(connection[0]) != helix) or
                            (int(connection[2]) >= 0 and
                             int(connection[2]) != helix)):
                        unsafe = True
                        break
                if unsafe:
                    break
            if unsafe:
                continue
            distance = (min(abs(index - old_three[1]),
                            size - abs(index - old_three[1]))
                        if helix == old_three[0] else size + helix)
            candidates.append((distance, helix, index, target))
    if not candidates:
        raise RuntimeError(
            "No safe longitudinal edge can receive the scaffold nick.")
    unused_distance, helix, source, target = min(candidates)

    old_three_entry = rows[old_three[0]]["scaf"][old_three[1]]
    old_five_entry = rows[old_five[0]]["scaf"][old_five[1]]
    old_three_entry[2:4] = [old_five[0], old_five[1]]
    old_five_entry[0:2] = [old_three[0], old_three[1]]
    source_entry = rows[helix]["scaf"][source]
    target_entry = rows[helix]["scaf"][target]
    source_entry[2:4] = [-1, -1]
    target_entry[0:2] = [-1, -1]
    return {"helix": helix, "source": source, "target": target}


def _heal_short_curved_staple_nicks(part):
    """Join a neighbouring nick when a periodic edge product is <21 nt.

    This never removes or moves a crossover.  It only restores a native
    same-helix phosphodiester connection at an existing staple nick, after
    which Autobreak can repartition the combined product normally.
    """
    from cadnano2.model.parts.part import _existingStapleNickBoundaries

    healed = []
    while True:
        short_oligos = set(
            oligo for oligo in part.oligos()
            if oligo.isStaple() and not oligo.isHybrid() and
            oligo.strand5p() is not None and oligo.actualLength() < 21)
        if not short_oligos:
            break
        candidates = []
        for helix, upper_index in _existingStapleNickBoundaries(part):
            vh = part.virtualHelix(helix)
            strand_set = vh.stapleStrandSet()
            lower = strand_set.getStrand(upper_index - 1)
            upper = strand_set.getStrand(upper_index)
            if lower is None or upper is None or lower is upper or \
                    not strand_set.strandsCanBeMerged(lower, upper):
                continue
            lower_oligo, upper_oligo = lower.oligo(), upper.oligo()
            involved = [oligo for oligo in (lower_oligo, upper_oligo)
                        if oligo in short_oligos]
            if not involved or lower_oligo is upper_oligo:
                continue
            combined = (lower_oligo.actualLength() +
                        upper_oligo.actualLength())
            # Prefer a direct 21--58 nt repair, then the shortest product
            # Autobreak can safely split further.
            score = (combined > 58, combined < 21,
                     abs(combined - 40), combined, helix, upper_index)
            priority = (lower if lower_oligo in short_oligos else upper)
            other = upper if priority is lower else lower
            candidates.append((score, strand_set, priority, other,
                               helix, upper_index, combined))
        if not candidates:
            break
        unused_score, strand_set, priority, other, helix, upper_index, \
            combined = min(candidates, key=lambda item: item[0])
        strand_set.mergeStrands(priority, other, useUndoStack=False)
        healed.append((helix, upper_index, combined))
    return healed


def _precondition_indel_blocked_staple_nicks(part, lattice):
    """Create hard-valid nicks whose only blocker is an indel endpoint."""
    from cadnano2.model.enum import StrandType
    from cadnano2.model.parts.part import (
        _bestStapleBreakPlan, _crossoverPositionsByHelix,
        _stapleOligoBaseRecords)

    if lattice == "square":
        settings = (7, 16, 8, 32)
    else:
        settings = (6, 14, 7, 28)
    staple_xovers = _crossoverPositionsByHelix(part, StrandType.Staple)
    scaffold_xovers = _crossoverPositionsByHelix(part, StrandType.Scaffold)
    planned = []
    for oligo in list(part.oligos()):
        if not oligo.isStaple() or oligo.isHybrid() or \
                oligo.strand5p() is None:
            continue
        records = _stapleOligoBaseRecords(oligo)
        deletion_dense = sum(
            int(record[3]) == 0 for record in records) >= 2
        plan_arguments = dict(
            preferredMinimum=30,
            preferredMaximum=50,
            targetLength=40,
            terminalMaximum=49,
            preferDeletionDense=deletion_dense)
        ordinary = _bestStapleBreakPlan(
            oligo, staple_xovers, scaffold_xovers,
            settings[0], settings[1], settings[2],
            settings[3],
            hardMaximum=(60 if deletion_dense else 57),
            **plan_arguments)
        if ordinary is not None:
            continue
        exceptional = _bestStapleBreakPlan(
            oligo, staple_xovers, scaffold_xovers,
            settings[0], settings[1], settings[2],
            settings[3],
            hardMaximum=64, **plan_arguments)
        if exceptional is not None or 58 <= oligo.actualLength() <= 64:
            continue
        fallback = _bestStapleBreakPlan(
            oligo, staple_xovers, scaffold_xovers,
            settings[0], settings[1], settings[2],
            settings[3],
            ignoreIndels=True,
            hardMaximum=(60 if deletion_dense else 57),
            **plan_arguments)
        if fallback is not None:
            planned.extend(fallback)

    applied = []
    for helix, upper_index in sorted(set(planned)):
        vh = part.virtualHelix(helix)
        if vh is None:
            continue
        strand_set = vh.stapleStrandSet()
        lower_index = upper_index - 1
        strand = strand_set.getStrand(lower_index)
        if strand is None or strand is not strand_set.getStrand(upper_index):
            continue
        split_index = (lower_index if strand_set.isDrawn5to3()
                       else upper_index)
        if strand_set.splitStrand(
                strand, split_index, useUndoStack=False):
            applied.append((helix, upper_index))
    return applied


def _normalize_curved_scaffold_crossovers(part, curvature_indels,
                                          layers=1,
                                          density_mode="periodic",
                                          density_multiple=1):
    """Replace geometry-only DNAxiS scaffold links with native lattice links.

    The imported route already supplies the intended single-scaffold
    topology.  Remove its unequal-index links, heal only their longitudinal
    cuts, and recreate the same directed links at nearest native cadnano
    phases.  Periodic designs then receive reciprocal blocks at the requested
    lattice-period multiple, retained only when the design remains one open
    scaffold oligo.  Minimum-density mode keeps only the normalized routing
    links already required for that topology.
    """
    from collections import defaultdict
    import contextlib
    import io
    from cadnano2.model.enum import StrandType
    from cadnano2.model.parts.part import (
        _filterAutoScaffoldCandidatesForPaths,
        _mergeRemovedScaffoldXoverBoundaries)

    old_records = [record for record in
                   part._existingScaffoldCrossoverRecords()
                   if record[0] != record[1]]
    if not old_records:
        raise RuntimeError("Curved scaffold contains no crossover route.")

    endpoints = set(
        (number, index)
        for first, second, first_index, second_index in old_records
        for number, index in ((first, first_index),
                              (second, second_index)))
    boundaries = set()
    for vh in part.getVirtualHelices():
        strands = sorted(list(vh.scaffoldStrandSet()),
                         key=lambda strand: strand.lowIdx())
        for left, right in zip(strands, strands[1:]):
            if left.highIdx() + 1 == right.lowIdx() and (
                    (vh.number(), left.highIdx()) in endpoints or
                    (vh.number(), right.lowIdx()) in endpoints):
                boundaries.add(
                    (vh.number(), left.highIdx(), right.lowIdx()))

    removed = 0
    for first, second, first_index, second_index in old_records:
        vh = part.virtualHelix(first)
        strand5p = (vh.scaffoldStrandSet().getStrand(first_index)
                    if vh is not None else None)
        strand3p = (strand5p.connection3p()
                    if strand5p is not None else None)
        if strand3p is None or \
                strand3p.virtualHelix().number() != second or \
                strand3p.idx5Prime() != second_index:
            continue
        part.removeXover(strand5p, strand3p, useUndoStack=False)
        removed += 1
    if removed != len(old_records):
        raise RuntimeError(
            "Could not remove every provisional scaffold crossover.")
    _mergeRemovedScaffoldXoverBoundaries(part, boundaries)

    nominal_sizes = {}
    for record in curvature_indels.get("rings", []):
        helix = int(record["helix"])
        nominal_sizes[helix] = int(record["nominal_bases"])
    density_mode = str(density_mode or "periodic").lower()
    if density_mode not in ("periodic", "minimum"):
        raise ValueError("Unknown Curved scaffold crossover density mode.")
    density_multiple = max(1, int(density_multiple or 1))
    requested_spacing = (part._step * density_multiple
                         if density_mode == "periodic" else None)

    def all_candidates():
        candidates = []
        seen = set()
        for vh in part.getVirtualHelices():
            strand_set = vh.scaffoldStrandSet()
            drawn_5_to_3 = strand_set.isDrawn5to3()
            for neighbor, index, strand_type, is_low in \
                    part.potentialCrossoverList(vh):
                # Scaffold topology has priority over provisional curvature
                # indels.  The final rebalance step relocates indels away
                # from every chosen scaffold/staple crossover; filtering a
                # native phase here would create false gaps in the periodic
                # route before that relocation can happen.
                if strand_type != StrandType.Scaffold:
                    continue
                from_is_5p = ((is_low and drawn_5_to_3) or
                              (not is_low and not drawn_5_to_3))
                if not from_is_5p:
                    continue
                record = (vh.number(), neighbor.number(), index, index)
                if record not in seen:
                    seen.add(record)
                    candidates.append(record)
        return candidates

    candidates = all_candidates()
    avoid_modulo = 7 if part._step == 21 else 8

    def record_is_avoided(record):
        return int(record[2] % avoid_modulo == 0 or
                   record[3] % avoid_modulo == 0)

    def block_is_avoided(first, second):
        return record_is_avoided(first) + record_is_avoided(second)

    candidates_by_direction = defaultdict(list)
    for record in candidates:
        candidates_by_direction[record[:2]].append(record)
    for records in candidates_by_direction.values():
        records.sort(key=lambda record: record[2])

    chosen = []
    used = set()
    old_by_pair = defaultdict(list)
    candidate_by_pair = defaultdict(list)
    for record in old_records:
        old_by_pair[tuple(sorted(record[:2]))].append(record)
    for record in candidates:
        candidate_by_pair[tuple(sorted(record[:2]))].append(record)

    unmatched_old = []
    for pair, records in sorted(old_by_pair.items()):
        if len(records) == 2 and records[0][:2] == records[1][:2][::-1]:
            blocks = []
            for first in candidate_by_pair.get(pair, ()):
                for second in candidate_by_pair.get(pair, ()):
                    if first[:2] == second[:2][::-1] and \
                            min(abs(first[2] - second[2]),
                                nominal_sizes.get(
                                    first[0], part.maxBaseIdx() + 1) -
                                abs(first[2] - second[2])) == 1:
                        match = []
                        for old in records:
                            replacement = (first if first[:2] == old[:2]
                                           else second)
                            size = nominal_sizes.get(
                                old[0], part.maxBaseIdx() + 1)
                            delta = abs(replacement[2] - old[2])
                            match.append(min(delta, size - delta))
                        blocks.append((block_is_avoided(first, second),
                                       sum(match), max(match),
                                       min(first[2], second[2]),
                                       first, second))
            if blocks:
                unused_avoided, unused_total, unused_peak, unused_index, \
                    first, second = min(blocks)
                for old in records:
                    replacement = (first if first[:2] == old[:2]
                                   else second)
                    chosen.append(replacement)
                    used.add(replacement)
                continue
        unmatched_old.extend(records)

    for old in sorted(unmatched_old,
                      key=lambda record: (record[2], record)):
        options = [record for record in
                   candidates_by_direction.get(old[:2], ())
                   if record not in used]
        if not options:
            raise RuntimeError(
                "No legal %s scaffold crossover can replace %d[%d]→%d[%d]." %
                ("Honeycomb" if part._step == 21 else "Square",
                 old[0], old[2], old[1], old[3]))
        size = nominal_sizes.get(old[0], part.maxBaseIdx() + 1)
        replacement = min(
            options,
            key=lambda record: (
                record_is_avoided(record),
                min(abs(record[2] - old[2]),
                    size - abs(record[2] - old[2])), record[2]))
        chosen.append(replacement)
        used.add(replacement)

    def create_record(record, undo=False):
        first, second, first_index, second_index = record
        first_vh = part.virtualHelix(first)
        second_vh = part.virtualHelix(second)
        strand5p = (first_vh.scaffoldStrandSet().getStrand(first_index)
                    if first_vh is not None else None)
        strand3p = (second_vh.scaffoldStrandSet().getStrand(second_index)
                    if second_vh is not None else None)
        if not part._canCreateScaffoldXover(
                strand5p, strand3p, first_index):
            return False
        part.createXover(
            strand5p, first_index, strand3p, second_index,
            useUndoStack=undo)
        return True

    # Legacy oligo bookkeeping can print a harmless duplicate-removal trace
    # while reciprocal crossovers are created.  Keep it out of the GUI and
    # preserve the explicit topology checks below as the source of truth.
    with contextlib.redirect_stdout(io.StringIO()), \
            contextlib.redirect_stderr(io.StringIO()):
        for record in chosen:
            if not create_record(record):
                raise RuntimeError(
                    "Legal scaffold crossover became unavailable at "
                    "%d[%d]." % (record[0], record[2]))

    def scaffold_is_single_open():
        oligos = [oligo for oligo in part.oligos()
                  if oligo.strand5p() is not None and
                  not oligo.isStaple()]
        return len(oligos) == 1 and not oligos[0].isLoop()

    if not scaffold_is_single_open():
        raise RuntimeError(
            "Legal scaffold phase replacement did not preserve one scaffold "
            "with one nick.")

    dense_added = 0
    normalized_part = part
    if density_mode != "minimum" and int(layers) >= 1:
        directions = dict(
            (vh.number(), vh.scaffoldStrandSet().isDrawn5to3())
            for vh in part.getVirtualHelices()
            if vh.number() in nominal_sizes)
        existing = set(
            record for record in
            part._existingScaffoldCrossoverRecords()
            if record[0] != record[1])
        baseline_outgoing = {}
        for helix, size in nominal_sizes.items():
            vh = part.virtualHelix(helix)
            strand_set = vh.scaffoldStrandSet()
            delta = 1 if directions[helix] else -1
            for index in range(size):
                strand = strand_set.getStrand(index)
                if strand is None:
                    raise RuntimeError(
                        "Curved scaffold is missing %d[%d]." %
                        (helix, index))
                if index != strand.idx3Prime():
                    baseline_outgoing[(helix, index)] = (
                        helix, (index + delta) % size)
                else:
                    connection = strand.connection3p()
                    baseline_outgoing[(helix, index)] = (
                        (connection.virtualHelix().number(),
                         connection.idx5Prime())
                        if connection is not None else None)
        # Heal only the longitudinal cuts created by the current cross-helix
        # route.  Same-helix circular boundary joins and the true scaffold
        # nick remain exactly as supplied by DNAxiS.
        for from_helix, to_helix, from_index, to_index in existing:
            from_delta = 1 if directions[from_helix] else -1
            baseline_outgoing[(from_helix, from_index)] = (
                from_helix,
                (from_index + from_delta) % nominal_sizes[from_helix])
            to_delta = 1 if directions[to_helix] else -1
            predecessor = (
                to_helix,
                (to_index - to_delta) % nominal_sizes[to_helix])
            baseline_outgoing[predecessor] = (to_helix, to_index)
        nick_sources = [node for node, target in baseline_outgoing.items()
                        if target is None]
        if len(nick_sources) != 1:
            raise RuntimeError(
                "Curved scaffold does not have exactly one source nick "
                "before dense routing.")
        all_nodes = [(helix, index)
                     for helix, size in nominal_sizes.items()
                     for index in range(size)]

        def outgoing_for_records(records):
            outgoing = dict(baseline_outgoing)
            for unused_from, to_helix, unused_from_index, to_index in records:
                delta = 1 if directions[to_helix] else -1
                predecessor = (
                    to_helix,
                    (to_index - delta) % nominal_sizes[to_helix])
                outgoing[predecessor] = None
            for from_helix, to_helix, from_index, to_index in records:
                outgoing[(from_helix, from_index)] = (
                    to_helix, to_index)
            return outgoing

        def topology_metrics(records):
            outgoing = outgoing_for_records(records)
            incoming = dict((node, 0) for node in all_nodes)
            for target in outgoing.values():
                if target is not None:
                    if target not in incoming:
                        return None
                    incoming[target] += 1
                    if incoming[target] > 1:
                        return None
            seen = set()
            open_components = 0
            cycles = 0
            for start in [node for node, count in incoming.items()
                          if count == 0]:
                if start in seen:
                    continue
                open_components += 1
                current = start
                while current is not None and current not in seen:
                    seen.add(current)
                    current = outgoing[current]
            for start in all_nodes:
                if start in seen:
                    continue
                cycles += 1
                current = start
                while current is not None and current not in seen:
                    seen.add(current)
                    current = outgoing[current]
            return {
                "components": open_components + cycles,
                "open_components": open_components,
                "cycles": cycles,
                "covered": len(seen)}

        route = [int(record["helix"]) for record in
                 curvature_indels.get("rings", [])]
        pool = sorted(set(all_candidates()).union(existing))
        dense_records = set(_filterAutoScaffoldCandidatesForPaths(
            pool, [route], part._step, part._step == 32,
            part.minBaseIdx(), part.maxBaseIdx(),
            dict((vh.number(), vh.coord())
                 for vh in part.getVirtualHelices()),
            densitySpacing=requested_spacing))
        # ``potentialCrossoverList`` suppresses a base that already carries
        # the provisional same-helix circular closure.  For Curved rings,
        # index 0/max is not a physical edge: it is one more member of the
        # selected periodic native register (e.g. Square 159/0).  Complete
        # every selected directed phase analytically before replacing the
        # topology, so exact-period ring lengths do not lose that block.
        selected_directions = defaultdict(set)
        for first, second, index, unused_to_index in dense_records:
            selected_directions[(first, second)].add(
                index % requested_spacing)
        for (first, second), phases in selected_directions.items():
            size = min(nominal_sizes[first], nominal_sizes[second])
            first_set = part.virtualHelix(first).scaffoldStrandSet()
            second_set = part.virtualHelix(second).scaffoldStrandSet()
            for phase in phases:
                for index in range(phase, size, requested_spacing):
                    if first_set.getStrand(index) is not None and \
                            second_set.getStrand(index) is not None:
                        dense_records.add(
                            (first, second, index, index))
        if not dense_records:
            raise RuntimeError(
                "No periodic 1/%d scaffold route was found." %
                requested_spacing)

        # The complete periodic register is the geometric template.  It may
        # consist of several small closed scaffold loops.  Remove the minimum
        # number of reciprocal blocks that merge those components, clustering
        # removals into one late internal seam and preserving both edge
        # registers.  With one pre-existing nick, the final component is one
        # open scaffold rather than a closed loop.
        physical_blocks = []
        used = set()
        for first in sorted(dense_records,
                            key=lambda item: (item[2], item[0], item[1])):
            if first in used:
                continue
            partners = [second for second in dense_records
                        if second not in used and
                        first[0] == second[1] and
                        first[1] == second[0] and
                        min(abs(first[2] - second[2]),
                            nominal_sizes.get(
                                first[0], part.maxBaseIdx() + 1) -
                            abs(first[2] - second[2])) == 1]
            if not partners:
                continue
            second = min(partners,
                         key=lambda item: (item[2], item[0], item[1]))
            block = tuple(sorted((first, second),
                                 key=lambda item: item[2]))
            physical_blocks.append(block)
            used.update(block)

        # Multilayer reference layouts use the maximum periodic register on
        # every consecutive helix pair, except the final (bottom) pair.  A
        # one-layer shell needs a different closure: start from the complete
        # register and remove only reciprocal blocks that strictly merge
        # scaffold components, stopping at one open scaffold.
        if len(route) < 2:
            raise RuntimeError(
                "A curved scaffold needs at least two helices.")
        if int(layers) > 1:
            bottom_pair = tuple(sorted((route[-2], route[-1])))
            bottom_blocks = [
                block for block in physical_blocks
                if tuple(sorted(block[0][:2])) == bottom_pair]
            if not bottom_blocks:
                raise RuntimeError(
                    "No legal scaffold seam exists on the bottom helix pair.")
            seam_block = min(
                bottom_blocks,
                key=lambda block: (
                    sum(record_is_avoided(record) for record in block),
                    min(record[2] for record in block)))
            retained_blocks = [
                block for block in physical_blocks
                if tuple(sorted(block[0][:2])) != bottom_pair]
            retained_blocks.append(seam_block)
            target_records = set(
                record for block in retained_blocks for record in block)
            metrics = topology_metrics(target_records)
        else:
            retained_blocks = list(physical_blocks)
            target_records = set(
                record for block in retained_blocks for record in block)
            metrics = topology_metrics(target_records)
            route_rank = dict((helix, index)
                              for index, helix in enumerate(route))
            while metrics is not None and not (
                    metrics["open_components"] == 1 and
                    metrics["cycles"] == 0):
                options = []
                # Removing different periodic blocks on the same helix pair
                # has the same component-merging role.  Test only the latest
                # remaining block per pair; this preserves the preferred late
                # seam while avoiding an O(blocks^2 * bases) global search.
                latest_by_pair = {}
                for block in retained_blocks:
                    pair = tuple(sorted(block[0][:2]))
                    previous = latest_by_pair.get(pair)
                    if previous is None or \
                            max(record[2] for record in block) > \
                            max(record[2] for record in previous):
                        latest_by_pair[pair] = block
                for block in latest_by_pair.values():
                    trial_records = target_records.difference(block)
                    trial_metrics = topology_metrics(trial_records)
                    if trial_metrics is None or \
                            trial_metrics["covered"] != len(all_nodes) or \
                            trial_metrics["components"] >= \
                            metrics["components"]:
                        continue
                    pair_rank = max(
                        route_rank[block[0][0]],
                        route_rank[block[0][1]])
                    options.append((
                        trial_metrics["components"],
                        trial_metrics["cycles"],
                        -pair_rank,
                        -max(record[2] for record in block),
                        block, trial_records, trial_metrics))
                if not options:
                    break
                unused_components, unused_cycles, unused_pair_rank, \
                    unused_index, removed_block, target_records, metrics = \
                    min(options, key=lambda option: option[:4])
                retained_blocks.remove(removed_block)
        if metrics is None or metrics["covered"] != len(all_nodes):
            raise RuntimeError(
                "The requested 1/%d scaffold template is topologically "
                "invalid." % requested_spacing)

        if metrics["open_components"] != 1 or metrics["cycles"] != 0:
            raise RuntimeError(
                "The bottom scaffold seam produced %d components instead "
                "of one scaffold with one nick." % metrics["components"])
        # Rebuild the legacy scaffold arrays directly from the verified graph.
        # The old batch crossover commands can leave stale longitudinal
        # boundaries while their record list still looks correct.  A complete
        # graph rewrite makes the saved JSON and the topology verifier use the
        # exact same connections.
        from cadnano2.model.document import Document
        from cadnano2.model.enum import LatticeType
        from cadnano2.model.io.legacydecoder import import_legacy_dict
        from cadnano2.model.io.legacyencoder import legacy_dict_from_part
        routed = legacy_dict_from_part(
            part, "curved-dense-routing", includeSequences=False)
        outgoing = outgoing_for_records(target_records)
        incoming = {}
        for source, target in outgoing.items():
            if target is None:
                continue
            if target in incoming:
                raise RuntimeError(
                    "Dense scaffold route has a branched target base.")
            incoming[target] = source
        rows = dict((int(row["num"]), row)
                    for row in routed.get("vstrands", []))
        for node in all_nodes:
            helix, index = node
            previous = incoming.get(node)
            following = outgoing.get(node)
            rows[helix]["scaf"][index] = [
                previous[0] if previous is not None else -1,
                previous[1] if previous is not None else -1,
                following[0] if following is not None else -1,
                following[1] if following is not None else -1]
        normalized_document = Document()
        normalized_part = import_legacy_dict(
            normalized_document, routed,
            LatticeType.Honeycomb if part._step == 21
            else LatticeType.Square,
            forceLatticeType=True)
        if normalized_part is None:
            raise RuntimeError(
                "Could not reload the verified dense scaffold graph.")
        normalized_scaffolds = [
            oligo for oligo in normalized_part.oligos()
            if oligo.strand5p() is not None and not oligo.isStaple()]
        if len(normalized_scaffolds) != 1 or \
                normalized_scaffolds[0].isLoop():
            raise RuntimeError(
                "Dense scaffold graph reload found %d components." %
                len(normalized_scaffolds))
        dense_added = max(0, len(target_records) - len(chosen))

    final_records = [record for record in
                     normalized_part._existingScaffoldCrossoverRecords()
                     if record[0] != record[1]]
    if any(first_index != second_index
           for unused_first, unused_second, first_index, second_index
           in final_records):
        raise RuntimeError("A curved scaffold crossover is still off-lattice.")
    spacing_by_direction = defaultdict(list)
    for record in final_records:
        spacing_by_direction[record[:2]].append(int(record[2]))
    actual_spacings = []
    for (first, unused_second), indices in spacing_by_direction.items():
        size = int(nominal_sizes.get(first, part.maxBaseIdx() + 1))
        ordered = sorted(set(index % size for index in indices))
        if len(ordered) == 1:
            actual_spacings.append(size)
        elif ordered:
            actual_spacings.extend(
                right - left for left, right in
                zip(ordered, ordered[1:] + [ordered[0] + size]))
    density_label = (
        "minimum" if density_mode == "minimum" else
        "maximum" if density_multiple == 1 else
        "periodic-%dx" % density_multiple)
    return ({
        "source": "cadnano native lattice phases",
        "density_mode": density_label,
        "density_multiple": (0 if density_mode == "minimum"
                             else density_multiple),
        "requested_spacing_bp": (0 if requested_spacing is None
                                 else requested_spacing),
        "native_spacing_minimum_bp": (
            min(actual_spacings) if actual_spacings else 0),
        "native_spacing_maximum_bp": (
            max(actual_spacings) if actual_spacings else 0),
        "replaced_crossovers": len(chosen),
        "added_dense_crossovers": dense_added,
        "final_crossovers": len(final_records),
        "scaffold_components": 1,
        "scaffold_nicks": 1}, normalized_part)


def _include_scaffold_pairs_in_staple_allowlist(allowed_pairs,
                                                 scaffold_records):
    """Allow AutoCS staple completion on every routed scaffold helix pair.

    Reinforced Curved designs start with a restricted set of inter-layer
    staple pairs.  Sparse scaffold routing (for example Honeycomb 1/42 or
    Square 1/64) must also expose each actual scaffold crossover pair to
    AutoCS so the legal native sites between scaffold crossovers can restore
    the normal total crossover density.  Pair direction is intentionally
    ignored, matching the scaffold/staple distance rule in ``part.py``.
    """
    if allowed_pairs is None:
        return None
    completed = set(tuple(sorted(pair)) for pair in allowed_pairs)
    for from_helix, to_helix, unused_from_index, unused_to_index in \
            scaffold_records:
        if from_helix != to_helix:
            completed.add(tuple(sorted((from_helix, to_helix))))
    return completed


def _apply_native_staple_rules(json_path, lattice):
    """Replace provisional DNAxiS staples with native lattice staples."""
    from cadnano2.model.document import Document
    from cadnano2.model.enum import LatticeType
    from cadnano2.model.io.legacydecoder import import_legacy_dict
    from cadnano2.model.io.legacyencoder import legacy_dict_from_part

    with open(json_path, "r", encoding="utf-8") as source:
        original = json.load(source)
    metadata = dict(original.get("curved_metadata") or {})
    lattice_type = (LatticeType.Honeycomb if lattice == "honeycomb"
                    else LatticeType.Square)
    document = Document()
    part = import_legacy_dict(
        document, original, lattice_type, forceLatticeType=True)
    if part is None:
        raise RuntimeError("Cannot decode the lattice-indel curved design.")
    curvature_indels = dict(original.get("curvature_indels") or {})
    if not curvature_indels.get("rings"):
        # Historical DNAxiS/Curved intermediates encoded the actual indels in
        # cadnano's loop/skip arrays but did not retain the corresponding
        # per-helix placement records.  The final native-AutoCS rebalance then
        # saw an empty plan and silently left the old distribution untouched.
        # Reconstruct the records losslessly from loop/skip before topology is
        # regenerated.  Counts, signs and nominal helix lengths are preserved
        # exactly; only their later safe coordinates may be redistributed.
        reconstructed = []
        for row in sorted(original.get("vstrands", []),
                          key=lambda item: int(item.get("num", 0))):
            loops = list(row.get("loop", []))
            skips = list(row.get("skip", []))
            insertions = []
            deletions = []
            for index, value in enumerate(loops):
                insertions.extend([index] * max(0, int(value)))
            for index, value in enumerate(skips):
                deletions.extend([index] * max(0, -int(value)))
            reconstructed.append({
                "helix": int(row["num"]),
                "nominal_bases": min(len(loops), len(skips)),
                "insertions": insertions,
                "deletions": deletions})
        curvature_indels["rings"] = reconstructed
        curvature_indels["record_source"] = (
            "lossless reconstruction from cadnano loop/skip arrays")

    def crossover_fingerprint(encoded):
        result = set()
        for row in encoded.get("vstrands", []):
            helix = int(row["num"])
            for strand_type in ("scaf", "stap"):
                for index, entry in enumerate(row.get(strand_type, [])):
                    for side, (helix_position, index_position) in enumerate(
                            ((0, 1), (2, 3))):
                        neighbour = int(entry[helix_position])
                        if neighbour >= 0 and neighbour != helix:
                            result.add((strand_type, helix, index, side,
                                        neighbour,
                                        int(entry[index_position])))
        return result
    scaffold_summary, part = _normalize_curved_scaffold_crossovers(
        part, curvature_indels, layers=int(metadata.get("layers", 1)),
        density_mode=metadata.get(
            "requested_scaffold_crossover_density_mode", "periodic"),
        density_multiple=int(metadata.get(
            "requested_scaffold_crossover_density_multiple", 1)))
    layers = int(metadata.get("layers", 1))
    allowed_staple_pairs = None
    if layers > 1:
        route = [int(record["helix"]) for record in
                 curvature_indels.get("rings", [])]
        if not route:
            route = sorted(vh.number() for vh in part.getVirtualHelices())
        rings_per_layer = int(metadata.get(
            "rings_per_layer", len(route) // layers))
        if rings_per_layer * layers == len(route):
            allowed_staple_pairs = set()

            def ring_at(layer, slice_index):
                offset = (slice_index if layer % 2 == 0 else
                          rings_per_layer - 1 - slice_index)
                return route[layer * rings_per_layer + offset]

            for layer in range(layers - 1):
                for slice_index in range(rings_per_layer):
                    allowed_staple_pairs.add(tuple(sorted((
                        ring_at(layer, slice_index),
                        ring_at(layer + 1, slice_index)))))
                # The serpentine scaffold already turns between these two
                # rings; the Square reference intentionally omits the staple
                # crossover family on that same turn pair.
                allowed_staple_pairs.discard(tuple(sorted((
                    route[(layer + 1) * rings_per_layer - 1],
                    route[(layer + 1) * rings_per_layer]))))
            # Retain the final reference seam family.
            allowed_staple_pairs.add(tuple(sorted((route[-2], route[-1]))))
            # The inter-layer allowlist above is only an extra reinforcement
            # set.  It must not hide the routed scaffold pairs themselves:
            # those pairs need their legal midpoint staple sites when the
            # selected scaffold density is sparser than the native AutoCS
            # staple density.
            allowed_staple_pairs = \
                _include_scaffold_pairs_in_staple_allowlist(
                    allowed_staple_pairs,
                    part._existingScaffoldCrossoverRecords())
    if not part.autoStaple(
            preservePeriodicCrossovers=True,
            allowedCrossoverPairs=allowed_staple_pairs):
        raise RuntimeError("Native AutoCS_staples could not create staples.")

    # AutoCS defines the final crossover topology.  Curved Design must never
    # add, delete, half-open or roll back any crossover after this point.
    # Only indels and subsequent staple nicks may move.
    preliminary = legacy_dict_from_part(
        part, os.path.basename(json_path), includeSequences=False)
    autocs_crossover_fingerprint = crossover_fingerprint(preliminary)
    closed_staple_boundaries = _close_circular_staple_boundaries(
        preliminary, curvature_indels)
    provisional = preliminary
    _rebalance_indels_against_final_crossovers(
        provisional, lattice, curvature_indels)

    # Autobreak must see the final indel positions because insertions and
    # deletions contribute to actual staple length and continuous-run tests.
    rebalanced_document = Document()
    rebalanced_part = import_legacy_dict(
        rebalanced_document, provisional, lattice_type,
        forceLatticeType=True)
    if rebalanced_part is None:
        raise RuntimeError(
            "Cannot reload the crossover-balanced curved design.")
    reloaded_scaffolds = [
        oligo for oligo in rebalanced_part.oligos()
        if oligo.strand5p() is not None and not oligo.isStaple()]
    if len(reloaded_scaffolds) != 1 or reloaded_scaffolds[0].isLoop():
        raise RuntimeError(
            "Final scaffold verification found %d components; Curved Design "
            "requires one scaffold with one nick." %
            len(reloaded_scaffolds))
    scaffold_summary["scaffold_components"] = 1
    scaffold_summary["scaffold_nicks"] = 1
    scaffold_summary["final_crossovers"] = len([
        record for record in
        rebalanced_part._existingScaffoldCrossoverRecords()
        if record[0] != record[1]])
    healed_short_staples = _heal_short_curved_staple_nicks(
        rebalanced_part)
    preconditioned_staple_nicks = \
        _precondition_indel_blocked_staple_nicks(
            rebalanced_part, lattice)
    pre_autobreak_indel_adjustment = {
        "relocated_conflicts": 0, "fused_insertions": 0}
    if preconditioned_staple_nicks:
        pre_autobreak_encoded = legacy_dict_from_part(
            rebalanced_part, os.path.basename(json_path),
            includeSequences=False)
        pre_autobreak_indel_adjustment = \
            _rebalance_indels_against_final_crossovers(
                pre_autobreak_encoded, lattice, curvature_indels,
                only_conflicts=True)
        pre_autobreak_document = Document()
        rebalanced_part = import_legacy_dict(
            pre_autobreak_document, pre_autobreak_encoded, lattice_type,
            forceLatticeType=True)
        if rebalanced_part is None:
            raise RuntimeError(
                "Cannot reload the indel-cleared Autobreak nick plan.")

    break_result = rebalanced_part.autoBreakStaples(
        preserveCrossovers=True, markUnbreakable=True,
        preferDeletionDense=True)
    encoded = legacy_dict_from_part(
        rebalanced_part, os.path.basename(json_path),
        includeSequences=False)
    # Autobreak has now fixed the final nick set.  Perform a second, local
    # conflict-only pass: indels on either nick endpoint are moved (or, for
    # insertions, fused into an existing site in the same crossover-free
    # interval) without globally redistributing the already balanced plan.
    post_nick_indel_adjustment = _rebalance_indels_against_final_crossovers(
        encoded, lattice, curvature_indels, only_conflicts=True)
    # The first pass above gives every individual helix the shared
    # equal-partition distribution.  Refine only the coordinates of those
    # already-budgeted indels against the final immutable AutoCS topology so
    # every curvature-bearing adjacent-helix pair bends in the requested
    # direction, remains inside its theoretical floor/ceiling interval, and
    # has the smallest practical axial curvature fluctuation.  Signed indel
    # totals, crossovers, nicks and strand topology are invariants here.
    # legacy_dict_from_part intentionally contains only caDNAno topology; add
    # the original geometry payload back before physical-neighbour analysis.
    # Without it the optimizer sees zero physical pairs and silently skips.
    encoded["curved_metadata"] = dict(metadata)
    pair_curvature_summary = optimize_curved_pair_curvature(
        encoded, curvature_indels, maximum_passes=3)
    pair_curvature_audit = dict(encoded.get("curved_metadata", {}).get(
        "pair_aware_indel_optimization") or {})
    final_indel_summary = _final_crossover_indel_statistics(
        encoded, lattice, curvature_indels)
    final_crossover_fingerprint = crossover_fingerprint(encoded)
    if final_crossover_fingerprint != autocs_crossover_fingerprint:
        added = final_crossover_fingerprint - autocs_crossover_fingerprint
        removed = autocs_crossover_fingerprint - final_crossover_fingerprint
        raise RuntimeError(
            "Curved Design changed the immutable AutoCS crossover set "
            "(added directed endpoints: %d; removed directed endpoints: %d)."
            % (len(added), len(removed)))
    if not final_indel_summary.get("domain_limit_feasible", False):
        raise RuntimeError(
            "Final Curved Design violates the +/-3 indel/domain hard limit.")
    metadata.update({
        "native_staple_rules_pending": False,
        "scaffold_crossover_source": scaffold_summary["source"],
        "scaffold_crossover_density_mode":
            scaffold_summary["density_mode"],
        "scaffold_crossover_density_multiple":
            scaffold_summary["density_multiple"],
        "scaffold_crossover_requested_spacing_bp":
            scaffold_summary["requested_spacing_bp"],
        "scaffold_crossover_native_spacing_minimum_bp":
            scaffold_summary["native_spacing_minimum_bp"],
        "scaffold_crossover_native_spacing_maximum_bp":
            scaffold_summary["native_spacing_maximum_bp"],
        "scaffold_crossover_count":
            scaffold_summary["final_crossovers"],
        "scaffold_crossover_replacements":
            scaffold_summary["replaced_crossovers"],
        "scaffold_crossover_dense_additions":
            scaffold_summary["added_dense_crossovers"],
        "scaffold_components": scaffold_summary["scaffold_components"],
        "scaffold_nicks": scaffold_summary["scaffold_nicks"],
        "staple_crossover_source": "cadnano AutoCS_staples",
        "staple_crossover_pairs": (
            [list(pair) for pair in sorted(allowed_staple_pairs)]
            if allowed_staple_pairs is not None else "all adjacent pairs"),
        "post_autocs_crossover_policy": "immutable",
        "post_autocs_crossovers_added": 0,
        "post_autocs_crossovers_removed": 0,
        "short_staple_nicks_healed_before_autobreak": [
            list(record) for record in healed_short_staples],
        "indel_blocked_nicks_preconditioned": [
            list(record) for record in preconditioned_staple_nicks],
        "pre_autobreak_indel_adjustment":
            pre_autobreak_indel_adjustment,
        "post_autobreak_indel_adjustment": post_nick_indel_adjustment,
        "pair_curvature_summary": pair_curvature_summary,
        "pair_aware_indel_optimization": pair_curvature_audit,
        "autobreak_lattice": lattice,
        "autobreak_continuous_minimum": (
            14 if lattice == "honeycomb" else 16),
        "autobreak_created_nicks": int(break_result.get("nicks", 0)),
        "autobreak_tolerated_59_64nt_staples": int(
            break_result.get("tolerated_long_staples", 0)),
        "closed_circular_staple_boundaries":
            int(closed_staple_boundaries),
        "red_staple_warning_count": int(break_result.get("skipped", 0)),
        "red_staple_warning": (
            "Red staples could not be partitioned while preserving every "
            "lattice crossover and require manual review."),
        "final_indel_summary": final_indel_summary})
    metadata.update({key: final_indel_summary[key] for key in (
        "domain_size_bp",
        "maximum_indel_per_domain_allowed",
        "maximum_insertion_per_domain",
        "maximum_deletion_per_domain",
        "maximum_insertion_per_crossover",
        "maximum_deletion_per_crossover",
        "effective_crossover_spacing_minimum",
        "effective_crossover_spacing_maximum")})
    encoded["curved_metadata"] = metadata
    if curvature_indels:
        curvature_indels["final_crossover_statistics"] = final_indel_summary
        encoded["curvature_indels"] = curvature_indels
    with open(json_path, "w", encoding="utf-8", newline="\n") as output:
        json.dump(encoded, output, separators=(",", ":"))
    return break_result, final_indel_summary


def _reoptimize_curved_staple_nicks(json_path, lattice):
    """Locally repartition short and deletion-dense final staples.

    The first Autobreak pass sees the complete oligo and can prefer a good
    *global* partition while still leaving an individual final product with
    two or more deletions at only 21--39 actual nucleotides.  Such a product
    has reduced uninterrupted duplex support and must be actively repaired,
    not merely receive a lower score.  Likewise, a legal ordinary 21--29 nt
    product must be reconsidered so that the local component can preferentially
    repartition into 30--50 nt products.  This pass therefore heals only the
    adjacent same-helix nick of a target product and reruns the native Autobreak
    planner on that local component.  AutoCS crossovers, scaffold topology and
    all indels are immutable.
    """
    from cadnano2.model.document import Document
    from cadnano2.model.enum import LatticeType
    from cadnano2.model.io.legacydecoder import import_legacy_dict
    from cadnano2.model.io.legacyencoder import legacy_dict_from_part
    from cadnano2.model.parts.part import (
        _existingStapleNickBoundaries, _stapleOligoBaseRecords)

    with open(json_path, "r", encoding="utf-8") as source:
        original = json.load(source)

    def crossover_fingerprint(design):
        result = set()
        for row in design.get("vstrands", ()):
            helix = int(row["num"])
            for strand_type in ("scaf", "stap"):
                for index, entry in enumerate(row.get(strand_type, ())):
                    for offset in (0, 2):
                        other, other_index = map(
                            int, entry[offset:offset + 2])
                        if other >= 0 and other != helix:
                            result.add((strand_type, helix, index, offset,
                                        other, other_index))
        return frozenset(result)

    before_crossovers = crossover_fingerprint(original)
    lattice_type = (LatticeType.Honeycomb if lattice == "honeycomb"
                    else LatticeType.Square)
    document = Document()
    part = import_legacy_dict(
        document, original, lattice_type, forceLatticeType=True)
    if part is None:
        raise RuntimeError(
            "Cannot reload Curved design for final staple optimization.")

    def staple_state():
        state = []
        for oligo in part.oligos():
            if not oligo.isStaple() or oligo.isHybrid() or \
                    oligo.strand5p() is None:
                continue
            records = _stapleOligoBaseRecords(oligo)
            state.append((
                oligo, int(oligo.actualLength()),
                sum(int(record[3]) == 0 for record in records)))
        return state

    initial_staple_state = staple_state()
    initial_normal_below_preferred = sum(
        deletion_count < 2 and 21 <= length < 30
        for unused_oligo, length, deletion_count in initial_staple_state)
    used_heal_boundaries = set()

    def heal_target_nicks():
        """Join one neighbour of each locally under-supported product.

        Ordinary staples use 30--50 nt as the preferred target.  Staples
        weakened by at least two deletions use 40--60 nt.  Joining is limited
        to an existing same-helix nick, after which native Autobreak performs
        the actual local repartition.
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
                boundary = (helix, upper_index)
                if boundary in used_heal_boundaries:
                    continue
                strand_set = part.virtualHelix(helix).stapleStrandSet()
                lower = strand_set.getStrand(upper_index - 1)
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
                target_strand = lower if lower_oligo in targets else upper
                other = upper if target_strand is lower else lower
                target_oligo = target_strand.oligo()
                unused_length, deletion_count = target_details[target_oligo]
                if deletion_count >= 2:
                    if 40 <= combined <= 60:
                        repair_class, deviation = 0, abs(combined - 50)
                    elif 80 <= combined <= 120:
                        repair_class, deviation = 1, abs(combined - 100)
                    elif 61 <= combined < 80:
                        repair_class, deviation = 2, abs(combined - 60)
                    elif combined < 40:
                        repair_class, deviation = 3, 40 - combined
                    else:
                        repair_class, deviation = 4, combined - 120
                else:
                    if 30 <= combined <= 50:
                        repair_class, deviation = 0, abs(combined - 40)
                    elif 60 <= combined <= 100:
                        repair_class, deviation = 1, abs(combined - 80)
                    elif 51 <= combined <= 57:
                        repair_class, deviation = 2, abs(combined - 50)
                    elif combined < 30:
                        repair_class, deviation = 3, 30 - combined
                    else:
                        repair_class, deviation = 4, abs(combined - 80)
                score = (repair_class, deviation, combined,
                         helix, upper_index)
                candidates.append((score, strand_set, target_strand, other,
                                   boundary))
            if not candidates:
                break
            unused_score, strand_set, target_strand, other, boundary = min(
                candidates, key=lambda item: item[0])
            strand_set.mergeStrands(
                target_strand, other, useUndoStack=False)
            used_heal_boundaries.add(boundary)
            healed += 1
        return healed

    total_created = total_healed = total_skipped = total_tolerated = 0
    pass_count = 0
    for pass_count in range(1, 21):
        part._autobreakStaplesApplied = False
        result = part.autoBreakStaples(
            preserveCrossovers=True, markUnbreakable=True,
            preferDeletionDense=True)
        created = int(result.get("nicks", 0))
        healed = heal_target_nicks()
        total_created += created
        total_healed += healed
        total_skipped += int(result.get("skipped", 0))
        total_tolerated += int(result.get("tolerated_long_staples", 0))
        if not created and not healed:
            break

    final_staple_state = staple_state()
    lengths = []
    dense_lengths = []
    exceptional_final = 0
    final_normal_below_preferred = 0
    for unused_oligo, length, deletion_count in final_staple_state:
        lengths.append(length)
        if deletion_count >= 2:
            dense_lengths.append(length)
            exceptional_final += int(length > 60)
        else:
            final_normal_below_preferred += int(21 <= length < 30)
            exceptional_final += int(length > 57)

    encoded = legacy_dict_from_part(
        part, os.path.basename(json_path), includeSequences=False)
    for key, value in original.items():
        if key not in ("name", "num_bases", "vstrands"):
            encoded[key] = value
    after_crossovers = crossover_fingerprint(encoded)
    if after_crossovers != before_crossovers:
        raise RuntimeError(
            "Curved staple optimization changed immutable AutoCS "
            "crossovers.")

    outside = [value for value in dense_lengths
               if not 40 <= value <= 60]
    audit = {
        "method": "iterative local short/deletion-dense nick repartition",
        "definition": (
            "ordinary 21--29 nt products and deletion-dense products "
            "below 40 nt"),
        "passes": int(pass_count),
        "healed_target_nicks": int(total_healed),
        "created_nicks": int(total_created),
        "optimizer_skipped_attempts": int(total_skipped),
        "exceptional_fallback_attempts": int(total_tolerated),
        "exceptional_58_64nt_components": int(exceptional_final),
        "minimum_staple_nt": min(lengths or [0]),
        "maximum_staple_nt": max(lengths or [0]),
        "normal_staple_range_nt": [21, 57],
        "normal_preferred_range_nt": [30, 50],
        "normal_below_preferred_before_repartition": int(
            initial_normal_below_preferred),
        "normal_below_preferred_after_repartition": int(
            final_normal_below_preferred),
        "exceptional_no_solution_range_nt": [58, 64],
        "deletion_dense_target_range_nt": [40, 60],
        "deletion_dense_staples": len(dense_lengths),
        "deletion_dense_in_40_60": sum(
            40 <= value <= 60 for value in dense_lengths),
        "deletion_dense_outside_40_60_nt": outside,
        "deletion_dense_lengths_nt": dense_lengths,
        "unresolved_deletion_dense_staples": len(outside),
        "preserved": ["all scaffold links", "all AutoCS crossovers",
                      "all indels"]}
    metadata = encoded.setdefault("curved_metadata", {})
    metadata["curved_staple_nick_optimization"] = audit
    metadata["autobreak_created_nicks"] = (
        int(metadata.get("autobreak_created_nicks", 0)) + total_created)
    metadata["autobreak_tolerated_59_64nt_staples"] = int(
        exceptional_final)
    metadata["autobreak_tolerated_58_64nt_staples"] = int(
        exceptional_final)
    encoded.setdefault("curvature_indels", {})[
        "curved_staple_nick_optimization"] = audit
    with open(json_path, "w", encoding="utf-8", newline="\n") as output:
        json.dump(encoded, output, separators=(",", ":"))
    return audit


def create_curved_project(spec, progress=None, cancelled=None):
    shape = str(spec["shape"]).lower()
    lattice = str(spec.get("lattice", "square")).lower()
    if lattice not in ("square", "honeycomb"):
        raise ValueError("Curved lattice must be square or honeycomb.")
    scaffold_density_mode = str(spec.get(
        "scaffold_crossover_density_mode", "periodic")).lower()
    if scaffold_density_mode not in ("periodic", "minimum"):
        raise ValueError("Unknown scaffold crossover density mode.")
    scaffold_density_multiple = max(1, int(spec.get(
        "scaffold_crossover_density_multiple", 1) or 1))
    scaffold_density_spacing = (
        (21 if lattice == "honeycomb" else 32) *
        scaffold_density_multiple
        if scaffold_density_mode == "periodic" else 0)
    requested_output_name = spec.get("output_name_override") or \
        curved_output_name(
            spec.get("name") or shape, shape, lattice,
            spec.get("layers", 1), spec["height_nm"],
            spec["maximum_diameter_nm"], spec["minimum_diameter_nm"])
    rings = build_rings(
        shape, spec["height_nm"], spec["maximum_diameter_nm"],
        spec["minimum_diameter_nm"], spec.get("layers", 1),
        lattice=lattice)
    rings, planned_indels = curved_indel_plan(rings, lattice)
    required_scaffold_bases = estimated_scaffold_bases(rings)
    if required_scaffold_bases > CURVED_SCAFFOLD_MAX_BASES:
        raise ValueError(
            "Curved Design requires %d scaffold bases; the supported maximum "
            "is %d bases." %
            (required_scaffold_bases, CURVED_SCAFFOLD_MAX_BASES))
    if not planned_indels.get("domain_limit_feasible", False):
        raise ValueError(
            "The requested geometry requires more than +/-3 bases in at "
            "least one native %d-bp domain. Reduce curvature, increase the "
            "axial span, or change the cross-section." %
            int(planned_indels["domain_size_bp"]))
    project_root, output_name, version_suffix, version_number = \
        unique_curved_project_target(
            spec["project_root"], requested_output_name)
    input_dir = os.path.join(project_root, "input")
    os.makedirs(input_dir, exist_ok=True)
    parameters = {
        "format": "cadnano-curved-project-v1",
        "name": safe_name(spec.get("name") or shape),
        "output_name": output_name,
        "requested_output_name": requested_output_name,
        "automatic_version_suffix": version_suffix,
        "automatic_version_number": version_number,
        "shape": shape,
        "lattice": lattice,
        "height_nm": float(spec["height_nm"]),
        "maximum_diameter_nm": float(spec["maximum_diameter_nm"]),
        "minimum_diameter_nm": float(spec["minimum_diameter_nm"]),
        "layers": int(spec.get("layers", 1)),
        "requested_scaffold_crossover_density_mode":
            scaffold_density_mode,
        "requested_scaffold_crossover_density_multiple": (
            scaffold_density_multiple
            if scaffold_density_mode == "periodic" else 0),
        "requested_scaffold_crossover_spacing_bp":
            scaffold_density_spacing,
        "ring_spacing_nm": RING_SPACING_NM,
        "dna_helix_radius_nm": DNA_HELIX_RADIUS_NM,
        "dimensions_reference": "outer-envelope",
        "profile_sampling": "meridional-arc-length",
        "layer_mode": ("honeycomb-brick-wall" if lattice == "honeycomb"
                       else "radial-offset"),
        "ring_count": len(rings),
        "rings_per_layer": len(rings) // int(spec.get("layers", 1)),
        "actual_meridian_spacing_nm": float(
            rings[0].get("meridian_spacing_nm", RING_SPACING_NM)),
        "actual_outer_height_nm": (
            max(float(ring["height_nm"]) for ring in rings) -
            min(float(ring["height_nm"]) for ring in rings) +
            2.0 * DNA_HELIX_RADIUS_NM),
        "actual_minimum_outer_diameter_nm": min(
            2.0 * (float(ring["radius_nm"]) + DNA_HELIX_RADIUS_NM)
            for ring in rings),
        "actual_maximum_outer_diameter_nm": max(
            2.0 * (float(ring["radius_nm"]) + DNA_HELIX_RADIUS_NM)
            for ring in rings),
        "minimum_ring_bp": min(int(ring["bp"]) for ring in rings),
        "maximum_ring_bp": max(int(ring["bp"]) for ring in rings),
        "large_diameter_warning": max(
            2.0 * (float(ring["radius_nm"]) + DNA_HELIX_RADIUS_NM)
            for ring in rings) > 60.0,
        "planned_indel_summary": planned_indels,
        "designed_scaffold_bases": required_scaffold_bases,
        "unrelaxed": True}
    stl_path = os.path.join(input_dir, output_name + "_shape.stl")
    settings_path = os.path.join(
        input_dir, output_name + "_design_settings.json")
    modules_path = os.path.join(
        input_dir, output_name + "_modules.csv")
    preview_path = os.path.join(
        input_dir, output_name + "_preview.png")
    with open(stl_path, "w", encoding="ascii", newline="\n") as output:
        output.write(_stl_text(rings))
    with open(settings_path, "w", encoding="utf-8", newline="\n") as output:
        json.dump(parameters, output, indent=2, sort_keys=True)
        output.write("\n")
    with open(modules_path, "w", encoding="utf-8", newline="\n") as output:
        output.write(
            "index,layer,slice,height_nm,radius_nm,requested_radius_nm,"
            "meridian_spacing_nm,bp,direction\n")
        for ring in rings:
            output.write(
                "{index},{layer},{slice},{height_nm:.6f},{radius_nm:.6f},"
                "{requested_radius_nm:.6f},{meridian_spacing_nm:.6f},"
                "{bp},{direction}\n".format(**ring))
    with open(preview_path, "wb") as output:
        output.write(_preview_png(rings))

    metadata = dict(parameters)
    metadata.update({
        "project_root": project_root,
        "input_stl": os.path.relpath(stl_path, project_root),
        "design_settings": os.path.relpath(settings_path, project_root),
        "modules_csv": os.path.relpath(modules_path, project_root),
        "preview_png": os.path.relpath(preview_path, project_root)})
    with tempfile.TemporaryDirectory(prefix="cadnano-dnaxis-") as temp_dir:
        runner_spec = {
            "output_dir": temp_dir,
            "project_root": project_root,
            "output_name": output_name,
            "lattice": lattice,
            "rings": rings,
            "metadata": metadata}
        spec_path = os.path.join(temp_dir, "runner_spec.json")
        result_path = os.path.join(temp_dir, "runner_result.json")
        with open(spec_path, "w", encoding="utf-8", newline="\n") as output:
            json.dump(runner_spec, output, separators=(",", ":"))
        _run([_dnaxis_python(), os.path.join(
            curved_root(), "headless_runner.py"), spec_path, result_path],
            progress=progress, cancelled=cancelled)
        with open(result_path, "r", encoding="utf-8") as source:
            result = json.load(source)
    break_result, final_indel_summary = _apply_native_staple_rules(
        result["json_path"], lattice)
    staple_audit = _reoptimize_curved_staple_nicks(
        result["json_path"], lattice)
    with open(result["json_path"], "r", encoding="utf-8") as source:
        final_metadata = dict(
            json.load(source).get("curved_metadata") or metadata)
    result["autobreak"] = break_result
    result["staple_nick_optimization"] = staple_audit
    result["indel_summary"] = final_indel_summary
    result.update({
        "project_root": project_root,
        "input_dir": input_dir,
        "metadata": final_metadata,
        "input_paths": [stl_path, settings_path, modules_path, preview_path]})
    return result
