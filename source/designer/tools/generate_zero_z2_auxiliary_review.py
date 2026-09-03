#!/usr/bin/env python3
"""Generate review files for Z2=0/8 using auxiliary SST channels.

The accepted Square/Kagome SST routers intentionally reject physical
occupation of the same (virtual helix, base) by the two layers.  For review of
the requested ultra-small spacings, this script moves every *complete* layer-2
polymer component that would collide with layer 1 to a second set of sixteen
virtual helices.  Components are never cropped or restarted, so the reviewed
absolute crossover phase and 32/48-nt chain topology remain unchanged.

This is deliberately a review-only representation.  It does not change the
production SST generator or claim that the auxiliary virtual helices are a
validated physical routing.
"""

from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = ROOT / "review_outputs" / \
    "Z2_0_8_auxiliary_channel_review_20260814"
SQUARE_SCRIPT = ROOT / "tools" / "generate_sst_absolute_phase_review.py"
KAGOME_SCRIPT = ROOT / "tools" / "generate_kagome_centered_sst_review.py"

CASES = ((64, 0, 64), (88, 8, 88))
EMPTY = [-1, -1, -1, -1]
PRIMARY = tuple(range(16))
AUXILIARY_SHIFT = 64
CAPTURE_COLOURS = {1: 0xFF2D8D, 2: 0xFF8A00}


def _module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SQUARE = _module("square_review", SQUARE_SCRIPT)
KAGOME = _module("kagome_review", KAGOME_SCRIPT)


def _components(field_rows):
    nodes = {
        (int(helix), int(base))
        for helix, records in field_rows.items()
        for base, record in enumerate(records)
        if record != EMPTY
    }
    output = []
    visited = set()
    for start in sorted(nodes):
        if start in visited:
            continue
        component = set()
        stack = [start]
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            helix, base = node
            record = field_rows[helix][base]
            for slot in (0, 2):
                partner = tuple(map(int, record[slot:slot + 2]))
                if partner in nodes and partner not in component:
                    stack.append(partner)
        visited.update(component)
        output.append(component)
    return output


def _occupied(layer):
    return {
        (helix, base)
        for field_rows in layer.values()
        for helix, records in field_rows.items()
        for base, record in enumerate(records)
        if record != EMPTY
    }


def _route_conflicts(layer1, layer2, array_length):
    """Detour only colliding layer-2 nodes through h64--79.

    The accepted visual convention is the one demonstrated manually on the
    first two helices: a polymer remains on its primary helix wherever no
    collision exists, enters the corresponding auxiliary helix for the
    colliding interval, then crosses directly back to the primary helix.  The
    chain therefore remains one connected caDNAno component and its relation
    to the original logical helix is visible without an external lookup.
    """
    occupied1 = _occupied(layer1)
    primary2 = {
        field: {helix: [EMPTY[:] for unused in range(array_length)]
                for helix in PRIMARY}
        for field in ("scaf", "stap")
    }
    auxiliary = {
        field: {
            helix + AUXILIARY_SHIFT: [EMPTY[:] for unused in range(array_length)]
            for helix in PRIMARY
        }
        for field in ("scaf", "stap")
    }
    destination = {}
    detours = []
    # Decide every destination before translating any connection.  A link
    # crossing the conflict boundary can then point from h64+n directly to
    # h+n (or the reverse), exactly like the user's reviewed first pair.
    for field in ("scaf", "stap"):
        for helix, records in layer2[field].items():
            for base, record in enumerate(records):
                if record == EMPTY:
                    continue
                destination[(field, helix, base)] = (
                    helix + AUXILIARY_SHIFT
                    if (helix, base) in occupied1 else helix,
                    base)

    for field in ("scaf", "stap"):
        for component_index, component in enumerate(_components(layer2[field])):
            collisions = sorted(component & occupied1)
            if collisions:
                boundary_links = []
                for helix, base in component:
                    record = layer2[field][helix][base]
                    here_aux = (helix, base) in occupied1
                    for slot in (0, 2):
                        partner = tuple(map(int, record[slot:slot + 2]))
                        if partner in component and \
                                here_aux != (partner in occupied1):
                            edge = tuple(sorted(((helix, base), partner)))
                            if edge not in boundary_links:
                                boundary_links.append(edge)
                detours.append({
                    "field": field,
                    "component_index": component_index,
                    "component_length_nt": len(component),
                    "detoured_node_count": len(collisions),
                    "detoured_nodes": [list(node) for node in collisions],
                    "detoured_helices": sorted({node[0] for node in collisions}),
                    "detoured_base_range": [
                        min(node[1] for node in collisions),
                        max(node[1] for node in collisions)],
                    "primary_auxiliary_boundary_links": [
                        [list(edge[0]), list(edge[1])]
                        for edge in sorted(boundary_links)],
                })

        for helix, records in layer2[field].items():
            for base, source_record in enumerate(records):
                if source_record == EMPTY:
                    continue
                target_helix, target_base = destination[(field, helix, base)]
                translated = copy.deepcopy(source_record)
                for slot in (0, 2):
                    partner, partner_base = map(
                        int, translated[slot:slot + 2])
                    if partner < 0:
                        continue
                    target_partner = destination.get(
                        (field, partner, partner_base))
                    if target_partner is None:
                        raise ValueError(
                            "layer-2 partner is missing from detour map")
                    translated[slot:slot + 2] = list(target_partner)
                if target_helix >= AUXILIARY_SHIFT:
                    auxiliary[field][target_helix][target_base] = translated
                else:
                    primary2[field][target_helix][target_base] = translated
    return primary2, auxiliary, destination, detours


def _merge(first, second):
    output = copy.deepcopy(first)
    for field in ("scaf", "stap"):
        for helix, records in second[field].items():
            for base, record in enumerate(records):
                if record == EMPTY:
                    continue
                if output[field][helix][base] != EMPTY:
                    raise ValueError("unresolved same-field collision")
                output[field][helix][base] = copy.deepcopy(record)
    return output


def _public_kagome_field(field_rows):
    output = {}
    for internal_helix, records in field_rows.items():
        public_helix = int(internal_helix) - 48
        output[public_helix] = copy.deepcopy(records)
        for record in output[public_helix]:
            for slot in (0, 2):
                partner = int(record[slot])
                if 48 <= partner <= 63:
                    record[slot] = partner - 48
    return output


def _translate_geometry(geometry, shift):
    output = copy.deepcopy(geometry)
    pair_keys = (
        "layer_ranges", "spacing_range", "scaffold_ranges",
        "complement_ranges", "capture_support_ranges", "target_envelope",
    )
    for key in pair_keys:
        value = output.get(key)
        if value is None:
            continue
        if value and isinstance(value[0], list):
            output[key] = [[int(low) + shift, int(high) + shift]
                           for low, high in value]
        else:
            output[key] = [int(item) + shift for item in value]
    output["sst_only_phase_adjustment_bp"] = int(shift)
    return output


def _square_layers(first, spacing, second, array_length):
    sst = SQUARE.load(SQUARE.SST_REFERENCE)
    source = {int(row["num"]): row for row in sst["vstrands"]
              if int(row["num"]) in SQUARE.SST_SOURCE_HELICES}
    requested, placement = SQUARE.desired_duplex_ranges(
        first, spacing, second)
    duplex, scaf_ranges, stap_ranges, placement, global_shift = \
        SQUARE.resolve_complete_sst_ranges(source, requested, placement)
    layers = []
    for index in (0, 1):
        layers.append({
            "scaf": SQUARE.build_field(
                source, "scaf", (scaf_ranges[index],), array_length),
            "stap": SQUARE.build_field(
                source, "stap", (stap_ranges[index],), array_length),
        })
    geometry = SQUARE.desired_duplex_ranges(first, spacing, second)[1]
    from moire_design_core.square_sst_geometry import \
        centered_square_sst_geometry
    geometry = centered_square_sst_geometry(first, spacing, second)
    return layers, tuple(duplex), tuple(scaf_ranges), tuple(stap_ranges), \
        geometry, global_shift


def _kagome_layers(first, spacing, second, array_length):
    unused_resource, source = KAGOME.kagome._source_rows()
    from moire_design_core.square_sst_geometry import \
        centered_square_sst_geometry
    base_geometry = centered_square_sst_geometry(first, spacing, second)
    original_duplex = tuple(tuple(item)
                            for item in base_geometry["layer_ranges"])
    phase_shift = None
    resolved = None
    for candidate in (0, 8, -8, 16, -16, 24, -24, 32, -32):
        candidate_duplex = tuple(
            (low + candidate, high + candidate)
            for low, high in original_duplex)
        try:
            candidate_resolved = [KAGOME.resolve_polymer_ranges(
                source, item, array_length) for item in candidate_duplex]
        except ValueError:
            continue
        phase_shift = candidate
        resolved = candidate_resolved
        duplex = candidate_duplex
        break
    if resolved is None:
        raise ValueError("unable to resolve complete Kagome SST phase")
    scaf_ranges = tuple(item[0] for item in resolved)
    stap_ranges = tuple(item[1] for item in resolved)
    layers = []
    for index in (0, 1):
        layers.append({
            "scaf": _public_kagome_field(KAGOME.build_field(
                source, "scaf", scaf_ranges[index], array_length)),
            "stap": _public_kagome_field(KAGOME.build_field(
                source, "stap", stap_ranges[index], array_length)),
        })
    geometry = _translate_geometry(base_geometry, phase_shift)
    # Seed capture columns stay on the fixed Seed template.  Only the SST is
    # shifted by the minimal 8-bp amount needed to select a legal Kagome edge
    # phase; the complete periodic topology itself is not altered.
    geometry["fixed_seed_capture_positions_by_layer"] = copy.deepcopy(
        base_geometry["seed_capture_positions_by_layer"])
    return layers, duplex, scaf_ranges, stap_ranges, geometry, 0


def _candidate_edges_square(layers, geometry):
    surface = {0, 1, 2, 3, 12, 13, 14, 15}
    output = []
    positions = geometry["seed_capture_positions_by_layer"]
    for layer_index, layer_positions in enumerate(positions, 1):
        rows = layers[layer_index - 1]["stap"]
        for position in layer_positions:
            for helix in sorted(surface):
                record = rows[helix][position]
                for slot in (0, 2):
                    partner, partner_base = map(int, record[slot:slot + 2])
                    if partner in surface and partner != helix and \
                            partner_base == position:
                        output.append({
                            "layer": layer_index, "position": position,
                            "helix": helix, "slot": slot,
                            "partner": partner,
                            "partner_base": partner_base,
                            "family": "square_u_shaped_16nt",
                        })
                        break
    return output


def _kagome_internal_layer_payload(layer, duplex_range, scaf_range,
                                    stap_range, array_length):
    rows = []
    for public in PRIMARY:
        row = {
            "num": public + 48,
            "scaf": copy.deepcopy(layer["scaf"][public]),
            "stap": copy.deepcopy(layer["stap"][public]),
            "loop": [0] * array_length,
            "skip": [0] * array_length,
            "stap_colors": [],
        }
        for field in ("scaf", "stap"):
            for record in row[field]:
                for slot in (0, 2):
                    partner = int(record[slot])
                    if 0 <= partner < 16:
                        record[slot] = partner + 48
        rows.append(row)
    return {
        "name": "kagome_layer",
        "vstrands": rows,
        "moire_structure_metadata": {
            "sst_duplex_ranges": [list(duplex_range)],
            "sst_scaffold_ranges": [list(scaf_range)],
            "sst_staple_ranges": [list(stap_range)],
        },
    }


def _candidate_edges_kagome(layers, duplex, scaf_ranges, stap_ranges,
                             geometry, array_length):
    fixed = geometry["fixed_seed_capture_positions_by_layer"]
    output = []
    for layer_index in (1, 2):
        payload = _kagome_internal_layer_payload(
            layers[layer_index - 1], duplex[layer_index - 1],
            scaf_ranges[layer_index - 1], stap_ranges[layer_index - 1],
            array_length)
        seed_positions = set(map(int, fixed[layer_index - 1]))
        for item in KAGOME.kagome.kagome_capture_anchor_candidates(payload):
            position = int(item["position"])
            helix = int(item["sst_helix"]) - 48
            if position not in seed_positions:
                continue
            if (layers[layer_index - 1]["scaf"][helix][position] == EMPTY or
                    layers[layer_index - 1]["stap"][helix][position] == EMPTY):
                continue
            output.append({
                "layer": layer_index,
                "position": position,
                "helix": helix,
                "slot": int(item["slot"]),
                "partner": (
                    None if int(item["original_partner"][0]) < 0 else
                    int(item["original_partner"][0]) - 48),
                "partner_base": (
                    None if int(item["original_partner"][0]) < 0 else
                    int(item["original_partner"][1])),
                "family": item["capture_family"],
                "origin_type": item["origin_type"],
            })
    return output


def _reciprocal_slot(record, node):
    matches = [slot for slot in (0, 2)
               if tuple(record[slot:slot + 2]) == tuple(node)]
    if len(matches) != 1:
        raise ValueError("capture edge is not uniquely reciprocal")
    return matches[0]


def _open_and_colour_capture(rows, candidates, destination):
    opened_edges = set()
    capture_nodes = []
    for item in candidates:
        layer = int(item["layer"])
        source_node = (int(item["helix"]), int(item["position"]))
        node = destination[layer][("stap",) + source_node]
        capture_nodes.append((layer, node, item))
        if item.get("partner") is None:
            # The reviewed Kagome topology already has a nick at this exact
            # endpoint.  It is a valid capture component and needs no cut.
            continue
        source_peer = (int(item["partner"]), int(item["partner_base"]))
        peer = destination[layer][("stap",) + source_peer]
        edge = tuple(sorted((node, peer)))
        if node == peer or edge in opened_edges:
            continue
        record = rows[node[0]]["stap"][node[1]]
        slot = _reciprocal_slot(record, peer)
        reverse = rows[peer[0]]["stap"][peer[1]]
        reverse_slot = _reciprocal_slot(reverse, node)
        record[slot:slot + 2] = [-1, -1]
        reverse[reverse_slot:reverse_slot + 2] = [-1, -1]
        opened_edges.add(edge)

    component_by_node = {}
    stap_rows = {number: row["stap"] for number, row in rows.items()
                 if number in PRIMARY or
                 AUXILIARY_SHIFT <= number < AUXILIARY_SHIFT + 16}
    for component in _components(stap_rows):
        for node in component:
            component_by_node[node] = component
    coloured = set()
    review = []
    for layer, node, item in capture_nodes:
        component = component_by_node.get(node)
        if not component:
            raise ValueError("capture node has no post-cut component")
        starts = [candidate for candidate in component
                  if rows[candidate[0]]["stap"][candidate[1]][0] == -1]
        if len(starts) != 1:
            raise ValueError("capture component does not have one 5' end")
        start = starts[0]
        if (layer, frozenset(component)) not in coloured:
            entries = rows[start[0]].setdefault("stap_colors", [])
            entries[:] = [entry for entry in entries
                          if int(entry[0]) != int(start[1])]
            entries.append([int(start[1]), CAPTURE_COLOURS[layer]])
            entries.sort(key=lambda entry: int(entry[0]))
            coloured.add((layer, frozenset(component)))
        review.append({
            **copy.deepcopy(item),
            "display_helix": int(node[0]),
            "display_base": int(node[1]),
            "capture_component_5prime": list(start),
            "capture_component_length_nt": len(component),
            "review_colour_hex": "#%06X" % CAPTURE_COLOURS[layer],
        })
    return review, len(opened_edges)


def _audit_reciprocity(rows):
    errors = []
    for field in ("scaf", "stap"):
        for helix, row in rows.items():
            for base, record in enumerate(row[field]):
                for slot in (0, 2):
                    partner, partner_base = map(int, record[slot:slot + 2])
                    if partner < 0:
                        continue
                    if partner not in rows or not 0 <= partner_base < len(
                            rows[partner][field]):
                        errors.append([field, helix, base, "missing"])
                        continue
                    reverse = rows[partner][field][partner_base]
                    if [helix, base] not in (reverse[:2], reverse[2:]):
                        errors.append([field, helix, base, "nonreciprocal"])
    return errors


def _make_rows(combined, auxiliary, seed, array_length, seed_shift):
    seed_geometry = {int(row["num"]): row for row in seed["vstrands"]}
    rows = []
    for helix in PRIMARY:
        geometry = seed_geometry[helix + 48]
        rows.append({
            "row": int(geometry["row"]), "col": int(geometry["col"]),
            "num": helix,
            "scaf": combined["scaf"][helix],
            "stap": combined["stap"][helix],
            "loop": [0] * array_length, "skip": [0] * array_length,
            "scafLoop": [], "stapLoop": [], "stap_colors": [],
        })
    for helix in PRIMARY:
        geometry = seed_geometry[helix + 48]
        number = helix + AUXILIARY_SHIFT
        rows.append({
            "row": int(geometry["row"]),
            "col": int(geometry["col"]) + 12,
            "num": number,
            "scaf": auxiliary["scaf"][number],
            "stap": auxiliary["stap"][number],
            "loop": [0] * array_length, "skip": [0] * array_length,
            "scafLoop": [], "stapLoop": [], "stap_colors": [],
        })
    rows.extend(SQUARE.seed_rows(seed, array_length, seed_shift))
    return sorted(rows, key=lambda row: int(row["num"]))


def build_case(lattice, first, spacing, second, capture_ready):
    array_length = 704
    seed = SQUARE.load(SQUARE.SEED_REFERENCE)
    if lattice == "square":
        layers, duplex, scaf_ranges, stap_ranges, geometry, seed_shift = \
            _square_layers(first, spacing, second, array_length)
        candidates = _candidate_edges_square(layers, geometry)
    elif lattice == "kagome":
        layers, duplex, scaf_ranges, stap_ranges, geometry, seed_shift = \
            _kagome_layers(first, spacing, second, array_length)
        candidates = _candidate_edges_kagome(
            layers, duplex, scaf_ranges, stap_ranges, geometry, array_length)
    else:
        raise ValueError(lattice)

    primary2, auxiliary, destination2, detours = _route_conflicts(
        layers[0], layers[1], array_length)
    combined = _merge(layers[0], primary2)
    rows_list = _make_rows(
        combined, auxiliary, seed, array_length, seed_shift)
    rows = {int(row["num"]): row for row in rows_list}
    destination = {
        1: {
            (field, helix, base): (helix, base)
            for field in ("scaf", "stap")
            for helix, records in layers[0][field].items()
            for base, record in enumerate(records) if record != EMPTY
        },
        2: destination2,
    }
    capture_review = []
    opened_edges = 0
    if capture_ready:
        capture_review, opened_edges = _open_and_colour_capture(
            rows, candidates, destination)

    reciprocity = _audit_reciprocity(rows)
    if reciprocity:
        raise ValueError("reciprocity audit failed: %s" % reciprocity[:5])
    primary_aux_nodes = {
        number: sum(record != EMPTY for field in ("scaf", "stap")
                    for record in rows[number][field])
        for number in range(AUXILIARY_SHIFT, AUXILIARY_SHIFT + 16)
    }
    if not any(primary_aux_nodes.values()):
        raise ValueError("requested auxiliary case did not move any polymer")
    mode = "capture_display" if capture_ready else "complete_sst"
    filename = (f"{lattice}_SST{first}_Z2_{spacing}_SST{second}_"
                f"seed_scaffold_aux16_{mode}_review.json")
    metadata = {
        "review_role": (
            "fixed Seed scaffold + complete SST with review-only auxiliary "
            "routing for Z2=0/8 occupation conflicts"),
        "lattice_type": lattice,
        "SST_Z2_SST": [first, spacing, second],
        "capture_displayed": bool(capture_ready),
        "capture_connection_stage": (
            "SST capture nicks opened and coloured; no Seed staple bridge "
            "drawn" if capture_ready else "complete SST; no capture nick"),
        "helix_numbering": (
            "primary SST 0-15; fixed Seed 16-63; auxiliary SST 64-79"),
        "sst_duplex_ranges": [list(item) for item in duplex],
        "sst_scaffold_ranges": [list(item) for item in scaf_ranges],
        "sst_complementary_chain_ranges": [list(item)
                                             for item in stap_ranges],
        "placement": geometry,
        "auxiliary_review_policy": {
            "production_rule_changed": False,
            "review_only": True,
            "physical_validation_pending": True,
            "collision_definition": (
                "same virtual helix/base occupied by different SST layers, "
                "checked across scaffold and complementary polymers"),
            "routing_unit": (
                "only occupied nodes are detoured; every polymer remains one "
                "connected component and returns to its primary helix"),
            "component_cropping_allowed": False,
            "absolute_crossover_phase_restarted": False,
            "auxiliary_helix_count": 16,
            "auxiliary_helix_numbers": list(range(64, 80)),
            "detoured_components": detours,
        },
        "capture_review_legend": ({
            "layer_1": "#%06X" % CAPTURE_COLOURS[1],
            "layer_2": "#%06X" % CAPTURE_COLOURS[2],
        } if capture_ready else {}),
        "capture_review_assignments": capture_review,
        "audit": {
            "reciprocity_error_count": 0,
            "detoured_component_count": len(detours),
            "auxiliary_nonempty_helix_count": sum(
                count > 0 for count in primary_aux_nodes.values()),
            "capture_candidate_endpoint_count": len(candidates),
            "opened_capture_edge_count": opened_edges,
            "seed_staple_record_count": sum(
                record != EMPTY for number in range(16, 64)
                for record in rows[number]["stap"]),
            "passed": True,
        },
    }
    payload = {
        "name": filename,
        "vstrands": sorted(rows.values(), key=lambda row: int(row["num"])),
        "num_bases": array_length,
        "lattice": "square",
        "moire_structure_metadata": metadata,
    }
    return filename, payload


def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    index = []
    for lattice in ("square", "kagome"):
        for case in CASES:
            for capture_ready in (False, True):
                filename, payload = build_case(
                    lattice, *case, capture_ready=capture_ready)
                target = OUTPUT_ROOT / filename
                target.write_text(json.dumps(payload, separators=(",", ":")),
                                  encoding="utf-8")
                index.append({
                    "file": filename,
                    "lattice": lattice,
                    "SST_Z2_SST": list(case),
                    "capture_displayed": capture_ready,
                    "audit": payload["moire_structure_metadata"]["audit"],
                })
    (OUTPUT_ROOT / "review_index.json").write_text(
        json.dumps(index, ensure_ascii=False, indent=2), encoding="utf-8")
    print(OUTPUT_ROOT)
    print("files", len(index), "all_passed",
          all(item["audit"]["passed"] for item in index))


if __name__ == "__main__":
    main()
