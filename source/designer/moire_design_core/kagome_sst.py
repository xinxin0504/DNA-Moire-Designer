"""Kagome SST-only routing and SST-side capture-anchor preparation.

This module deliberately knows nothing about Seed geometry.  It generates the
validated 12-of-16 Kagome SST cross-section and exposes the legal SST endpoints
that a future Seed-specific contact solver may use.  Square routing remains in
``structure.py`` and must not import Kagome assumptions.
"""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict, deque
from copy import deepcopy
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from .square_sst_geometry import (
    centered_square_sst_geometry,
    refresh_seed_overlap_metadata,
    seed_template_capture_columns,
)
from .sst_auxiliary_routing import (
    AUXILIARY_INTERNAL,
    actual_helix as auxiliary_actual_helix,
    route_layer2_conflicts,
)

INTERNAL_HELICES = tuple(range(48, 64))
ACTIVE_HELICES = (48, 49, 50, 51, 52, 54, 56, 57, 58, 59, 60, 62)
HOLE_HELICES = (53, 55, 61, 63)
LINE_HELICES = {
    "scaf": (52, 54, 60, 62),
    "stap": (48, 50, 56, 58),
}
SOURCE_PHASE_START = {
    "scaf": {16: 16, 0: 32},
    "stap": {8: 8, 24: 24},
}
CAPTURE_ANCHOR_SIDES = {
    # slot 0 is prev/5', slot 2 is next/3' in legacy cadnano records.
    "crossover": ((49, 0), (51, 0), (60, 2), (62, 2)),
    "linear": ((48, 2), (50, 2)),
}
LEFT_COLOR = 5614080
RIGHT_COLOR = 11141375


def _reference_path() -> Path:
    candidates = (
        Path(__file__).with_name("resources") / "kagome_resource_128.json",
        Path.home() / "Desktop" / "kagome_resource_128.json",
    )
    return next((path for path in candidates if path.is_file()), candidates[0])


KAGOME_SST_REFERENCE = _reference_path()


def _empty_record() -> List[int]:
    return [-1, -1, -1, -1]


def _is_empty(record: Sequence[int]) -> bool:
    return list(record) == [-1, -1, -1, -1]


def kagome_active_ranges(z1_bp: int, z2_bp: int, z3_bp: int,
                         global_shift_bp: int = 0
                         ) -> Dict[str, List[List[int]]]:
    """Return validated outward-growing ranges for two Kagome SST layers."""
    z1_bp, z2_bp, z3_bp = map(int, (z1_bp, z2_bp, z3_bp))
    if z1_bp < 64 or z3_bp < 64:
        raise ValueError("Kagome SST 1st layer和2nd layer至少需要64 bp。")
    if any(value % 8 for value in (z1_bp, z2_bp, z3_bp)):
        raise ValueError("Kagome SST长度和spacing必须位于8-bp网格。")
    if z2_bp < 0:
        raise ValueError("Kagome SST spacing不能小于0 bp。")
    shift = int(global_shift_bp)
    if shift % 32:
        raise ValueError("Kagome SST全局平移必须是32 bp整数倍。")
    k1 = (128 - z1_bp) // 8
    k3 = (128 - z3_bp) // 8
    spacing_delta = z2_bp - 32
    return {
        "scaf": [
            [16 + 16 * math.floor(k1 / 2) + shift, 143 + shift],
            [176 + spacing_delta + shift,
             303 + spacing_delta - 16 * math.floor(k3 / 2) + shift],
        ],
        "stap": [
            [8 + 16 * math.ceil(k1 / 2) + shift, 151 + shift],
            [168 + spacing_delta + shift,
             311 + spacing_delta - 16 * math.ceil(k3 / 2) + shift],
        ],
    }


def required_global_shift(z1_bp: int, z2_bp: int = 32,
                          z3_bp: Optional[int] = None,
                          minimum_duplex_start: int = 32) -> int:
    """Reserve left-side Seed-routing room without changing SST phase."""
    z3_bp = int(z1_bp if z3_bp is None else z3_bp)
    ranges = kagome_active_ranges(z1_bp, z2_bp, z3_bp, 0)
    minimum_active = min(low for values in ranges.values()
                         for low, unused_high in values)
    left_duplex_start = max(ranges["scaf"][0][0],
                            ranges["stap"][0][0])
    required = max(0, -minimum_active,
                   int(minimum_duplex_start) - left_duplex_start)
    return 32 * math.ceil(required / 32)


def _line_intervals(low: int, high: int,
                    origin: int) -> List[Tuple[int, int]]:
    """Clip the immutable 32-nt nick lattice at the physical SST boundary.

    ``origin`` is the accepted template nick/capture phase.  It is deliberately
    *not* re-anchored to ``low``: doing that may make more 32-mers, but it also
    moves the legal capture sites.  A clipped 16-nt outer remnant is therefore
    joined to its adjacent 32-nt component, exactly as in the validated Kagome
    fixtures.  The one exceptional 64-nt boundary state is ``16+32+16``:
    only the left remnant is joined to the middle component, yielding the
    template-phased ``48+16`` layout.  Moving every nick by 16 bp to turn
    this into ``32+32`` is forbidden because it changes capture phase.
    """
    low, high, origin = map(int, (low, high, origin))
    total = high-low+1
    if total < 32 or total % 16:
        raise ValueError(
            "Kagome线型SST活动区间必须是至少32 nt的16-nt整数倍。")
    cuts = [low]
    first_nick = origin + 32 * math.ceil((low-origin) / 32)
    for nick in range(first_nick, high+1, 32):
        if low < nick <= high:
            cuts.append(nick)
    cuts.append(high+1)
    intervals = [(cuts[index], cuts[index+1]-1)
                 for index in range(len(cuts)-1)]
    if [end-start+1 for start, end in intervals] == [16, 32, 16]:
        intervals[0:2] = [(intervals[0][0], intervals[1][1])]
        return intervals
    if len(intervals) > 1 and intervals[0][1]-intervals[0][0]+1 < 32:
        intervals[0:2] = [(intervals[0][0], intervals[1][1])]
    if len(intervals) > 1 and intervals[-1][1]-intervals[-1][0]+1 < 32:
        intervals[-2:] = [(intervals[-2][0], intervals[-1][1])]
    if any(end-start+1 not in (32, 48) for start, end in intervals):
        raise ValueError("Kagome线型SST边界无法形成合法32/48-nt组件。")
    return intervals


def _write_linear(records, helix: int, low: int, high: int,
                  field: str) -> None:
    for base in range(low, high + 1):
        if field == "scaf":
            records[helix][base] = [
                helix if base > low else -1,
                base - 1 if base > low else -1,
                helix if base < high else -1,
                base + 1 if base < high else -1,
            ]
        else:
            records[helix][base] = [
                helix if base < high else -1,
                base + 1 if base < high else -1,
                helix if base > low else -1,
                base - 1 if base > low else -1,
            ]


def _components(rows: Dict[int, Dict[str, Any]], field: str):
    graph = defaultdict(set)
    nodes = set()
    for helix, row in rows.items():
        for base, record in enumerate(row.get(field, [])):
            if _is_empty(record):
                continue
            node = (helix, base)
            nodes.add(node)
            for slot in (0, 2):
                partner, partner_base = map(int, record[slot:slot + 2])
                if partner >= 0:
                    graph[node].add((partner, partner_base))
                    graph[(partner, partner_base)].add(node)
    output = []
    unseen = set(nodes)
    while unseen:
        stack = [next(iter(unseen))]
        component = set()
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            stack.extend(graph[node] - component)
        unseen -= component
        output.append(component)
    return output


def _rebuild_colors(rows: Dict[int, Dict[str, Any]], left_high: int) -> None:
    for row in rows.values():
        row["stap_colors"] = []
    for component in _components(rows, "stap"):
        starts = [(helix, base) for helix, base in component
                  if rows[helix]["stap"][base][0] == -1]
        if len(starts) != 1:
            raise ValueError("Kagome SST组件必须恰有一个5′端。")
        helix, base = starts[0]
        color = (LEFT_COLOR if min(index for unused, index in component)
                 <= left_high else RIGHT_COLOR)
        rows[helix]["stap_colors"].append([base, color])
    for row in rows.values():
        row["stap_colors"].sort(key=lambda item: item[0])


@lru_cache(maxsize=1)
def _source_rows() -> Tuple[Dict[str, Any], Dict[int, Dict[str, Any]]]:
    if not KAGOME_SST_REFERENCE.is_file():
        raise FileNotFoundError("找不到kagome_resource_128.json。")
    payload = json.loads(KAGOME_SST_REFERENCE.read_text(encoding="utf-8"))
    by_number = {int(row["num"]): row for row in payload["vstrands"]}
    missing = [number for number in range(16, 32) if number not in by_number]
    if missing:
        raise ValueError("Kagome SST模板缺少helix：%s。" % missing)
    return payload, {number: by_number[number] for number in range(16, 32)}


def _translate_geometry(geometry: Dict[str, Any], shift: int
                        ) -> Dict[str, Any]:
    output = deepcopy(geometry)
    if not shift:
        return output
    pair_keys = (
        "layer_ranges", "spacing_range", "scaffold_ranges",
        "complement_ranges", "seed_layer_ranges", "capture_support_ranges",
        "reference_envelope", "seed_partition_ranges", "target_envelope",
        "optimized_overlap_ranges",
    )
    for key in pair_keys:
        values = output.get(key)
        if values is None:
            continue
        if values and isinstance(values[0], list):
            output[key] = [[int(low)+shift, int(high)+shift]
                           for low, high in values]
        else:
            output[key] = [int(value)+shift for value in values]
    for key in ("seed_capture_positions_by_layer",
                "theoretical_capture_positions_by_layer"):
        output[key] = [[int(value)+shift for value in values]
                       for values in output.get(key, [])]
    output["capture_phase_reference_origin"] = int(
        output["capture_phase_reference_origin"]) + shift
    output["coordinate_shift_bp"] = int(
        output.get("coordinate_shift_bp", 0)) + shift
    output["kagome_additional_seed_routing_shift_bp"] = shift
    return output


def _range_field(source_rows, field: str, active_range, array_length: int):
    low, high = map(int, active_range)
    result = {helix: [_empty_record() for unused in range(array_length)]
              for helix in INTERNAL_HELICES}
    straight = set(LINE_HELICES[field])
    active = set(ACTIVE_HELICES)
    for target_start in range(low, high + 1, 16):
        source_start = SOURCE_PHASE_START[field].get(target_start % 32)
        if source_start is None:
            raise ValueError("Kagome %s绝对相位非法：%d。" %
                             (field, target_start))
        base_shift = target_start - source_start
        for helix in sorted(active - straight):
            source_helix = helix - 32
            for offset in range(16):
                result[helix][target_start + offset] = _translated_record(
                    source_rows[source_helix][field][source_start + offset],
                    base_shift)
    # This origin is immutable: it is shared by the routing template and the
    # capture catalogue.  Never move it merely to obtain prettier lengths.
    origin = 16 if field == "scaf" else 8
    for helix in straight:
        for component_low, component_high in _line_intervals(
                low, high, origin):
            _write_linear(result, helix, component_low, component_high, field)
    return result


def _field_component_lengths(records, field: str):
    rows = {helix: {field: values} for helix, values in records.items()}
    return {len(component) for component in _components(rows, field)}


def _resolve_polymer_ranges(source_rows, duplex_range, array_length: int):
    """Select the unique 32/48-nt Kagome boundary state for one layer."""
    state = _resolve_polymer_state(source_rows, duplex_range, array_length)
    return state["scaffold_range"], state["staple_range"]


def _resolve_polymer_state(source_rows, duplex_range, array_length: int):
    """Select ranges while preserving the immutable template nick phase."""
    low, high = map(int, duplex_range)
    candidates = []
    for scaf_left in (0, 8, 16):
        for scaf_right in (0, 8, 16):
            for stap_left in (0, 8, 16):
                for stap_right in (0, 8, 16):
                    if scaf_left + scaf_right + stap_left + stap_right != 16:
                        continue
                    scaf_range = (low-scaf_left, high+scaf_right)
                    stap_range = (low-stap_left, high+stap_right)
                    if min(scaf_range + stap_range) < 0:
                        continue
                    if (max(scaf_range[0], stap_range[0]) != low or
                            min(scaf_range[1], stap_range[1]) != high):
                        continue
                    try:
                        scaf = _range_field(
                            source_rows, "scaf", scaf_range, array_length)
                        stap = _range_field(
                            source_rows, "stap", stap_range, array_length)
                    except ValueError:
                        continue
                    scaf_lengths = _field_component_lengths(scaf, "scaf")
                    stap_lengths = _field_component_lengths(stap, "stap")
                    # A 16-nt complete-SST component is not a generally
                    # available Kagome edge length.  It is the single
                    # template-phased exception produced when the *entire*
                    # non-scaffold range is exactly 64 nt and clips as
                    # 16+32+16 (later represented as 48+16).  In particular,
                    # a 72-bp duplex uses an 80-nt non-scaffold range and must
                    # remain 32+48, never 48+16 or any 24-nt surrogate.
                    allowed_stap_lengths = {32, 48}
                    if stap_range[1] - stap_range[0] + 1 == 64:
                        allowed_stap_lengths.add(16)
                    if (scaf_lengths <= {32, 48} and
                            stap_lengths <= allowed_stap_lengths):
                        candidates.append({
                            "scaffold_range": scaf_range,
                            "staple_range": stap_range,
                            "scaffold_line_phase_shift_bp": 0,
                            "staple_line_phase_shift_bp": 0,
                        })
    if not candidates:
        raise ValueError("Kagome边缘状态不唯一：%s -> %s。" %
                         (duplex_range, candidates))
    if len(candidates) != 1:
        raise ValueError("Kagome边缘状态不唯一：%s -> %s。" %
                         (duplex_range, candidates))
    return candidates[0]


def build_kagome_layer_fields(duplex_range, array_length: int):
    """Build one phase-correct Kagome layer on an explicit duplex range.

    The mixed Square--Kagome generator uses this public, topology-only
    adapter so layer 1 remains a genuine Square SST while layer 2 follows
    the same 32/48-nt Kagome phase solver as a two-Kagome design.
    """
    unused_resource, source_rows = _source_rows()
    scaffold_range, staple_range = _resolve_polymer_ranges(
        source_rows, duplex_range, int(array_length))
    return {
        "scaf": _range_field(
            source_rows, "scaf", scaffold_range, int(array_length)),
        "stap": _range_field(
            source_rows, "stap", staple_range, int(array_length)),
        "scaffold_range": list(map(int, scaffold_range)),
        "staple_range": list(map(int, staple_range)),
    }


def kagome_layer_capture_catalogue(layer: int, staple_records,
                                    duplex_range):
    """Return the immutable Kagome capture candidates for one layer."""
    duplex_low, duplex_high = map(int, duplex_range)
    output = []
    first_position = duplex_low + ((8-duplex_low) % 16)
    for position in range(first_position, duplex_high + 1, 16):
        phase = position % 32
        family = "crossover" if phase == 24 else "linear"
        for helix, slot in CAPTURE_ANCHOR_SIDES[family]:
            record = list(staple_records[helix][position])
            category = _classify_anchor(record, helix, position, slot)
            if category == "illegal":
                raise ValueError(
                    "Kagome SST非法capture候选 h%d:%d slot%d=%s。" %
                    (helix, position, slot, record))
            output.append({
                "id": "L%d-h%d-b%d-s%d" %
                      (int(layer), helix, position, slot),
                "layer": int(layer),
                "position": int(position),
                "capture_family": (
                    "u_shaped_16nt" if family == "crossover" else
                    "linear_32nt_or_right_edge_16nt"),
                "template_phase_mod32": int(phase),
                "sst_helix": int(helix),
                "logical_sst_helix": int(helix),
                "slot": int(slot),
                "side": "prev" if slot == 0 else "next",
                "origin_type": category,
                "original_partner": list(record[slot:slot + 2]),
                "seed_helix": None,
                "seed_mapping_pending": True,
            })
    return output


def _translated_record(record: Iterable[int], base_shift: int) -> List[int]:
    output = []
    values = list(record)
    for slot in (0, 2):
        partner, partner_base = map(int, values[slot:slot + 2])
        if partner < 0:
            output.extend((-1, -1))
        elif 16 <= partner < 32:
            output.extend((partner + 32, partner_base + base_shift))
        else:
            raise ValueError("Kagome SST模板包含跨出SST截面的连接。")
    return output


def build_kagome_sst_payload(name: str, z1_bp: int = 128,
                             z2_bp: int = 32, z3_bp: int = 128,
                             array_length: Optional[int] = None,
                             layers_design_sequence_identical: Optional[
                                 bool] = None
                             ) -> Dict[str, Any]:
    """Build complete two-layer Kagome SST routing without Seed/capture."""
    z1_bp, z2_bp, z3_bp = map(int, (z1_bp, z2_bp, z3_bp))
    if z1_bp < 64 or z3_bp < 64 or z2_bp < 0 or z2_bp > 160 or \
            any(value % 8 for value in (z1_bp, z2_bp, z3_bp)):
        raise ValueError(
            "Kagome SST长度至少64 bp、spacing为0-160 bp，且均须为8 bp整数倍。")
    resource, source_rows = _source_rows()
    geometry = centered_square_sst_geometry(z1_bp, z2_bp, z3_bp)
    # Keep the shared centre placement, then reserve only the minimum whole
    # 32-bp canvas translation required by the accepted Seed view.
    extra_shift = 0
    first_low = int(geometry["layer_ranges"][0][0])
    while first_low + extra_shift < 32:
        extra_shift += 32
    geometry = _translate_geometry(geometry, extra_shift)
    fixed_seed_columns = deepcopy(geometry["seed_capture_positions_by_layer"])
    base_geometry = deepcopy(geometry)
    resolved = None
    resolved_states = None
    phase_adjustment = None
    for candidate in (0, 8, -8, 16, -16, 24, -24, 32, -32):
        candidate_geometry = deepcopy(base_geometry)
        for key in ("layer_ranges", "capture_support_ranges"):
            candidate_geometry[key] = [
                [int(low)+candidate, int(high)+candidate]
                for low, high in candidate_geometry[key]]
        candidate_geometry["spacing_range"] = [
            int(value)+candidate
            for value in candidate_geometry["spacing_range"]]
        candidate_geometry["target_envelope"] = [
            int(value)+candidate
            for value in candidate_geometry["target_envelope"]]
        candidate_geometry["theoretical_capture_positions_by_layer"] = [
            [int(value)+candidate for value in values]
            for values in candidate_geometry[
                "theoretical_capture_positions_by_layer"]]
        refresh_seed_overlap_metadata(candidate_geometry)
        candidate_ranges = [list(map(int, item))
                            for item in candidate_geometry["layer_ranges"]]
        provisional_length = max(
            544, int(array_length or 0),
            32 * math.ceil((int(candidate_ranges[-1][1]) + 97) / 32))
        try:
            candidate_states = [_resolve_polymer_state(
                source_rows, item, provisional_length)
                for item in candidate_ranges]
        except ValueError:
            continue
        geometry = candidate_geometry
        duplex_ranges = candidate_ranges
        resolved_states = candidate_states
        resolved = [(item["scaffold_range"], item["staple_range"])
                    for item in candidate_states]
        phase_adjustment = candidate
        break
    if resolved is None:
        raise ValueError(
            "Kagome SST无法在范本nick相位上形成合法边界组件。")
    geometry["seed_capture_positions_by_layer"] = fixed_seed_columns
    geometry["sst_only_phase_adjustment_bp"] = int(phase_adjustment)
    ranges = {
        "scaf": [list(item[0]) for item in resolved],
        "stap": [list(item[1]) for item in resolved],
    }
    linear_phase_shifts = [{
        "scaf": int(item["scaffold_line_phase_shift_bp"]),
        "stap": int(item["staple_line_phase_shift_bp"]),
    } for item in resolved_states]
    shift = int(geometry.get("coordinate_shift_bp", 0))
    maximum_base = max(high for values in ranges.values()
                       for unused_low, high in values)
    length = max(544, int(array_length or 0),
                 32 * math.ceil((maximum_base + 65) / 32))
    rows = {}
    layer_records = [{field: {} for field in ("scaf", "stap")}
                     for unused in range(2)]
    for source_number in range(16, 32):
        source = source_rows[source_number]
        number = source_number + 32
        row = deepcopy(source)
        row["num"] = number
        # Same 4x4 cross-section placement used by the validated capture file.
        row["row"] = int(source["row"]) + 5
        row["col"] = int(source["col"]) + 7
        row["scaf"] = [_empty_record() for unused in range(length)]
        row["stap"] = [_empty_record() for unused in range(length)]
        row["loop"] = [0] * length
        row["skip"] = [0] * length
        row["scafLoop"] = []
        row["stapLoop"] = []
        row["stap_colors"] = []
        rows[number] = row
        for layer in layer_records:
            layer["scaf"][number] = [
                _empty_record() for unused in range(length)]
            layer["stap"][number] = [
                _empty_record() for unused in range(length)]

    for field in ("scaf", "stap"):
        for layer_index, active_range in enumerate(ranges[field]):
            layer_records[layer_index][field] = _range_field(
                source_rows, field, active_range, length)

    combined, auxiliary, unused_destination, auxiliary_metadata = \
        route_layer2_conflicts(layer_records[0], layer_records[1], length)
    for number in INTERNAL_HELICES:
        rows[number]["scaf"] = combined["scaf"][number]
        rows[number]["stap"] = combined["stap"][number]
    if auxiliary_metadata["enabled"]:
        for number in AUXILIARY_INTERNAL:
            logical = number - 16
            source = rows[logical]
            row = deepcopy(source)
            row["num"] = number
            row["col"] = int(source["col"]) + 12
            row["scaf"] = auxiliary["scaf"][number]
            row["stap"] = auxiliary["stap"][number]
            row["loop"] = [0] * length
            row["skip"] = [0] * length
            row["stap_colors"] = []
            rows[number] = row

    _rebuild_colors(rows, duplex_ranges[0][1])
    payload = {
        "name": name,
        "vstrands": [rows[number] for number in sorted(rows)],
        "num_bases": length,
        "lattice": resource.get("lattice", "square"),
        "scaffold_colors": resource.get("scaffold_colors", []),
        "moire_structure_metadata": {
            "stage": "kagome_sst_complete_only",
            "lattice_type": "kagome",
            "seed_design_status": "not generated",
            "seed_capture_mapping_status": "pending Seed-specific geometry",
            "z1_bp": int(z1_bp),
            "z2_bp": int(z2_bp),
            "z3_bp": int(z3_bp),
            "sst_scaffold_ranges": deepcopy(ranges["scaf"]),
            "sst_staple_ranges": deepcopy(ranges["stap"]),
            "sst_duplex_ranges": deepcopy(duplex_ranges),
            "layer_ranges": deepcopy(duplex_ranges),
            "global_base_shift_bp": shift,
            "minimum_left_duplex_start": 32,
            "left_layer_growth_direction": "left/outward",
            "right_layer_growth_direction": "right/outward",
            "active_sst_helices_internal": list(ACTIVE_HELICES),
            "hole_sst_helices_internal": list(HOLE_HELICES),
            "u_sst_policy": "add/remove complete 32-nt U-shaped SST",
            "linear_sst_policy": (
                "trim by 8-bp duplex increments; merge a 16-nt edge "
                "remainder with its adjacent 32-nt linear SST to make 48 nt; "
                "preserve a terminal 16-nt component only for the immutable "
                "template-phased 48+16 boundary case"),
            "linear_nick_phase_shifts_by_layer": deepcopy(
                linear_phase_shifts),
            "layers_design_sequence_identical": (
                None if layers_design_sequence_identical is None else
                bool(layers_design_sequence_identical)),
            "capture_gaps_reserved": False,
            "auxiliary_sst_routing": deepcopy(auxiliary_metadata),
            "source_reference": "kagome_resource_128.json",
            "centered_geometry": deepcopy(geometry),
        },
    }
    metadata = payload["moire_structure_metadata"]
    candidates = []
    for layer, (staple_range, duplex_range) in enumerate(
            zip(ranges["stap"], duplex_ranges), 1):
        unused_staple_low, unused_staple_high = staple_range
        duplex_low, duplex_high = duplex_range
        first_position = int(duplex_low) + ((8-int(duplex_low)) % 16)
        for position in range(first_position, int(duplex_high) + 1, 16):
            phase = position % 32
            family = "crossover" if phase == 24 else "linear"
            for helix, slot in CAPTURE_ANCHOR_SIDES[family]:
                record = list(layer_records[layer-1]["stap"][helix][position])
                category = _classify_anchor(record, helix, position, slot)
                if category == "illegal":
                    raise ValueError(
                        "Kagome SST非法capture候选 h%d:%d slot%d=%s。" %
                        (helix, position, slot, record))
                candidates.append({
                    "id": "L%d-h%d-b%d-s%d" %
                          (layer, helix, position, slot),
                    "layer": layer,
                    "position": position,
                    "capture_family": (
                        "u_shaped_16nt" if family == "crossover" else
                        "linear_32nt_or_right_edge_16nt"),
                    "template_phase_mod32": phase,
                    "sst_helix": helix,
                    "logical_sst_helix": helix,
                    "slot": slot,
                    "side": "prev" if slot == 0 else "next",
                    "origin_type": category,
                    "original_partner": list(record[slot:slot + 2]),
                    "seed_helix": None,
                    "seed_mapping_pending": True,
                })
    metadata["kagome_capture_anchor_catalogue_complete"] = deepcopy(
        candidates)
    by_layer = defaultdict(set)
    for item in candidates:
        by_layer[int(item["layer"])].add(int(item["position"]))
    metadata["kagome_theoretical_capture_positions_by_layer"] = [
        sorted(by_layer[index]) for index in (1, 2)]
    # The accepted Seed reference already lives in the public design canvas;
    # ``shift`` here is the Kagome SST router's internal left-room allowance,
    # not an additional Seed translation.  Use the physical template contact
    # columns in the same absolute coordinates as ``duplex_ranges``.
    actual_overlap_ranges = deepcopy(geometry["optimized_overlap_ranges"])
    actual_columns = [
        [int(position) for position in positions
         if int(actual_overlap_ranges[layer][0]) <= int(position) <=
         int(actual_overlap_ranges[layer][1])]
        for layer, positions in enumerate(
            geometry["seed_capture_positions_by_layer"])]
    metadata["kagome_theoretical_capture_anchor_count"] = len(candidates)
    metadata["capture_count_semantics"] = (
        "theoretical SST candidates only; physical captures require the "
        "fixed Seed-template/SST-duplex overlap intersection")
    metadata["variable_length_layout"] = {
        "lattice_type": "kagome",
        "z1_bp": int(z1_bp), "z2_bp": int(z2_bp), "z3_bp": int(z3_bp),
        "array_length": length,
        # Seed placement is part of the shared centred geometry.  Keep it at
        # the top level as well as in ``square_centered_geometry`` so stage 2
        # cannot silently fall back to the unshifted frozen Seed template.
        "coordinate_shift_bp": int(
            geometry.get("coordinate_shift_bp", 0)),
        "seed_layer_ranges": deepcopy(geometry["seed_layer_ranges"]),
        "seed_partition_ranges": deepcopy(
            geometry["seed_partition_ranges"]),
        "seed_partition_lengths_bp": deepcopy(
            geometry["seed_partition_lengths_bp"]),
        "layer_ranges": deepcopy(duplex_ranges),
        "scaffold_ranges": deepcopy(ranges["scaf"]),
        "staple_ranges": deepcopy(ranges["stap"]),
        "theoretical_capture_positions_by_layer": [
            sorted(by_layer[index]) for index in (1, 2)],
        "seed_capture_positions_by_layer": deepcopy(actual_columns),
        "capture_support_ranges": deepcopy(actual_overlap_ranges),
        "capture_phase_reference_origin": 56,
        "overlap_ranges": deepcopy(
            geometry["optimized_overlap_ranges"]),
        "seed_sst_overlap_bp": deepcopy(
            geometry["optimized_seed_overlap_bp"]),
        "actual_capture_positions_by_layer": [[], []],
        "actual_capture_positions": [],
        "actual_capture_count_pending_seed_overlap": True,
        "capture_site_assignments": [],
        "capture_export_site_assignments": [],
        "seed_mapping_pending": True,
        "auxiliary_sst_routing": deepcopy(auxiliary_metadata),
        "square_centered_geometry": deepcopy(geometry),
        "spacing_range": deepcopy(geometry["spacing_range"]),
        "seed_z2_range": deepcopy(geometry["spacing_range"]),
    }
    metadata["variable_length_layout"]["auxiliary_sst_routing"] = deepcopy(
        auxiliary_metadata)
    return payload


def _classify_anchor(record: Sequence[int], helix: int,
                     position: int, slot: int) -> str:
    partner, partner_base = map(int, record[slot:slot + 2])
    if partner < 0:
        return "preexisting_nick"
    if partner == helix and abs(partner_base - position) == 1:
        return "linear_continuous"
    if partner in ACTIVE_HELICES and partner != helix and \
            partner_base == position:
        return "sst_crossover"
    return "illegal"


def kagome_capture_anchor_candidates(payload: Dict[str, Any]
                                      ) -> List[Dict[str, Any]]:
    """Return *theoretical SST-side* Kagome capture endpoints.

    These candidates describe only the immutable SST periodic topology.  They
    are deliberately not the physical Seed--SST captures: an actual capture
    must additionally coincide with a capture coordinate from the fixed Seed
    template and lie inside the current Seed/SST duplex overlap.  Use
    :func:`kagome_actual_capture_anchor_candidates` for that intersection.
    In particular, neither a 128-bp SST nor any other SST length implies a
    fixed number of physical captures or a fixed family for the first one.
    """
    rows = {int(row["num"]): row for row in payload.get("vstrands", [])}
    if set(rows) == set(range(16)):
        raise ValueError("请先将SST-first payload转换为内部helix编号。")
    metadata = payload.get("moire_structure_metadata", {})
    catalogue = metadata.get("kagome_capture_anchor_catalogue_complete")
    if catalogue:
        layout = metadata.get("variable_length_layout", {})
        output = []
        for source in catalogue:
            item = deepcopy(source)
            logical = int(item.get("logical_sst_helix",
                                   item["sst_helix"]))
            actual = auxiliary_actual_helix(
                layout, int(item["layer"]), "stap", logical,
                int(item["position"]))
            row = rows.get(actual)
            if row is None:
                raise ValueError(
                    "Kagome辅助capture helix不存在：h%d。" % actual)
            record = list(row["stap"][int(item["position"])])
            slot = int(item["slot"])
            item.update({
                "sst_helix": actual,
                "logical_sst_helix": logical,
                "original_partner": list(record[slot:slot + 2]),
                "auxiliary_detour": actual != logical,
            })
            output.append(item)
        return output
    staple_ranges = metadata.get("sst_staple_ranges")
    if not staple_ranges:
        raise ValueError("Kagome SST缺少staple活动区间元数据。")
    duplex_ranges = metadata.get("sst_duplex_ranges")
    if not duplex_ranges or len(duplex_ranges) != len(staple_ranges):
        raise ValueError("Kagome SST缺少与staple层对应的duplex区间。")
    output = []
    for layer, ((unused_staple_low, unused_staple_high),
                (duplex_low, duplex_high)) in enumerate(
                    zip(staple_ranges, duplex_ranges), 1):
        # Capture columns come from the immutable 2L Seed template grid, not
        # from ``staple_low+16``.  The first column inside a clipped SST layer
        # may therefore be either a four-endpoint U/crossover column or a
        # two-endpoint linear column.  56/72 in the 3L Kagome Seed reference
        # establish the absolute phases 24/8 (mod 32), respectively.
        first_position = int(duplex_low) + ((8-int(duplex_low)) % 16)
        for position in range(first_position, int(duplex_high) + 1, 16):
            phase = position % 32
            # The fixed Seed template establishes two absolute capture
            # phases: base 24 (mod 32) opens the four U-shaped endpoints and
            # base 8 (mod 32) opens the two linear endpoints.  The SST length
            # never restarts this cycle; the real Seed overlap merely decides
            # which of these absolute columns are physically used first/last.
            family = "crossover" if phase == 24 else "linear"
            for helix, slot in CAPTURE_ANCHOR_SIDES[family]:
                record = list(rows[helix]["stap"][position])
                category = _classify_anchor(record, helix, position, slot)
                if category == "illegal":
                    raise ValueError(
                        "Kagome SST非法capture候选 h%d:%d slot%d=%s。" %
                        (helix, position, slot, record))
                output.append({
                    "id": "L%d-h%d-b%d-s%d" %
                          (layer, helix, position, slot),
                    "layer": layer,
                    "position": position,
                    "capture_family": ("u_shaped_16nt" if
                                       family == "crossover" else
                                       "linear_32nt_or_right_edge_16nt"),
                    "template_phase_mod32": phase,
                    "sst_helix": helix,
                    "slot": slot,
                    "side": "prev" if slot == 0 else "next",
                    "origin_type": category,
                    "original_partner": list(record[slot:slot + 2]),
                    "seed_helix": None,
                    "seed_mapping_pending": True,
                })
    return output


def kagome_actual_capture_anchor_candidates(
        payload: Dict[str, Any],
        seed_capture_positions_by_layer: Iterable[Iterable[int]],
        seed_layer_ranges: Optional[Iterable[Iterable[int]]] = None
        ) -> List[Dict[str, Any]]:
    """Intersect theoretical Kagome endpoints with the real Seed overlap.

    The order of operations is intentionally Seed-first:

    1. read the capture columns from the fixed Seed template after its current
       common coordinate translation;
    2. intersect those columns with the actual SST duplex and optional Seed
       support interval for the same layer;
    3. retain the already classified SST endpoint family at that exact base.

    Consequently the first physical column can be either the two-endpoint
    linear family or the four-endpoint U-shaped family.  No decision is made
    from SST length, local column number, or a hard-coded endpoint count.
    """
    metadata = payload.get("moire_structure_metadata", {})
    duplex_ranges = [tuple(map(int, values)) for values in
                     metadata.get("sst_duplex_ranges", [])]
    seed_columns = [set(map(int, values)) for values in
                    seed_capture_positions_by_layer]
    if len(duplex_ranges) != 2 or len(seed_columns) != 2:
        raise ValueError("Kagome capture重叠必须按两层提供SST和Seed坐标。")
    seed_ranges = (None if seed_layer_ranges is None else
                   [tuple(map(int, values)) for values in seed_layer_ranges])
    if seed_ranges is not None and len(seed_ranges) != 2:
        raise ValueError("Kagome Seed支撑区必须恰好包含两层。")

    rows = {int(row["num"]): row for row in payload.get("vstrands", [])}
    output = []
    for item in kagome_capture_anchor_candidates(payload):
        layer_index = int(item["layer"]) - 1
        position = int(item["position"])
        low, high = duplex_ranges[layer_index]
        if seed_ranges is not None:
            low = max(low, seed_ranges[layer_index][0])
            high = min(high, seed_ranges[layer_index][1])
        if position not in seed_columns[layer_index] or not low <= position <= high:
            continue
        # A coordinate at a polymer overhang is not a capture site even when
        # it shares the correct absolute 16-bp phase.  Require this exact SST
        # endpoint to be double stranded.  This preserves the template rule:
        # a single-stranded start is skipped, while a genuinely duplex first
        # column (such as layer-2 base 200) is retained.
        row = rows.get(int(item["sst_helix"]))
        if row is None or position >= len(row.get("scaf", [])) or \
                position >= len(row.get("stap", [])):
            continue
        if (row["scaf"][position] == [-1, -1, -1, -1] or
                row["stap"][position] == [-1, -1, -1, -1]):
            continue
        physical = deepcopy(item)
        physical.update({
            "candidate_scope": "actual_seed_sst_duplex_overlap",
            "seed_overlap_low": low,
            "seed_overlap_high": high,
        })
        output.append(physical)
    return output


def _disconnect(rows: Dict[int, Dict[str, Any]], helix: int, position: int,
                slot: int) -> int:
    record = rows[helix]["stap"][position]
    partner, partner_base = map(int, record[slot:slot + 2])
    if partner < 0:
        return 0
    reverse = rows[partner]["stap"][partner_base]
    reverse_slots = [other for other in (0, 2)
                     if reverse[other:other + 2] == [helix, position]]
    if len(reverse_slots) != 1:
        raise ValueError("Kagome SST待切连接不是唯一互反连接。")
    record[slot:slot + 2] = [-1, -1]
    other = reverse_slots[0]
    reverse[other:other + 2] = [-1, -1]
    return 2


def prepare_kagome_capture_sites(
        source: Dict[str, Any], name: Optional[str] = None,
        positions: Optional[Iterable[int]] = None,
        seed_capture_positions_by_layer: Optional[
            Iterable[Iterable[int]]] = None,
        seed_layer_ranges: Optional[Iterable[Iterable[int]]] = None
        ) -> Dict[str, Any]:
    """Open selected legal SST endpoints without assigning any Seed helix.

    All anchor types are classified before any mutation.  This prevents a cut
    made for one capture from turning a later linear-continuous site into an
    apparent pre-existing nick.
    """
    payload = deepcopy(source)
    if name:
        payload["name"] = name
    baseline = (kagome_capture_anchor_candidates(source)
                if seed_capture_positions_by_layer is None else
                kagome_actual_capture_anchor_candidates(
                    source, seed_capture_positions_by_layer,
                    seed_layer_ranges))
    selected_positions = (None if positions is None else
                          {int(value) for value in positions})
    selected = [item for item in baseline
                if selected_positions is None or
                int(item["position"]) in selected_positions]
    rows = {int(row["num"]): row for row in payload["vstrands"]}
    changed_slots = 0
    for item in selected:
        if item["origin_type"] == "preexisting_nick":
            continue
        changed_slots += _disconnect(
            rows, int(item["sst_helix"]), int(item["position"]),
            int(item["slot"]))

    component_by_node = {}
    for component in _components(rows, "stap"):
        for node in component:
            component_by_node[node] = component
    prepared = []
    for item in selected:
        item = deepcopy(item)
        component = component_by_node.get(
            (int(item["sst_helix"]), int(item["position"])), set())
        item["capture_extension_nt"] = len(component)
        item["extension_supported_without_seed_exception"] = \
            len(component) in (16, 32, 48)
        prepared.append(item)
    unsupported = [item for item in prepared
                   if not item["extension_supported_without_seed_exception"]]
    if unsupported:
        raise ValueError("Kagome capture延伸产生非16/32/48-nt组件：%s。" %
                         [(item["id"], item["capture_extension_nt"])
                          for item in unsupported])

    metadata = payload.setdefault("moire_structure_metadata", {})
    counts = Counter(item["origin_type"] for item in prepared)
    extensions = Counter(int(item["capture_extension_nt"])
                         for item in prepared)
    selected_positions_by_layer = {1: set(), 2: set()}
    for item in prepared:
        selected_positions_by_layer[int(item["layer"])].add(
            int(item["position"]))
    selected_positions_by_layer = [
        sorted(selected_positions_by_layer[layer]) for layer in (1, 2)]
    metadata.update({
        "stage": "kagome_sst_capture_ready",
        "capture_gaps_reserved": True,
        "capture_sites_prepared": True,
        "capture_seed_mapping_pending": True,
        "capture_anchor_classification_baseline":
            "immutable complete Kagome SST routing",
        "capture_candidate_scope": (
            "theoretical_sst_candidates" if
            seed_capture_positions_by_layer is None else
            "actual_fixed_seed_sst_duplex_overlap"),
        "actual_capture_positions_by_layer": (
            selected_positions_by_layer if
            seed_capture_positions_by_layer is not None else [[], []]),
        "actual_capture_positions": (
            [position for values in selected_positions_by_layer
             for position in values] if
            seed_capture_positions_by_layer is not None else []),
        "actual_capture_anchor_count": (
            len(prepared) if seed_capture_positions_by_layer is not None
            else 0),
        "actual_capture_count_pending_seed_overlap": (
            seed_capture_positions_by_layer is None),
        "kagome_capture_anchor_assignments_sst_only": prepared,
        "kagome_capture_anchor_type_counts": dict(counts),
        "kagome_capture_extension_counts": {
            str(length): count for length, count in sorted(extensions.items())},
        "capture_modified_record_slot_count": changed_slots,
        "complete_sst_source_unchanged": True,
    })
    variable_layout = metadata.get("variable_length_layout")
    if isinstance(variable_layout, dict):
        variable_layout["actual_capture_positions_by_layer"] = deepcopy(
            metadata["actual_capture_positions_by_layer"])
        variable_layout["actual_capture_positions"] = list(
            metadata["actual_capture_positions"])
        variable_layout["actual_capture_count_pending_seed_overlap"] = bool(
            metadata["actual_capture_count_pending_seed_overlap"])
    return payload


def validate_kagome_sst_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Validate Kagome cross-section, reciprocity, components and metadata."""
    rows = {int(row["num"]): row for row in payload.get("vstrands", [])}
    metadata = payload.get("moire_structure_metadata", {})
    errors = []
    warnings = []
    auxiliary_enabled = bool(payload.get(
        "moire_structure_metadata", {}).get(
            "auxiliary_sst_routing", {}).get("enabled"))
    expected_rows = set(INTERNAL_HELICES) | (
        set(AUXILIARY_INTERNAL) if auxiliary_enabled else set())
    if set(rows) != expected_rows or len(payload.get("vstrands", [])) != len(
            expected_rows):
        errors.append(
            "Kagome SST helix集合错误：应为%s。" % sorted(expected_rows))
    for helix in HOLE_HELICES:
        for actual in (helix, helix + 16):
            row = rows.get(actual)
            if row and any(not _is_empty(record)
                           for field in ("scaf", "stap")
                           for record in row.get(field, [])):
                errors.append("Kagome空位helix %d不得含有strand。" % actual)
    component_report = {}
    capture_ready = bool(payload.get("moire_structure_metadata", {}).get(
        "capture_sites_prepared"))
    for field in ("scaf", "stap"):
        active_helices = set()
        for helix, row in rows.items():
            for base, record in enumerate(row.get(field, [])):
                if _is_empty(record):
                    continue
                active_helices.add(helix)
                for slot in (0, 2):
                    partner, partner_base = map(int, record[slot:slot + 2])
                    if partner < 0:
                        continue
                    if partner not in rows or not (0 <= partner_base <
                                                   len(rows[partner][field])):
                        errors.append("%s越界连接 h%d:%d。" %
                                      (field, helix, base))
                        continue
                    reverse = rows[partner][field][partner_base]
                    if ([helix, base] not in
                            (reverse[0:2], reverse[2:4])):
                        errors.append("%s非互反连接 h%d:%d。" %
                                      (field, helix, base))
        logical_active = {
            helix - 16 if helix in AUXILIARY_INTERNAL else helix
            for helix in active_helices}
        if logical_active != set(ACTIVE_HELICES):
            errors.append("%s逻辑活动helix错误：%s。" %
                          (field, sorted(logical_active)))
        components = _components(rows, field)
        lengths = Counter(len(component) for component in components)
        illegal = []
        staple_ranges = [tuple(map(int, bounds)) for bounds in
                          metadata.get("sst_staple_ranges", [])]
        short_staple_ranges = [bounds for bounds in staple_ranges
                                if bounds[1] - bounds[0] + 1 == 64]
        for component in components:
            length = len(component)
            if length in (32, 48):
                continue
            if length == 16 and field == "stap" and capture_ready:
                continue
            if length == 16 and field == "stap":
                logical_helices = {
                    helix-16 if helix in AUXILIARY_INTERNAL else helix
                    for helix, unused_base in component}
                bases = {base for unused_helix, base in component}
                if (logical_helices <= set(LINE_HELICES["stap"]) and
                        any(low <= min(bases) <= max(bases) <= high and
                            (min(bases) == low or max(bases) == high)
                            for low, high in short_staple_ranges)):
                    continue
            illegal.append(length)
        illegal = sorted(set(illegal))
        if illegal:
            errors.append("%s含非法SST组件长度：%s。" % (field, illegal))
        component_report[field] = {
            "component_count": sum(lengths.values()),
            "length_counts": dict(sorted(lengths.items())),
        }
    if metadata.get("lattice_type") != "kagome":
        errors.append("Kagome SST缺少lattice_type=kagome元数据。")
    duplex_ranges = metadata.get("sst_duplex_ranges", [])
    expected_lengths = [metadata.get("z1_bp"), metadata.get("z3_bp")]
    if len(duplex_ranges) != 2:
        errors.append("Kagome SST必须包含两层duplex range。")
    else:
        for bounds, expected in zip(duplex_ranges, expected_lengths):
            if int(bounds[1]) - int(bounds[0]) + 1 != int(expected):
                errors.append("Kagome SST duplex长度与参数不一致。")
        if int(duplex_ranges[0][0]) < 32:
            errors.append("Kagome SST起始双链位置必须大于等于32。")
        spacing = int(duplex_ranges[1][0]) - int(duplex_ranges[0][1]) - 1
        if spacing != int(metadata.get("z2_bp", spacing)):
            errors.append("Kagome SST层间距与Z2不一致。")
    return {
        "valid": not errors,
        "errors": errors,
        "warnings": warnings,
        "lattice_type": "kagome",
        "sst_ranges": [
            {"layer": index + 1, "range": list(bounds)}
            for index, bounds in enumerate(duplex_ranges)],
        "active_helices": list(ACTIVE_HELICES),
        "hole_helices": list(HOLE_HELICES),
        "components": component_report,
        "capture_anchor_type_counts": metadata.get(
            "kagome_capture_anchor_type_counts", {}),
        "capture_extension_counts": metadata.get(
            "kagome_capture_extension_counts", {}),
    }
