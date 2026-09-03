#!/usr/bin/env python3
"""Generate centred complete-Kagome-SST review designs.

This is deliberately a review-only generator.  It reuses the accepted
Kagome 12-of-16 topology, places SST1/Z2/SST2 with the same centred geometry
used by the Square review, and keeps the accepted two-layer Seed scaffold as
an immutable positional reference.  Capture gaps are not opened here.
"""

from __future__ import annotations

import copy
import json
import shutil
from collections import Counter
from pathlib import Path

from moire_design_core import kagome_sst as kagome
from moire_design_core.square_sst_geometry import (
    centered_square_sst_geometry,
)


ROOT = Path(__file__).resolve().parents[1]
RESOURCE_ROOT = ROOT / "moire_design_core" / "resources"
SEED_REFERENCE = RESOURCE_ROOT / "Square_Seed_2L_newtemplate.json"
OUTPUT_ROOT = ROOT / "review_outputs" / \
    "Z2_position_preview_complete_SST_final_kagome_20260813"

CASES = (
    ("linked", 128, 32, 128, True),
    ("linked", 120, 40, 120, True),
    ("linked", 112, 48, 112, True),
    ("independent", 128, 40, 96, False),
    ("independent", 96, 56, 128, False),
    ("linked", 136, 56, 136, True),
    ("linked", 152, 72, 152, True),
    ("linked", 160, 64, 160, True),
    ("linked", 96, 64, 96, True),
    ("linked", 80, 48, 80, True),
    ("independent", 136, 32, 136, False),
    ("independent", 80, 32, 80, False),
    ("linked", 160, 96, 160, True),
    ("independent", 128, 16, 128, False),
    ("independent", 128, 24, 128, False),
)

EMPTY = [-1, -1, -1, -1]
SST_INTERNAL_HELICES = tuple(range(48, 64))
SST_PUBLIC_HELICES = tuple(range(16))
SEED_SOURCE_HELICES = tuple(range(48))
SEED_PUBLIC_SHIFT = 16


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def translated_record(record, base_shift):
    return kagome._translated_record(record, base_shift)


def build_field(source_rows, field, active_range, array_length):
    """Build one absolute-phase Kagome polymer range."""
    low, high = map(int, active_range)
    if low < 0 or (high-low+1) % 16:
        raise ValueError("Kagome polymer range must contain 16-bp domains")
    result = {
        helix: [EMPTY[:] for unused in range(array_length)]
        for helix in SST_INTERNAL_HELICES
    }
    straight = set(kagome.LINE_HELICES[field])
    active = set(kagome.ACTIVE_HELICES)
    for target_start in range(low, high+1, 16):
        source_start = kagome.SOURCE_PHASE_START[field].get(
            target_start % 32)
        if source_start is None:
            raise ValueError(
                "Unknown Kagome %s absolute phase at base %d" %
                (field, target_start))
        base_shift = target_start-source_start
        for helix in sorted(active-straight):
            source_helix = helix-32
            for offset in range(16):
                result[helix][target_start+offset] = translated_record(
                    source_rows[source_helix][field][source_start+offset],
                    base_shift)
    line_origin = 16 if field == "scaf" else 8
    for helix in sorted(straight):
        for component_low, component_high in kagome._line_intervals(
                low, high, line_origin):
            kagome._write_linear(
                result, helix, component_low, component_high, field)
    return result


def field_component_lengths(field_rows, field):
    rows = {
        helix: {field: records}
        for helix, records in field_rows.items()
    }
    return sorted(len(component)
                  for component in kagome._components(rows, field))


def resolve_polymer_ranges(source_rows, duplex_range, array_length):
    """Find the unique legal blue/grey Kagome boundary state."""
    low, high = map(int, duplex_range)
    candidates = []
    for scaf_left in (0, 8, 16):
        for scaf_right in (0, 8, 16):
            for stap_left in (0, 8, 16):
                for stap_right in (0, 8, 16):
                    if (scaf_left+scaf_right+stap_left+stap_right) != 16:
                        continue
                    scaf_range = (low-scaf_left, high+scaf_right)
                    stap_range = (low-stap_left, high+stap_right)
                    if min(scaf_range+stap_range) < 0:
                        continue
                    if (max(scaf_range[0], stap_range[0]) != low or
                            min(scaf_range[1], stap_range[1]) != high):
                        continue
                    try:
                        scaf = build_field(
                            source_rows, "scaf", scaf_range, array_length)
                        stap = build_field(
                            source_rows, "stap", stap_range, array_length)
                    except ValueError:
                        continue
                    scaf_lengths = set(field_component_lengths(scaf, "scaf"))
                    stap_lengths = set(field_component_lengths(stap, "stap"))
                    if scaf_lengths <= {32, 48} and stap_lengths <= {32, 48}:
                        candidates.append((scaf_range, stap_range))
    if len(candidates) != 1:
        raise ValueError(
            "Kagome boundary is not uniquely resolvable for duplex %s: %s" %
            (duplex_range, candidates))
    return candidates[0]


def merge_fields(first, second):
    result = copy.deepcopy(first)
    for helix in result:
        for base, record in enumerate(second[helix]):
            if record == EMPTY:
                continue
            if result[helix][base] != EMPTY:
                raise ValueError(
                    "same-polymer layer overlap h%d:%d" % (helix, base))
            result[helix][base] = copy.deepcopy(record)
    return result


def shifted_seed_rows(seed, array_length, base_shift):
    rows = []
    for source in seed["vstrands"]:
        if int(source["num"]) not in SEED_SOURCE_HELICES:
            continue
        row = copy.deepcopy(source)
        row["num"] = int(source["num"])+SEED_PUBLIC_SHIFT
        row["scaf"] = [EMPTY[:] for unused in range(array_length)]
        row["stap"] = [EMPTY[:] for unused in range(array_length)]
        row["loop"] = [0]*array_length
        row["skip"] = [0]*array_length
        row["stap_colors"] = []
        for source_base, record in enumerate(source.get("scaf", [])):
            target_base = source_base+base_shift
            if not 0 <= target_base < array_length:
                continue
            translated = copy.deepcopy(record)
            for slot in (0, 2):
                if translated[slot] >= 0:
                    translated[slot] += SEED_PUBLIC_SHIFT
                    translated[slot+1] += base_shift
            row["scaf"][target_base] = translated
        for field in ("loop", "skip"):
            for source_base, value in enumerate(source.get(field, [])):
                target_base = source_base+base_shift
                if 0 <= target_base < array_length:
                    row[field][target_base] = value
        rows.append(row)
    return rows


def internal_sst_payload(name, source_rows, duplex_ranges,
                         scaf_ranges, stap_ranges, array_length):
    scaf = merge_fields(
        build_field(source_rows, "scaf", scaf_ranges[0], array_length),
        build_field(source_rows, "scaf", scaf_ranges[1], array_length))
    stap = merge_fields(
        build_field(source_rows, "stap", stap_ranges[0], array_length),
        build_field(source_rows, "stap", stap_ranges[1], array_length))
    rows = []
    for helix in SST_INTERNAL_HELICES:
        source = source_rows[helix-32]
        row = {
            "row": int(source["row"])+5,
            "col": int(source["col"])+7,
            "num": helix,
            "scaf": scaf[helix],
            "stap": stap[helix],
            "loop": [0]*array_length,
            "skip": [0]*array_length,
            "scafLoop": [],
            "stapLoop": [],
            "stap_colors": [],
        }
        rows.append(row)
    kagome._rebuild_colors(
        {int(row["num"]): row for row in rows}, duplex_ranges[0][1])
    return {
        "name": name,
        "vstrands": rows,
        "num_bases": array_length,
        "lattice": "square",
        "moire_structure_metadata": {
            "stage": "kagome_sst_complete_only",
            "lattice_type": "kagome",
            "z1_bp": duplex_ranges[0][1]-duplex_ranges[0][0]+1,
            "z2_bp": duplex_ranges[1][0]-duplex_ranges[0][1]-1,
            "z3_bp": duplex_ranges[1][1]-duplex_ranges[1][0]+1,
            "sst_duplex_ranges": [list(item) for item in duplex_ranges],
            "sst_scaffold_ranges": [list(item) for item in scaf_ranges],
            "sst_staple_ranges": [list(item) for item in stap_ranges],
            "capture_gaps_reserved": False,
        },
    }


def to_public_sst_rows(payload, seed_geometry):
    output = []
    for source in payload["vstrands"]:
        old = int(source["num"])
        row = copy.deepcopy(source)
        row["num"] = old-48
        geometry = seed_geometry[old]
        row["row"] = int(geometry["row"])
        row["col"] = int(geometry["col"])
        for field in ("scaf", "stap"):
            for record in row[field]:
                for slot in (0, 2):
                    if 48 <= int(record[slot]) <= 63:
                        record[slot] -= 48
        output.append(row)
    return output


def interval_intersection(first, second):
    low = max(first[0], second[0])
    high = min(first[1], second[1])
    return [low, high] if low <= high else None


def reserve_kagome_seed_routing_space(geometry):
    """Translate Seed and both SSTs together until duplex begins at base 32."""
    output = copy.deepcopy(geometry)
    first_low = int(output["layer_ranges"][0][0])
    extra_shift = 0
    while first_low+extra_shift < 32:
        extra_shift += 32
    if not extra_shift:
        return output
    pair_keys = (
        "layer_ranges", "spacing_range", "scaffold_ranges",
        "complement_ranges", "seed_layer_ranges", "capture_support_ranges",
        "reference_envelope",
        "seed_partition_ranges", "target_envelope",
    )
    for key in pair_keys:
        values = output.get(key)
        if values is None:
            continue
        if values and isinstance(values[0], list):
            output[key] = [
                [int(low)+extra_shift, int(high)+extra_shift]
                for low, high in values]
        else:
            output[key] = [int(value)+extra_shift for value in values]
    output["seed_capture_positions_by_layer"] = [
        [int(value)+extra_shift for value in values]
        for values in output["seed_capture_positions_by_layer"]]
    output["theoretical_capture_positions_by_layer"] = [
        [int(value)+extra_shift for value in values]
        for values in output["theoretical_capture_positions_by_layer"]]
    output["capture_phase_reference_origin"] = int(
        output["capture_phase_reference_origin"])+extra_shift
    output["coordinate_shift_bp"] = (
        int(output.get("coordinate_shift_bp", 0))+extra_shift)
    output["kagome_additional_seed_routing_shift_bp"] = extra_shift
    return output


def build_case(kind, first_length, spacing, second_length, linked,
               seed, source_rows):
    if linked and (first_length != second_length or
                   (first_length+spacing) % 32):
        raise ValueError("linked Kagome layers require identical length and 32-bp translation")
    geometry = reserve_kagome_seed_routing_space(
        centered_square_sst_geometry(first_length, spacing, second_length))
    duplex = tuple(tuple(item) for item in geometry["layer_ranges"])
    provisional_length = max(704, duplex[-1][1]+128)
    resolved = [resolve_polymer_ranges(
        source_rows, item, provisional_length) for item in duplex]
    scaf_ranges = tuple(item[0] for item in resolved)
    stap_ranges = tuple(item[1] for item in resolved)
    array_length = max(
        704,
        32*((max(scaf_ranges[-1][1], stap_ranges[-1][1])+97+31)//32))
    filename = (f"{kind}_SST{first_length}_Z2_{spacing}_SST"
                f"{second_length}_complete_kagome_review.json")
    internal = internal_sst_payload(
        filename, source_rows, duplex, scaf_ranges, stap_ranges, array_length)
    validation = kagome.validate_kagome_sst_payload(internal)
    if not validation["valid"]:
        raise AssertionError(validation["errors"])
    seed_geometry = {
        int(row["num"]): row for row in seed["vstrands"]
        if 48 <= int(row["num"]) <= 63
    }
    rows = to_public_sst_rows(internal, seed_geometry)
    rows.extend(shifted_seed_rows(
        seed, array_length, int(geometry["coordinate_shift_bp"])))
    rows.sort(key=lambda row: int(row["num"]))
    spacing_range = (duplex[0][1]+1, duplex[1][0]-1)
    scaf_counts = Counter(validation["components"]["scaf"]["length_counts"])
    stap_counts = Counter(validation["components"]["stap"]["length_counts"])
    metadata = {
        "review_role": "fixed Seed scaffold + complete Kagome SST scaffold/complementary chain",
        "helix_numbering": "SST 0-15; fixed Seed scaffold 16-63",
        "seed_reference": "Square_Seed_2L_newtemplate.json; scaffold topology unchanged",
        "sst_topology_reference": "Kagome_SST_original_128.json helices 16-31",
        "SST_Z2_SST": [first_length, spacing, second_length],
        "linked": bool(linked),
        "sst_duplex_ranges": [list(item) for item in duplex],
        "sst_scaffold_ranges": [list(item) for item in scaf_ranges],
        "sst_complementary_chain_ranges": [list(item) for item in stap_ranges],
        "spacing_range": list(spacing_range),
        "spacing_definition": "duplex-to-duplex coordinate gap; Kagome linear/U-shaped single strands may occupy part of it",
        "capture_gaps_reserved": False,
        "capture_connection_stage": "not applied",
        "active_kagome_sst_helices_public": [value-48 for value in kagome.ACTIVE_HELICES],
        "hole_kagome_sst_helices_public": [value-48 for value in kagome.HOLE_HELICES],
        "single_strands_inside_Z2": {
            "first_SST_blue_scaffold": interval_intersection(scaf_ranges[0], spacing_range),
            "first_SST_grey_complement": interval_intersection(stap_ranges[0], spacing_range),
            "second_SST_blue_scaffold": interval_intersection(scaf_ranges[1], spacing_range),
            "second_SST_grey_complement": interval_intersection(stap_ranges[1], spacing_range),
        },
        "placement": geometry,
        "audit": {
            "valid_kagome_topology": True,
            "hole_helices_empty": True,
            "scaffold_component_lengths": validation["components"]["scaf"]["length_counts"],
            "complement_component_lengths": validation["components"]["stap"]["length_counts"],
            "same_polymer_layer_overlap_count": 0,
            "capture_nick_count": 0,
            "passed": True,
        },
    }
    payload = {
        "name": filename,
        "vstrands": rows,
        "num_bases": array_length,
        "lattice": "square",
        "moire_structure_metadata": metadata,
    }
    return filename, payload, metadata


def main():
    seed = load(SEED_REFERENCE)
    unused_resource, source_rows = kagome._source_rows()
    if OUTPUT_ROOT.exists():
        shutil.rmtree(OUTPUT_ROOT)
    OUTPUT_ROOT.mkdir(parents=True)
    index = []
    for case in CASES:
        filename, payload, metadata = build_case(
            *case, seed, source_rows)
        (OUTPUT_ROOT / filename).write_text(
            json.dumps(payload, separators=(",", ":")), encoding="utf-8")
        index.append({"file": filename, **metadata})
    (OUTPUT_ROOT / "review_index.json").write_text(
        json.dumps(index, ensure_ascii=False, indent=2), encoding="utf-8")
    archive = shutil.make_archive(str(OUTPUT_ROOT), "zip", root_dir=OUTPUT_ROOT)
    print(OUTPUT_ROOT)
    print(archive)
    print("files", len(index), "all_passed", all(
        item["audit"]["passed"] for item in index))


if __name__ == "__main__":
    main()
