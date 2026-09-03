#!/usr/bin/env python3
"""Generate review-only Seed-scaffold + complete-SST caDNAno JSON files.

This script deliberately does not call the variable-length SST generator.  It
uses the reviewed SST helices 16--31 as a single absolute 32-bp phase table.
Changing SST length crops or extends that table; it never restarts the phase at
a local layer boundary.
"""

from __future__ import annotations

import copy
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESOURCE_ROOT = ROOT / "moire_design_core" / "resources"
SST_REFERENCE = RESOURCE_ROOT / "Square_SST_original_128.json"
SEED_REFERENCE = RESOURCE_ROOT / "Square_Seed_2L_newtemplate.json"
OUTPUT_ROOT = ROOT / "review_outputs" / \
    "Z2_position_preview_complete_SST_final_square_20260813"

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
)

EMPTY = [-1, -1, -1, -1]
SST_SOURCE_HELICES = tuple(range(16, 32))
SST_TARGET_HELICES = tuple(range(16))
SEED_SOURCE_HELICES = tuple(range(48))
SEED_TARGET_SHIFT = 16
PHASE_ANCHOR_BASE = 48
PHASE_REPEAT_BP = 32
CANONICAL_SOURCE_BASE = 64
REFERENCE_ENVELOPE = (48, 335)


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def pad(values, length, fill):
    result = copy.deepcopy(list(values[:length]))
    result.extend(copy.deepcopy(fill) for unused in range(length-len(result)))
    return result


def desired_duplex_ranges(first_length: int, spacing: int,
                          second_length: int):
    for label, value in (("first SST", first_length),
                         ("Z2 spacing", spacing),
                         ("second SST", second_length)):
        if value < 0 or value % 8:
            raise ValueError(f"{label} must be a non-negative multiple of 8 bp")
    if first_length == 0 or second_length == 0:
        raise ValueError("SST duplex lengths must be greater than zero")
    if spacing > 160:
        raise ValueError("Z2 spacing must be within 0-160 bp")
    # Position the complete SST1/Z2/SST2 envelope around the centre of the
    # reviewed 128/32/128 reference.  Boundaries remain on the 8-bp domain
    # grid.  An odd number of spare 8-bp domains is assigned to the right;
    # the next step restores exact symmetry.  This is the same deterministic
    # left/right alternation used by the side-view partition preview.
    reference_length = REFERENCE_ENVELOPE[1]-REFERENCE_ENVELOPE[0]+1
    target_length = first_length+spacing+second_length
    spare_domains = (reference_length-target_length)//8
    left_domains = spare_domains//2
    right_domains = spare_domains-left_domains
    ideal_start = REFERENCE_ENVELOPE[0] + left_domains*8
    first_start = ideal_start
    first = (first_start, first_start + first_length - 1)
    second = (first[1] + 1 + spacing,
              first[1] + spacing + second_length)
    placement = {
        "reference_envelope": list(REFERENCE_ENVELOPE),
        "target_envelope": [first[0], second[1]],
        "ideal_start_before_base0_limit": ideal_start,
        "global_32bp_canvas_shift": 0,
        "left_margin_bp": first[0]-REFERENCE_ENVELOPE[0],
        "right_margin_bp": REFERENCE_ENVELOPE[1]-second[1],
        "envelope_center_offset_bp": (
            (first[0]+second[1])/2.0 -
            sum(REFERENCE_ENVELOPE)/2.0),
        "z2_center_offset_bp": (
            ((first[1]+1)+(second[0]-1))/2.0 -
            ((175+1)+(208-1))/2.0),
        "placement_step_bp": 8,
        "placement_policy": "reference-centred; odd 8-bp remainder alternates and differs by at most 8 bp",
    }
    return (first, second), placement


def canonical_source_base(target_base: int) -> int:
    return CANONICAL_SOURCE_BASE + (
        (target_base-CANONICAL_SOURCE_BASE) % PHASE_REPEAT_BP)


def translated_record(source_record, source_base, target_base,
                      active_range):
    output = []
    low, high = active_range
    for side in (0, 2):
        partner, partner_base = map(int, source_record[side:side+2])
        if partner < 0:
            output.extend((-1, -1))
            continue
        partner_target = target_base + (partner_base-source_base)
        partner_target_helix = partner-16
        if not (0 <= partner_target_helix < 16 and
                low <= partner_target <= high):
            output.extend((-1, -1))
        else:
            output.extend((partner_target_helix, partner_target))
    return output


def build_field(source_rows, field, ranges, array_length):
    result = {number: [EMPTY[:] for unused in range(array_length)]
              for number in SST_TARGET_HELICES}
    for active_range in ranges:
        for target_base in range(active_range[0], active_range[1]+1):
            source_base = canonical_source_base(target_base)
            for target_helix in SST_TARGET_HELICES:
                source_record = source_rows[target_helix+16][field][source_base]
                result[target_helix][target_base] = translated_record(
                    source_record, source_base, target_base, active_range)
    return result


def component_lengths_for_field(field_rows, field):
    rows = [{"num": helix, field: records}
            for helix, records in field_rows.items()]
    return audit_component_lengths(rows, field)


def legal_polymer_ranges(source_rows, duplex_range, array_length):
    """Resolve a complete Square SST boundary without cutting a 32-nt U.

    The blue scaffold and grey complementary polymer may each extend 8 bp
    beyond the duplex into Z2.  Their intersection is exactly the requested
    duplex.  The two polymers contain precisely 16 bp of combined overhang,
    which is the reviewed four-state 8-bp boundary cycle.
    """
    low, high = duplex_range
    candidates = []
    offsets = (0, 8, 16)
    for scaf_left in offsets:
        for scaf_right in offsets:
            for stap_left in offsets:
                for stap_right in offsets:
                    if (scaf_left + scaf_right + stap_left + stap_right != 16):
                        continue
                    scaf_range = (low-scaf_left, high+scaf_right)
                    stap_range = (low-stap_left, high+stap_right)
                    if min(scaf_range + stap_range) < 0:
                        continue
                    if (max(scaf_range[0], stap_range[0]) != low or
                            min(scaf_range[1], stap_range[1]) != high):
                        continue
                    scaf = build_field(source_rows, "scaf", (scaf_range,),
                                       array_length)
                    stap = build_field(source_rows, "stap", (stap_range,),
                                       array_length)
                    scaf_lengths = component_lengths_for_field(scaf, "scaf")
                    stap_lengths = component_lengths_for_field(stap, "stap")
                    if (scaf_lengths and stap_lengths and
                            set(scaf_lengths) == {32} and
                            set(stap_lengths) == {32}):
                        candidates.append((scaf_range, stap_range))
    if len(candidates) != 1:
        raise ValueError(
            "Square SST boundary is not uniquely resolvable for duplex "
            f"{duplex_range}: {candidates}")
    return candidates[0]


def resolve_complete_sst_ranges(source_rows, duplex, placement):
    """Apply only a phase-preserving 32*N canvas shift when base 0 requires it."""
    global_shift = 0
    while global_shift <= 4096:
        shifted = tuple((low+global_shift, high+global_shift)
                        for low, high in duplex)
        array_length = max(640, shifted[-1][1]+96)
        try:
            resolved = [legal_polymer_ranges(
                source_rows, item, array_length) for item in shifted]
        except ValueError:
            global_shift += PHASE_REPEAT_BP
            continue
        scaffold_ranges = tuple(item[0] for item in resolved)
        complement_ranges = tuple(item[1] for item in resolved)
        actual_placement = copy.deepcopy(placement)
        actual_placement["global_32bp_canvas_shift"] = global_shift
        actual_placement["reference_envelope"] = [
            value+global_shift for value in REFERENCE_ENVELOPE]
        actual_placement["target_envelope"] = [
            shifted[0][0], shifted[-1][1]]
        # Relative centring is unchanged because Seed and both SST layers are
        # translated together by exactly 32*N bp.
        return (shifted, scaffold_ranges, complement_ranges,
                actual_placement, global_shift)
    raise ValueError("unable to place complete Square SST on caDNAno canvas")


def seed_rows(seed, array_length, base_shift=0):
    rows = []
    for source in seed["vstrands"]:
        if int(source["num"]) not in SEED_SOURCE_HELICES:
            continue
        row = copy.deepcopy(source)
        row["num"] = int(row["num"]) + SEED_TARGET_SHIFT
        shifted_scaf = [EMPTY[:] for unused in range(array_length)]
        for source_base, source_record in enumerate(source.get("scaf", [])):
            target_base = source_base + base_shift
            if not (0 <= target_base < array_length):
                continue
            record = copy.deepcopy(source_record)
            for side in (0, 2):
                if record[side] >= 0:
                    record[side] += SEED_TARGET_SHIFT
                    record[side+1] += base_shift
            shifted_scaf[target_base] = record
        row["scaf"] = shifted_scaf
        row["stap"] = [EMPTY[:] for unused in range(array_length)]
        for field in ("loop", "skip"):
            shifted = [0]*array_length
            for source_base, value in enumerate(source.get(field, [])):
                target_base = source_base+base_shift
                if 0 <= target_base < array_length:
                    shifted[target_base] = value
            row[field] = shifted
        row["stap_colors"] = []
        rows.append(row)
    return rows


def audit_reciprocity(rows, field, helices):
    by_num = {int(row["num"]): row for row in rows}
    errors = []
    for helix in helices:
        records = by_num[helix][field]
        for base, record in enumerate(records):
            for side in (0, 2):
                partner, partner_base = record[side:side+2]
                if partner < 0:
                    continue
                if partner not in by_num or not (0 <= partner_base < len(
                        by_num[partner][field])):
                    errors.append([helix, base, partner, partner_base,
                                   "missing partner"])
                    continue
                peer = by_num[partner][field][partner_base]
                if [helix, base] not in (peer[:2], peer[2:]):
                    errors.append([helix, base, partner, partner_base,
                                   "not reciprocal"])
    return errors


def audit_phase(rows, source_rows, field, ranges):
    by_num = {int(row["num"]): row for row in rows}
    errors = []
    crossover_residues = set()
    for active_range in ranges:
        for base in range(active_range[0], active_range[1]+1):
            source_base = canonical_source_base(base)
            for helix in SST_TARGET_HELICES:
                expected = translated_record(
                    source_rows[helix+16][field][source_base],
                    source_base, base, active_range)
                actual = by_num[helix][field][base]
                if actual != expected:
                    errors.append([helix, base, actual, expected])
                if ((actual[0] >= 0 and actual[0] != helix) or
                        (actual[2] >= 0 and actual[2] != helix)):
                    crossover_residues.add(base % PHASE_REPEAT_BP)
    return errors, sorted(crossover_residues)


def audit_duplex(rows, duplex_ranges):
    by_num = {int(row["num"]): row for row in rows}
    missing = []
    for low, high in duplex_ranges:
        for helix in SST_TARGET_HELICES:
            for base in range(low, high+1):
                if (by_num[helix]["scaf"][base] == EMPTY or
                        by_num[helix]["stap"][base] == EMPTY):
                    missing.append([helix, base])
    return missing


def audit_component_lengths(rows, field):
    """Return component lengths for the 16 SST helices only."""
    by_num = {int(row["num"]): row for row in rows
              if int(row["num"]) in SST_TARGET_HELICES}
    nodes = {
        (helix, base)
        for helix, row in by_num.items()
        for base, record in enumerate(row[field]) if record != EMPTY}
    visited = set()
    lengths = []
    for start in sorted(nodes):
        if start in visited:
            continue
        stack = [start]
        component = set()
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            helix, base = node
            record = by_num[helix][field][base]
            for side in (0, 2):
                partner = tuple(record[side:side+2])
                if partner in nodes and partner not in component:
                    stack.append(partner)
        visited.update(component)
        lengths.append(len(component))
    return sorted(lengths)


def interval_intersection(first, second):
    low = max(first[0], second[0])
    high = min(first[1], second[1])
    return [low, high] if low <= high else None


def build_case(kind, first_length, spacing, second_length, linked,
               sst, seed):
    if linked and (first_length != second_length):
        raise ValueError("linked/identical SST layers must have equal lengths")
    if linked and ((first_length+spacing) % PHASE_REPEAT_BP):
        raise ValueError(
            "linked/identical SST layers require a 32-bp phase translation")
    source_rows = {int(row["num"]): row for row in sst["vstrands"]
                   if int(row["num"]) in SST_SOURCE_HELICES}
    requested_duplex, placement = desired_duplex_ranges(
        first_length, spacing, second_length)
    (duplex, scaffold_ranges, complement_ranges, placement,
     global_shift) = resolve_complete_sst_ranges(
        source_rows, requested_duplex, placement)
    highest = max(seed.get("num_bases", 0)+global_shift,
                  duplex[-1][1]+64, complement_ranges[-1][1]+32)
    array_length = max(640, ((highest+31)//32)*32)
    scaffold = build_field(source_rows, "scaf", scaffold_ranges,
                           array_length)
    complement = build_field(source_rows, "stap", complement_ranges,
                             array_length)
    seed_geometry = {int(row["num"]): row for row in seed["vstrands"]}
    output_rows = []
    for target_helix in SST_TARGET_HELICES:
        geometry = seed_geometry[target_helix+48]
        row = {
            "row": int(geometry["row"]),
            "col": int(geometry["col"]),
            "num": target_helix,
            "scaf": scaffold[target_helix],
            "stap": complement[target_helix],
            "loop": [0]*array_length,
            "skip": [0]*array_length,
            "scafLoop": [],
            "stapLoop": [],
            "stap_colors": [],
        }
        output_rows.append(row)
    output_rows.extend(seed_rows(seed, array_length, global_shift))
    output_rows.sort(key=lambda row: int(row["num"]))

    phase_scaf, scaf_residues = audit_phase(
        output_rows, source_rows, "scaf", scaffold_ranges)
    phase_stap, stap_residues = audit_phase(
        output_rows, source_rows, "stap", complement_ranges)
    reciprocal_scaf = audit_reciprocity(
        output_rows, "scaf", SST_TARGET_HELICES)
    reciprocal_stap = audit_reciprocity(
        output_rows, "stap", SST_TARGET_HELICES)
    missing_duplex = audit_duplex(output_rows, duplex)
    scaf_component_lengths = audit_component_lengths(output_rows, "scaf")
    stap_component_lengths = audit_component_lengths(output_rows, "stap")
    invalid_components = [
        ["scaf", length] for length in scaf_component_lengths
        if length != 32] + [
        ["stap", length] for length in stap_component_lengths
        if length != 32]
    spacing_range = (duplex[0][1]+1, duplex[1][0]-1)
    same_polymer_overlaps = []
    for field, ranges in (("scaf", scaffold_ranges),
                          ("stap", complement_ranges)):
        overlap = interval_intersection(ranges[0], ranges[1])
        if overlap is not None:
            same_polymer_overlaps.append([field, overlap])
    duplex_in_spacing = interval_intersection(
        (max(scaffold_ranges[0][0], complement_ranges[0][0]),
         min(scaffold_ranges[0][1], complement_ranges[0][1])),
        spacing_range)
    duplex_in_spacing_2 = interval_intersection(
        (max(scaffold_ranges[1][0], complement_ranges[1][0]),
         min(scaffold_ranges[1][1], complement_ranges[1][1])),
        spacing_range)
    linked_phase = ((duplex[1][0]-duplex[0][0]) % 32 == 0)
    if linked and not linked_phase:
        raise AssertionError("linked SST layers are not a 32-bp translation")
    all_errors = (phase_scaf + phase_stap + reciprocal_scaf +
                  reciprocal_stap + missing_duplex + invalid_components +
                  same_polymer_overlaps)
    if duplex_in_spacing is not None or duplex_in_spacing_2 is not None:
        all_errors.append(["unintended duplex inside Z2",
                           duplex_in_spacing, duplex_in_spacing_2])
    if all_errors:
        raise AssertionError("audit failed: %r" % all_errors[:10])

    filename = (f"{kind}_SST{first_length}_Z2_{spacing}_SST"
                f"{second_length}_complete_square_final.json")
    metadata = {
        "review_role": "Seed scaffold + complete SST scaffold/complementary chain",
        "helix_numbering": "SST 0-15; fixed Seed scaffold 16-63",
        "seed_reference": "Square_Seed_2L_newtemplate.json helices 0-47; scaffold topology unchanged",
        "sst_geometry_reference": "Square_Seed_2L_newtemplate.json helices 48-63",
        "sst_topology_reference": "Square_SST_original_128.json helices 16-31",
        "absolute_phase_anchor_base": PHASE_ANCHOR_BASE,
        "absolute_phase_repeat_bp": PHASE_REPEAT_BP,
        "SST_Z2_SST": [first_length, spacing, second_length],
        "linked": bool(linked),
        "sst_duplex_ranges": [list(item) for item in duplex],
        "sst_scaffold_ranges": [list(item) for item in scaffold_ranges],
        "sst_complementary_chain_ranges": [list(item) for item in complement_ranges],
        "spacing_range": list(spacing_range),
        "spacing_definition": "duplex-to-duplex coordinate gap; SST single-strand overhangs may occupy part of it",
        "square_boundary_rule": "complete reciprocal 32-nt U units; blue and grey single strands may enter Z2",
        "capture_gaps_reserved": False,
        "capture_connection_stage": "not applied; capture nicks are opened only in the later Seed-connection step",
        "single_strands_inside_Z2": {
            "first_SST_blue_scaffold": interval_intersection(
                scaffold_ranges[0], spacing_range),
            "first_SST_grey_complement": interval_intersection(
                complement_ranges[0], spacing_range),
            "second_SST_blue_scaffold": interval_intersection(
                scaffold_ranges[1], spacing_range),
            "second_SST_grey_complement": interval_intersection(
                complement_ranges[1], spacing_range),
        },
        "placement": placement,
        "audit": {
            "absolute_phase_errors": 0,
            "scaffold_reciprocity_errors": 0,
            "complement_reciprocity_errors": 0,
            "missing_duplex_bases": 0,
            "scaffold_32nt_U_component_count": len(scaf_component_lengths),
            "complement_32nt_U_component_count": len(stap_component_lengths),
            "non_32nt_component_count": 0,
            "same_polymer_layer_overlap_count": 0,
            "unintended_duplex_inside_Z2_count": 0,
            "linked_layer_translation_mod_32": (
                (duplex[1][0]-duplex[0][0]) % 32),
            "scaffold_crossover_base_residues_mod_32": scaf_residues,
            "complement_crossover_base_residues_mod_32": stap_residues,
            "passed": True,
        },
    }
    payload = {
        "name": filename,
        "vstrands": output_rows,
        "num_bases": array_length,
        "lattice": "square",
        "moire_structure_metadata": metadata,
    }
    return filename, payload, metadata


def main():
    sst = load(SST_REFERENCE)
    seed = load(SEED_REFERENCE)
    if OUTPUT_ROOT.exists():
        shutil.rmtree(OUTPUT_ROOT)
    OUTPUT_ROOT.mkdir(parents=True)
    index = []
    for case in CASES:
        filename, payload, metadata = build_case(*case, sst, seed)
        (OUTPUT_ROOT / filename).write_text(
            json.dumps(payload, separators=(",", ":")), encoding="utf-8")
        index.append({"file": filename, **metadata})
    (OUTPUT_ROOT / "review_index.json").write_text(
        json.dumps(index, ensure_ascii=False, indent=2), encoding="utf-8")
    archive = shutil.make_archive(str(OUTPUT_ROOT), "zip",
                                  root_dir=OUTPUT_ROOT)
    print(OUTPUT_ROOT)
    print(archive)
    print("files", len(index), "all_passed", all(
        item["audit"]["passed"] for item in index))


if __name__ == "__main__":
    main()
