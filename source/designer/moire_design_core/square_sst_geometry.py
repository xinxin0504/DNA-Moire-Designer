"""Shared placement geometry for complete Square SST superlattices.

This module is intentionally topology-free.  Both the caDNAno generator and
the 1B side-view preview consume the exact same centred layer/range result.
"""

from __future__ import annotations

from typing import Dict, Iterable, Tuple


DOMAIN_BP = 8
REPEAT_BP = 32
REFERENCE_ENVELOPE = (48, 335)
REFERENCE_SEED_LAYER_RANGES = ((48, 175), (208, 335))
REFERENCE_CAPTURE_POSITIONS_BY_LAYER = (
    (56, 72, 88, 104, 120, 136, 152, 168),
    (216, 232, 248, 264, 280, 296, 312, 328),
)
REFERENCE_CAPTURE_CONTACT_RANGE = (56, 328)


def _capture_columns_in_duplex(layer_range: Iterable[int]) -> Tuple[int, ...]:
    """Return absolute 16-bp capture phases inside one real duplex range.

    The fixed reference lists above describe only the canonical 128/32/128
    drawing.  They must never clip a moved SST layer: physical support starts
    at the current SST duplex start and ends at its duplex end.  A legal
    absolute phase outside that double-stranded interval is not connectable.
    """
    low, high = map(int, layer_range)
    first = low + ((8-low) % 16)
    return tuple(range(first, high+1, 16))


def seed_template_capture_columns(layer_range: Iterable[int],
                                  canvas_shift: int = 0
                                  ) -> Tuple[int, ...]:
    """Return real fixed-Seed contact columns inside one SST duplex.

    The accepted 2L Seed contains an unbroken 16-bp contact phase from base
    56 through 328.  Bases 184 and 200 are the two open boundary nicks and
    are therefore valid when a moved SST layer is genuinely duplex there.
    Coordinates beyond this translated template contact range are theoretical
    SST candidates only, not physical Seed captures.
    """
    low, high = map(int, layer_range)
    first = REFERENCE_CAPTURE_CONTACT_RANGE[0] + int(canvas_shift)
    last = REFERENCE_CAPTURE_CONTACT_RANGE[1] + int(canvas_shift)
    return tuple(position for position in range(first, last+1, 16)
                 if low <= position <= high)


def _validate(first_length: int, spacing: int, second_length: int) -> None:
    values = (int(first_length), int(spacing), int(second_length))
    if values[0] < 64 or values[2] < 64:
        raise ValueError("Square SST两层长度均须至少为64 bp。")
    if not 0 <= values[1] <= 160:
        raise ValueError("Square SST spacing必须位于0-160 bp。")
    if any(value % DOMAIN_BP for value in values):
        raise ValueError("Square SST长度与spacing必须采用8 bp步长。")


def _complete_u_polymer_ranges(layer_range: Iterable[int]):
    """Return blue/grey ranges without cutting a reciprocal 32-nt U.

    The four legal boundary states are selected by the absolute 8-bp phase.
    Consequently this function never restarts a local phase at either SST
    boundary.  The intersection of the two polymers is exactly the duplex.
    """
    low, high = map(int, layer_range)
    length = high-low+1
    if low % DOMAIN_BP or length % DOMAIN_BP:
        raise ValueError("Square SST边界必须位于8 bp网格。")
    low_phase = low % 16
    length_phase = length % 16
    if low_phase == 0 and length_phase == 0:
        scaffold = (low, high)
        complement = (low-8, high+8)
    elif low_phase == 0 and length_phase == 8:
        scaffold = (low, high+8)
        complement = (low-8, high)
    elif low_phase == 8 and length_phase == 0:
        scaffold = (low-8, high+8)
        complement = (low, high)
    elif low_phase == 8 and length_phase == 8:
        scaffold = (low-8, high)
        complement = (low, high+8)
    else:  # pragma: no cover - guarded by 8-bp validation
        raise ValueError("Square SST边界相位不可解析。")
    return scaffold, complement


def complete_square_polymer_ranges(layer_range: Iterable[int]):
    """Public adapter for one phase-correct complete-U Square layer."""
    scaffold, complement = _complete_u_polymer_ranges(layer_range)
    return list(scaffold), list(complement)


def refresh_seed_overlap_metadata(geometry: Dict[str, object]
                                  ) -> Dict[str, object]:
    """Recompute three-segment support/overlap after a rigid translation.

    Lattice-specific routers may apply one additional 8-bp phase translation
    after the shared placement is selected.  Keeping this calculation in one
    function prevents preview, generated JSON and scaffold review from
    reporting different Z1/Z3 overlaps.
    """
    first, second = [tuple(map(int, item)) for item in
                     geometry["layer_ranges"]]
    spacing_range = tuple(map(int, geometry["spacing_range"]))
    reference_low, reference_high = map(
        int, geometry["reference_envelope"])
    seed_partition_ranges = (
        (reference_low, spacing_range[0]-1),
        spacing_range,
        (spacing_range[1]+1, reference_high),
    )
    seed_partition_lengths = tuple(
        max(0, high-low+1) for low, high in seed_partition_ranges)
    if sum(seed_partition_lengths) != 288:
        raise ValueError(
            "当前SST位置无法在固定288-bp Seed内定义完整Z1/Z2/Z3分区。")
    overlap_ranges = []
    overlap_bp = []
    for layer, support in zip(
            (first, second),
            (seed_partition_ranges[0], seed_partition_ranges[2])):
        low = max(layer[0], support[0])
        high = min(layer[1], support[1])
        overlap_ranges.append([low, high])
        overlap_bp.append(max(0, high-low+1))
    geometry["seed_partition_ranges"] = [
        list(item) for item in seed_partition_ranges]
    geometry["seed_partition_lengths_bp"] = list(seed_partition_lengths)
    geometry["optimized_overlap_ranges"] = overlap_ranges
    geometry["optimized_seed_overlap_bp"] = overlap_bp
    geometry["optimized_seed_support_bp"] = [
        seed_partition_lengths[0], seed_partition_lengths[2]]
    geometry["maximin_overlap_bp"] = min(overlap_bp)
    geometry["total_overlap_bp"] = sum(overlap_bp)
    geometry["left_spare_domains"] = int(
        (first[0]-reference_low)//DOMAIN_BP)
    geometry["right_spare_domains"] = int(
        (reference_high-second[1])//DOMAIN_BP)
    geometry["envelope_center_offset_bp"] = (
        (first[0]+second[1])/2.0-
        (reference_low+reference_high)/2.0)
    geometry["spacing_center_offset_bp"] = (
        (spacing_range[0]+spacing_range[1])/2.0-
        (reference_low+reference_high)/2.0)
    return geometry


def _validate_template_phase(geometry: Dict[str, object]) -> None:
    """Reject a placement that no longer follows the reviewed template.

    Relative placement is allowed to move in 8-bp domain steps, but the SST
    graph itself is tiled on one immutable absolute 32-bp phase.  Only a
    common 32*N canvas translation may move that absolute phase together
    with the Seed.  This audit is deliberately kept beside the placement
    solver so a future centring/overlap change cannot silently generate a
    geometrically attractive but phase-incompatible design.
    """
    canvas_shift = int(geometry.get("coordinate_shift_bp", 0))
    if canvas_shift % REPEAT_BP:
        raise ValueError("Square SST公共画布平移必须保持32-bp范本相位。")
    layers = [tuple(map(int, item)) for item in geometry["layer_ranges"]]
    expected_scaffold = []
    expected_complement = []
    for layer in layers:
        scaffold, complement = _complete_u_polymer_ranges(layer)
        expected_scaffold.append(list(scaffold))
        expected_complement.append(list(complement))
    if expected_scaffold != geometry["scaffold_ranges"] or \
            expected_complement != geometry["complement_ranges"]:
        raise ValueError("Square SST边界与完整32-nt U型范本相位不一致。")
    origin = REFERENCE_CAPTURE_CONTACT_RANGE[0] + canvas_shift
    if int(geometry["capture_phase_reference_origin"]) != origin:
        raise ValueError("Square capture绝对相位原点发生偏移。")
    for layer, positions in zip(
            layers, geometry["seed_capture_positions_by_layer"]):
        expected = list(seed_template_capture_columns(layer, canvas_shift))
        if list(map(int, positions)) != expected:
            raise ValueError("Square capture列未保持固定Seed范本相位。")
        if any((int(position)-origin) % 16 for position in positions):
            raise ValueError("Square capture列不在固定16-bp接触相位上。")


def centered_square_sst_geometry(first_length: int, spacing: int,
                                 second_length: int) -> Dict[str, object]:
    """Place one rigid SST1/Z2/SST2 assembly against the fixed Seed.

    The assembly is translated only in integral 8-bp steps; its two SST
    lengths and their spacing never change.  Candidate translations are
    ranked lexicographically by the fixed-Seed three-segment policy:

    1. maximize the smaller of the two actual Seed/SST overlaps;
    2. maximize the total actual overlap, so support that a short SST cannot
       use is reassigned to the longer SST rather than wasted;
    3. among equal overlap solutions, centre the actual Z2/spacing interval
       in the 288-bp Seed as closely as the 8-bp grid permits; and
    4. use the lower coordinate as the deterministic half-domain tie-break.

    Thus Z2 is centred for equal or sufficiently long layers, but may shift
    when one SST is too short to use its nominal half of the available Seed.
    If the selected polymer position would cross base 0, Seed and both SST
    layers receive one common 32*N canvas translation; relative overlap and
    crossover phase are unchanged.
    """
    first_length, spacing, second_length = map(
        int, (first_length, spacing, second_length))
    _validate(first_length, spacing, second_length)
    reference_low, reference_high = REFERENCE_ENVELOPE
    reference_length = reference_high-reference_low+1
    target_length = first_length+spacing+second_length

    # Keep the full spacing interval inside the 288-bp Seed envelope so the
    # three displayed segments always sum exactly to 288 bp.  Every endpoint
    # is on the same 8-bp grid as the accepted reference.
    first_candidate = reference_low-first_length
    last_candidate = reference_high-first_length-spacing+1
    candidates = range(first_candidate, last_candidate+1, DOMAIN_BP)
    reference_center = (reference_low+reference_high)/2.0

    def candidate_metrics(start):
        first_end = start+first_length-1
        second_start = first_end+1+spacing
        z1_support = first_end-reference_low+1
        z3_support = reference_high-second_start+1
        overlap_first = min(first_length, z1_support)
        overlap_second = min(second_length, z3_support)
        envelope_center = start+(target_length-1)/2.0
        spacing_center = first_end + (spacing+1)/2.0
        score = (
            min(overlap_first, overlap_second),
            overlap_first+overlap_second,
            -abs(spacing_center-reference_center),
            -start,
        )
        return score, (z1_support, z3_support), (
            overlap_first, overlap_second), envelope_center-reference_center

    ranked = [(candidate_metrics(start), start) for start in candidates]
    if not ranked:  # pragma: no cover - spacing validation guarantees this
        raise ValueError("固定288-bp Seed内无法放置当前SST spacing。")
    (unused_score, unshifted_support, optimized_overlap,
     unshifted_center_offset), ideal_start = max(ranked)
    canvas_shift = 0
    # Every legal boundary state may contain an 8-bp single-strand overhang.
    while ideal_start+canvas_shift < DOMAIN_BP:
        canvas_shift += REPEAT_BP
    first_start = ideal_start+canvas_shift
    first = (first_start, first_start+first_length-1)
    second = (first[1]+1+spacing,
              first[1]+spacing+second_length)
    scaffold_ranges = []
    complement_ranges = []
    for layer in (first, second):
        scaffold, complement = _complete_u_polymer_ranges(layer)
        scaffold_ranges.append(scaffold)
        complement_ranges.append(complement)
    seed_ranges = tuple(
        (low+canvas_shift, high+canvas_shift)
        for low, high in REFERENCE_SEED_LAYER_RANGES)
    theoretical_capture_positions = tuple(
        _capture_columns_in_duplex(layer) for layer in (first, second))
    capture_positions = tuple(
        seed_template_capture_columns(layer, canvas_shift)
        for layer in (first, second))
    reference_envelope = tuple(
        value+canvas_shift for value in REFERENCE_ENVELOPE)
    spacing_range = (first[1]+1, second[0]-1)
    seed_partition_ranges = (
        (reference_envelope[0], spacing_range[0]-1),
        spacing_range,
        (spacing_range[1]+1, reference_envelope[1]),
    )
    seed_partition_lengths = tuple(
        max(0, high-low+1) for low, high in seed_partition_ranges)
    if sum(seed_partition_lengths) != 288:
        raise ValueError(
            "当前SST位置无法在固定288-bp Seed内定义完整Z1/Z2/Z3分区。")
    result = {
        "layer_ranges": [list(first), list(second)],
        # Capture support follows the current physical SST duplex.  Keep it
        # separate from the immutable Seed reference envelope so downstream
        # code cannot accidentally clip a moved layer with the canonical
        # 128/32/128 coordinates.
        "capture_support_ranges": [list(first), list(second)],
        "spacing_range": list(spacing_range),
        "scaffold_ranges": [list(item) for item in scaffold_ranges],
        "complement_ranges": [list(item) for item in complement_ranges],
        "seed_layer_ranges": [list(item) for item in seed_ranges],
        "seed_capture_positions_by_layer": [
            list(item) for item in capture_positions],
        "theoretical_capture_positions_by_layer": [
            list(item) for item in theoretical_capture_positions],
        "capture_phase_reference_origin": (
            REFERENCE_CAPTURE_CONTACT_RANGE[0] + canvas_shift),
        "reference_envelope": list(reference_envelope),
        "seed_partition_ranges": [
            list(item) for item in seed_partition_ranges],
        "seed_partition_lengths_bp": list(seed_partition_lengths),
        "target_envelope": [first[0], second[1]],
        "ideal_start_before_canvas_shift": ideal_start,
        "coordinate_shift_bp": canvas_shift,
        "left_spare_domains": int(
            (first_start-reference_envelope[0])//DOMAIN_BP),
        "right_spare_domains": int(
            (reference_envelope[1]-second[1])//DOMAIN_BP),
        "envelope_center_offset_bp": (
            (first[0]+second[1])/2.0-sum(reference_envelope)/2.0),
        "optimized_seed_overlap_bp": list(map(int, optimized_overlap)),
        "optimized_seed_support_bp": list(map(int, unshifted_support)),
        "maximin_overlap_bp": int(min(optimized_overlap)),
        "total_overlap_bp": int(sum(optimized_overlap)),
        "unshifted_envelope_center_offset_bp": float(
            unshifted_center_offset),
        "placement_step_bp": DOMAIN_BP,
        "placement_policy": (
            "rigid 8-bp translation; maximize minimum actual overlap, then "
            "total actual overlap, then centre Z2 among tied solutions; "
            "complete reciprocal 32-nt U on the immutable absolute template "
            "phase; only a shared 32*N canvas translation may avoid base 0; "
            "capture support is each current SST duplex range, not the "
            "canonical reference Seed interval"),
    }
    refresh_seed_overlap_metadata(result)
    _validate_template_phase(result)
    result["template_phase_preserved"] = True
    result["template_phase_period_bp"] = REPEAT_BP
    result["template_phase_anchor_bp"] = 64 + canvas_shift
    result["capture_phase_period_bp"] = 16
    return result
