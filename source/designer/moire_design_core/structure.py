"""Deterministic staged structure generation for the Square moire seed.

The public functions in this module deliberately keep caDNAno model imports
out of the GUI process.  The model worker runs in a short-lived subprocess,
then this module validates the legacy JSON before the UI exposes it.
"""

from __future__ import annotations

import json
import subprocess
import sys
from collections import defaultdict
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from moire_runtime import worker_command

from .kagome_sst import (
    ACTIVE_HELICES as KAGOME_ACTIVE_HELICES,
    LINE_HELICES as KAGOME_LINE_HELICES,
    build_kagome_layer_fields,
    build_kagome_sst_payload,
    kagome_layer_capture_catalogue,
    prepare_kagome_capture_sites,
    validate_kagome_sst_payload,
)
from .structure_rules import (
    KAGOME_CAPTURE_TOPOLOGY_REFERENCE,
    SEED_ROUTING_REFERENCE,
    SQUARE_CAPTURE_REFERENCE,
)
from .square_sst_geometry import (
    complete_square_polymer_ranges,
    centered_square_sst_geometry,
    refresh_seed_overlap_metadata,
    seed_template_capture_columns,
)
from .sst_auxiliary_routing import (
    AUXILIARY_INTERNAL,
    PRIMARY_INTERNAL,
    actual_helix as auxiliary_actual_helix,
    route_layer2_conflicts,
)


def _bundled_reference(filename):
    """Resolve a validated design resource bundled with this release."""
    return Path(__file__).with_name("resources") / filename


SST_REFERENCE = _bundled_reference("Square_SST_original_128.json")
SEED_CAPTURE_REFERENCE = _bundled_reference(
    "Square_Seed_2L_newtemplate.json")
# Keep these names separate even though Square currently uses the same 2L
# file for routing and capture.  In particular, the Kagome 3L file is an
# endpoint/topology catalogue only and must never become a Seed router.
SEED_ROUTING_REFERENCE = Path(SEED_ROUTING_REFERENCE)
SQUARE_CAPTURE_REFERENCE = Path(SQUARE_CAPTURE_REFERENCE)
KAGOME_CAPTURE_TOPOLOGY_REFERENCE = Path(
    KAGOME_CAPTURE_TOPOLOGY_REFERENCE)
SST_BASE_SHIFT = 32
SST_LAYER_RANGES = ((48, 175), (208, 335))
SST_STAPLE_RANGES = ((40, 183), (200, 343))
SST_ARRAY_LENGTH = 544
CAPTURE_SEED_HELICES = tuple(range(0, 8)) + tuple(range(24, 32))
CAPTURE_OUTPUT_HELICES = tuple(range(48, 64))
# The Seed is an immutable physical object.  SST lengths and SST spacing may
# change, but they may only change the geometric overlap with this accepted
# two-layer reference; they must never crop, extend or reroute the Seed.
FIXED_SEED_LAYER_RANGES = ((48, 175), (208, 335))
FIXED_SEED_CAPTURE_POSITIONS_BY_LAYER = (
    (56, 72, 88, 104, 120, 136, 152, 168),
    (216, 232, 248, 264, 280, 296, 312, 328),
)
FIXED_SEED_NOMINAL_SUPPORT_BP = (128, 128)
CAPTURE_FACE_DEFINITIONS = (
    {
        "id": "face1", "label": "Face 1 · upper Seed edge",
        "internal_seed_helices": tuple(range(0, 8)),
        "sst_first_seed_helices": tuple(range(16, 24)),
        "physical_internal_seed_helices": (0, 1, 2, 3),
        "physical_sst_first_seed_helices": (16, 17, 18, 19),
        "export_only_internal_seed_helices": (4, 5, 6, 7),
        "export_only_sst_first_seed_helices": (20, 21, 22, 23),
        "color": "#7b61b8",
    },
    {
        "id": "face2", "label": "Face 2 · lower Seed edge",
        "internal_seed_helices": tuple(range(24, 32)),
        "sst_first_seed_helices": tuple(range(40, 48)),
        "physical_internal_seed_helices": (31, 30, 29, 28),
        "physical_sst_first_seed_helices": (47, 46, 45, 44),
        "export_only_internal_seed_helices": (27, 26, 25, 24),
        "export_only_sst_first_seed_helices": (43, 42, 41, 40),
        "color": "#2a9d8f",
    },
)
CAPTURE_PHASE_CYCLE = ("A0", "B0")
CAPTURE_EXPORT_PHASE_CYCLE = ("A0", "B0", "A1", "B1")
# Each 32-nt U unit supplies one A and one B column.  Only the origin four
# helices on each face are connected to SST in the structure JSON.  The
# translated four helices are generated as alternative capture sequences at
# export time, so every face helix still receives one site per 32 bp without
# drawing eight simultaneous Seed–SST bridges at each 16-bp column.
CAPTURE_PHASE_MAPPINGS = {
    "A0": ((48, 31), (50, 29), (61, 2), (63, 0)),
    "B0": ((49, 30), (51, 28), (60, 3), (62, 1)),
    "A1": ((48, 27), (50, 25), (61, 6), (63, 4)),
    "B1": ((49, 26), (51, 24), (60, 7), (62, 5)),
}
CAPTURE_DIRECT_POSITIONS = (
    56, 72, 88, 104, 120, 136, 152, 168,
    216, 232, 248, 264, 280, 296, 312, 328,
)
CAPTURE_PAIR_COLORS = (
    0xff9600, 0x0433ff, 0xaa75c2, 0x00fdff,
    0x00f900, 0xffadff, 0x942192, 0xff2100,
)
# Capture ``pair`` is a physical/cooperative routing concept: two adjacent
# 16-bp columns can cooperate as one pair.  It must not be reused as the
# display/export identity.  The immutable 2L Seed has eighteen real-space
# Capture columns, and every column needs its own colour even when it is not
# reached by the current SST.  The two Z2 reserve columns (184 and 200) did
# not have distinct colours in the historical template, so two additional
# colours are assigned here while the other sixteen retain their reference
# appearance.
CAPTURE_TEMPLATE_COLUMNS = tuple(range(56, 329, 16))
CAPTURE_COLUMN_COLORS = (
    0x0433ff, 0xff9300, 0x00fdff, 0x942192,
    0x00f900, 0xff2600, 0xff40ff, 0xfffb00,
    0x6a5acd, 0x8b4513,
    0xe8af1a, 0xe64e68, 0x82ddbc, 0xdfd9a0,
    0xb853e4, 0xe32175, 0x51d9cb, 0x4563db,
)
# Historical reference colours identify every non-Z2 Capture column.  This
# map is used only to recover column identity from the frozen template; final
# JSON/XLSX colours always come from ``CAPTURE_COLUMN_COLORS`` above.
CAPTURE_REFERENCE_COLUMN_BY_COLOR = {
    0x0433ff: 56,
    0xff9300: 72,
    0x00fdff: 88,
    0x942192: 104,
    0x00f900: 120,
    0xff2600: 136,
    0xff40ff: 152,
    0xfffb00: 168,
    0xe8af1a: 216,
    0xe64e68: 232,
    0x82ddbc: 248,
    0xdfd9a0: 264,
    0xb853e4: 280,
    0xe32175: 296,
    0x51d9cb: 312,
    0x4563db: 328,
}
CAPTURE_EXTENSION_NT = 16
# The Square reference uses a 16-nt capture extension.  Other lattices may
# use 32 nt.  This describes the exported extension only; the Seed-side
# capture core is copied verbatim from the accepted template and has no
# runtime length policy.
NORMAL_STAPLE_MIN_NT = 21
NORMAL_STAPLE_MAX_NT = 58


def capture_extension_nt(layout: Dict[str, Any], layer: Optional[int] = None,
                         assignment: Optional[Dict[str, Any]] = None) -> int:
    """Return the capture extension length for one assignment.

    A future Kagome layout may set either ``capture_extension_nt`` globally
    or ``capture_extension_nt_by_layer``.  Existing Square JSON files omit
    both keys and therefore retain the validated 16-nt behavior.
    """
    if assignment is not None and assignment.get("capture_extension_nt"):
        value = assignment["capture_extension_nt"]
    else:
        layer = (assignment or {}).get("layer", layer)
        by_layer = layout.get("capture_extension_nt_by_layer")
        if by_layer and layer is not None and 0 < int(layer) <= len(by_layer):
            value = by_layer[int(layer) - 1]
        else:
            value = layout.get("capture_extension_nt", CAPTURE_EXTENSION_NT)
    value = int(value)
    if value not in (16, 32):
        raise ValueError("capture延伸长度目前只支持16 nt或32 nt。")
    return value


SCAFFOLD_CAPACITY_P8064 = 8064
SCAFFOLD_CAPACITY_ORTHOGONAL = 7557
SCAFFOLD_MAX_COUNT = 3


def scaffold_capacity_plan(total_nt: int) -> Dict[str, int]:
    """Return the hard Seed-scaffold allocation for ``total_nt`` bases.

    A single Seed scaffold may use P8064.  Multi-scaffold designs use the
    mutually orthogonal 7557-nt scaffold family, so two and three molecules
    have exact aggregate capacities of ``2 * 7557`` and ``3 * 7557``.
    The scaffold count is fixed by these thresholds; routing must balance the
    selected count and must never silently add another molecule merely to make
    a difficult partition easier.
    """
    total_nt = int(total_nt)
    if total_nt < 1:
        raise ValueError("Seed scaffold总长度必须大于0 nt。")
    if total_nt <= SCAFFOLD_CAPACITY_P8064:
        return {
            "count": 1,
            "per_scaffold_capacity_nt": SCAFFOLD_CAPACITY_P8064,
            "total_capacity_nt": SCAFFOLD_CAPACITY_P8064,
        }
    if total_nt <= 2 * SCAFFOLD_CAPACITY_ORTHOGONAL:
        return {
            "count": 2,
            "per_scaffold_capacity_nt": SCAFFOLD_CAPACITY_ORTHOGONAL,
            "total_capacity_nt": 2 * SCAFFOLD_CAPACITY_ORTHOGONAL,
        }
    if total_nt <= SCAFFOLD_MAX_COUNT * SCAFFOLD_CAPACITY_ORTHOGONAL:
        return {
            "count": 3,
            "per_scaffold_capacity_nt": SCAFFOLD_CAPACITY_ORTHOGONAL,
            "total_capacity_nt": (
                SCAFFOLD_MAX_COUNT * SCAFFOLD_CAPACITY_ORTHOGONAL),
        }
    raise ValueError(
        "Seed scaffold总长度%d nt超过3条正交scaffold的总容量%d nt；"
        "请减小Seed长度。" % (
            total_nt,
            SCAFFOLD_MAX_COUNT * SCAFFOLD_CAPACITY_ORTHOGONAL))
SQUARE_SCAF_LOW = (
    (4, 26, 15), (18, 28, 7), (10, 20, 31), (2, 12, 23))
SQUARE_SCAF_HIGH = (
    (5, 27, 16), (19, 29, 8), (11, 21, 0), (3, 13, 24))
SST_FIRST_HELICES = tuple(range(16))
SEED_AFTER_SST_HELICES = tuple(range(16, 64))
INTERNAL_TO_SST_FIRST = {
    **{number: number + 16 for number in range(48)},
    **{number: number - 48 for number in range(48, 64)},
}
SST_FIRST_TO_INTERNAL = {
    value: key for key, value in INTERNAL_TO_SST_FIRST.items()}


def fixed_seed_overlap_layout(
        sst_layer_ranges: Iterable[Iterable[int]],
        legal_positions_by_layer: Optional[
            Iterable[Iterable[int]]] = None,
        lattice_type: str = "square",
        seed_layer_ranges: Optional[Iterable[Iterable[int]]] = None,
        seed_capture_positions_by_layer: Optional[
            Iterable[Iterable[int]]] = None) -> Dict[str, Any]:
    """Intersect the current SST with the immutable two-layer Seed.

    This is the only Seed/SST overlap calculation used by the Moire
    Designer.  It is deliberately read-only: no Seed base, scaffold,
    staple, nick, seam or capture core is generated or resized here.
    """
    sst_ranges = [tuple(map(int, item)) for item in sst_layer_ranges]
    if len(sst_ranges) != 2:
        raise ValueError("SST必须恰好包含两层实际双链范围。")
    legal = (None if legal_positions_by_layer is None else
             [set(map(int, values)) for values in
              legal_positions_by_layer])
    if legal is not None and len(legal) != 2:
        raise ValueError("capture合法位点必须按两层提供。")
    # "Support" means the actual duplex start-to-end interval of each current
    # SST layer.  The old fixed [48..175]/[208..335] clipping was valid only
    # for the canonical 128/32/128 drawing and incorrectly rejected a legal
    # first duplex column after a centred SST movement (for example base 200
    # in the second layer of 96/56/128).  A caller may still supply a narrower
    # physical contact mask, but absence of one must never reinstate the old
    # canonical clipping.
    seed_ranges = tuple(tuple(map(int, item)) for item in
                        (seed_layer_ranges or sst_ranges))
    capture_grid = tuple(tuple(map(int, item)) for item in
                         (seed_capture_positions_by_layer or (
                             tuple(range(56, 329, 16)),
                             tuple(range(56, 329, 16)))))
    if len(seed_ranges) != 2 or len(capture_grid) != 2:
        raise ValueError("Seed支撑区与capture网格必须恰好包含两层。")
    overlap_ranges = []
    positions_by_layer = []
    overlap_bp = []
    for layer, ((sst_low, sst_high), (seed_low, seed_high),
                fixed_positions) in enumerate(zip(
                    sst_ranges, seed_ranges, capture_grid)):
        low, high = max(sst_low, seed_low), min(sst_high, seed_high)
        overlap_ranges.append([low, high])
        overlap_bp.append(max(0, high - low + 1))
        allowed = None if legal is None else legal[layer]
        positions_by_layer.append([
            position for position in fixed_positions
            if low <= position <= high and
            (allowed is None or position in allowed)])
    capture_columns = [len(values) for values in positions_by_layer]
    normalized_lattice = str(lattice_type).strip().lower().replace("-", "_")
    if normalized_lattice in ("square_kagome", "square_kagome_mixed"):
        lattice_by_layer = ("square", "kagome")
    elif normalized_lattice.startswith("kagome"):
        lattice_by_layer = ("kagome", "kagome")
    else:
        lattice_by_layer = ("square", "square")
    if "kagome" in lattice_by_layer:
        # The frozen Kagome catalogue alternates a B column (four physical
        # origin bridges plus four translated export bridges) with an A
        # column (two plus two): 8 + 4 captures form one cooperative pair.
        # Derive the phase from the immutable absolute 16-bp grid so cropping
        # an odd column from either edge does not accidentally swap 8 and 4.
        sites_per_column = [
            ([8 for unused in values] if layer_lattice == "square" else
             [8 if ((int(position) - 56) // 16) % 2 == 0 else 4
              for position in values])
            for values, layer_lattice in zip(
                positions_by_layer, lattice_by_layer)]
        count_semantics = (
            "Square columns contain eight final capture sites; Kagome "
            "columns follow the immutable 8+4 capture cycle; two adjacent "
            "columns form one cooperative capture pair")
    else:
        # Square has eight capture helices on every retained column.  Two
        # adjacent 8-site columns form one cooperative pair.
        sites_per_column = [[8 for unused in values]
                            for values in positions_by_layer]
        count_semantics = (
            "Square columns each contain eight final capture sites; two "
            "adjacent columns form one cooperative capture pair")
    # ``pair_counts`` is the number of palette/routing groups required by
    # the implementation, so a final singleton column still needs one group.
    # It is not the physical number of complete column pairs.  Preserve it
    # for routing compatibility and expose the exact half-integer equivalent
    # separately for UI/reporting (7 columns = 3.5 pair-equivalents).
    pair_counts = [(count + 1) // 2 for count in capture_columns]
    pair_equivalents = [count / 2.0 for count in capture_columns]
    return {
        "lattice_type": str(lattice_type),
        "seed_geometry_policy": "immutable_2L_reference",
        "seed_length_adjustment_enabled": False,
        "seed_layer_ranges": [list(item)
                              for item in seed_ranges],
        "sst_layer_ranges": [list(item) for item in sst_ranges],
        "overlap_ranges": overlap_ranges,
        "seed_sst_overlap_bp": overlap_bp,
        "capture_positions_by_layer": positions_by_layer,
        "actual_capture_positions_by_layer": positions_by_layer,
        "capture_positions": [position
                              for values in positions_by_layer
                              for position in values],
        "actual_capture_positions": [position
                                     for values in positions_by_layer
                                     for position in values],
        "capture_columns_by_layer": capture_columns,
        "capture_sites_per_column_by_layer": sites_per_column,
        "capture_sites_by_layer": [sum(values)
                                   for values in sites_per_column],
        "capture_pair_equivalents_by_layer": pair_equivalents,
        "capture_pair_equivalents": sum(pair_equivalents),
        "pair_count_by_layer": pair_counts,
        "pair_count": sum(pair_counts),
        "minimum_capture_pairs_per_layer": 2,
        "capture_pairs_valid": all(count >= 2 for count in pair_counts),
        "capture_count_semantics": (
            "actual physical captures are template-valid columns inside "
            "each current SST layer's real duplex start-to-end support; "
            "single-stranded boundary coordinates are excluded; " +
            count_semantics),
    }


def structure_layout(z1_bp: int = 128, z2_bp: int = 32,
                     z3_bp: int = 128, seed_z1_bp: Optional[int] = None,
                     seed_z3_bp: Optional[int] = None,
                     capture_extension_length_nt: int = 16) -> Dict[str, Any]:
    """Return dynamic SST ranges against the immutable two-layer Seed.

    ``seed_z1_bp`` and ``seed_z3_bp`` are accepted only so old project files
    continue to load.  Their values no longer affect geometry.  The physical
    Seed, its scaffold, ordinary staples, nicks, seams and capture candidates
    always come from ``Square_Seed_2L_newtemplate.json``.  Only the real-space
    intersection with the current SST is recomputed.
    """
    z1_bp, z2_bp, z3_bp = map(int, (z1_bp, z2_bp, z3_bp))
    ignored_seed_z1_bp = None if seed_z1_bp is None else int(seed_z1_bp)
    ignored_seed_z3_bp = None if seed_z3_bp is None else int(seed_z3_bp)
    seed_z1_bp, seed_z3_bp = FIXED_SEED_NOMINAL_SUPPORT_BP
    capture_extension_length_nt = int(capture_extension_length_nt)
    # Capture cores are immutable template topology.  Their length is not a
    # design input and is never re-cut, rejected or repaired by Moire
    # Designer.  The extension value is retained only for sequence export.
    if capture_extension_length_nt not in (16, 32):
        raise ValueError("capture延伸长度目前只支持16 nt或32 nt。")
    if z1_bp < 64 or z3_bp < 64:
        raise ValueError("SST 1st layer和2nd layer至少需要64 bp。")
    if any(value % 8 for value in (z1_bp, z2_bp, z3_bp)):
        raise ValueError("SST 1st layer、spacing和2nd layer必须是8 bp整数倍。")
    geometry = centered_square_sst_geometry(z1_bp, z2_bp, z3_bp)
    coordinate_shift = int(geometry["coordinate_shift_bp"])
    layer_ranges = tuple(tuple(item) for item in geometry["layer_ranges"])
    spacing_range = tuple(geometry["spacing_range"])
    scaffold_ranges = tuple(tuple(item)
                            for item in geometry["scaffold_ranges"])
    staple_ranges = tuple(tuple(item)
                          for item in geometry["complement_ranges"])
    seed_layer_ranges = tuple(tuple(item)
                              for item in geometry["seed_layer_ranges"])
    seed_capture_positions = tuple(tuple(item) for item in
                                   geometry[
                                       "seed_capture_positions_by_layer"])
    # Capture is possible only where the current SST duplex intersects the
    # current Seed Z1/Z3 partitions.  ``seed_layer_ranges`` describes the
    # frozen reference's nominal 128-bp supports and is therefore not the
    # right mask when Z2 is shorter/longer than 32 bp.  In particular, a
    # 16-bp Z2 makes base 200 the first real Z3 support coordinate.
    seed_partition_support = (
        tuple(map(int, geometry["seed_partition_ranges"][0])),
        tuple(map(int, geometry["seed_partition_ranges"][2])),
    )
    overlap = fixed_seed_overlap_layout(
        layer_ranges, seed_layer_ranges=seed_partition_support,
        seed_capture_positions_by_layer=seed_capture_positions)
    overlap_ranges = tuple(map(tuple, geometry[
        "optimized_overlap_ranges"]))
    capture_support_ranges = tuple(
        map(tuple, overlap["overlap_ranges"]))
    # Capture candidates are properties of the fixed Seed.  An SST may use a
    # candidate only when that exact base remains inside the actual duplex
    # overlap.  No new capture coordinate is extrapolated when SST moves.
    # A/B phase belongs to the immutable Seed template, not to whichever
    # column happens to be first after an SST move.  Base 56 is B and base 72
    # is A; a shared canvas translation preserves that absolute cycle.
    capture_grid_origin = int(geometry.get(
        "capture_phase_reference_origin",
        seed_capture_positions[0][0] if seed_capture_positions and
        seed_capture_positions[0] else 56 + coordinate_shift))
    capture_positions_by_layer = overlap["capture_positions_by_layer"]
    # Only fixed-template capture columns inside the real duplex overlap are
    # active.  The Seed is never extended or rerouted to create another one.
    omitted_edge_capture_positions = []
    capture_positions = tuple(
        position for positions in capture_positions_by_layer
        for position in positions)
    pair_count_by_layer = overlap["pair_count_by_layer"]
    if any(count < 2 for count in pair_count_by_layer):
        raise ValueError(
            "每层至少需要2列capture pair；当前为%s。固定2L Seed"
            "不会缩进、增长或改变routing，请减小SST spacing。" %
            "/".join(map(str, pair_count_by_layer)))
    # Array allocation may grow to contain a long SST, but all Seed records
    # remain the untouched two-layer reference topology.
    seed_capture_helix_range = (
        seed_layer_ranges[0][0], seed_layer_ranges[1][1])
    minimum_length = max(layer_ranges[-1][1],
                         seed_capture_helix_range[1]) + 64
    array_length = max(
        SST_ARRAY_LENGTH,
        ((minimum_length + 31) // 32) * 32)
    result = {
        "z1_bp": z1_bp,
        "z2_bp": z2_bp,
        "z3_bp": z3_bp,
        "sst_z1_bp": z1_bp,
        "sst_spacing_bp": z2_bp,
        "sst_z3_bp": z3_bp,
        "seed_z1_requested_bp": seed_z1_bp,
        "seed_z3_requested_bp": seed_z3_bp,
        "ignored_legacy_seed_z1_input_bp": ignored_seed_z1_bp,
        "ignored_legacy_seed_z3_input_bp": ignored_seed_z3_bp,
        "seed_z1_actual_bp": FIXED_SEED_NOMINAL_SUPPORT_BP[0],
        "seed_z3_actual_bp": FIXED_SEED_NOMINAL_SUPPORT_BP[1],
        "seed_geometry_policy": "immutable_2L_reference",
        "seed_length_adjustment_enabled": False,
        "coordinate_shift_bp": coordinate_shift,
        "square_centered_geometry": geometry,
        "layer_ranges": [list(item) for item in layer_ranges],
        "spacing_range": list(spacing_range),
        "seed_z2_range": list(spacing_range),
        "spacing_seed_z2_coincident": True,
        "seed_partition_ranges": deepcopy(
            geometry["seed_partition_ranges"]),
        "seed_partition_lengths_bp": list(map(
            int, geometry["seed_partition_lengths_bp"])),
        "staple_ranges": [list(item) for item in staple_ranges],
        "sst_scaffold_ranges": [list(item) for item in scaffold_ranges],
        "sst_complementary_chain_ranges": [
            list(item) for item in staple_ranges],
        "seed_layer_ranges": [list(item) for item in seed_layer_ranges],
        "capture_support_ranges": [list(item)
                                   for item in capture_support_ranges],
        "overlap_ranges": [list(item) for item in overlap_ranges],
        "seed_sst_overlap_bp": list(
            geometry["optimized_seed_overlap_bp"]),
        "full_coverage": [
            max(0, high - low + 1) == sst_high - sst_low + 1
            for (low, high), (sst_low, sst_high) in zip(
                overlap_ranges, layer_ranges)],
        "seed_capture_helix_range": list(seed_capture_helix_range),
        "capture_positions_by_layer": capture_positions_by_layer,
        "actual_capture_positions_by_layer": capture_positions_by_layer,
        "capture_grid_origin": capture_grid_origin,
        "capture_positions": list(capture_positions),
        "actual_capture_positions": list(capture_positions),
        "omitted_edge_capture_positions":
            list(omitted_edge_capture_positions),
        "expected_capture_bridges": 4 * len(capture_positions),
        "expected_capture_export_sequences": 8 * len(capture_positions),
        "capture_column_count": len(capture_positions),
        "actual_capture_column_count": len(capture_positions),
        "capture_columns_by_layer": overlap["capture_columns_by_layer"],
        "capture_sites_per_column_by_layer": overlap[
            "capture_sites_per_column_by_layer"],
        "capture_sites_by_layer": overlap["capture_sites_by_layer"],
        "capture_pair_equivalents_by_layer": overlap[
            "capture_pair_equivalents_by_layer"],
        "capture_pair_equivalents": overlap[
            "capture_pair_equivalents"],
        "capture_count_semantics": (
            "actual fixed-Seed/SST-duplex intersection; no fixed count is "
            "inferred from SST length; each column contains eight final "
            "capture sites and equals one-half cooperative pair"),
        "capture_extension_nt": capture_extension_length_nt,
        "pair_count_by_layer": pair_count_by_layer,
        "pair_count": sum(pair_count_by_layer),
        "minimum_capture_pairs_per_layer": 2,
        "capture_pair_count_is_hard_constraint": True,
        "outer_capture_policy": (
            "fixed 2L Seed capture candidates intersected with the actual "
            "SST duplex; Seed topology is never changed"),
        "array_length": array_length,
    }
    result["capture_site_assignments"] = capture_site_assignments(result)
    result["capture_export_site_assignments"] = \
        capture_export_site_assignments(result)
    return result


def _layout_from_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    metadata = payload.get("moire_structure_metadata", {})
    configured = metadata.get("variable_length_layout")
    if configured:
        if configured.get("lattice_type") in (
                "kagome", "square_kagome", "square_kagome_mixed"):
            # Kagome and mixed Square--Kagome designs use phase-shifted
            # scaffold/staple windows and a layer-specific SST-side capture
            # catalogue.  Do not rebuild either through the all-Square
            # layout constructor: doing so silently changes the second-layer
            # phase and the expected bridge count after a saved JSON is
            # reopened.
            layout = deepcopy(configured)
            metadata_assignments = metadata.get(
                "kagome_capture_anchor_assignments_sst_only", [])
            if metadata_assignments and configured.get(
                    "lattice_type") in ("kagome", "square_kagome",
                                        "square_kagome_mixed"):
                layout["kagome_capture_anchor_assignments"] = deepcopy(
                    metadata_assignments)
            return layout
        rebuilt = structure_layout(
            configured.get("z1_bp", 128), configured.get("z2_bp", 32),
            configured.get("z3_bp", 128),
            configured.get("seed_z1_requested_bp"),
            configured.get("seed_z3_requested_bp"),
            configured.get("capture_extension_nt", CAPTURE_EXTENSION_NT))
        # Recompute every overlap/capture assignment from the immutable Seed.
        # Do not revive retired Seed growth, crop, edge-seam or AutoCS fields
        # merely because they exist in an older project JSON.
        for key in (
                "seed_cross_section_preset", "lattice_type",
                "layers_design_sequence_identical",
                "capture_grid_origins_by_layer",
                "sst_layer_fixture_lengths_bp",
                "sst_layer_translation_bp",
                "seed_routing_is_frozen_reference",
                "seed_routing_reference", "mean_indel_per_helix",
                "mean_indel_per_helix_requested",
                "mean_indel_per_helix_actual",
                "actual_z2_spacing_bp",
                "seed_z2_indel_range",
                "seed_z2_indel_placements",
                "seed_z2_indel_distribution",
                "auxiliary_sst_routing"):
            if key in configured:
                rebuilt[key] = configured[key]
        rebuilt["capture_site_assignments"] = capture_site_assignments(
            rebuilt)
        rebuilt["capture_export_site_assignments"] = \
            capture_export_site_assignments(rebuilt)
        rebuilt["expected_capture_bridges"] = sum(
            len(item["bridges"])
            for item in rebuilt["capture_site_assignments"])
        rebuilt["expected_capture_export_sequences"] = sum(
            len(item["bridges"])
            for item in rebuilt["capture_export_site_assignments"])
        return rebuilt
    return structure_layout()


def capture_pair_index(base_index: int, layout: Dict[str, Any]):
    """Return the pair index for a legal capture column, including a singleton."""
    offset = 0
    by_layer = layout.get("capture_positions_by_layer")
    if by_layer is None:
        by_layer = [
            [position for position in layout.get("capture_positions", [])
             if low <= position <= high + 8]
            for low, high in layout["layer_ranges"]]
    for positions in by_layer:
        if not positions:
            continue
        nearest = min(
            range(len(positions)),
            key=lambda index: abs(int(positions[index]) - int(base_index)))
        if abs(int(positions[nearest]) - int(base_index)) <= 8:
            return offset + nearest // 2
        offset += (len(positions) + 1) // 2
    return None


def capture_column_index(base_index: int, layout: Dict[str, Any]):
    """Return the immutable Seed Capture-column index for an exact base.

    A long SST can translate the complete frozen Seed by a whole 32-bp
    canvas period.  Column identity follows that translation, but it is not
    affected by which columns currently overlap the SST or by cooperative
    pair grouping.
    """
    shift = int(layout.get("coordinate_shift_bp", 0) or 0)
    logical_base = int(base_index) - shift
    try:
        return CAPTURE_TEMPLATE_COLUMNS.index(logical_base)
    except ValueError:
        return None


def capture_column_color(base_index: int, layout: Dict[str, Any]):
    """Return the unique RGB integer assigned to one Capture column."""
    index = capture_column_index(base_index, layout)
    if index is None:
        return None
    return CAPTURE_COLUMN_COLORS[index]


def _load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sst_reference_rows(resource: Optional[Dict[str, Any]] = None):
    """Return the validated 4x4 SST block normalized to helix 0--15.

    The original archived design stores Seed on 0--15 and the reusable SST
    on 16--31.  A cropped standalone resource stores that same SST block on
    0--15.  Prefer the archived 16--31 block when present, then normalize its
    helix ids so the remaining tiling code has one invariant input format.
    """
    resource = _load_json(SST_REFERENCE) if resource is None else resource
    by_number = {
        int(row["num"]): row for row in resource.get("vstrands", [])}
    source_start = 16 if all(number in by_number
                             for number in range(16, 32)) else 0
    source_numbers = range(source_start, source_start + 16)
    missing = [number for number in source_numbers if number not in by_number]
    if missing:
        raise ValueError(
            "SST模板缺少helix：" + ", ".join(map(str, missing)))
    rows = [deepcopy(by_number[number]) for number in source_numbers]
    if source_start:
        for row in rows:
            row["num"] = int(row["num"]) - source_start
            for field in ("scaf", "stap"):
                for record in row.get(field, []):
                    for offset in (0, 2):
                        partner = int(record[offset])
                        if source_start <= partner < source_start + 16:
                            record[offset] = partner - source_start
    return rows


def _remap_payload_helices(payload, mapping, numbering):
    """Return a deep copy with row, crossover and sequence helix ids remapped."""
    result = deepcopy(payload)
    for row in result.get("vstrands", []):
        old_number = int(row["num"])
        if old_number not in mapping:
            continue
        row["num"] = mapping[old_number]
        for field in ("scaf", "stap"):
            for record in row.get(field, []):
                for offset in (0, 2):
                    partner = int(record[offset])
                    if partner in mapping:
                        record[offset] = mapping[partner]
        row["stap_colors"] = [list(item)
                              for item in row.get("stap_colors", [])]
    result["vstrands"] = sorted(
        result.get("vstrands", []), key=lambda row: int(row["num"]))
    for entry in result.get("scaffold_sequences", []):
        start = int(entry.get("start_vh", -1))
        if start in mapping:
            entry["start_vh"] = mapping[start]
    metadata = result.setdefault("moire_structure_metadata", {})
    metadata["helix_numbering"] = numbering
    routing = metadata.get("auxiliary_sst_routing", {})
    auxiliary = (list(range(64, 80)) if isinstance(routing, dict) and
                 routing.get("enabled") else [])
    if numbering == "sst_first":
        metadata.update({
            "sst_helix_numbers": list(SST_FIRST_HELICES) + auxiliary,
            "seed_helix_numbers": list(SEED_AFTER_SST_HELICES),
            "numbering_policy": (
                "SST 0-15; Seed 16-63; auxiliary SST 64-79 when present"),
        })
    else:
        metadata.update({
            "sst_helix_numbers": list(CAPTURE_OUTPUT_HELICES) + auxiliary,
            "seed_helix_numbers": list(range(48)),
            "numbering_policy": (
                "internal Seed 0-47; SST 48-63; auxiliary SST 64-79"),
        })
    return result


def payload_to_internal_numbering(payload):
    """Normalize a public SST-first payload for the legacy routing engine."""
    metadata = payload.get("moire_structure_metadata", {})
    rows = {int(row["num"]) for row in payload.get("vstrands", [])}
    if metadata.get("helix_numbering") == "sst_first":
        mapping = dict(SST_FIRST_TO_INTERNAL)
        mapping.update({number: number for number in range(64, 80)})
        return _remap_payload_helices(
            payload, mapping, "internal_seed_first")
    if rows == set(SST_FIRST_HELICES) and metadata.get("stage", "").startswith(
            "sst_"):
        return _remap_payload_helices(
            payload, {number: number + 48 for number in SST_FIRST_HELICES},
            "internal_seed_first")
    return deepcopy(payload)


def payload_to_sst_first_numbering(payload):
    """Expose stable SST ids 0-15, followed by Seed ids 16-63."""
    metadata = payload.get("moire_structure_metadata", {})
    if metadata.get("helix_numbering") == "sst_first":
        return deepcopy(payload)
    rows = {int(row["num"]) for row in payload.get("vstrands", [])}
    if any(number < 48 for number in rows):
        mapping = dict(INTERNAL_TO_SST_FIRST)
        mapping.update({number: number for number in range(64, 80)})
    else:
        mapping = {number: number - 48 for number in CAPTURE_OUTPUT_HELICES}
        mapping.update({number: number for number in range(64, 80)})
    return _remap_payload_helices(payload, mapping, "sst_first")


def _empty_record() -> List[int]:
    return [-1, -1, -1, -1]


def _translate_record(record: Iterable[int], helix_shift: int,
                      base_shift: int) -> List[int]:
    values = list(record)
    translated: List[int] = []
    for offset in (0, 2):
        helix, base = int(values[offset]), int(values[offset + 1])
        translated.extend(
            (-1, -1) if helix < 0 else
            (helix + helix_shift, base + base_shift))
    return translated


def _periodic_source_index(source_low, source_high, offset, target_length,
                           interior_start=32, repeat_bp=32):
    """Use native ends and a validated interior repeat.

    Scaffold crossover phases repeat every 32 bases, but the accepted Seed
    staple/nick pattern is a complete 128-base routing program.  Treating the
    staple program as a 32-base repeat creates a nick at every capture column
    and leaves an internal 16-base Seed core.  Callers therefore select the
    repeat appropriate to the field they are copying.
    """
    if offset < 8:
        return source_low + offset
    if offset >= target_length - 8:
        return source_high - (target_length - 1 - offset)
    return (source_low + int(interior_start) +
            ((offset - 8) % int(repeat_bp)))


def _copy_periodic_segment(source_records, source_range, target_records,
                           target_range, helix_shift=0,
                           allowed_external=None,
                           minimum_internal_helix=48,
                           interior_start=32, repeat_bp=32):
    source_low, source_high = source_range
    target_low, target_high = target_range
    target_length = target_high - target_low + 1
    for offset, target_index in enumerate(range(target_low, target_high + 1)):
        source_index = _periodic_source_index(
            source_low, source_high, offset, target_length,
            interior_start=interior_start, repeat_bp=repeat_bp)
        record = list(source_records[source_index])
        translated = []
        for side in (0, 2):
            partner, partner_base = map(int, record[side:side + 2])
            if partner < 0:
                translated.extend((-1, -1))
                continue
            partner_target = target_index + (partner_base - source_index)
            partner_shifted = partner + helix_shift
            internal = (partner_shifted >= minimum_internal_helix and
                        source_low <= partner_base <= source_high and
                        target_low <= partner_target <= target_high)
            external = (allowed_external is not None and
                        allowed_external(partner_shifted, partner_target))
            translated.extend(
                (partner_shifted, partner_target)
                if internal or external else (-1, -1))
        target_records[target_index] = translated


def _sst_active_ranges(layout):
    """Return the shared centred complete-U Square polymer ranges."""
    scaffold = layout.get("sst_scaffold_ranges")
    complement = layout.get("sst_complementary_chain_ranges")
    if scaffold and complement:
        return deepcopy(scaffold), deepcopy(complement)
    geometry = centered_square_sst_geometry(
        int(layout["z1_bp"]), int(layout["z2_bp"]),
        int(layout["z3_bp"]))
    return (deepcopy(geometry["scaffold_ranges"]),
            deepcopy(geometry["complement_ranges"]))


def capture_site_assignments(layout, translation="origin"):
    """Return physical or translated Seed/SST assignments.

    Physical structure files use only ``origin`` (A0/B0): face1 helices 0–3
    and face2 helices 31–28.  ``translated`` (A1/B1) is reserved for sequence
    export and supplies the other four helices on each face.
    """
    if translation not in ("origin", "translated"):
        raise ValueError("capture translation必须是origin或translated。")
    unused_scaffold, staple_ranges = _sst_active_ranges(layout)
    assignments = []
    face_by_helix = {
        int(number): face["id"]
        for face in CAPTURE_FACE_DEFINITIONS
        for number in face["internal_seed_helices"]}
    if layout.get("lattice_type") in (
            "square_kagome", "square_kagome_mixed"):
        square_layout = deepcopy(layout)
        square_layout["lattice_type"] = "square"
        square_layout["capture_positions_by_layer"] = [
            list(layout.get("capture_positions_by_layer", [[], []])[0]),
            [],
        ]
        kagome_layout = deepcopy(layout)
        kagome_layout["lattice_type"] = "kagome"
        kagome_layout["capture_positions_by_layer"] = [
            [], list(layout.get("capture_positions_by_layer", [[], []])[1]),
        ]
        return sorted(
            capture_site_assignments(square_layout, translation) +
            capture_site_assignments(kagome_layout, translation),
            key=lambda item: (int(item["layer"]),
                              int(item["position"])))
    if layout.get("lattice_type") == "kagome":
        seed_by_sst_origin = {
            48: 31, 50: 29,
            49: 30, 51: 28, 60: 3, 62: 1,
        }
        translated_seed = {31: 27, 29: 25, 30: 26, 28: 24,
                           3: 7, 1: 5}
        grouped = defaultdict(list)
        allowed_positions = {
            (layer_index + 1, int(position))
            for layer_index, positions in enumerate(
                layout.get("capture_positions_by_layer", []))
            for position in positions}
        for item in layout.get("kagome_capture_anchor_assignments", []):
            if (int(item["layer"]), int(item["position"])) not in \
                    allowed_positions:
                continue
            grouped[(int(item["layer"]), int(item["position"]),
                     int(item["capture_extension_nt"]))].append(item)
        for (layer, position, extension_nt), items in sorted(grouped.items()):
            bridges = []
            for item in items:
                logical_sst_helix = int(item.get(
                    "logical_sst_helix", item["sst_helix"]))
                sst_helix = auxiliary_actual_helix(
                    layout, layer, "stap", logical_sst_helix, position)
                seed_helix = seed_by_sst_origin[logical_sst_helix]
                if translation == "translated":
                    seed_helix = translated_seed[seed_helix]
                bridges.append({
                    "sst_helix": sst_helix,
                    "logical_sst_helix": logical_sst_helix,
                    "seed_helix": seed_helix,
                    "sst_slot": int(item["slot"]),
                    "face": face_by_helix[seed_helix],
                })
            assignments.append({
                "layer": layer,
                "position": position,
                "phase": "K",
                "translation": translation,
                "export_only": translation == "translated",
                "cycle": "K0" if translation == "origin" else "K1",
                "capture_extension_nt": extension_nt,
                "bridges": bridges,
            })
        return assignments
    origins = layout.get("capture_grid_origins_by_layer")
    for layer_index, positions in enumerate(
            layout["capture_positions_by_layer"]):
        capture_grid_origin = int(
            origins[layer_index] if origins and layer_index < len(origins)
            else layout.get(
                "capture_grid_origin", CAPTURE_DIRECT_POSITIONS[0]))
        for position in positions:
            # The validated reference starts with B at base 56 and A at 72.
            # Determine phase from that fixed Seed grid; using the variable
            # SST staple range here re-phased capture columns onto AutoCS
            # crossover bases for non-128-bp supports.
            unit_index = (int(position) - capture_grid_origin) // 16
            phase = "B" if unit_index % 2 == 0 else "A"
            cycle = phase + ("0" if translation == "origin" else "1")
            assignments.append({
                "layer": layer_index + 1,
                "position": int(position),
                "phase": cycle[0],
                "translation": translation,
                "export_only": translation == "translated",
                "cycle": cycle,
                "capture_extension_nt": capture_extension_nt(
                    layout, layer=layer_index + 1),
                "bridges": [
                    {
                        "sst_helix": auxiliary_actual_helix(
                            layout, layer_index + 1, "stap",
                            int(sst_helix), int(position)),
                        "logical_sst_helix": int(sst_helix),
                        "seed_helix": int(seed_helix),
                        "face": face_by_helix[int(seed_helix)],
                    }
                    for sst_helix, seed_helix in
                    CAPTURE_PHASE_MAPPINGS[cycle]
                ],
            })
    return assignments


def capture_export_site_assignments(layout):
    """Return origin plus translated sites required in sequence exports."""
    assignments = (
        capture_site_assignments(layout, "origin") +
        capture_site_assignments(layout, "translated"))
    return sorted(assignments, key=lambda item: (
        int(item["layer"]), int(item["position"]),
        bool(item["export_only"])))


def _tile_complete_sst_field(resource_rows, field, target_ranges,
                             array_length):
    """Copy the reviewed U topology using one global absolute 32-bp phase."""
    source_by_number = {int(row["num"]): row for row in resource_rows}
    output = {
        int(source["num"]) + 48:
        [_empty_record() for unused in range(array_length)]
        for source in resource_rows}
    for low, high in target_ranges:
        if low < 0 or high >= array_length:
            raise ValueError("SST活动区间超出caDNAno坐标范围。")
        for target_index in range(low, high + 1):
            source_index = 64 + ((target_index - 64) % 32)
            for source_number, source in source_by_number.items():
                target_number = source_number + 48
                source_record = source[field][source_index]
                translated = []
                for offset in (0, 2):
                    partner, partner_base = map(
                        int, source_record[offset:offset + 2])
                    if partner < 0:
                        translated.extend((-1, -1))
                        continue
                    partner_target = target_index + (
                        partner_base-source_index)
                    if not (low <= partner_target <= high):
                        translated.extend((-1, -1))
                    else:
                        translated.extend((partner+48, partner_target))
                output[target_number][target_index] = translated
    return output


def _open_capture_sites(rows, layout):
    """Open only the U crossovers later replaced by Seed capture bridges."""
    opened = 0
    opened_edges = set()
    opened_nodes = set()
    # A physical Square capture column uses four SST nucleotides, i.e. the
    # two U-crossover edges named by the A0/B0 bridge assignment.  Earlier
    # code opened every crossover on all eight surface helices.  That left
    # four unrelated dangling SST edges in every realized capture column.
    # Sequence-only A1/B1 translations are deliberately excluded here.
    for assignment in capture_site_assignments(layout, "origin"):
        position = int(assignment["position"])
        for bridge in assignment["bridges"]:
            helix = int(bridge["sst_helix"])
            if (helix, position) in opened_nodes:
                continue
            record = rows[helix]["stap"][position]
            candidates = []
            for offset in (0, 2):
                partner, partner_base = map(
                    int, record[offset:offset + 2])
                if (partner >= 48 and partner != helix and
                        partner_base == position):
                    edge = tuple(sorted(
                        ((helix, position), (partner, partner_base))))
                    candidates.append((offset, partner, partner_base, edge))
            if len(candidates) != 1:
                raise ValueError(
                    "Square SST capture端点没有唯一的U型crossover边。")
            offset, partner, partner_base, edge = candidates[0]
            if edge in opened_edges:
                continue
            reverse = rows[partner]["stap"][partner_base]
            reverse_slots = [
                slot for slot in (0, 2)
                if reverse[slot:slot + 2] == [helix, position]]
            if len(reverse_slots) != 1:
                raise ValueError(
                    "Square SST capture待切边不是唯一互反边。")
            record[offset:offset + 2] = [-1, -1]
            reverse_slot = reverse_slots[0]
            reverse[reverse_slot:reverse_slot + 2] = [-1, -1]
            opened_edges.add(edge)
            opened_nodes.update(edge)
            opened += 2
    return opened


def _square_layer_fields(resource_rows, scaffold_range, staple_range,
                         array_length):
    return {
        "scaf": _tile_complete_sst_field(
            resource_rows, "scaf", (scaffold_range,), array_length),
        "stap": _tile_complete_sst_field(
            resource_rows, "stap", (staple_range,), array_length),
    }


def _square_component_colours(rows, resource_rows):
    """Rebuild colours at the actual 5' end after an auxiliary detour."""
    source = {int(row["num"]) + 48: {
        int(index): int(colour)
        for index, colour in row.get("stap_colors", [])}
        for row in resource_rows}
    for row in rows.values():
        row["stap_colors"] = []
    nodes = {(number, base)
             for number, row in rows.items()
             for base, record in enumerate(row["stap"])
             if record != _empty_record()}
    unseen = set(nodes)
    while unseen:
        component = set()
        stack = [next(iter(unseen))]
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            number, base = node
            record = rows[number]["stap"][base]
            for slot in (0, 2):
                partner = tuple(map(int, record[slot:slot + 2]))
                if partner in nodes and partner not in component:
                    stack.append(partner)
        unseen -= component
        starts = [(number, base) for number, base in component
                  if int(rows[number]["stap"][base][0]) < 0]
        if len(starts) != 1:
            raise ValueError("Square SST组件必须恰有一个5′端。")
        actual, base = starts[0]
        logical = actual - 16 if 64 <= actual < 80 else actual
        source_index = 64 + ((base - 64) % 32)
        colour = source.get(logical, {}).get(source_index, 0x0066CC)
        rows[actual]["stap_colors"].append([base, colour])
    for row in rows.values():
        row["stap_colors"].sort(key=lambda item: int(item[0]))


def _variable_sst_payload(name, layout, reserve_capture_gaps):
    resource = _load_json(SST_REFERENCE)
    resource_rows = _sst_reference_rows(resource)
    full_rows = {}
    if SEED_CAPTURE_REFERENCE.is_file():
        full_rows = {int(row["num"]): row for row in
                     _load_json(SEED_CAPTURE_REFERENCE)["vstrands"]}
    array_length = layout["array_length"]
    scaffold_ranges, staple_ranges = _sst_active_ranges(layout)
    layers = [_square_layer_fields(
        resource_rows, scaffold_ranges[index], staple_ranges[index],
        array_length) for index in (0, 1)]
    combined, auxiliary, unused_destination, auxiliary_metadata = \
        route_layer2_conflicts(layers[0], layers[1], array_length)
    layout["auxiliary_sst_routing"] = deepcopy(auxiliary_metadata)
    output = []
    for source in resource_rows:
        number = int(source["num"]) + 48
        reference_row = full_rows.get(number)
        row = dict(source)
        row["num"] = number
        if reference_row is not None:
            row["row"], row["col"] = (
                int(reference_row["row"]), int(reference_row["col"]))
        row["scaf"] = combined["scaf"][number]
        row["stap"] = combined["stap"][number]
        row["loop"] = [0] * array_length
        row["skip"] = [0] * array_length
        row["stap_colors"] = []
        output.append(row)
    if auxiliary_metadata["enabled"]:
        for source in resource_rows:
            logical = int(source["num"]) + 48
            number = logical + 16
            reference_row = full_rows.get(logical)
            row = dict(source)
            row["num"] = number
            if reference_row is not None:
                row["row"] = int(reference_row["row"])
                row["col"] = int(reference_row["col"]) + 12
            row["scaf"] = auxiliary["scaf"][number]
            row["stap"] = auxiliary["stap"][number]
            row["loop"] = [0] * array_length
            row["skip"] = [0] * array_length
            row["stap_colors"] = []
            output.append(row)
    rows = {int(row["num"]): row for row in output}
    _square_component_colours(rows, resource_rows)
    opened = _open_capture_sites(rows, layout) if reserve_capture_gaps else 0
    return {
        "name": name,
        "vstrands": output,
        "num_bases": array_length,
        "lattice": "square",
        "scaffold_colors": resource.get("scaffold_colors", []),
        "moire_structure_metadata": {
            "stage": ("sst_capture_ready" if reserve_capture_gaps
                      else "sst_complete_only"),
            "base_shift_bp": SST_BASE_SHIFT,
            "sst_scaffold_ranges": scaffold_ranges,
            "sst_staple_ranges": staple_ranges,
            "sst_duplex_ranges": layout["layer_ranges"],
            "capture_anchor_positions": layout["capture_positions"],
            "variable_length_layout": layout,
            "capture_gaps_reserved": bool(reserve_capture_gaps),
            "capture_gap_endpoint_count": opened,
            "auxiliary_sst_routing": deepcopy(auxiliary_metadata),
            "sst_unit_policy": (
                "shared centred geometry; complete reciprocal 32-nt U; "
                "one global absolute crossover phase"),
            "sequence_assignment": "pending step 3",
        },
    }


def _mixed_square_kagome_payload(name, z1_bp, z2_bp, z3_bp,
                                 reserve_capture_gaps=False,
                                 array_length=None):
    """Build Square layer 1 plus Kagome layer 2 on one centred canvas."""
    geometry = centered_square_sst_geometry(z1_bp, z2_bp, z3_bp)
    duplex_ranges = [list(map(int, item))
                     for item in geometry["layer_ranges"]]
    maximum = int(duplex_ranges[-1][1]) + 97
    length = max(544, int(array_length or 0),
                 32 * ((maximum + 31) // 32))
    square_resource = _load_json(SST_REFERENCE)
    square_rows = _sst_reference_rows(square_resource)
    # A mixed design still has one exact user spacing.  Some centred boundary
    # phases are legal for only one of the two polymer families (notably
    # 64/0/64).  Translate both layers together in 8-bp steps, preserving
    # their lengths and spacing, until Square layer 1 and Kagome layer 2 are
    # simultaneously legal.  This is a whole-design phase translation, not a
    # change in Z2 and not a movement of one layer relative to the other.
    square_layer = None
    square_duplex = None
    square_ranges = None
    kagome_layer = None
    kagome_duplex = None
    phase_adjustment = None
    for adjustment in (0, 8, -8, 16, -16, 24, -24, 32, -32):
        first_candidate = [duplex_ranges[0][0] + adjustment,
                           duplex_ranges[0][1] + adjustment]
        second_candidate = [duplex_ranges[1][0] + adjustment,
                            duplex_ranges[1][1] + adjustment]
        candidate_geometry = deepcopy(geometry)
        candidate_geometry["layer_ranges"] = [
            first_candidate, second_candidate]
        candidate_geometry["capture_support_ranges"] = [
            first_candidate, second_candidate]
        candidate_geometry["spacing_range"] = [
            int(value)+adjustment
            for value in geometry["spacing_range"]]
        candidate_geometry["target_envelope"] = [
            int(value)+adjustment
            for value in geometry["target_envelope"]]
        candidate_geometry["theoretical_capture_positions_by_layer"] = [
            [int(value)+adjustment for value in values]
            for values in geometry[
                "theoretical_capture_positions_by_layer"]]
        try:
            refresh_seed_overlap_metadata(candidate_geometry)
            candidate_geometry["capture_support_ranges"] = deepcopy(
                candidate_geometry["optimized_overlap_ranges"])
            scaffold_range, complement_range = \
                complete_square_polymer_ranges(first_candidate)
            candidate_square = _square_layer_fields(
                square_rows, scaffold_range, complement_range, length)
            candidate_kagome = build_kagome_layer_fields(
                second_candidate, length)
        except ValueError:
            continue
        square_layer = candidate_square
        square_duplex = first_candidate
        square_ranges = (scaffold_range, complement_range)
        kagome_layer = candidate_kagome
        kagome_duplex = second_candidate
        phase_adjustment = int(adjustment)
        geometry = candidate_geometry
        break
    if kagome_layer is None:
        raise ValueError(
            "Square–Kagome第二层无法在合法Kagome相位形成32/48-nt组件。")
    layers = [square_layer, {
        "scaf": kagome_layer["scaf"],
        "stap": kagome_layer["stap"],
    }]
    combined, auxiliary, unused_destination, auxiliary_metadata = \
        route_layer2_conflicts(layers[0], layers[1], length)

    full_rows = {}
    if SEED_CAPTURE_REFERENCE.is_file():
        full_rows = {int(row["num"]): row for row in
                     _load_json(SEED_CAPTURE_REFERENCE)["vstrands"]}
    output = []
    kagome_active = set(KAGOME_ACTIVE_HELICES)
    for source in square_rows:
        logical = int(source["num"]) + 48
        reference_row = full_rows.get(logical)
        row = dict(source)
        row["num"] = logical
        if reference_row is not None:
            row["row"], row["col"] = (
                int(reference_row["row"]), int(reference_row["col"]))
        row["scaf"] = combined["scaf"][logical]
        row["stap"] = combined["stap"][logical]
        row["loop"] = [0] * length
        row["skip"] = [0] * length
        row["stap_colors"] = []
        output.append(row)
    if auxiliary_metadata["enabled"]:
        for source in square_rows:
            logical = int(source["num"]) + 48
            number = logical + 16
            reference_row = full_rows.get(logical)
            row = dict(source)
            row["num"] = number
            if reference_row is not None:
                row["row"] = int(reference_row["row"])
                row["col"] = int(reference_row["col"]) + 12
            row["scaf"] = auxiliary["scaf"][number]
            row["stap"] = auxiliary["stap"][number]
            row["loop"] = [0] * length
            row["skip"] = [0] * length
            row["stap_colors"] = []
            output.append(row)
    rows = {int(row["num"]): row for row in output}
    _square_component_colours(rows, square_rows)

    # Capture support is the actual duplex intersection.  Square layer 1
    # follows the immutable Square phase; Kagome layer 2 uses its own endpoint
    # catalogue and active-mask/slot rules.
    actual_overlap_ranges = [list(map(int, item)) for item in
                             geometry["optimized_overlap_ranges"]]
    square_columns = list(seed_template_capture_columns(
        actual_overlap_ranges[0],
        int(geometry.get("coordinate_shift_bp", 0))))
    kagome_seed_columns = list(seed_template_capture_columns(
        actual_overlap_ranges[1],
        int(geometry.get("coordinate_shift_bp", 0))))
    kagome_candidates = kagome_layer_capture_catalogue(
        2, kagome_layer["stap"], kagome_duplex)
    allowed = set(kagome_seed_columns)
    kagome_candidates = [item for item in kagome_candidates
                         if int(item["position"]) in allowed]
    for item in kagome_candidates:
        item["capture_extension_nt"] = (
            16 if item["capture_family"] == "u_shaped_16nt" else 32)
    layout = {
        "lattice_type": "square_kagome",
        "lattice_by_layer": ["square", "kagome"],
        "z1_bp": int(z1_bp), "z2_bp": int(z2_bp),
        "z3_bp": int(z3_bp), "array_length": int(length),
        # The frozen Seed follows the shared centred placement, not either
        # SST duplex range.  Promote these fields so the scaffold subprocess
        # receives the same coordinate frame used by the preview and SST.
        "coordinate_shift_bp": int(
            geometry.get("coordinate_shift_bp", 0)),
        "seed_layer_ranges": deepcopy(geometry["seed_layer_ranges"]),
        "layer_ranges": [square_duplex, kagome_duplex],
        "capture_support_ranges": deepcopy(actual_overlap_ranges),
        "spacing_range": deepcopy(geometry["spacing_range"]),
        "seed_z2_range": deepcopy(geometry["spacing_range"]),
        "seed_partition_ranges": deepcopy(
            geometry["seed_partition_ranges"]),
        "seed_partition_lengths_bp": deepcopy(
            geometry["seed_partition_lengths_bp"]),
        "overlap_ranges": deepcopy(
            geometry["optimized_overlap_ranges"]),
        "seed_sst_overlap_bp": deepcopy(
            geometry["optimized_seed_overlap_bp"]),
        "scaffold_ranges": [list(square_ranges[0]),
                            list(kagome_layer["scaffold_range"])],
        "staple_ranges": [list(square_ranges[1]),
                          list(kagome_layer["staple_range"])],
        "sst_scaffold_ranges": [list(square_ranges[0]),
                                list(kagome_layer["scaffold_range"])],
        "sst_complementary_chain_ranges": [
            list(square_ranges[1]),
            list(kagome_layer["staple_range"])],
        "seed_capture_positions_by_layer": [
            square_columns, kagome_seed_columns],
        "capture_positions_by_layer": [
            square_columns, kagome_seed_columns],
        "theoretical_capture_positions_by_layer": [
            list(range(int(square_duplex[0]) +
                       ((8-int(square_duplex[0])) % 16),
                       int(square_duplex[1]) + 1, 16)),
            sorted({int(item["position"]) for item in
                    kagome_layer_capture_catalogue(
                        2, kagome_layer["stap"], kagome_duplex)})],
        "kagome_capture_anchor_assignments": deepcopy(kagome_candidates),
        "capture_extension_nt_by_layer": [16, 32],
        "auxiliary_sst_routing": deepcopy(auxiliary_metadata),
        "mixed_kagome_phase_adjustment_bp": phase_adjustment,
        "mixed_global_phase_translation_bp": phase_adjustment,
        "requested_spacing_bp": int(z2_bp),
        "actual_spacing_bp": int(kagome_duplex[0]) -
                             int(square_duplex[1]) - 1,
        "centered_geometry": deepcopy(geometry),
        "square_centered_geometry": deepcopy(geometry),
    }
    layout["capture_positions"] = [
        value for values in layout["capture_positions_by_layer"]
        for value in values]
    layout["capture_column_count"] = len(layout["capture_positions"])
    layout["pair_count_by_layer"] = [
        (len(values) + 1) // 2
        for values in layout["capture_positions_by_layer"]]
    layout["pair_count"] = sum(layout["pair_count_by_layer"])
    layout["capture_site_assignments"] = capture_site_assignments(layout)
    layout["capture_export_site_assignments"] = \
        capture_export_site_assignments(layout)

    if reserve_capture_gaps:
        # Square openings are the ordinary U crossovers.  Kagome candidates
        # carry their explicit slot, including actual auxiliary helices.
        square_open_layout = deepcopy(layout)
        square_open_layout["capture_positions_by_layer"] = [
            square_columns, []]
        _open_capture_sites(rows, square_open_layout)
        for item in kagome_candidates:
            logical = int(item["logical_sst_helix"])
            actual = auxiliary_actual_helix(
                layout, 2, "stap", logical, int(item["position"]))
            record = rows[actual]["stap"][int(item["position"])]
            slot = int(item["slot"])
            partner, partner_base = map(int, record[slot:slot + 2])
            if partner < 0:
                continue
            reverse = rows[partner]["stap"][partner_base]
            reverse_slots = [value for value in (0, 2)
                             if reverse[value:value + 2] ==
                             [actual, int(item["position"])]]
            if len(reverse_slots) != 1:
                raise ValueError("Square–Kagome capture待切边不是唯一互反边。")
            record[slot:slot + 2] = [-1, -1]
            other = reverse_slots[0]
            reverse[other:other + 2] = [-1, -1]
    return {
        "name": name,
        "vstrands": [rows[number] for number in sorted(rows)],
        "num_bases": length,
        "lattice": "square",
        "scaffold_colors": square_resource.get("scaffold_colors", []),
        "moire_structure_metadata": {
            "stage": ("sst_capture_ready" if reserve_capture_gaps else
                      "sst_complete_only"),
            "lattice_type": "square_kagome",
            "lattice_by_layer": ["square", "kagome"],
            "sst_scaffold_ranges": deepcopy(layout["scaffold_ranges"]),
            "sst_staple_ranges": deepcopy(layout["staple_ranges"]),
            "sst_duplex_ranges": deepcopy(layout["layer_ranges"]),
            "capture_gaps_reserved": bool(reserve_capture_gaps),
            "kagome_capture_anchor_assignments_sst_only": deepcopy(
                kagome_candidates),
            "kagome_theoretical_capture_positions_by_layer": [
                [], sorted({int(item["position"])
                            for item in kagome_layer_capture_catalogue(
                                2, kagome_layer["stap"], kagome_duplex)})],
            "auxiliary_sst_routing": deepcopy(auxiliary_metadata),
            "variable_length_layout": layout,
            "sequence_assignment": "pending step 3",
        },
    }


def build_shifted_sst_payload(name: str, reserve_capture_gaps: bool = False,
                              array_length: Optional[int] = None,
                              z1_bp: int = 128, z2_bp: int = 32,
                              z3_bp: int = 128,
                              seed_z1_bp: Optional[int] = None,
                              seed_z3_bp: Optional[int] = None,
                              capture_extension_length_nt: int = 16,
                              lattice_type: str = "square",
                              layers_design_sequence_identical: Optional[
                                  bool] = None
                              ) -> Dict[str, Any]:
    """Return complete SST-only layers, or an explicit capture-ready copy."""
    seed_z1_bp = int(z1_bp if seed_z1_bp is None else seed_z1_bp)
    seed_z3_bp = int(z3_bp if seed_z3_bp is None else seed_z3_bp)
    normalized_lattice = str(lattice_type).strip().lower().replace("-", "_")
    if normalized_lattice in ("square_kagome", "square_kagome_mixed"):
        payload = _mixed_square_kagome_payload(
            name, int(z1_bp), int(z2_bp), int(z3_bp),
            reserve_capture_gaps=reserve_capture_gaps,
            array_length=array_length)
        return payload_to_sst_first_numbering(payload)
    if normalized_lattice in ("kagome", "kagome_kagome"):
        from .calculations import phase_is_compatible
        if (layers_design_sequence_identical is True and
                not phase_is_compatible(
                    int(z1_bp), int(z2_bp), int(z3_bp))):
            raise ValueError(
                "Kagome双层设计与序列一致时，Z1/Z2/Z3必须满足"
                "32-bp相位联动规则。")
        payload = build_kagome_sst_payload(
            name, z1_bp=z1_bp, z2_bp=z2_bp, z3_bp=z3_bp,
            array_length=array_length,
            layers_design_sequence_identical=
                layers_design_sequence_identical)
        if reserve_capture_gaps:
            # Production capture-ready SSTs open only physical sites in the
            # fixed Seed/SST duplex intersection.  The complete theoretical
            # catalogue remains metadata/audit information and is never cut
            # merely because it exists in a 128-bp or other SST period.
            metadata = payload["moire_structure_metadata"]
            overlap = fixed_seed_overlap_layout(
                metadata["sst_duplex_ranges"],
                metadata["kagome_theoretical_capture_positions_by_layer"],
                lattice_type="kagome",
                seed_layer_ranges=metadata["sst_duplex_ranges"],
                seed_capture_positions_by_layer=metadata[
                    "variable_length_layout"][
                        "seed_capture_positions_by_layer"])
            payload = prepare_kagome_capture_sites(
                payload, name,
                seed_capture_positions_by_layer=overlap[
                    "actual_capture_positions_by_layer"],
                seed_layer_ranges=overlap["seed_layer_ranges"])
        return payload_to_sst_first_numbering(payload)
    if normalized_lattice not in ("square", "square_square",
                                   "square_square_c4"):
        raise ValueError("未知SST点阵类型：%s。" % lattice_type)
    from .calculations import phase_is_compatible
    if layers_design_sequence_identical is True and not phase_is_compatible(
            int(z1_bp), int(z2_bp), int(z3_bp)):
        raise ValueError(
            "双层设计与序列一致时，Z1/Z2/Z3必须满足"
            "Square 32-bp相位联动规则。")
    layout = structure_layout(
        int(z1_bp), int(z2_bp), int(z3_bp), seed_z1_bp, seed_z3_bp,
        int(capture_extension_length_nt))
    layout["seed_cross_section_preset"] = "s8_r4x4"
    layout["lattice_type"] = "square"
    layout["layers_design_sequence_identical"] = (
        None if layers_design_sequence_identical is None else
        bool(layers_design_sequence_identical))
    payload = _variable_sst_payload(
        name, layout, reserve_capture_gaps)
    return payload_to_sst_first_numbering(payload)


def build_output_sst_snapshot_payload(source: Dict[str, Any],
                                      name: str) -> Dict[str, Any]:
    """Close SST capture gaps in a copy without mutating the live design."""
    if not SST_REFERENCE.is_file():
        raise FileNotFoundError("找不到 Square_origamiseed-resource.json。")
    source_is_sst_first = source.get("moire_structure_metadata", {}).get(
        "helix_numbering") == "sst_first"
    payload = payload_to_internal_numbering(source)
    source_metadata = payload.get("moire_structure_metadata", {})
    lattice_type = source_metadata.get("lattice_type", "square")
    if lattice_type == "kagome":
        layout = deepcopy(source_metadata.get("variable_length_layout", {}))
        if not layout:
            raise ValueError("Kagome SST缺少variable_length_layout元数据。")
    else:
        layout = _layout_from_payload(payload)
    layer_ranges = tuple(tuple(item) for item in layout["layer_ranges"])
    payload["name"] = name
    rows = {int(row["num"]): row for row in payload.get("vstrands", [])}
    # Disconnect the Seed-facing half of every capture bridge.
    for number, row in rows.items():
        if number >= 48:
            continue
        for record in row.get("stap", []):
            for offset in (0, 2):
                partner, base = int(record[offset]), int(record[offset + 1])
                if partner >= 48:
                    record[offset:offset + 2] = [-1, -1]

    # At the same visible two-layer coordinates, replace only the SST staple
    # graph with the intact supplied template. Scaffold and applied scaffold
    # sequences remain unchanged, so complementary output sequences can be
    # regenerated deterministically.
    intact = payload_to_internal_numbering(build_shifted_sst_payload(
        name, reserve_capture_gaps=False,
        z1_bp=layout["z1_bp"], z2_bp=layout["z2_bp"],
        z3_bp=layout["z3_bp"],
        seed_z1_bp=layout.get("seed_z1_requested_bp"),
        seed_z3_bp=layout.get("seed_z3_requested_bp"),
        capture_extension_length_nt=layout.get(
            "capture_extension_nt", CAPTURE_EXTENSION_NT),
        lattice_type=lattice_type,
        layers_design_sequence_identical=layout.get(
            "layers_design_sequence_identical")))
    intact_rows = {int(row["num"]): row for row in intact["vstrands"]}
    actual_sst_helices = sorted(
        number for number in rows if number >= 48)
    for number in actual_sst_helices:
        row = rows.get(number)
        source_row = intact_rows.get(number)
        if row is None or source_row is None:
            continue
        for low, high in layout["staple_ranges"]:
            for index in range(low, min(high + 1, len(row.get("stap", [])))):
                row["stap"][index] = _empty_record()
        for index, record in enumerate(source_row.get("stap", [])):
            if index >= len(row.get("stap", [])):
                break
            if any(low <= index <= high
                   for low, high in layout["staple_ranges"]):
                row["stap"][index] = list(record)
        colors = {
            int(index): int(color)
            for index, color in row.get("stap_colors", [])
            if not any(low <= int(index) <= high
                       for low, high in layout["staple_ranges"])}
        colors.update({
            int(index): int(color)
            for index, color in source_row.get("stap_colors", [])
            if int(index) < len(row.get("stap", []))})
        row["stap_colors"] = [list(item) for item in sorted(colors.items())]

    metadata = payload.setdefault("moire_structure_metadata", {})
    metadata.update({
        "export_role": "output",
        "capture_connections": "disconnected in export snapshot only",
        "sst_capture_gaps": "closed from original SST template",
        "source_design_unchanged": True,
    })
    return (payload_to_sst_first_numbering(payload)
            if source_is_sst_first else payload)


def build_capture_ready_sst_payload(source: Dict[str, Any],
                                    name: str) -> Dict[str, Any]:
    """Open capture sites in a copy; never mutate the SST-only source."""
    source_is_sst_first = source.get("moire_structure_metadata", {}).get(
        "helix_numbering") == "sst_first"
    payload = payload_to_internal_numbering(source)
    lattice_type = payload.get("moire_structure_metadata", {}).get(
        "lattice_type")
    if lattice_type in ("square_kagome", "square_kagome_mixed"):
        source_metadata = payload["moire_structure_metadata"]
        layout = source_metadata["variable_length_layout"]
        prepared = _mixed_square_kagome_payload(
            name, layout["z1_bp"], layout["z2_bp"], layout["z3_bp"],
            reserve_capture_gaps=True,
            array_length=layout.get("array_length"))
        prepared_metadata = prepared["moire_structure_metadata"]
        prepared_layout = prepared_metadata["variable_length_layout"]
        # The mixed-lattice builder reconstructs the SST topology. Preserve
        # accepted project inputs that are not derived topology fields;
        # otherwise scaffold review silently loses the requested Seed Z2
        # indel and the final design reverts to zero insertion/deletion.
        for key, value in layout.items():
            prepared_layout.setdefault(key, deepcopy(value))
        for key, value in source_metadata.items():
            prepared_metadata.setdefault(key, deepcopy(value))
        return payload_to_sst_first_numbering(prepared)
    if lattice_type == "kagome":
        metadata = payload["moire_structure_metadata"]
        theoretical = (metadata.get(
            "kagome_theoretical_capture_positions_by_layer") or
            metadata.get("variable_length_layout", {}).get(
                "theoretical_capture_positions_by_layer"))
        if not theoretical:
            raise ValueError("Kagome SST缺少理论capture候选目录。")
        overlap = fixed_seed_overlap_layout(
            metadata.get("sst_duplex_ranges",
                         metadata["variable_length_layout"]["layer_ranges"]),
            theoretical, lattice_type="kagome",
            seed_layer_ranges=metadata.get(
                "sst_duplex_ranges",
                metadata["variable_length_layout"]["layer_ranges"]),
            seed_capture_positions_by_layer=metadata[
                "variable_length_layout"][
                    "seed_capture_positions_by_layer"])
        prepared = prepare_kagome_capture_sites(
            payload, name,
            seed_capture_positions_by_layer=overlap[
                "actual_capture_positions_by_layer"],
            seed_layer_ranges=overlap["seed_layer_ranges"])
        return (payload_to_sst_first_numbering(prepared)
                if source_is_sst_first else prepared)
    payload["name"] = name
    layout = _layout_from_payload(payload)
    rows = {int(row["num"]): row for row in payload.get("vstrands", [])}
    opened = _open_capture_sites(rows, layout)
    metadata = payload.setdefault("moire_structure_metadata", {})
    metadata.update({
        "stage": "sst_capture_ready",
        "capture_gaps_reserved": True,
        "capture_gap_endpoint_count": opened,
        "derived_from_complete_sst": True,
        "complete_sst_source_unchanged": True,
    })
    return (payload_to_sst_first_numbering(payload)
            if source_is_sst_first else payload)


def build_complete_sst_only_payload(source: Dict[str, Any],
                                    name: str) -> Dict[str, Any]:
    """Return a sequence-bearing, Seed-free snapshot of complete SST units."""
    payload = build_output_sst_snapshot_payload(
        payload_to_internal_numbering(source), name)
    kept = []
    for row in payload.get("vstrands", []):
        if int(row["num"]) < 48:
            continue
        row = deepcopy(row)
        for field in ("scaf", "stap"):
            for record in row.get(field, []):
                for offset in (0, 2):
                    partner = int(record[offset])
                    if 0 <= partner < 48:
                        record[offset:offset + 2] = [-1, -1]
        kept.append(row)
    payload["vstrands"] = kept
    payload["scaffold_sequences"] = [
        deepcopy(item)
        for item in payload.get("scaffold_sequences", [])
        if int(item.get("start_vh", -1)) in set(
            payload.get("moire_structure_metadata", {}).get(
                "sst_helix_numbers", CAPTURE_OUTPUT_HELICES))]
    payload = payload_to_sst_first_numbering(payload)
    metadata = payload.setdefault("moire_structure_metadata", {})
    metadata.update({
        "stage": "sst_complete_only_with_sequences",
        "export_role": "complete_sst",
        "seed_helices_removed": True,
        "sst_helix_export_map": {
            str(old): old - 48 for old in CAPTURE_OUTPUT_HELICES},
        "capture_connections": "not present",
        "capture_gaps_reserved": False,
        "sst_unit_policy": "complete 32-nt U units",
    })
    return payload


def write_shifted_sst(filename: str, z1_bp: int = 128,
                      z2_bp: int = 32, z3_bp: int = 128,
                      seed_z1_bp: Optional[int] = None,
                      seed_z3_bp: Optional[int] = None,
                      capture_extension_length_nt: int = 16,
                      lattice_type: str = "square",
                      seed_cross_section_preset: str = "s8_r4x4",
                      layers_design_sequence_identical: Optional[
                          bool] = None,
                      mean_indel_per_helix: float = 0.0) -> Path:
    target = Path(filename).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = build_shifted_sst_payload(
        target.name, reserve_capture_gaps=False,
        z1_bp=z1_bp, z2_bp=z2_bp, z3_bp=z3_bp,
        seed_z1_bp=seed_z1_bp, seed_z3_bp=seed_z3_bp,
        capture_extension_length_nt=capture_extension_length_nt,
        lattice_type=lattice_type,
        layers_design_sequence_identical=
            layers_design_sequence_identical)
    metadata = payload.setdefault("moire_structure_metadata", {})
    metadata["seed_cross_section_preset"] = str(seed_cross_section_preset)
    layout = metadata.setdefault("variable_length_layout", {})
    layout.update({
        "seed_cross_section_preset": str(seed_cross_section_preset),
        # Legacy callers may still pass Seed Z1/Z3 values.  They are ignored:
        # the physical Seed is the immutable 2L template in every project.
        "seed_z1_requested_bp": FIXED_SEED_NOMINAL_SUPPORT_BP[0],
        "seed_z3_requested_bp": FIXED_SEED_NOMINAL_SUPPORT_BP[1],
        "capture_extension_nt": int(capture_extension_length_nt),
        "layers_design_sequence_identical": (
            None if layers_design_sequence_identical is None else
            bool(layers_design_sequence_identical)),
        # Placement is performed later against the immutable Seed template.
        # Keeping the requested mean in the SST metadata makes scaffold
        # review and finalization use the same deterministic Z2 allocation.
        "mean_indel_per_helix": float(mean_indel_per_helix),
    })
    target.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                      encoding="utf-8")
    return target


def _worker(mode: str, output: str, input_json: Optional[str] = None) -> Dict[str, Any]:
    command = worker_command("structure", mode, str(output))
    if input_json:
        command.extend(["--input", str(input_json)])
    completed = subprocess.run(
        command, check=False, text=True, capture_output=True)
    if completed.returncode:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError("结构生成失败：%s" % detail)
    try:
        return json.loads(completed.stdout)
    except Exception as error:
        raise RuntimeError("结构生成器没有返回有效报告。") from error


def generate_scaffold_review(filename: str,
                             sst_json: Optional[str] = None) -> Dict[str, Any]:
    return _worker("scaffold", filename, sst_json)


def estimate_scaffold_capacity(sst_json: str) -> Dict[str, Any]:
    """Return exact routing-aware capacity before the design is accepted."""
    return _worker("capacity", str(Path(sst_json).with_suffix(".capacity")),
                   sst_json)


def finalize_structure(scaffold_json: str, filename: str) -> Dict[str, Any]:
    return _worker("finalize", filename, scaffold_json)


def export_sequence_variants(structure_json: str, output_directory: str,
                             base_name: Optional[str] = None) -> Dict[str, Any]:
    """Export capture and complete-SST JSON, XLSX and SVG snapshots."""
    command = worker_command(
        "sequence-export", str(structure_json), str(output_directory))
    if base_name:
        command.extend(["--name", str(base_name)])
    completed = subprocess.run(
        command, check=False, text=True, capture_output=True)
    if completed.returncode:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError("序列与SVG导出失败：%s" % detail)
    try:
        return json.loads(completed.stdout)
    except Exception as error:
        raise RuntimeError("导出器没有返回有效报告。") from error


def validate_sst(filename: str) -> Dict[str, Any]:
    """Validate an SST-only expert file without requiring Seed helices."""
    path = Path(filename).expanduser().resolve()
    payload = payload_to_internal_numbering(_load_json(path))
    payload_lattice = payload.get("moire_structure_metadata", {}).get(
        "lattice_type")
    if payload_lattice == "kagome":
        return validate_kagome_sst_payload(payload)
    raw_rows = payload.get("vstrands", [])
    rows = {int(row["num"]): row for row in raw_rows}
    errors = []
    warnings = []

    layout = _layout_from_payload(payload)
    layer_ranges = tuple(tuple(item) for item in layout["layer_ranges"])
    metadata = payload.get("moire_structure_metadata", {})
    scaffold_ranges = tuple(tuple(item) for item in metadata.get(
        "sst_scaffold_ranges", layer_ranges))
    staple_ranges = tuple(tuple(item) for item in metadata.get(
        "sst_staple_ranges", layout["staple_ranges"]))
    routing = metadata.get("auxiliary_sst_routing", {})
    auxiliary_enabled = bool(
        isinstance(routing, dict) and routing.get("enabled"))
    auxiliary_numbers = set(AUXILIARY_INTERNAL) if auxiliary_enabled else set()
    expected_numbers = set(CAPTURE_OUTPUT_HELICES) | auxiliary_numbers
    if set(rows) != expected_numbers or len(raw_rows) != len(expected_numbers):
        expected_label = ("helix 48–63及辅助helix 64–79"
                          if auxiliary_enabled else "helix 48–63")
        errors.append("SST文件必须只包含%s，且不能有重复编号。" %
                      expected_label)
    reference_coordinates = {}
    if SEED_CAPTURE_REFERENCE.is_file():
        reference = _load_json(SEED_CAPTURE_REFERENCE)
        reference_coordinates = {
            int(row["num"]): (int(row["row"]), int(row["col"]))
            for row in reference["vstrands"]
            if int(row["num"]) in expected_numbers}
    changed_coordinates = []
    for number, expected in reference_coordinates.items():
        row = rows.get(number)
        if row is not None and (int(row["row"]), int(row["col"])) != expected:
            changed_coordinates.append(number)
    if changed_coordinates:
        errors.append("SST截面不可修改；坐标变化的helix：%s。" %
                      ", ".join(map(str, changed_coordinates)))
    protected_gap_count = 0
    full_rows = {}
    if SEED_CAPTURE_REFERENCE.is_file():
        full_rows = {int(row["num"]): row for row in
                     _load_json(SEED_CAPTURE_REFERENCE)["vstrands"]}
    lattice_by_layer = list(layout.get("lattice_by_layer") or [])
    if not lattice_by_layer:
        lattice_by_layer = [payload_lattice] * len(scaffold_ranges)

    def _logical_expected_ranges(number, ranges):
        """Return only ranges occupied by this helix in each real layer."""
        expected_ranges = []
        for range_index, bounds in enumerate(ranges):
            layer_lattice = (lattice_by_layer[range_index]
                             if range_index < len(lattice_by_layer)
                             else payload_lattice)
            if layer_lattice == "kagome" and number not in \
                    KAGOME_ACTIVE_HELICES:
                continue
            expected_ranges.append(tuple(map(int, bounds)))
        return expected_ranges

    # h64–79 are sequence-topology detours for h48–63, not additional
    # physical SST helices.  Validate each logical helix using the disjoint
    # union of its primary and auxiliary records.
    for number in CAPTURE_OUTPUT_HELICES:
        row = rows.get(number)
        auxiliary_row = rows.get(number + 16) if auxiliary_enabled else None
        if row is None:
            continue
        scaffold = row.get("scaf", [])
        required_length = max(high for unused_low, high in
                              scaffold_ranges + staple_ranges) + 1
        if len(scaffold) < required_length:
            errors.append("helix %d 的base范围不足%d。" %
                          (number, required_length))
            continue
        occupied_primary = {
            index for index, record in enumerate(scaffold)
            if record != _empty_record()}
        occupied_auxiliary = ({
            index for index, record in enumerate(
                auxiliary_row.get("scaf", []))
            if record != _empty_record()
        } if auxiliary_row is not None else set())
        # Logical primary/auxiliary overlap is expected: it is precisely the
        # caDNAno display collision that the auxiliary channel resolves.
        # The actual rows are different, so only their union is audited.
        occupied = occupied_primary | occupied_auxiliary
        expected_scaffold_ranges = _logical_expected_ranges(
            number, scaffold_ranges)
        expected = set().union(*(
            set(range(low, high + 1))
            for low, high in expected_scaffold_ranges)) \
            if expected_scaffold_ranges else set()
        if occupied != expected:
            errors.append(
                "helix %d 的SST scaffold必须恰好覆盖本层活动区间：%s。" %
                (number, "、".join("%d–%d" % item
                                    for item in expected_scaffold_ranges)))
        staple_primary = {
            index for index, record in enumerate(row.get("stap", []))
            if record != _empty_record()}
        staple_auxiliary = ({
            index for index, record in enumerate(
                auxiliary_row.get("stap", []))
            if record != _empty_record()
        } if auxiliary_row is not None else set())
        staple_occupied = staple_primary | staple_auxiliary
        expected_staple_ranges = _logical_expected_ranges(
            number, staple_ranges)
        expected_staples = set().union(*(
            set(range(low, high + 1))
            for low, high in expected_staple_ranges)) \
            if expected_staple_ranges else set()
        if staple_occupied != expected_staples:
            errors.append(
                "helix %d 的完整SST链必须恰好覆盖：%s。" %
                (number, "、".join("%d–%d" % item
                                    for item in expected_staple_ranges)))
        reference_row = full_rows.get(number)
        if reference_row is not None:
            for actual_row in filter(None, (row, auxiliary_row)):
                for low, high in staple_ranges:
                    for index in range(low, min(
                            high + 1, len(actual_row.get("stap", [])))):
                        current = actual_row["stap"][index]
                        for offset in (0, 2):
                            if 0 <= int(current[offset]) < 48:
                                errors.append(
                                    "SST-only文件不能提前连接Seed：helix %d base %d。" %
                                    (int(actual_row["num"]), index))
                            if current[offset:offset + 2] == [-1, -1]:
                                protected_gap_count += 1
    scaffold_components = _components(payload, "scaf", expected_numbers)
    staple_components = _components(payload, "stap", expected_numbers)
    allowed_scaffold_lengths = ({32, 48}
                                if payload_lattice == "square_kagome"
                                else {32})
    capture_ready = bool(metadata.get("capture_gaps_reserved"))
    allowed_staple_lengths = ({16, 32, 48}
                              if payload_lattice == "square_kagome" and
                              capture_ready else allowed_scaffold_lengths)
    kagome_short_staple_ranges = [
        tuple(map(int, bounds))
        for index, bounds in enumerate(staple_ranges)
        if index < len(lattice_by_layer) and
        lattice_by_layer[index] == "kagome" and
        int(bounds[1]) - int(bounds[0]) + 1 == 64
    ]

    def _is_template_phased_kagome_16(component):
        """Allow only the complete-SST 64-nt 48+16 edge exception.

        Mixed validation must mirror the Kagome layer solver without
        weakening Square validation or treating every 16-mer as legal.  The
        component must be a reciprocal non-scaffold line component, lie
        wholly inside the Kagome layer's exactly-64-nt staple range, and
        touch that physical range boundary.
        """
        if payload_lattice != "square_kagome" or capture_ready or \
                component["length"] != 16 or not component["reciprocal"]:
            return False
        logical_helices = {
            helix - 16 if helix in AUXILIARY_INTERNAL else helix
            for helix in component["helices"]
        }
        if not logical_helices <= set(KAGOME_LINE_HELICES["stap"]):
            return False
        low = int(component["base_min"])
        high = int(component["base_max"])
        return any(range_low <= low <= high <= range_high and
                   (low == range_low or high == range_high)
                   for range_low, range_high in
                   kagome_short_staple_ranges)

    invalid_scaffold = [item for item in scaffold_components
                        if item["length"] not in allowed_scaffold_lengths or
                        not item["reciprocal"]]
    invalid_sst = [item for item in staple_components
                   if (item["length"] not in allowed_staple_lengths and
                       not _is_template_phased_kagome_16(item)) or
                   not item["reciprocal"]]
    if invalid_scaffold:
        errors.append("SST-only中的scaffold必须是本点阵合法、互反的32/48-nt单元。")
    if invalid_sst:
        errors.append(
            "SST-only中的SST必须是本点阵合法、互反的32/48-nt单元；"
            "Kagome 16-nt仅允许范本相位的64-nt边界48+16特例。")
    return {
        "path": str(path),
        "valid": not errors,
        "errors": errors,
        "warnings": warnings,
        "helix_count": len(rows),
        "sst_ranges": [
            {"range": list(item), "complete": not errors}
            for item in layer_ranges],
        "sst_first_base": layer_ranges[0][0],
        "variable_length_layout": layout,
        "protected_capture_gap_endpoints": protected_gap_count,
        "scaffold_32nt_component_count": len(scaffold_components),
        "sst_32nt_component_count": len(staple_components),
        "complete_32nt_units": not invalid_scaffold and not invalid_sst,
    }


def _components(payload: Dict[str, Any], field: str,
                helix_filter: Optional[set] = None) -> List[Dict[str, Any]]:
    rows = {int(row["num"]): row for row in payload.get("vstrands", [])}
    color_map = ({
        (number, int(index)): int(color)
        for number, row in rows.items()
        for index, color in row.get("stap_colors", [])
    } if field == "stap" else {})
    nodes = set()
    for number, row in rows.items():
        if helix_filter is not None and number not in helix_filter:
            continue
        for index, record in enumerate(row.get(field, [])):
            if record != _empty_record():
                nodes.add((number, index))
    neighbors = {node: [] for node in nodes}
    reciprocal = True
    for node in nodes:
        number, index = node
        record = rows[number][field][index]
        for offset in (0, 2):
            other = (int(record[offset]), int(record[offset + 1]))
            if other[0] < 0:
                continue
            if other not in nodes:
                reciprocal = False
                continue
            neighbors[node].append(other)
            other_record = rows[other[0]][field][other[1]]
            if list(node) not in (other_record[0:2], other_record[2:4]):
                reciprocal = False
    result = []
    visited = set()
    for node in sorted(nodes):
        if node in visited:
            continue
        stack = [node]
        visited.add(node)
        component = set()
        while stack:
            current = stack.pop()
            component.add(current)
            for other in neighbors[current]:
                if other not in visited:
                    visited.add(other)
                    stack.append(other)
        degrees = [len(neighbors[item]) for item in component]
        actual_length = sum(
            1 + int(rows[number].get("loop", [])[index]) +
            int(rows[number].get("skip", [])[index])
            for number, index in component)
        result.append({
            "length": len(component),
            "actual_length": actual_length,
            "is_loop": bool(component) and all(value == 2 for value in degrees),
            "end_count": sum(value == 1 for value in degrees),
            "helices": sorted({item[0] for item in component}),
            "reciprocal": reciprocal,
            "colors": sorted({color_map[item] for item in component
                              if item in color_map}),
            "base_min": min(item[1] for item in component),
            "base_max": max(item[1] for item in component),
            "output_base_min": min(
                (item[1] for item in component if item[0] >= 48),
                default=None),
        })
    return result


def _staple_component_details(payload: Dict[str, Any]):
    """Return node labels, actual lengths and color markers for staples."""
    rows = {int(row["num"]): row for row in payload.get("vstrands", [])}
    nodes = {
        (number, index)
        for number, row in rows.items()
        for index, record in enumerate(row.get("stap", []))
        if record != _empty_record()}
    adjacency = {node: set() for node in nodes}
    for number, index in nodes:
        record = rows[number]["stap"][index]
        for offset in (0, 2):
            other = (int(record[offset]), int(record[offset + 1]))
            if other in nodes:
                adjacency[(number, index)].add(other)
    color_map = {
        (number, int(index)): int(color)
        for number, row in rows.items()
        for index, color in row.get("stap_colors", [])}
    components = []
    labels = {}
    visited = set()
    for first in sorted(nodes):
        if first in visited:
            continue
        component = {first}
        stack = [first]
        visited.add(first)
        while stack:
            node = stack.pop()
            for other in adjacency[node]:
                if other not in visited:
                    visited.add(other)
                    component.add(other)
                    stack.append(other)
        component_index = len(components)
        for node in component:
            labels[node] = component_index
        actual_length = sum(
            1 + int(rows[number].get("loop", [])[index]) +
            int(rows[number].get("skip", [])[index])
            for number, index in component)
        components.append({
            "nodes": component,
            "actual_length": actual_length,
            "colors": {
                color_map[node] for node in component if node in color_map},
        })
    return components, labels




def _short_staple_audit(payload: Dict[str, Any], minimum_length: int = 21):
    """Classify residual short staples after automatic merge attempts."""
    rows = {int(row["num"]): row for row in payload.get("vstrands", [])}
    nodes = {
        (number, index)
        for number, row in rows.items()
        for index, record in enumerate(row.get("stap", []))
        if record != _empty_record()}
    neighbors = {node: set() for node in nodes}
    for number, index in nodes:
        record = rows[number]["stap"][index]
        for offset in (0, 2):
            other = (int(record[offset]), int(record[offset + 1]))
            if other in nodes:
                neighbors[(number, index)].add(other)
    components = []
    labels = {}
    visited = set()
    for first in sorted(nodes):
        if first in visited:
            continue
        component = set([first])
        stack = [first]
        visited.add(first)
        while stack:
            node = stack.pop()
            for other in neighbors[node]:
                if other not in visited:
                    visited.add(other)
                    component.add(other)
                    stack.append(other)
        component_index = len(components)
        for node in component:
            labels[node] = component_index
        components.append(component)
    layout = payload.get("moire_structure_metadata", {}).get(
        "variable_length_layout", structure_layout())
    capture_targets = {
        (int(bridge["seed_helix"]), int(assignment["position"]))
        for assignment in capture_export_site_assignments(layout)
        for bridge in assignment["bridges"]}
    capture_flags = [
        ((any(number < 48 for number, unused in component) and
          any(number >= 48 for number, unused in component)) or
         any(node in capture_targets for node in component))
        for component in components]
    result = []
    for component_index, component in enumerate(components):
        # SST-only strands are the designed SST output topology rather than
        # Seed staples.  Its native 16-nt domains must not be judged by the
        # Seed ordinary-staple 21--58-nt rule.
        if not any(number < 48 for number, unused in component):
            continue
        # Capture cores come directly from the accepted Seed template.  They
        # are not ordinary staples and must never enter a length audit.
        if capture_flags[component_index]:
            continue
        if len(component) >= minimum_length:
            continue
        physical_edge = False
        blocked_by_capture = False
        for number, index in component:
            scaffold = rows[number].get("scaf", [])
            for neighbor_index in (index - 1, index + 1):
                if (neighbor_index < 0 or neighbor_index >= len(scaffold) or
                        scaffold[neighbor_index] == _empty_record()):
                    physical_edge = True
                    continue
                other_component = labels.get((number, neighbor_index))
                if other_component is not None and \
                        other_component != component_index and \
                        capture_flags[other_component]:
                    blocked_by_capture = True
        if blocked_by_capture:
            reason = "adjacent_nick_reserved_for_capture_extension"
        elif physical_edge:
            reason = "physical_structure_edge_without_capture_protection"
        else:
            reason = "unprotected_internal_short_staple"
        result.append({
            "length": len(component),
            "helices": sorted({number for number, unused in component}),
            "base_min": min(index for unused, index in component),
            "base_max": max(index for unused, index in component),
            "protected": bool(blocked_by_capture or
                              reason ==
                              "physical_structure_edge_without_capture_protection"),
            "reason": reason,
        })
    return result


def _seed_scaffold_crossover_audit(
        rows: Dict[int, Dict[str, Any]], seed_preset: str = "s8_r4x4"):
    """Audit Seed scaffold crossovers against native Square lattice sites."""
    categories = {
        "non_neighbor": 0,
        "mismatched_index": 0,
        "illegal_phase": 0,
        "nonreciprocal": 0,
        "different_direction_clearance": 0,
    }
    examples = []
    total = 0
    invalid_edges = 0
    directed_events = defaultdict(list)
    for number in range(48):
        row = rows.get(number)
        if row is None:
            continue
        source_coord = (int(row["row"]), int(row["col"]))
        even = source_coord[0] % 2 == source_coord[1] % 2
        neighbor_coords = (
            ((source_coord[0], source_coord[1] + 1),
             (source_coord[0] + 1, source_coord[1]),
             (source_coord[0], source_coord[1] - 1),
             (source_coord[0] - 1, source_coord[1]))
            if even else
            ((source_coord[0], source_coord[1] - 1),
             (source_coord[0] - 1, source_coord[1]),
             (source_coord[0], source_coord[1] + 1),
             (source_coord[0] + 1, source_coord[1])))
        for index, record in enumerate(row.get("scaf", [])):
            partner, partner_index = map(int, record[2:4])
            if partner < 0 or partner == number:
                continue
            directed_events[number].append((index, partner, partner_index))
            total += 1
            failed = []
            partner_row = rows.get(partner)
            if partner_row is None or partner not in range(48):
                failed.append("non_neighbor")
            else:
                partner_coord = (int(partner_row["row"]),
                                 int(partner_row["col"]))
                if partner_coord not in neighbor_coords:
                    failed.append("non_neighbor")
                else:
                    direction = neighbor_coords.index(partner_coord)
                    residues = set(SQUARE_SCAF_LOW[direction]) | \
                        set(SQUARE_SCAF_HIGH[direction])
                    if index % 32 not in residues:
                        failed.append("illegal_phase")
                if partner_index != index:
                    failed.append("mismatched_index")
                reciprocal = False
                if 0 <= partner_index < len(partner_row.get("scaf", [])):
                    partner_record = partner_row["scaf"][partner_index]
                    reciprocal = any(
                        list(map(int, partner_record[offset:offset + 2])) ==
                        [number, index]
                        for offset in (0, 2))
                if not reciprocal:
                    failed.append("nonreciprocal")
            if failed:
                invalid_edges += 1
                for category in set(failed):
                    categories[category] += 1
                if len(examples) < 8:
                    examples.append({
                        "from": [number, index],
                        "to": [partner, partner_index],
                        "reasons": sorted(set(failed)),
                    })
    clearance_conflicts = 0
    minimum_clearance = None
    for number, events in directed_events.items():
        ordered = sorted(events)
        for left, right in zip(ordered, ordered[1:]):
            left_index, left_partner, unused_left_partner_index = left
            right_index, right_partner, unused_right_partner_index = right
            if left_partner == right_partner:
                continue
            index_gap = right_index - left_index
            minimum_clearance = (index_gap if minimum_clearance is None
                                 else min(minimum_clearance, index_gap))
            # Inclusive distance >=8 bp is index difference >=7.  This is a
            # hard cadnano scaffold rule, not an optimization preference.
            if index_gap >= 7:
                continue
            clearance_conflicts += 1
            categories["different_direction_clearance"] += 1
            if len(examples) < 8:
                examples.append({
                    "helix": number,
                    "from": [left_partner, left_index],
                    "to": [right_partner, right_index],
                    "index_gap": index_gap,
                    "inclusive_distance_bp": index_gap + 1,
                    "reasons": ["different_direction_clearance"],
                })
    return {
        "total": total,
        "invalid": invalid_edges + clearance_conflicts,
        "invalid_crossover_edges": invalid_edges,
        "different_direction_clearance_conflicts": clearance_conflicts,
        "minimum_different_direction_index_gap": minimum_clearance,
        "categories": categories,
        "examples": examples,
    }


def validate_structure(filename: str, require_staples: bool = False) -> Dict[str, Any]:
    path = Path(filename).expanduser().resolve()
    payload = payload_to_internal_numbering(_load_json(path))
    raw_rows = payload.get("vstrands", [])
    rows = {int(row["num"]): row for row in raw_rows}
    layout = _layout_from_payload(payload)
    layer_ranges = tuple(tuple(item) for item in layout["layer_ranges"])
    metadata = payload.get("moire_structure_metadata", {})
    routing = metadata.get("auxiliary_sst_routing", {})
    auxiliary_enabled = bool(
        isinstance(routing, dict) and routing.get("enabled"))
    expected = set(range(64)) | (
        set(AUXILIARY_INTERNAL) if auxiliary_enabled else set())
    errors = []
    warnings = []

    # Workflow Z1/Z3 are the real Seed--SST duplex overlaps, not the wider
    # scaffold edge envelope.  Workflow Z2 is a physical length and therefore
    # includes every insertion/deletion stored in the frozen Seed rows.
    overlap_values = (layout.get("seed_sst_overlap_bp") or
                      layout.get("seed_partition_lengths_bp") or [])
    seed_z1_overlap = (
        float(overlap_values[0]) if len(overlap_values) >= 1 else None)
    seed_z3_overlap = (
        float(overlap_values[1]) if len(overlap_values) == 2 else
        float(overlap_values[2]) if len(overlap_values) >= 3 else None)
    z2_range = (layout.get("seed_z2_indel_range") or
                layout.get("seed_z2_range") or
                layout.get("spacing_range"))
    physical_z2_values = []
    if isinstance(z2_range, (list, tuple)) and len(z2_range) == 2:
        z2_low, z2_high = map(int, z2_range)
        for number in range(48):
            row = rows.get(number, {})
            loops = row.get("loop", [])
            skips = row.get("skip", [])
            if z2_low < 0 or z2_high >= min(len(loops), len(skips)):
                continue
            physical_z2_values.append(sum(
                1 + int(loops[index]) + int(skips[index])
                for index in range(z2_low, z2_high + 1)))
    if physical_z2_values:
        seed_z2_physical = (
            sum(physical_z2_values) / float(len(physical_z2_values)))
    else:
        seed_z2_physical = float(layout.get(
            "actual_z2_spacing_bp",
            float(layout.get("z2_bp", 0)) + float(layout.get(
                "mean_indel_per_helix_actual", 0.0) or 0.0)))
    if set(rows) != expected:
        errors.append(
            "helix编号必须为0–63%s。" %
            ("，并在启用辅助绕行时包含64–79"
             if auxiliary_enabled else ""))
    if len(raw_rows) != len(rows):
        errors.append("存在重复的helix编号。")
    if len({(row.get("row"), row.get("col")) for row in rows.values()}) != len(rows):
        errors.append("存在重复的左视图helix坐标。")
    seed_preset = layout.get("seed_cross_section_preset", "s8_r4x4")
    coordinate_reference = SEED_CAPTURE_REFERENCE
    if coordinate_reference.is_file() and set(range(64)).issubset(rows):
        reference_rows = {
            int(row["num"]): row for row in
            _load_json(coordinate_reference)["vstrands"]
            if int(row["num"]) < 64}
        changed = [
            number for number, row in rows.items()
            if number in reference_rows and
            (int(row["row"]), int(row["col"])) !=
            (int(reference_rows[number]["row"]),
             int(reference_rows[number]["col"]))]
        if changed:
            errors.append("结构截面不可修改；坐标变化的helix：%s。" %
                          ", ".join(map(str, changed)))
    frozen_seed = layout.get("seed_geometry_policy") == \
        "immutable_2L_reference"
    crossover_audit = (_seed_scaffold_crossover_audit(rows, seed_preset)
                       if not frozen_seed else {
                           "total": 0, "invalid": 0,
                           "invalid_crossover_edges": 0,
                           "different_direction_clearance_conflicts": 0,
                           "minimum_different_direction_index_gap": None,
                           "categories": {
                               "non_neighbor": 0,
                               "mismatched_index": 0,
                               "illegal_phase": 0,
                               "nonreciprocal": 0,
                               "different_direction_clearance": 0,
                           },
                           "examples": [],
                           "skipped_for_immutable_template": True,
                       })
    if crossover_audit["invalid"]:
        categories = crossover_audit["categories"]
        errors.append(
            "Seed scaffold存在%d个非法crossover（非相邻%d、"
            "端点索引不一致%d、Square相位非法%d、非互反%d、"
            "同helix不同向间距小于8 bp %d）。" %
            (crossover_audit["invalid"], categories["non_neighbor"],
             categories["mismatched_index"],
             categories["illegal_phase"],
             categories["nonreciprocal"],
             categories["different_direction_clearance"]))
    seed_edge_ranges = {}
    for number in range(48):
        row = rows.get(number, {})
        occupied = [
            index for index, record in enumerate(row.get("scaf", []))
            if record != _empty_record()]
        if occupied:
            seed_edge_ranges[number] = (min(occupied), max(occupied))
    selected_edge_limit = layout.get("seed_edge_stagger_limit_used_bp")
    preferred_11_count = layout.get(
        "seed_edge_preferred_11_scaffold_count")
    selected_edge_count = layout.get(
        "seed_edge_selected_scaffold_count")
    relaxed_reduced_count = bool(layout.get(
        "seed_edge_21_relaxation_reduced_scaffold_count", False))
    capture_low = capture_high = None
    frozen_two_layer_reference = bool(layout.get(
        "seed_routing_is_frozen_reference", False))
    reviewed_edge_growth = bool(layout.get(
        "seed_edge_growth_uses_reviewed_envelope", False))
    if len(seed_edge_ranges) == 48:
        low_edge_stagger = (
            max(value[0] for value in seed_edge_ranges.values()) -
            min(value[0] for value in seed_edge_ranges.values()))
        high_edge_stagger = (
            max(value[1] for value in seed_edge_ranges.values()) -
            min(value[1] for value in seed_edge_ranges.values()))
        if max(low_edge_stagger, high_edge_stagger) > 21:
            errors.append(
                "Seed同侧边缘错位必须不超过21 bp；当前左侧%d bp、"
                "右侧%d bp。" % (low_edge_stagger, high_edge_stagger))
        elif max(low_edge_stagger, high_edge_stagger) > 11:
            if frozen_two_layer_reference or reviewed_edge_growth:
                warnings.append(
                    "冻结2L范本边缘错位为%d bp；保留范本routing，"
                    "仍满足21-bp硬上限。" % high_edge_stagger)
            elif selected_edge_limit != 21 or \
                    preferred_11_count is None or \
                    selected_edge_count is None or \
                    selected_edge_count >= preferred_11_count or \
                    not relaxed_reduced_count:
                errors.append(
                    "Seed只有在21-bp候选能减少scaffold数量时才允许"
                    "同侧错位超过11 bp；当前左侧%d bp、右侧%d bp，"
                    "但没有记录有效的减链收益。" %
                    (low_edge_stagger, high_edge_stagger))
            else:
                warnings.append(
                    "Seed为减少scaffold数量而使用21-bp合法放宽；"
                    "当前左侧%d bp、右侧%d bp。" %
                    (low_edge_stagger, high_edge_stagger))
        elif selected_edge_limit == 21:
            errors.append(
                "Seed已能在11-bp优选范围内完成，禁止保留无收益的"
                "21-bp放宽状态。")
        capture_helices = set(range(0, 8)) | set(range(24, 32))
        capture_low = min(seed_edge_ranges[number][0]
                          for number in capture_helices)
        capture_high = max(seed_edge_ranges[number][1]
                           for number in capture_helices)
        if not frozen_two_layer_reference and not reviewed_edge_growth and (
                capture_low != min(value[0]
                                   for value in seed_edge_ranges.values()) or
                capture_high != max(
                    value[1] for value in seed_edge_ranges.values())):
            errors.append(
                "Seed的最大实际Z1/Z3长度必须由capture helix定义；"
                "检测到普通helix伸得更远。")
        first_layer_end = int(layer_ranges[0][1])
        second_layer_start = int(layer_ranges[1][0])
        actual_seed_z1 = max(
            first_layer_end - low + 1
            for low, unused_high in seed_edge_ranges.values())
        actual_seed_z3 = max(
            high - second_layer_start + 1
            for unused_low, high in seed_edge_ranges.values())
        requested_seed_z1 = int(layout.get(
            "seed_z1_requested_bp", actual_seed_z1))
        requested_seed_z3 = int(layout.get(
            "seed_z3_requested_bp", actual_seed_z3))
        seed_z1_growth = max(0, actual_seed_z1 - requested_seed_z1)
        seed_z3_growth = max(0, actual_seed_z3 - requested_seed_z3)
        pair_counts = list(layout.get("pair_count_by_layer", [0, 0]))
        capture_growth_exception = bool(layout.get(
            "phase_quantized_edge_growth_exception_used",
            layout.get("capture_pair_growth_exception_used", False)))
        if max(seed_z1_growth, seed_z3_growth) > 11 and \
                not frozen_two_layer_reference and \
                not capture_growth_exception:
            errors.append(
                "Seed Z1/Z3超过11 bp但没有记录8-bp输入参数经合法"
                "scaffold相位量化产生的边缘增长例外。")
    else:
        low_edge_stagger = high_edge_stagger = None
        actual_seed_z1 = actual_seed_z3 = None
        requested_seed_z1 = requested_seed_z3 = None
        seed_z1_growth = seed_z3_growth = None
        capture_growth_exception = False
    seed_components = _components(payload, "scaf", set(range(48)))
    seed_loops = [item for item in seed_components if item["is_loop"]]
    single_nick_seed = [item for item in seed_components
                        if not item["is_loop"] and item["end_count"] == 2]
    invalid_seed = [item for item in seed_components
                    if not item["is_loop"] and item["end_count"] != 2]
    expected_scaffolds = int(payload.get("moire_structure_metadata", {}).get(
        "seed_scaffold_count", 2))
    seed_scaffold_total_nt = sum(
        int(item["actual_length"]) for item in seed_components)
    try:
        capacity_plan = scaffold_capacity_plan(seed_scaffold_total_nt)
    except ValueError as error:
        capacity_plan = None
        errors.append(str(error))
    if capacity_plan is not None and \
            expected_scaffolds != capacity_plan["count"]:
        errors.append(
            "Seed scaffold元数据分段数%d与总长度%d nt规定的%d条不一致。" %
            (expected_scaffolds, seed_scaffold_total_nt,
             capacity_plan["count"]))
    if len(seed_components) != expected_scaffolds:
        errors.append("Seed scaffold应按容量分为%d条；当前为%d条。" %
                      (expected_scaffolds, len(seed_components)))
        short_components = [item["length"] for item in seed_components
                            if item["length"] <= 64]
        if short_components:
            errors.append(
                "检测到未通过seam并入主scaffold的边缘32-bp小模块：%s。" %
                short_components)
    if invalid_seed:
        errors.append("Seed存在分支或多个nick的异常scaffold组件。")
    scaffold_capacity = (
        capacity_plan["per_scaffold_capacity_nt"]
        if capacity_plan is not None else SCAFFOLD_CAPACITY_ORTHOGONAL)
    oversized = [item["actual_length"] for item in seed_components
                 if item["actual_length"] > scaffold_capacity]
    if oversized:
        errors.append("scaffold超过当前分段允许的%d nt：%s" %
                      (scaffold_capacity, oversized))
    if capacity_plan is not None and \
            len(seed_components) != capacity_plan["count"]:
        errors.append(
            "Seed scaffold总长度%d nt必须使用%d条；当前为%d条。" %
            (seed_scaffold_total_nt, capacity_plan["count"],
             len(seed_components)))
    # Validate only the legal capture columns inside the actual support
    # overlap.  SST overhang is intentional and must not be treated as a
    # missing Seed scaffold error.
    missing_anchors = []
    for number in CAPTURE_SEED_HELICES:
        row = rows.get(number)
        if row is None:
            continue
        for position in layout["capture_positions"]:
            if len(row.get("scaf", [])) <= position or \
                    row["scaf"][position] == _empty_record():
                missing_anchors.append({"helix": number,
                                        "position": position})
    if missing_anchors:
        errors.append("Seed scaffold没有覆盖全部合法capture列。")
    sst_ranges = []
    # Square occupies all sixteen SST sites.  Kagome deliberately uses the
    # validated 12-of-16 cross-section, so its four geometric holes must not
    # be reported as missing scaffold.  Capture finalization never changes
    # the SST scaffold graph; this check therefore audits every active Kagome
    # helix across the accepted scaffold windows and leaves the holes empty.
    lattice_by_layer = layout.get("lattice_by_layer")
    metadata_scaffold_ranges = payload.get(
        "moire_structure_metadata", {}).get("sst_scaffold_ranges")
    checked_sst_ranges = tuple(
        tuple(map(int, item)) for item in
        (metadata_scaffold_ranges or layer_ranges))
    for range_index, (low, high) in enumerate(checked_sst_ranges):
        layer_lattice = (
            lattice_by_layer[range_index]
            if lattice_by_layer and range_index < len(lattice_by_layer)
            else layout.get("lattice_type"))
        active_output_helices = (
            tuple(KAGOME_ACTIVE_HELICES)
            if layer_lattice == "kagome" else CAPTURE_OUTPUT_HELICES)
        occupied = all(
            all(
                (rows[number]["scaf"][index] != _empty_record()) or
                (auxiliary_enabled and
                 rows[number + 16]["scaf"][index] != _empty_record())
                for index in range(low, high + 1))
            for number in active_output_helices if number in rows)
        sst_ranges.append({"range": [low, high], "complete": occupied})
        if not occupied:
            errors.append("capture SST scaffold区间%d–%d不完整。" %
                          (low, high))
    staple_components = _components(payload, "stap") if require_staples else []
    capture_components = [
        item for item in staple_components
        if item["helices"] and min(item["helices"]) < 48 and
        max(item["helices"]) >= 48]
    detailed_staples, unused_staple_labels = \
        _staple_component_details(payload) if require_staples else ([], {})
    normal_staple_lengths = [
        int(item["actual_length"]) for item in detailed_staples
        if not (any(number < 48 for number, unused in item["nodes"]) and
                any(number >= 48 for number, unused in item["nodes"]))]
    normal_staple_components = [
        item for item in detailed_staples
        if not (any(number < 48 for number, unused in item["nodes"]) and
                any(number >= 48 for number, unused in item["nodes"]))]
    normal_staple_length_histogram = {}
    if normal_staple_lengths:
        for lower in range(10, max(normal_staple_lengths) + 1, 10):
            count = sum(lower <= length <= lower + 9
                        for length in normal_staple_lengths)
            if count:
                normal_staple_length_histogram[
                    "%d–%d" % (lower, lower + 9)] = round(
                        100.0 * count / len(normal_staple_lengths), 1)

    def _has_continuous_16(component):
        by_helix = {}
        for number, index in component["nodes"]:
            by_helix.setdefault(int(number), set()).add(int(index))
        for indices in by_helix.values():
            longest = current = 0
            previous = None
            for index in sorted(indices):
                current = current + 1 if previous is not None and \
                    index == previous + 1 else 1
                longest = max(longest, current)
                previous = index
            if longest >= 16:
                return True
        return False

    continuous_16_count = sum(
        _has_continuous_16(item) for item in normal_staple_components)
    continuous_16_percentage = (
        100.0 * continuous_16_count / len(normal_staple_components)
        if normal_staple_components else 0.0)
    oversized_normal_staples = [
        length for length in normal_staple_lengths
        if length > NORMAL_STAPLE_MAX_NT]
    if require_staples and oversized_normal_staples:
        message = "范本中有%d条超过%d nt的普通staple；最长为%d nt。" % (
            len(oversized_normal_staples), NORMAL_STAPLE_MAX_NT,
            max(oversized_normal_staples))
        if layout.get("seed_geometry_policy") == "immutable_2L_reference":
            warnings.append(message + " 固定2L Seed不重新break。")
        else:
            errors.append(message)
    capture_palette = {
        color for item in capture_components for color in item["colors"]}
    if require_staples and not staple_components:
        errors.append("尚未生成staple/capture链。")
    seed_staple_missing = []
    if require_staples:
        staple_required_range = layout.get(
            "seed_staple_required_coverage_range")
        if staple_required_range is None:
            staple_required_range = [
                int(layer_ranges[0][0]), int(layer_ranges[-1][1])]
        staple_low, staple_high = map(int, staple_required_range)
        for number in range(48):
            row = rows.get(number, {})
            scaffold = row.get("scaf", [])
            staples = row.get("stap", [])
            for index, record in enumerate(scaffold):
                # Frozen-short Seed routing intentionally retains scaffold-
                # only support outside the requested physical duplex.  Only
                # the exact designed staple band must be paired.
                if staple_low <= index <= staple_high and \
                        record != _empty_record() and (
                        index >= len(staples) or
                        staples[index] == _empty_record()):
                    seed_staple_missing.append([number, index])
        if seed_staple_missing:
            errors.append(
                "Seed scaffold上有%d个碱基未被正常staple覆盖。" %
                len(seed_staple_missing))
    short_staples = (_short_staple_audit(payload)
                     if require_staples else [])
    protected_short_staples = [
        item for item in short_staples if item["protected"]]
    invalid_short_staples = [
        item for item in short_staples if not item["protected"]]
    if require_staples and protected_short_staples:
        warnings.append(
            "保留%d条无法安全并入的Seed边缘短staple（相邻nick用于"
            "capture延伸）；最短为%d nt。" % (
                len(protected_short_staples),
                min(item["length"] for item in protected_short_staples)))
    if require_staples and invalid_short_staples:
        errors.append(
            "存在%d条不受capture/物理边缘保护的内部短staple。" %
            len(invalid_short_staples))
    expected_bridges = int(layout["expected_capture_bridges"])
    if require_staples and len(capture_components) != expected_bridges:
        errors.append("当前长度设计应有%d条Seed–SST capture桥；当前为%d条。" %
                      (expected_bridges, len(capture_components)))
    assignments = capture_site_assignments(layout)
    expected_capture_edges = {
        (int(bridge["seed_helix"]), int(bridge["sst_helix"]),
         int(assignment["position"]))
        for assignment in assignments
        for bridge in assignment["bridges"]}
    actual_capture_edges = set()
    capture_bases_sharing_seed_crossovers = set()
    if require_staples:
        for number in range(48):
            row = rows.get(number, {})
            for position, record in enumerate(row.get("stap", [])):
                has_capture = any(
                    int(record[offset]) >= 48 and
                    int(record[offset + 1]) == position
                    for offset in (0, 2))
                has_seed_crossover = any(
                    0 <= int(record[offset]) < 48 and
                    int(record[offset]) != number
                    for offset in (0, 2))
                if has_capture and has_seed_crossover:
                    capture_bases_sharing_seed_crossovers.add(
                        (number, position))
                for offset in (0, 2):
                    partner = int(record[offset])
                    partner_position = int(record[offset + 1])
                    if partner >= 48 and partner_position == position:
                        actual_capture_edges.add(
                            (number, partner, position))
    if require_staples and capture_bases_sharing_seed_crossovers:
        errors.append(
            "capture base不得与Seed内部staple crossover共位；发现%s。" %
            sorted(capture_bases_sharing_seed_crossovers))
    missing_capture_edges = sorted(
        expected_capture_edges - actual_capture_edges)
    unexpected_capture_edges = sorted(
        actual_capture_edges - expected_capture_edges)
    if require_staples and missing_capture_edges:
        errors.append(
            "Seed capture face缺少%d条A0/B0物理映射连接。" %
            len(missing_capture_edges))
    if require_staples and unexpected_capture_edges:
        errors.append(
            "Seed capture face存在%d条不属于A0/B0物理映射的连接。" %
            len(unexpected_capture_edges))
    capture_face_coverage = {}
    expected_seed_helices_by_face = defaultdict(set)
    for assignment in assignments:
        for bridge in assignment["bridges"]:
            expected_seed_helices_by_face[str(bridge["face"])].add(
                int(bridge["seed_helix"]))
    for face in CAPTURE_FACE_DEFINITIONS:
        all_face = set(map(int, face["internal_seed_helices"]))
        # Derive coverage from the accepted lattice-specific assignment.
        # This remains the same four helices for Square, while Kagome uses
        # only the legal anchors present in its 12-of-16 SST cross-section.
        expected_face = expected_seed_helices_by_face[face["id"]]
        export_only_face = set(map(
            int, face["export_only_internal_seed_helices"]))
        actual_face = {
            seed for seed, unused_sst, unused_position in
            actual_capture_edges if seed in expected_face}
        capture_face_coverage[face["id"]] = {
            "all_internal_helices": sorted(all_face),
            "expected_internal_helices": sorted(expected_face),
            "export_only_internal_helices": sorted(export_only_face),
            "actual_internal_helices": sorted(actual_face),
            "missing_internal_helices": sorted(expected_face - actual_face),
            "complete": actual_face == expected_face,
        }
        if require_staples and actual_face != expected_face:
            errors.append(
                "%s未覆盖应与%s SST直接连接的%d根Seed helix；"
                "缺少%s。" %
                (face["id"], layout.get("lattice_type", "square"),
                 len(expected_face), sorted(expected_face - actual_face)))
    capture_pair_colors = {}
    capture_colors_by_layer = {}
    if require_staples:
        grouped_colors = {index: set()
                          for index in range(layout["pair_count"])}
        for item in capture_components:
            base = (item.get("output_base_min")
                    if item.get("output_base_min") is not None
                    else item["base_min"])
            pair_index = capture_pair_index(base, layout)
            if pair_index is None:
                errors.append("发现不属于两层capture范围的连接链。")
                continue
            if len(item["colors"]) != 1:
                errors.append(
                    "capture pair %d 中有链未使用唯一颜色。" %
                    (pair_index + 1))
            grouped_colors[pair_index].update(item["colors"])
        capture_pair_colors = {
            index + 1: sorted(colors)
            for index, colors in grouped_colors.items()}
        layer_colors = {
            layer_index + 1: set()
            for layer_index in range(len(
                layout["capture_positions_by_layer"]))}
        # Capture cores are immutable catalogue/template topology.  Validate
        # their presence only; never calculate or enforce a core length.
        for assignment in capture_export_site_assignments(layout):
            for bridge in assignment["bridges"]:
                node = (int(bridge["seed_helix"]),
                        int(assignment["position"]))
                component_index = unused_staple_labels.get(node)
                if component_index is None:
                    errors.append(
                        "合法capture位点helix %d/base %d没有staple。" % node)
                    continue
                component = detailed_staples[component_index]
                colors = set(component["colors"])
                layer_colors[int(assignment["layer"])].update(colors)
        capture_colors_by_layer = {
            layer: sorted(colors) for layer, colors in layer_colors.items()}
    return {
        "path": str(path),
        "valid": not errors,
        "errors": errors,
        "warnings": warnings,
        "seed_scaffold_components": seed_components,
        "seed_scaffold_lengths": sorted(
            item["actual_length"] for item in seed_components),
        "seed_scaffold_total_nt": seed_scaffold_total_nt,
        "seed_scaffold_capacity_plan": capacity_plan,
        "seed_scaffold_loop_count": len(seed_loops),
        "seed_scaffold_single_nick_count": len(single_nick_seed),
        "seed_scaffold_crossover_audit": crossover_audit,
        "seed_scaffold_edge_ranges": seed_edge_ranges,
        "seed_scaffold_low_edge_stagger_bp": low_edge_stagger,
        "seed_scaffold_high_edge_stagger_bp": high_edge_stagger,
        "seed_edge_stagger_limit_used_bp": selected_edge_limit,
        "seed_edge_preferred_11_scaffold_count": preferred_11_count,
        "seed_edge_selected_scaffold_count": selected_edge_count,
        "seed_edge_21_relaxation_reduced_scaffold_count":
            relaxed_reduced_count,
        "capture_helices_define_maximum_actual_length": (
            len(seed_edge_ranges) == 48 and
            capture_low == min(value[0]
                               for value in seed_edge_ranges.values()) and
            capture_high == max(value[1]
                                for value in seed_edge_ranges.values())),
        "seed_z1_requested_bp": requested_seed_z1,
        "seed_z3_requested_bp": requested_seed_z3,
        "seed_z1_actual_bp": actual_seed_z1,
        "seed_z3_actual_bp": actual_seed_z3,
        "seed_z1_overlap_bp": seed_z1_overlap,
        "seed_z2_actual_bp": seed_z2_physical,
        "seed_z3_overlap_bp": seed_z3_overlap,
        "seed_z1_edge_growth_bp": seed_z1_growth,
        "seed_z3_edge_growth_bp": seed_z3_growth,
        "capture_pair_growth_exception_used": capture_growth_exception,
        "phase_quantized_edge_growth_exception_used":
            capture_growth_exception,
        "capture_anchor_missing": missing_anchors,
        "sst_ranges": sst_ranges,
        "staple_component_count": len(staple_components),
        "seed_staple_missing_base_count": len(seed_staple_missing),
        "seed_staple_missing_examples": seed_staple_missing[:16],
        "short_staple_count": len(short_staples),
        "protected_short_staple_count": len(protected_short_staples),
        "invalid_short_staple_count": len(invalid_short_staples),
        "short_staple_audit": short_staples,
        "minimum_staple_length": min(
            (item["length"] for item in staple_components), default=0),
        "minimum_normal_staple_length": min(
            normal_staple_lengths, default=0),
        "maximum_normal_staple_length": max(
            normal_staple_lengths, default=0),
        "normal_staple_minimum_nt": NORMAL_STAPLE_MIN_NT,
        "normal_staple_maximum_nt": NORMAL_STAPLE_MAX_NT,
        "normal_staple_length_histogram":
            normal_staple_length_histogram,
        "continuous_16_base_count": continuous_16_count,
        "continuous_16_base_percentage": round(
            continuous_16_percentage, 1),
        "capture_bridge_component_count": len(capture_components),
        "capture_color_count": len(capture_palette),
        "capture_pair_colors": capture_pair_colors,
        "capture_colors_by_layer": capture_colors_by_layer,
        "capture_extension_nt": capture_extension_nt(layout),
        "capture_core_source": "immutable_template",
        "capture_core_length_evaluated": False,
        "capture_face_coverage": capture_face_coverage,
        "capture_mapping_missing": [list(item)
                                    for item in missing_capture_edges],
        "capture_mapping_unexpected": [list(item)
                                       for item in unexpected_capture_edges],
        "expected_seed_scaffold_count": expected_scaffolds,
        "expected_capture_bridge_count": expected_bridges,
        "variable_length_layout": layout,
    }
