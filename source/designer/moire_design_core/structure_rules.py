"""Authoritative rules for Moiré Designer stage-2 structure generation.

This module is deliberately data-only.  It must not import caDNAno and it
must not inspect metadata cached by an older generated JSON.  Stage-2 code
may derive geometry from the immutable fixtures, but every policy decision
comes from :data:`RULES` below.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple


RESOURCE_ROOT = Path(__file__).with_name("resources")


@dataclass(frozen=True)
class ScaffoldCapacityRule:
    p8064: int = 8064
    orthogonal: int = 7557
    maximum_count: int = 3

    def count_for(self, total_nt: int) -> int:
        total_nt = int(total_nt)
        if total_nt <= self.p8064:
            return 1
        if total_nt <= 2 * self.orthogonal:
            return 2
        if total_nt <= 3 * self.orthogonal:
            return 3
        raise ValueError(
            "Seed scaffold总长度%d nt超过3条正交scaffold的容量%d nt；"
            "请减小Seed长度。" %
            (total_nt, 3 * self.orthogonal))


@dataclass(frozen=True)
class CrossoverRule:
    square_period_bp: int = 32
    square_inclusive_clearance_bp: int = 8
    scaffold_staple_same_pair_exclusion_bp: int = 10
    square_avoid_coordinate_multiple: int = 8

    @property
    def square_minimum_index_gap(self) -> int:
        return self.square_inclusive_clearance_bp - 1


@dataclass(frozen=True)
class Stage2Rules:
    capacities: ScaffoldCapacityRule = ScaffoldCapacityRule()
    crossovers: CrossoverRule = CrossoverRule()
    minimum_capture_pairs_per_layer: int = 2
    allowed_seed_cross_section: str = "8x8_minus_4x4_pore"
    sst_length_step_bp: int = 8
    sst_translation_period_bp: int = 32
    validated_lengths_bp: Tuple[int, ...] = (96, 104, 112, 120, 128)


RULES = Stage2Rules()


FIXTURES = {
    "square": {
        length: RESOURCE_ROOT / ("square_sst_%dbp_fixture.json" % length)
        for length in RULES.validated_lengths_bp
    },
    "kagome": {
        length: RESOURCE_ROOT / ("kagome_sst_%dbp_fixture.json" % length)
        for length in RULES.validated_lengths_bp
    },
}

# Seed routing and capture topology intentionally have different sources.
#
# The two-layer Square file is the *only* routing baseline used by stage 2.
# The old three-layer files must never be used to generate, crop, repeat or
# partition a two-layer Seed.  Kagome's historical three-layer design is kept
# solely as an endpoint catalogue for the Kagome-specific capture projection;
# it is not a Seed scaffold/staple routing template.
SEED_ROUTING_REFERENCE = (
    RESOURCE_ROOT / "Square_Seed_2L_newtemplate.json")
SQUARE_CAPTURE_REFERENCE = SEED_ROUTING_REFERENCE
KAGOME_CAPTURE_TOPOLOGY_REFERENCE = (
    RESOURCE_ROOT / "Kagome_Seed_Ka-seed-pore_3L.json")

# Compatibility alias for readers that only need the selected two-layer
# routing reference.  Do not add a Kagome routing entry here.
SEED_REFERENCES = {
    "square": SEED_ROUTING_REFERENCE,
    "kagome": SEED_ROUTING_REFERENCE,
}


def validated_sst_fixture(lattice: str, length_bp: int) -> Path:
    lattice = str(lattice).strip().lower()
    length_bp = int(length_bp)
    try:
        path = FIXTURES[lattice][length_bp]
    except KeyError as exc:
        raise ValueError(
            "没有%s %d-bp的冻结SST范本。" %
            (lattice, length_bp)) from exc
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def assert_supported_seed_preset(preset: str) -> None:
    if str(preset) not in ("s8_r4x4", "8x8_minus_4x4_pore"):
        raise ValueError(
            "第2步只支持8×8 Seed减4×4 pore；不支持其他Seed截面。")


def assert_capture_base_is_free(record, seed_helix_limit: int = 48) -> None:
    """Reject a capture base that contains a Seed-internal crossover."""
    for offset in (0, 2):
        partner = int(record[offset])
        if 0 <= partner < int(seed_helix_limit):
            raise ValueError(
                "capture base与Seed内部staple crossover共位；"
                "这表示SST/capture的32-bp整体平移或相位错误。"
                "不得删除或移动AutoCS crossover。")
