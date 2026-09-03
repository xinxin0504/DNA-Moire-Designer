"""Calibrated calculations for the first Square moire-bilayer prototype."""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Dict, Iterable, List

from .models import MoireProject, SquareBilayerSettings
from .missing_helix_calibration import (
    S8_R4X4C_BRANCH_RMSE, S8_R4X4C_CALIBRATION_VERSION,
    S8_R4X4C_FAILED_INDEL_POINTS, S8_R4X4C_TWIST_POINTS,
    S8_R4X4C_VALIDATED_INDEL_RANGE, calibration_profile,
    default_seed_cells, reference_profile)
from .square_sst_geometry import centered_square_sst_geometry


DNA_RISE_NM = 0.34
SQUARE_SST_HEIGHT_BP = 128
SQUARE_DOMAIN_BP = 8
SQUARE_REPEAT_BP = 32
SQUARE_NATIVE_PITCH_BP_PER_TURN = 10.67
SQUARE_ELASTIC_PITCH_LOW = 10.33
SQUARE_ELASTIC_PITCH_HIGH = 10.67
S8_R4X4C_CALIBRATION_SPAN_BP = 96.0
MAX_SEED_INSERTION_PER_HELIX = 10.0
MAX_SEED_DELETIONS_PER_DOMAIN = 3
FIXED_SEED_TOTAL_BP = 288
REFERENCE_SEED_Z2_BP = 32

# Direct SNUPI calibration of the S8-R4x4C cross-section.  The extended
# curve contains all nine converged observations from -10 through +6
# indels/helix.  The failed +8/+10 simulations are metadata only.
S8_R4X4C_NATIVE_TWIST_DEG_PER_BASE = next(
    twist for indel, twist in S8_R4X4C_TWIST_POINTS if indel == 0.0)
CAPTURE_BASES = 16
CAPTURE_BASELINE_DEG = (
    CAPTURE_BASES * S8_R4X4C_NATIVE_TWIST_DEG_PER_BASE)
# Kube et al. elastic-model anchors used by cadnano Twist/Bend.  The Square
# 6x6/8x8 ordering follows the physically monotonic correction documented in
# cadnano: increasing cross-section stiffness lowers the native twist.
SQUARE_ELASTIC_J_NM4 = (31.0, 527.0, 2695.0, 8545.0)
SQUARE_ELASTIC_TWIST_1067 = (.9375, .5625, .15625, .078125)
SQUARE_ELASTIC_TWIST_1033 = (.265625, 0.0, .046875, 0.0)
LATTICE_LAYERS = {
    "square_square_c4": ("square", "square"),
    "kagome_kagome": ("kagome", "kagome"),
    "square_kagome": ("square", "kagome"),
}


def moire_period_from_angle(angle_deg: float, lattice_constant_nm: float) -> float:
    """Return the moire period for two identical, unstrained 2D lattices."""
    angle = abs(float(angle_deg))
    if angle <= 1e-12:
        return math.inf
    return float(lattice_constant_nm) / (
        2.0 * math.sin(math.radians(angle) / 2.0))


def angle_from_moire_period(period_nm: float, lattice_constant_nm: float) -> float:
    """Return the small relative angle for identical unstrained lattices."""
    period = float(period_nm)
    lattice = float(lattice_constant_nm)
    if period <= 0:
        raise ValueError("Moiré period must be positive.")
    ratio = lattice / (2.0 * period)
    if ratio > 1.0:
        raise ValueError("Moiré period is smaller than the geometric limit.")
    return math.degrees(2.0 * math.asin(ratio))


def _piecewise_linear(value: float, points: Iterable[tuple]) -> float:
    """Interpolate a monotonic calibration and linearly extrapolate its ends."""
    ordered = tuple(sorted((float(x), float(y)) for x, y in points))
    value = float(value)
    pair = ordered[:2] if value <= ordered[0][0] else ordered[-2:]
    for left, right in zip(ordered, ordered[1:]):
        if left[0] <= value <= right[0]:
            pair = (left, right)
            break
    (x0, y0), (x1, y1) = pair
    return y0 + (value-x0)*(y1-y0)/(x1-x0)


def calibrated_twist_per_base(indel_per_helix: float) -> float:
    """Interpolate the direct 96-bp S8-R4x4C SNUPI observations.

    This function is retained for calibration provenance.  Design prediction
    uses :func:`elastic_calibrated_twist_for_cross_section`, which first
    converts indel density to effective pitch and an elastic twist estimate.
    """
    return _piecewise_linear(indel_per_helix, S8_R4X4C_TWIST_POINTS)


def calibrated_indel_for_twist_per_base(twist_per_base: float) -> float:
    """Inverse of :func:`calibrated_twist_per_base`."""
    inverse = tuple((twist, indel) for indel, twist in S8_R4X4C_TWIST_POINTS)
    return _piecewise_linear(twist_per_base, inverse)


def calibrated_twist_for_cross_section(indel_per_helix: float, cells,
                                       size: int = 8) -> float:
    """Return topology-aware twist while preserving the direct S8 curve."""
    info = calibration_profile(cells, size)
    selected = info["profile"]
    reference = reference_profile()
    reference_twist = calibrated_twist_per_base(indel_per_helix)
    side = 2 if float(indel_per_helix) < 0 else 3
    sensitivity = selected[side] / max(reference[side], 1e-12)
    return selected[1] + sensitivity * (reference_twist - reference[1])


def calibrated_indel_for_cross_section(twist_per_base: float, cells,
                                       size: int = 8) -> float:
    """Invert the topology-aware calibration, including linear extrapolation.

    The previous bounded bisection silently saturated at +/-64 indels per
    helix.  That changed a user-entered large Twist into a smaller angle
    before the one-sided +10 insertion feasibility check could report it.
    Each branch is an affine transform of the frozen reference curve, so its
    inverse can be evaluated directly without imposing an artificial bound.
    """
    target = float(twist_per_base)
    selected = calibration_profile(cells, size)["profile"]
    reference = reference_profile()
    branch = 2 if target < selected[1] else 3
    sensitivity = selected[branch] / max(reference[branch], 1e-12)
    reference_target = (
        reference[1] + (target-selected[1]) / sensitivity)
    return calibrated_indel_for_twist_per_base(reference_target)


def _log_interpolate(value: float, xs, ys) -> float:
    """Interpolate an elastic anchor against log(J), clamping its ends."""
    value = max(1e-12, float(value))
    if value <= xs[0]:
        return float(ys[0])
    if value >= xs[-1]:
        return float(ys[-1])
    log_value = math.log(value)
    for index in range(len(xs)-1):
        if xs[index] <= value <= xs[index+1]:
            fraction = ((log_value-math.log(xs[index])) /
                        (math.log(xs[index+1])-math.log(xs[index])))
            return float(ys[index] +
                         fraction*(ys[index+1]-ys[index]))
    return float(ys[-1])


def _square_polar_moment_nm4(cells) -> float:
    """Return the 2-nm mechanical polar moment used by the calibration."""
    points = [(float(row)*2.0, float(column)*2.0)
              for row, column in cells]
    if not points:
        return 0.0
    center_x = sum(point[0] for point in points) / len(points)
    center_y = sum(point[1] for point in points) / len(points)
    intrinsic = math.pi / 2.0
    area = math.pi
    return sum(intrinsic + area*((x-center_x)**2 + (y-center_y)**2)
               for x, y in points)


def _square_elastic_twist_per_base(effective_pitch: float, cells) -> float:
    """Raw Square elastic twist density at one effective helical pitch."""
    polar_moment = _square_polar_moment_nm4(cells)
    at_low = _log_interpolate(
        polar_moment, SQUARE_ELASTIC_J_NM4, SQUARE_ELASTIC_TWIST_1033)
    at_high = _log_interpolate(
        polar_moment, SQUARE_ELASTIC_J_NM4, SQUARE_ELASTIC_TWIST_1067)
    fraction = ((float(effective_pitch)-SQUARE_ELASTIC_PITCH_LOW) /
                (SQUARE_ELASTIC_PITCH_HIGH-SQUARE_ELASTIC_PITCH_LOW))
    return at_low + fraction*(at_high-at_low)


def elastic_calibrated_twist_for_cross_section(
        indel_per_helix: float, nominal_span_bp: float, cells,
        size: int = 8, return_metadata: bool = False):
    """Predict twist via effective pitch, elasticity and SNUPI calibration.

    ``indel_per_helix`` is the total edit in ``nominal_span_bp``.  Therefore
    +9 bases over 96 bp and +9 bases over 32 bp no longer share an erroneous
    lookup value: their effective pitches, and hence their twist densities,
    are different.
    """
    nominal_span = float(nominal_span_bp)
    if nominal_span <= 0.0:
        raise ValueError("Twist预测的名义选区长度必须大于0 bp。")
    cells = tuple(cells)
    selected = calibration_profile(cells, size)["profile"]
    x0, t0, k_minus, k_plus = selected
    effective_pitch = (SQUARE_NATIVE_PITCH_BP_PER_TURN *
                       (nominal_span+float(indel_per_helix)) /
                       nominal_span)
    native_elastic = _square_elastic_twist_per_base(
        SQUARE_NATIVE_PITCH_BP_PER_TURN, cells)
    pitch_elastic = _square_elastic_twist_per_base(effective_pitch, cells)
    # Anchor the literature elastic response at the frozen calibration x0.
    # This preserves the directly simulated native T0 exactly while retaining
    # the physical pitch response away from native.
    raw_twist = x0 + (pitch_elastic-native_elastic)
    coefficient = k_minus if raw_twist < x0 else k_plus
    calibrated = t0 + coefficient*(raw_twist-x0)
    if not return_metadata:
        return calibrated
    return calibrated, {
        "effective_pitch_bp_per_turn": effective_pitch,
        "polar_moment_nm4": _square_polar_moment_nm4(cells),
        "uncalibrated_twist_deg_per_base": raw_twist,
        "calibration_x0": x0,
        "calibration_T0": t0,
        "calibration_k_deletion": k_minus,
        "calibration_k_insertion": k_plus,
        "calibration_equivalent_indel_per_96bp": (
            float(indel_per_helix)*S8_R4X4C_CALIBRATION_SPAN_BP /
            nominal_span),
    }


def elastic_calibrated_indel_for_cross_section(
        twist_per_base: float, nominal_span_bp: float, cells,
        size: int = 8) -> float:
    """Invert the effective-pitch elastic plus SNUPI calibration chain."""
    nominal_span = float(nominal_span_bp)
    if nominal_span <= 0.0:
        raise ValueError("Twist反算的名义选区长度必须大于0 bp。")
    cells = tuple(cells)
    x0, t0, k_minus, k_plus = calibration_profile(
        cells, size)["profile"]
    coefficient = k_minus if float(twist_per_base) < t0 else k_plus
    raw_target = x0 + (float(twist_per_base)-t0) / max(coefficient, 1e-12)
    native_elastic = _square_elastic_twist_per_base(
        SQUARE_NATIVE_PITCH_BP_PER_TURN, cells)
    elastic_target = native_elastic + (raw_target-x0)
    polar_moment = _square_polar_moment_nm4(cells)
    at_low = _log_interpolate(
        polar_moment, SQUARE_ELASTIC_J_NM4, SQUARE_ELASTIC_TWIST_1033)
    at_high = _log_interpolate(
        polar_moment, SQUARE_ELASTIC_J_NM4, SQUARE_ELASTIC_TWIST_1067)
    response = at_high-at_low
    if abs(response) <= 1e-12:
        raise ValueError("当前Seed截面的弹性Twist响应接近零，无法反算indel。")
    effective_pitch = (SQUARE_ELASTIC_PITCH_LOW +
                       (elastic_target-at_low) /
                       response*(SQUARE_ELASTIC_PITCH_HIGH-
                                 SQUARE_ELASTIC_PITCH_LOW))
    return nominal_span*(effective_pitch /
                         SQUARE_NATIVE_PITCH_BP_PER_TURN-1.0)


def lattice_layers(symmetry: str):
    return LATTICE_LAYERS.get(str(symmetry), LATTICE_LAYERS["square_square_c4"])


def phase_residue_for_z2(growth_bp: int) -> int:
    """Required Z2 residue for an identical-layer Square SST design."""
    return (-int(growth_bp)) % SQUARE_REPEAT_BP


def compatible_z2_values(growth_bp: int, maximum: int = 400) -> List[int]:
    residue = phase_residue_for_z2(growth_bp)
    return list(range(residue, int(maximum)+1, SQUARE_REPEAT_BP))


def compatible_growth_values(spacer_bp: int, minimum: int = 64,
                             maximum: int = 400) -> List[int]:
    residue = (-int(spacer_bp)) % SQUARE_REPEAT_BP
    first = int(minimum)
    while first % SQUARE_REPEAT_BP != residue:
        first += SQUARE_DOMAIN_BP
    return list(range(first, int(maximum)+1, SQUARE_REPEAT_BP))


def preview_seed_partition(spacer_bp: int, sst_z1_bp: int,
                           sst_z3_bp: int, linked: bool,
                           total_bp: int = FIXED_SEED_TOTAL_BP) -> Dict[str, object]:
    """Partition the fixed 288-bp Seed for the side-view preview.

    Z2 always equals the nominal SST spacing.  Its two boundaries move in
    integral 8-bp domain steps.  When an 8-bp step cannot be split equally,
    linked layers use a deterministic alternating direction around the
    128/32/128 reference; independent layers keep the extra 8-bp support on
    the side that gives the better SST overlap.

    This is preview geometry only and deliberately does not modify the fixed
    Seed reference or the downstream SST/Scaffold generator.
    """
    spacing = int(spacer_bp)
    total = int(total_bp)
    first_sst = int(sst_z1_bp)
    second_sst = int(sst_z3_bp)
    if total != FIXED_SEED_TOTAL_BP:
        raise ValueError("当前预览只支持固定288-bp Seed。")
    geometry = centered_square_sst_geometry(
        first_sst, spacing, second_sst)
    z1, z2, z3 = map(int, geometry["seed_partition_lengths_bp"])
    phase_compatible = (
        not linked or
        (first_sst == second_sst and
         (first_sst+spacing) % SQUARE_REPEAT_BP == 0))
    return {
        "total_bp": total,
        "z1_bp": z1,
        "z2_bp": z2,
        "z3_bp": z3,
        "center_offset_bp": float(geometry[
            "envelope_center_offset_bp"]),
        "linked": bool(linked),
        "phase_compatible": bool(phase_compatible),
        "selection": (
            "shared Square overlap-optimized, Z2-centred-tiebreak "
            "complete-U geometry"),
        "sst_overlap_z1_bp": min(z1, first_sst),
        "sst_overlap_z3_bp": min(z3, second_sst),
        "sst_layer_ranges": geometry["layer_ranges"],
        "sst_scaffold_ranges": geometry["scaffold_ranges"],
        "sst_complementary_chain_ranges": geometry[
            "complement_ranges"],
        "seed_partition_ranges": geometry["seed_partition_ranges"],
        "coordinate_shift_bp": geometry["coordinate_shift_bp"],
    }


def phase_is_compatible(z1_bp: int, z2_bp: int, z3_bp: int) -> bool:
    return (
        z1_bp % SQUARE_DOMAIN_BP == 0 and
        z2_bp % SQUARE_DOMAIN_BP == 0 and
        z3_bp % SQUARE_DOMAIN_BP == 0 and
        z1_bp % SQUARE_REPEAT_BP == z3_bp % SQUARE_REPEAT_BP and
        z2_bp % SQUARE_REPEAT_BP == phase_residue_for_z2(z1_bp)
    )


def actual_z2_spacing_length(nominal_bp: int,
                             mean_indel_per_helix: float) -> float:
    """Return mean physical Z2 length after distributed insertions/deletions.

    Individual helices receive integral loop/skip edits; the calibrated input
    is their mean, so the shared SST-spacing/Seed-Z2 readout may be fractional.
    """
    return max(0.0, float(nominal_bp) + float(mean_indel_per_helix))


def seed_indel_limits_for_spacing(nominal_spacing_bp: int) -> tuple:
    """Return the feasible mean indel interval for one Seed helix.

    Z2 is divided into complete 8-bp domains and each domain may carry at
    most three insertions or deletions on a helix. Insertions retain the
    established global +10 structural cap; deletions are limited only by the
    number of available domains.
    """
    spacing = int(nominal_spacing_bp)
    if spacing < 0 or spacing % SQUARE_DOMAIN_BP:
        raise ValueError(
            "Nominal spacing must be a non-negative multiple of 8 bp.")
    domain_count = spacing // SQUARE_DOMAIN_BP
    domain_limit = float(MAX_SEED_DELETIONS_PER_DOMAIN * domain_count)
    return (-domain_limit, min(domain_limit, MAX_SEED_INSERTION_PER_HELIX))


def minimum_seed_deletion_per_helix(nominal_spacing_bp: int) -> float:
    """Return the spacing-dependent Seed deletion limit per helix.

    Z2/spacing is divided into complete 8-bp domains.  At most three bases
    may be deleted from each domain, and the structure generator distributes
    those deletions across all domains rather than concentrating them in one
    local segment.
    """
    return seed_indel_limits_for_spacing(nominal_spacing_bp)[0]


def maximum_seed_insertion_per_helix(nominal_spacing_bp: int) -> float:
    """Return the spacing-dependent Seed insertion limit per helix."""
    return seed_indel_limits_for_spacing(nominal_spacing_bp)[1]


def _status(level: str, title: str, detail: str) -> Dict[str, str]:
    return {"level": level, "title": title, "detail": detail}


def solve_square_bilayer(settings: SquareBilayerSettings) -> MoireProject:
    """Build a calibrated prediction for the selected bilayer and Square Seed.

    Nominal Z lengths remain independent design inputs.  Angle/period solves
    the distributed mean indel, which changes the mean physical Z2/spacing
    length while leaving its shared SST/Seed crossover coordinates intact.
    """
    cfg = replace(settings)
    layer_types = lattice_layers(cfg.lattice_symmetry)
    period_available = layer_types[0] == layer_types[1]
    # ``lattice_constant_nm`` is the legacy single-layer field consumed by
    # the geometric period formula. Keep it synchronized with the selected
    # first layer for direct core callers and restored projects as well.
    cfg.lattice_constant_nm = float(cfg.layer1_lattice_constant_nm)
    seed_size = int(getattr(cfg, "seed_cross_section_size", 8))
    seed_cells = tuple(tuple(map(int, cell)) for cell in
                       (getattr(cfg, "seed_cross_section_cells", None) or
                        default_seed_cells(seed_size)))
    section_calibration = calibration_profile(seed_cells, seed_size)
    sst_lengths = (
        cfg.sst_growth_bp_z1, cfg.spacer_bp_z2, cfg.sst_growth_bp_z3)
    if cfg.sst_growth_bp_z1 < 64 or cfg.sst_growth_bp_z3 < 64:
        raise ValueError("SST 1st layer和2nd layer至少需要64 bp。")
    if any(int(value) % SQUARE_DOMAIN_BP for value in sst_lengths):
        raise ValueError("SST 1st layer、spacing和2nd layer必须是8 bp的整数倍。")
    # Legacy Seed Z1/Z3 settings are ignored.  The physical two-layer Seed is
    # an immutable accepted template; only its overlap with the SST changes.
    cfg.growth_bp_z1 = 128
    cfg.growth_bp_z3 = 128
    if cfg.layers_design_sequence_identical and not phase_is_compatible(
            cfg.sst_growth_bp_z1, cfg.spacer_bp_z2,
            cfg.sst_growth_bp_z3):
        raise ValueError(
            "SST 1st layer、spacing、2nd layer不满足Square SST相位："
            "两层余数必须相同，且1st layer + spacing必须为32 bp整数倍。")

    baseline = CAPTURE_BASES * elastic_calibrated_twist_for_cross_section(
        0.0, CAPTURE_BASES, seed_cells, seed_size)
    zero_spacing = cfg.spacer_bp_z2 == 0
    if zero_spacing:
        # With no nominal Z2 domain there is nowhere to place an insertion or
        # deletion. Therefore the only physically valid coupled solution is
        # zero local Twist and zero mean indel; do not preserve a stale angle
        # or period from a previously selected non-zero spacing.
        cfg.target_mode = "angle"
        cfg.target_angle_deg = 0.0
        cfg.mean_indel_per_helix = 0.0
    if cfg.target_mode == "period" and period_available:
        cfg.target_angle_deg = angle_from_moire_period(
            cfg.target_period_nm, cfg.lattice_constant_nm)
    elif cfg.target_mode == "period":
        cfg.target_mode = "angle"

    if cfg.target_mode in ("angle", "period"):
        requested_local = float(cfg.target_angle_deg)
        if cfg.target_definition == "experimental_total":
            requested_local -= baseline
        if cfg.spacer_bp_z2 > 0:
            requested_rate = requested_local / cfg.spacer_bp_z2
            cfg.mean_indel_per_helix = (
                elastic_calibrated_indel_for_cross_section(
                    requested_rate, cfg.spacer_bp_z2,
                    seed_cells, seed_size))
            twist_per_base, twist_model = (
                elastic_calibrated_twist_for_cross_section(
                    cfg.mean_indel_per_helix, cfg.spacer_bp_z2,
                    seed_cells, seed_size, return_metadata=True))
            predicted_local = cfg.spacer_bp_z2 * twist_per_base
        else:
            cfg.mean_indel_per_helix = 0.0
            twist_per_base = 0.0
            twist_model = {
                "effective_pitch_bp_per_turn":
                    SQUARE_NATIVE_PITCH_BP_PER_TURN,
                "polar_moment_nm4": _square_polar_moment_nm4(seed_cells),
                "uncalibrated_twist_deg_per_base": 0.0,
                "calibration_equivalent_indel_per_96bp": 0.0,
            }
            predicted_local = 0.0
    else:
        if cfg.spacer_bp_z2 <= 0:
            twist_per_base = 0.0
            twist_model = {
                "effective_pitch_bp_per_turn":
                    SQUARE_NATIVE_PITCH_BP_PER_TURN,
                "polar_moment_nm4": _square_polar_moment_nm4(seed_cells),
                "uncalibrated_twist_deg_per_base": 0.0,
                "calibration_equivalent_indel_per_96bp": 0.0,
            }
        else:
            twist_per_base, twist_model = (
                elastic_calibrated_twist_for_cross_section(
                    cfg.mean_indel_per_helix, cfg.spacer_bp_z2,
                    seed_cells, seed_size, return_metadata=True))
        predicted_local = cfg.spacer_bp_z2 * twist_per_base
    predicted_total = (0.0 if zero_spacing else predicted_local + baseline)
    reported_angle = (0.0 if zero_spacing else
        predicted_total if cfg.target_definition == "experimental_total"
        else predicted_local)
    cfg.target_angle_deg = reported_angle
    cfg.target_period_nm = (moire_period_from_angle(
        reported_angle, cfg.lattice_constant_nm) if period_available else 0.0)
    actual_z2_bp = actual_z2_spacing_length(
        cfg.spacer_bp_z2, cfg.mean_indel_per_helix)
    preview_partition = preview_seed_partition(
        cfg.spacer_bp_z2, cfg.sst_growth_bp_z1,
        cfg.sst_growth_bp_z3, cfg.layers_design_sequence_identical)

    validation: List[Dict[str, str]] = []
    symmetry_label = {
        "square_square_c4": "Square–Square",
        "kagome_kagome": "Kagome–Kagome",
        "square_kagome": "Square–Kagome",
    }.get(cfg.lattice_symmetry, cfg.lattice_symmetry)
    calibration_neighbors = section_calibration["neighbors"]
    validation.append(_status(
        "pass" if section_calibration["exact"] else "info",
        "Square Seed截面校准",
        "%d×%d Square网格，%d根helix；%s。" % (
            seed_size, seed_size, len(seed_cells),
            ("匹配冻结校准节点%s" % calibration_neighbors[0][0]
             if section_calibration["exact"] else
             "按相邻冻结校准节点进行拓扑插值"))))
    if section_calibration["extrapolated"]:
        validation.append(_status(
            "warning", "Seed截面校准范围",
            "当前截面与冻结校准节点距离较大；可继续比较趋势，但数值应在"
            "补充力学模拟后用于最终实验设计。"))
    validation.append(_status(
        "pass" if cfg.lattice_symmetry == "square_square_c4" else "info",
        "双层点阵对称性",
        "%s；1st a=%.1f nm，2nd a=%.1f nm。%s" % (
            symmetry_label, cfg.layer1_lattice_constant_nm,
            cfg.layer2_lattice_constant_nm,
            ("两层点阵不同，Moiré period不定义，仅使用Twist。"
             if not period_available else "使用同型点阵公式计算period。"))))
    if cfg.layers_design_sequence_identical:
        validation.append(_status(
            "pass", "Z1/Z2/Z3相位",
            "双层设计与序列一致；1st=%d bp、spacing/Seed Z2=%d bp、"
            "2nd=%d bp；满足8 bp domain与32 bp重复相位。" %
            (cfg.sst_growth_bp_z1, cfg.spacer_bp_z2,
             cfg.sst_growth_bp_z3)))
    else:
        validation.append(_status(
            "info", "独立双层设计",
            "双层设计或序列不一致；SST三项仅限制为8 bp整数倍，"
            "不施加彼此的32 bp相位关系。"))
    seed_actual_z1 = seed_actual_z3 = 128
    validation.append(_status(
        "info", "Seed支撑与routing",
        "Seed使用固定2L范本；scaffold、staple、capture core、nick和seam"
        "均不修改，只计算与当前SST的实际重叠。"))
    validation.append(_status(
        "pass", "校准角度",
        "当前Seed截面native %.3f°/base；16 bp capture基线=%.2f°。" %
        (elastic_calibrated_twist_for_cross_section(
            0.0, CAPTURE_BASES, seed_cells, seed_size), baseline)))
    calibration_low, calibration_high = S8_R4X4C_VALIDATED_INDEL_RANGE
    equivalent_indel = float(twist_model.get(
        "calibration_equivalent_indel_per_96bp", 0.0))
    calibration_domain_exceeded = not (
        calibration_low - 1e-9 <= equivalent_indel <=
        calibration_high + 1e-9)
    if calibration_domain_exceeded:
        validation.append(_status(
            "warning", "校准外推",
            "当前indel密度折算为96-bp校准区后为%.2f indel/helix，超出"
            "已收敛模拟的−10到+6范围；当前弹性模型后的SNUPI校正属于"
            "校准域外推，+8和+10 insertion模拟未收敛且未参与拟合。" %
            equivalent_indel))
    minimum_deletion, maximum_insertion = seed_indel_limits_for_spacing(
        cfg.spacer_bp_z2)
    insertion_limit_exceeded = (
        cfg.mean_indel_per_helix > maximum_insertion + 1e-9)
    deletion_limit_exceeded = (
        cfg.mean_indel_per_helix < minimum_deletion - 1e-9)
    indel_limit_exceeded = (
        insertion_limit_exceeded or deletion_limit_exceeded)
    # Keep feasibility in the live prediction without raising while Twist
    # (1.1) and spacing (1.2) are still being edited.  The final combination
    # is blocked when the user accepts 1.2.
    if zero_spacing:
        validation.append(_status(
            "pass", "Z2 = 0 parameter linkage",
            "A 0-bp spacing contains no 8-bp domain, so Twist and mean "
            "insertion/deletion are both fixed at 0."))
    if cfg.capture_mode != "fully_cooperative":
        validation.append(_status(
            "error", "Capture模式",
            "第一版只将fully cooperative capture作为可实验设计。"))
    else:
        validation.append(_status(
            "pass", "Capture协同性",
            "capture-0/capture-1双domain配对沿z方向连续。"))
    validation.append(_status(
        "info", "参数联动",
        "名义Z2不由目标角度反算；Twist所需增删使实际Z2/spacing为%.1f bp，"
        "%s" % (actual_z2_bp,
                  ("Moiré period由同型点阵公式独立计算。"
                   if period_available else
                   "Square–Kagome不提供period。"))))

    prediction = {
        "requested_angle_deg": float(cfg.target_angle_deg),
        "requested_period_nm": float(cfg.target_period_nm),
        "predicted_local_surface_angle_deg": predicted_local,
        "predicted_experimental_total_deg": predicted_total,
        "reported_angle_deg": reported_angle,
        "predicted_moire_period_nm": (
            cfg.target_period_nm if period_available else None),
        "period_available": period_available,
        "layer_lattice_types": list(layer_types),
        "layer_lattice_constants_nm": [
            cfg.layer1_lattice_constant_nm,
            cfg.layer2_lattice_constant_nm],
        "twist_deg_per_base": twist_per_base,
        "uncalibrated_twist_deg_per_base": twist_model.get(
            "uncalibrated_twist_deg_per_base"),
        "effective_pitch_bp_per_turn": twist_model.get(
            "effective_pitch_bp_per_turn"),
        "polar_moment_nm4": twist_model.get("polar_moment_nm4"),
        "calibration_equivalent_indel_per_96bp": equivalent_indel,
        "twist_prediction_model": "effective-pitch-elastic-then-SNUPI",
        "mean_indel_per_helix": cfg.mean_indel_per_helix,
        "maximum_seed_insertion_per_helix":
            maximum_insertion,
        "maximum_seed_deletions_per_domain":
            MAX_SEED_DELETIONS_PER_DOMAIN,
        "minimum_seed_deletion_per_helix": minimum_deletion,
        "seed_insertion_limit_exceeded": insertion_limit_exceeded,
        "seed_deletion_limit_exceeded": deletion_limit_exceeded,
        "seed_indel_limit_exceeded": indel_limit_exceeded,
        "seed_twist_calibration_domain_exceeded":
            calibration_domain_exceeded,
        "nominal_z2_spacing_bp": cfg.spacer_bp_z2,
        "actual_z2_spacing_bp": actual_z2_bp,
        "nominal_interlayer_spacing_nm": cfg.spacer_bp_z2 * DNA_RISE_NM,
        "actual_interlayer_spacing_nm": actual_z2_bp * DNA_RISE_NM,
        "preview_seed_partition": preview_partition,
        "phase": {
            "z1_z3_residue_mod32": (
                cfg.sst_growth_bp_z1 % SQUARE_REPEAT_BP),
            "z2_residue_mod32": cfg.spacer_bp_z2 % SQUARE_REPEAT_BP,
            "z2_options": compatible_z2_values(
                cfg.sst_growth_bp_z1, maximum=400),
        },
        "calibration": {
            "version": S8_R4X4C_CALIBRATION_VERSION,
            "section": calibration_neighbors[0][0],
            "seed_cross_section_cells": [list(cell) for cell in seed_cells],
            "seed_cross_section_size": seed_size,
            "exact_node": section_calibration["exact"],
            "neighbors": [list(item) for item in calibration_neighbors],
            "descriptor": list(section_calibration["descriptor"]),
            "points": [list(point) for point in S8_R4X4C_TWIST_POINTS],
            "validated_indel_per_helix_range":
                list(S8_R4X4C_VALIDATED_INDEL_RANGE),
            "failed_indel_per_helix_points":
                list(S8_R4X4C_FAILED_INDEL_POINTS),
            "branch_rmse_deg_per_base": dict(S8_R4X4C_BRANCH_RMSE),
            "failed_points_excluded_from_fit": True,
            "source_calibration_span_bp":
                S8_R4X4C_CALIBRATION_SPAN_BP,
            "prediction_chain":
                "indel density -> effective pitch -> elastic twist -> "
                "topology-aware SNUPI calibration",
            "native_twist_deg_per_base": S8_R4X4C_NATIVE_TWIST_DEG_PER_BASE,
            "native_32bp_angle_deg": 32*S8_R4X4C_NATIVE_TWIST_DEG_PER_BASE,
            "capture_16bp_angle_deg": baseline,
        },
        "confidence": (
            "direct calibrated node" if section_calibration["exact"] else
            "topology interpolation; outside calibrated section range" if
            section_calibration["extrapolated"] else
            "topology interpolation within neighboring calibrated nodes"),
    }
    seed_plan = {
        "template": cfg.seed_template,
        "lattice": "square",
        "outer_cross_section": [seed_size, seed_size],
        "cross_section_cells": [list(cell) for cell in seed_cells],
        "occupied_helices": len(seed_cells),
        "segments": [
            {"name": "Seed Z1", "role": "fixed layer 1 support",
             "bp": cfg.growth_bp_z1, "routing_bp": seed_actual_z1},
            {"name": "Z2", "role": "spacing + relative twist",
             "bp": cfg.spacer_bp_z2, "actual_bp": actual_z2_bp},
            {"name": "Seed Z3", "role": "fixed layer 2 support",
             "bp": cfg.growth_bp_z3, "routing_bp": seed_actual_z3},
        ],
        "sst_segments": [
            {"name": "SST 1st layer", "bp": cfg.sst_growth_bp_z1},
            {"name": "SST spacing / Seed Z2", "bp": cfg.spacer_bp_z2,
             "actual_bp": actual_z2_bp},
            {"name": "SST 2nd layer", "bp": cfg.sst_growth_bp_z3},
        ],
        "routing_status": (
            "staged SST/scaffold/staple generator available" if
            cfg.lattice_symmetry == "square_square_c4" and
            set(seed_cells) == set(default_seed_cells(8)) else
            "prediction only; downstream structure generator not enabled"),
    }
    capture_plan = {
        "mode": cfg.capture_mode,
        "surfaces": ["north", "south"],
        "inactive_surfaces": ["east", "west"],
        "pair": ["capture-0", "capture-1"],
        "capture_length_nt": CAPTURE_BASES,
        "axial_continuity_required": True,
        "domain_length_nt": SQUARE_DOMAIN_BP,
        "sst_length_nt": SQUARE_REPEAT_BP,
        "sst_sets": ["SST-a", "SST-a*"],
    }
    return MoireProject(
        settings=cfg,
        prediction=prediction,
        validation=validation,
        capture_plan=capture_plan,
        seed_plan=seed_plan,
    )
