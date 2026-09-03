"""Serializable data models used by both the standalone and cadnano UI."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields
from datetime import datetime, timezone
from typing import Any, Dict, List


SCHEMA_VERSION = "1.2"
CORE_VERSION = "0.4.0-sst-seed-support-split"


def _default_seed_cross_section():
    return [[row, col] for row in range(8) for col in range(8)
            if not (2 <= row <= 5 and 2 <= col <= 5)]


@dataclass
class SquareBilayerSettings:
    project_name: str = "square_moire_bilayer"
    interface_language: str = "en"
    target_mode: str = "angle"
    target_definition: str = "local_surface"
    target_angle_deg: float = 3.2967555036483183
    target_period_nm: float = 48.669158335514105
    lattice_constant_nm: float = 2.8
    lattice_context: str = "solution_cryo"
    lattice_symmetry: str = "square_square_c4"
    layer1_lattice_constant_nm: float = 2.8
    layer2_lattice_constant_nm: float = 2.8
    seed_cross_section_size: int = 8
    seed_cross_section_cells: List[List[int]] = field(
        default_factory=_default_seed_cross_section)
    # Compatibility fields for legacy project JSON.  The physical Seed is a
    # fixed 128/spacing/128 accepted template; callers must not use these to
    # resize or reroute it.
    growth_bp_z1: int = 128
    spacer_bp_z2: int = 32
    growth_bp_z3: int = 128
    # SST layer lengths retain the validated 8-bp domain/32-bp phase rules.
    # spacer_bp_z2 is the single shared SST spacing / Seed Z2 value.
    sst_growth_bp_z1: int = 128
    sst_growth_bp_z3: int = 128
    capture_spacing_bp: int = 16
    capture_baseline_deg: float = 1.6483777518241591
    snap_spacer_to_bp: int = 8
    auto_solve_spacer: bool = False
    mean_indel_per_helix: float = 0.0
    layers_design_sequence_identical: bool = True
    seed_template: str = "Seed S(F), 48-helix square frame"
    capture_mode: str = "fully_cooperative"


@dataclass
class MoireProject:
    settings: SquareBilayerSettings
    prediction: Dict[str, Any]
    validation: List[Dict[str, str]]
    capture_plan: Dict[str, Any]
    seed_plan: Dict[str, Any]
    measurements: List[Dict[str, Any]] = field(default_factory=list)
    schema_version: str = SCHEMA_VERSION
    core_version: str = CORE_VERSION
    created_at: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat())
    modified_at: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, payload: Dict[str, Any]) -> "MoireProject":
        raw_settings = dict(payload["settings"])
        # The Windows product is English-only.  Ignore language preferences
        # retained by projects created with an older bilingual release.
        raw_settings["interface_language"] = "en"
        # Version 1.1 used growth_bp_z1/z3 for both Seed and SST.  Initializing
        # the new SST fields from those values keeps every legacy project
        # geometrically unchanged on first load.
        raw_settings.setdefault(
            "sst_growth_bp_z1", raw_settings.get("growth_bp_z1", 128))
        raw_settings.setdefault(
            "sst_growth_bp_z3", raw_settings.get("growth_bp_z3", 128))
        raw_settings.setdefault(
            "layer1_lattice_constant_nm",
            raw_settings.get("lattice_constant_nm", 2.8))
        raw_settings.setdefault(
            "layer2_lattice_constant_nm",
            raw_settings.get("lattice_constant_nm", 2.8))
        allowed = {item.name for item in fields(SquareBilayerSettings)}
        settings = SquareBilayerSettings(**{
            key: value for key, value in raw_settings.items()
            if key in allowed})
        known = dict(payload)
        known.pop("settings", None)
        return cls(settings=settings, **known)
