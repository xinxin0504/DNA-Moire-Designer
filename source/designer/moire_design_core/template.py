"""Reference Seed-S template handling for the prototype."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from .models import MoireProject


REFERENCE_SEED_CANDIDATES = (
    Path(__file__).with_name("resources") / "Square_Seed_2L_newtemplate.json",
    Path(__file__).with_name("resources") / "Square_Seed_2L_original.json",
)


def reference_seed_path() -> Optional[Path]:
    for candidate in REFERENCE_SEED_CANDIDATES:
        if candidate.is_file():
            return candidate
    return None


def export_reference_seed(project: MoireProject, filename: str) -> Path:
    source = reference_seed_path()
    if source is None:
        raise FileNotFoundError("找不到论文 Seed S 的参考 cadnano JSON。")
    payload = json.loads(source.read_text(encoding="utf-8"))
    target = Path(filename).expanduser().resolve()
    payload["name"] = target.name
    payload["moire_metadata"] = {
        "schema_version": project.schema_version,
        "project_name": project.settings.project_name,
        "target_angle_deg": project.settings.target_angle_deg,
        "predicted_local_surface_angle_deg": project.prediction[
            "predicted_local_surface_angle_deg"],
        "z2_bp": project.settings.spacer_bp_z2,
        "prototype_note": (
            "Reference paper Seed-S routing. Z2 metadata does not yet rewrite "
            "strand indices in this first prototype."),
    }
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    return target
