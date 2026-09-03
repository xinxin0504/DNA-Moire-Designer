"""Project persistence and deterministic export helpers."""

from __future__ import annotations

import csv
import json
import math
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from .models import MoireProject


def save_project(project: MoireProject, filename: str) -> Path:
    target = Path(filename).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    project = replace(
        project, modified_at=datetime.now(timezone.utc).isoformat())
    target.write_text(
        json.dumps(project.to_dict(), ensure_ascii=False, indent=2),
        encoding="utf-8")
    return target


def load_project(filename: str) -> MoireProject:
    payload = json.loads(Path(filename).expanduser().read_text(encoding="utf-8"))
    return MoireProject.from_dict(payload)


def export_capture_map(project: MoireProject, filename: str) -> Path:
    target = Path(filename).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "layer", "growth_segment", "surface", "capture_group",
            "capture_length_nt", "sst_set", "status"])
        for layer, segment in ((1, "Z1"), (2, "Z3")):
            for surface in project.capture_plan["surfaces"]:
                writer.writerow([
                    layer, segment, surface, "capture-0 + capture-1",
                    project.capture_plan["capture_length_nt"],
                    "SST-a / SST-a*", "sequence assignment pending"])
    try:
        from moire_designer.i18n import localize_csv
        localize_csv(target, getattr(
            project.settings, "interface_language", "en"))
    except Exception:
        pass
    return target


def add_measurement(project: MoireProject, angle_deg: float,
                    period_nm: Optional[float] = None,
                    source: str = "manual") -> None:
    predicted = float(project.prediction["reported_angle_deg"])
    project.measurements.append({
        "source": source,
        "angle_deg": float(angle_deg),
        "period_nm": (None if period_nm is None or math.isnan(period_nm)
                      else float(period_nm)),
        "prediction_error_deg": float(angle_deg)-predicted,
        "recorded_at": datetime.now(timezone.utc).isoformat(),
    })
