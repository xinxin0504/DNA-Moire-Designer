"""Generate Kagome Seed shrink/growth JSONs for visual review."""

from __future__ import annotations

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from moire_design_core.structure import (
    finalize_structure,
    generate_scaffold_review,
    validate_structure,
    write_shifted_sst,
)


SHRINK_CASES = ((64, 64), (96, 96), (104, 104), (112, 112),
                (120, 120), (96, 120), (120, 96))
GROWTH_CASES = (
    (136, 136), (136, 144), (136, 152), (136, 160),
    (144, 136), (144, 144), (144, 152), (144, 160),
    (152, 136), (152, 144), (152, 152), (152, 160),
    (160, 136), (160, 144), (160, 152), (160, 160),
)


def generate(root: Path):
    root.mkdir(parents=True, exist_ok=True)
    rows = []
    for group, cases in (("shrink", SHRINK_CASES),
                         ("growth", GROWTH_CASES)):
        for z1, z3 in cases:
            label = f"Kagome_Seed_Z1_{z1}_Z2_32_Z3_{z3}"
            folder = root / group / label
            folder.mkdir(parents=True, exist_ok=True)
            sst = folder / f"{label}_sst.json"
            scaffold = folder / f"{label}_scaffold.json"
            complete = folder / f"{label}_complete.json"
            write_shifted_sst(
                str(sst), 128, 32, 128, z1, z3, 32, "kagome")
            review = generate_scaffold_review(str(scaffold), str(sst))
            finalize_structure(str(scaffold), str(complete))
            validation = validate_structure(
                str(complete), require_staples=True)
            payload = json.loads(complete.read_text(encoding="utf-8"))
            layout = payload["moire_structure_metadata"][
                "variable_length_layout"]
            audit = validation.get("seed_scaffold_crossover_audit") or {}
            rows.append({
                "group": group,
                "z1_bp": z1, "z2_bp": 32, "z3_bp": z3,
                "valid": bool(validation["valid"]),
                "errors": list(validation["errors"]),
                "complete_json": str(complete),
                "seed_scaffold_lengths": review["seed_scaffold_lengths"],
                "seed_staple_physical_range": layout.get(
                    "seed_staple_physical_range"),
                "capture_pair_count_by_layer": layout.get(
                    "pair_count_by_layer"),
                "capture_bridge_component_count": validation.get(
                    "capture_bridge_component_count"),
                "omitted_nonoverlap_capture_positions": layout.get(
                    "omitted_nonoverlap_kagome_capture_positions", []),
                "invalid_scaffold_crossovers": audit.get("invalid", 0),
                "invalid_short_staples": validation.get(
                    "invalid_short_staple_count"),
                "protected_short_staples": validation.get(
                    "protected_short_staple_count"),
                "normal_staple_length_range_nt": [
                    validation.get("minimum_normal_staple_length"),
                    validation.get("maximum_normal_staple_length")],
            })
            print(label, validation["valid"], flush=True)
    report = {
        "scope": ("Kagome SST 128/32/128; Square 8x8 Seed; "
                  "Moiré Designer only"),
        "shrink_policy": (
            "freeze Seed scaffold and interior staples; trim Seed edge "
            "staples and omit only Kagome anchors outside actual overlap"),
        "growth_policy": (
            "freeze reviewed interior; use legal 10/11-bp Seed scaffold "
            "edge growth and capacity-minimal 2/3-scaffold partition"),
        "case_count": len(rows),
        "valid_count": sum(item["valid"] for item in rows),
        "cases": rows,
    }
    (root / "review_index.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return report


if __name__ == "__main__":
    destination = (Path(__file__).resolve().parents[1] / "review_outputs" /
                   "Kagome_Seed_shrink_growth_review")
    print(json.dumps(generate(destination), ensure_ascii=False, indent=2))
