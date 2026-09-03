"""Re-audit the saved Z2 review JSONs and write a compact review index."""

from __future__ import annotations

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from moire_design_core.structure import validate_structure


def main() -> None:
    root = (Path(__file__).resolve().parents[1] / "review_outputs" /
            "Z2_identical_vs_independent")
    summary_path = root / "validation_summary.json"
    report = json.loads(summary_path.read_text(encoding="utf-8"))
    compact_cases = []
    for row in report["cases"]:
        z1 = int(row["z1_bp"])
        z2 = int(row["z2_bp"])
        z3 = int(row["z3_bp"])
        group = ("identical_linked" if
                 row["layers_design_sequence_identical"] else
                 "different_independent")
        label = f"Z1_{z1}_Z2_{z2}_Z3_{z3}"
        final_path = root / group / label / f"{label}_complete.json"
        item = {
            "case": label,
            "mode": group,
            "z1_bp": z1,
            "z2_bp": z2,
            "z3_bp": z3,
            "identical_layer_phase_compatible": bool(
                row["identical_layer_phase_compatible"]),
            "status": row["status"],
            "seed_scaffold_lengths": row.get("seed_scaffold_lengths", []),
        }
        if final_path.exists():
            validation = validate_structure(
                str(final_path), require_staples=True)
            audit = validation.get("seed_scaffold_crossover_audit") or {}
            item.update({
                "complete_json": str(final_path),
                "valid": bool(validation["valid"]),
                "errors": list(validation["errors"]),
                "capture_bridge_component_count": validation.get(
                    "capture_bridge_component_count"),
                "normal_staple_length_range_nt": [
                    validation.get("minimum_normal_staple_length"),
                    validation.get("maximum_normal_staple_length")],
                "invalid_scaffold_crossovers": audit.get("invalid", 0),
                "minimum_different_direction_index_gap": audit.get(
                    "minimum_different_direction_index_gap"),
            })
            row["capture_bridge_component_count"] = item[
                "capture_bridge_component_count"]
            row["final_valid"] = item["valid"]
            row["final_errors"] = item["errors"]
        else:
            error = str(row.get("error", ""))
            item["rejection_reason"] = next(
                (line.strip() for line in reversed(error.splitlines())
                 if line.strip().startswith(("RuntimeError:",
                                             "ValueError:"))),
                error[-400:])
        compact_cases.append(item)
    report["complete_count"] = sum(
        row["status"] == "complete" for row in report["cases"])
    report["rejected_count"] = sum(
        row["status"] == "rejected" for row in report["cases"])
    summary_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    (root / "review_index.json").write_text(json.dumps({
        "scope": report["scope"],
        "identical_rule": report["identical_rule"],
        "different_rule": report["different_rule"],
        "complete_count": report["complete_count"],
        "rejected_count": report["rejected_count"],
        "cases": compact_cases,
    }, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
