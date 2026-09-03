"""Generate audited Square Moiré Z2 review designs.

This is an explicit review utility, not a routing implementation.  It drives
the same public staged API used by the application and keeps rejected cases
separate from complete, strictly validated structures.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

from moire_design_core.calculations import phase_is_compatible
from moire_design_core.structure import (
    finalize_structure,
    generate_scaffold_review,
    validate_sst,
    validate_structure,
    write_shifted_sst,
)


IDENTICAL_CASES = (
    (96, 32, 96),
    (96, 64, 96),
    (112, 48, 112),
    (120, 40, 120),
    (128, 64, 128),
)

INDEPENDENT_CASES = (
    (96, 40, 112),
    (96, 48, 128),
    (104, 32, 120),
    (112, 24, 128),
    (120, 48, 104),
)


def _case_summary(case, identical, sst_validation=None,
                  scaffold_report=None, scaffold_validation=None,
                  final_validation=None, error=None):
    z1, z2, z3 = case
    summary = {
        "z1_bp": z1,
        "z2_bp": z2,
        "z3_bp": z3,
        "layers_design_sequence_identical": bool(identical),
        "identical_layer_phase_compatible": phase_is_compatible(
            z1, z2, z3),
        "expected_phase_policy": (
            "linked_32bp_phase" if identical else
            "independent_8bp_domains"),
        "status": "rejected" if error else "complete",
    }
    if error:
        summary["error"] = str(error)
    if sst_validation:
        summary["sst_valid"] = bool(sst_validation["valid"])
        summary["sst_errors"] = list(sst_validation["errors"])
    if scaffold_report:
        summary["seed_scaffold_lengths"] = list(
            scaffold_report["seed_scaffold_lengths"])
        summary["seed_scaffold_count"] = len(
            scaffold_report["seed_scaffold_lengths"])
    if scaffold_validation:
        summary["scaffold_valid"] = bool(scaffold_validation["valid"])
        summary["scaffold_errors"] = list(scaffold_validation["errors"])
        summary["scaffold_crossover_audit"] = scaffold_validation.get(
            "seed_scaffold_crossover_audit")
    if final_validation:
        summary["final_valid"] = bool(final_validation["valid"])
        summary["final_errors"] = list(final_validation["errors"])
        summary["capture_bridge_component_count"] = final_validation.get(
            "capture_bridge_component_count")
        summary["normal_staple_length_range_nt"] = [
            final_validation.get("minimum_normal_staple_length"),
            final_validation.get("maximum_normal_staple_length")]
        summary["capture_core_source"] = final_validation.get(
            "capture_core_source", "immutable_template")
        summary["capture_core_length_evaluated"] = False
    return summary


def generate(root: Path):
    root.mkdir(parents=True, exist_ok=True)
    rejected = root / "rejected"
    rejected.mkdir(exist_ok=True)
    summaries = []
    for identical, cases, group_name in (
            (True, IDENTICAL_CASES, "identical_linked"),
            (False, INDEPENDENT_CASES, "different_independent")):
        group = root / group_name
        group.mkdir(exist_ok=True)
        for case in cases:
            z1, z2, z3 = case
            label = "Z1_%d_Z2_%d_Z3_%d" % case
            case_dir = group / label
            case_dir.mkdir(exist_ok=True)
            sst_path = case_dir / (label + "_sst.json")
            scaffold_path = case_dir / (label + "_scaffold.json")
            final_path = case_dir / (label + "_complete.json")
            sst_validation = scaffold_report = scaffold_validation = None
            final_validation = None
            error = None
            try:
                write_shifted_sst(
                    str(sst_path), z1, z2, z3, z1, z3, 16,
                    "square", "s8_r4x4", identical)
                sst_validation = validate_sst(str(sst_path))
                if not sst_validation["valid"]:
                    raise RuntimeError("SST validation: %s" %
                                       sst_validation["errors"])
                scaffold_report = generate_scaffold_review(
                    str(scaffold_path), str(sst_path))
                scaffold_validation = validate_structure(
                    str(scaffold_path))
                if not scaffold_validation["valid"]:
                    raise RuntimeError("Scaffold validation: %s" %
                                       scaffold_validation["errors"])
                finalize_structure(str(scaffold_path), str(final_path))
                final_validation = validate_structure(
                    str(final_path), require_staples=True)
                if not final_validation["valid"]:
                    raise RuntimeError("Final validation: %s" %
                                       final_validation["errors"])
            except Exception as caught:
                error = caught
                target = rejected / (group_name + "_" + label)
                if target.exists():
                    shutil.rmtree(target)
                shutil.move(str(case_dir), str(target))
            summary = _case_summary(
                case, identical, sst_validation, scaffold_report,
                scaffold_validation, final_validation, error)
            summaries.append(summary)
            print(json.dumps(summary, ensure_ascii=False), flush=True)
    report = {
        "scope": "Moiré Designer Square 8x8 Seed only; caDNAno unchanged",
        "identical_rule": (
            "Z1%32 == Z3%32 and Z2%32 == (-Z1)%32"),
        "different_rule": (
            "Z1/Z2/Z3 are independent 8-bp domains; no shared 32-bp phase"),
        "complete_count": sum(item["status"] == "complete"
                              for item in summaries),
        "rejected_count": sum(item["status"] == "rejected"
                              for item in summaries),
        "cases": summaries,
    }
    (root / "validation_summary.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return report


if __name__ == "__main__":
    destination = (Path(__file__).resolve().parents[1] / "review_outputs" /
                   "Z2_identical_vs_independent")
    report = generate(destination)
    print(json.dumps({
        "output": str(destination),
        "complete_count": report["complete_count"],
        "rejected_count": report["rejected_count"],
    }, ensure_ascii=False, indent=2))
