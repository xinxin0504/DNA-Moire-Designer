"""Canonical 8x8/4x4-pore Kagome stage-2 regression."""

import json
import tempfile
import unittest
from collections import Counter
from pathlib import Path

from moire_design_core.structure import (
    build_shifted_sst_payload,
    finalize_structure,
    generate_scaffold_review,
    validate_structure,
)
from moire_design_core.sequence_workflow_worker import (
    _sequence_sheets,
    analyze,
)


class KagomeFrozenFinalTests(unittest.TestCase):

    def test_two_layer_kagome_uses_2l_seed_and_capture_only_catalogue(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            sst_path = folder / "kagome_sst.json"
            scaffold_path = folder / "kagome_scaffold.json"
            final_path = folder / "kagome_final.json"
            sst = build_shifted_sst_payload(
                sst_path.name, reserve_capture_gaps=True,
                z1_bp=128, z2_bp=32, z3_bp=128,
                capture_extension_length_nt=32,
                lattice_type="kagome",
                layers_design_sequence_identical=True)
            sst_path.write_text(json.dumps(sst), encoding="utf-8")
            scaffold_report = generate_scaffold_review(
                str(scaffold_path), str(sst_path))
            self.assertEqual(
                scaffold_report["seed_scaffold_lengths"], [7300, 7336])
            self.assertEqual(
                scaffold_report["seed_routing_lengths"], [128, 128])
            finalize_report = finalize_structure(
                str(scaffold_path), str(final_path))
            self.assertEqual(
                finalize_report["capture_half_crossover_fallback_positions"],
                [])

            payload = json.loads(final_path.read_text())
            metadata = payload["moire_structure_metadata"]
            self.assertEqual(
                metadata["seed_routing_source"],
                "Square_Seed_2L_newtemplate.json")
            self.assertIn(
                "first two layers only",
                metadata["kagome_capture_topology_source"])
            self.assertEqual(payload["num_bases"], 544)
            self.assertEqual({
                (len(row["scaf"]), len(row["stap"]),
                 len(row["loop"]), len(row["skip"]))
                for row in payload["vstrands"]}, {(544, 544, 544, 544)})
            sequence_report = analyze(str(final_path))
            self.assertGreaterEqual(
                sequence_report["summary"]["seed_scaffold"]["count"], 1)
            # Sequence export is intentionally strict: an unassigned stage-2
            # structure must not fabricate SST input bases.  Assignment and
            # capture-manifest coverage are exercised by the sequence-workflow
            # tests after accepted inputs have been added.
            with self.assertRaisesRegex(
                    ValueError, "lacks an accepted SST sublattice input base"):
                _sequence_sheets(payload)
            report = validate_structure(str(final_path), require_staples=True)
            self.assertTrue(report["valid"], report["errors"])
            self.assertEqual(report["seed_scaffold_lengths"], [7300, 7336])
            self.assertEqual(report["capture_bridge_component_count"], 48)
            self.assertGreaterEqual(report["capture_color_count"], 8)
            self.assertGreaterEqual(report["minimum_normal_staple_length"], 21)
            self.assertEqual(report["maximum_normal_staple_length"], 61)
            self.assertTrue(any("固定2L Seed不重新break" in warning
                                for warning in report["warnings"]))


if __name__ == "__main__":
    unittest.main()
