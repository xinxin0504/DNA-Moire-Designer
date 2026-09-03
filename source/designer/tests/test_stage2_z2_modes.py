import json
from collections import Counter
from pathlib import Path
import tempfile
import unittest

from moire_design_core.calculations import phase_is_compatible
from moire_design_core.structure import (
    build_output_sst_snapshot_payload,
    generate_scaffold_review,
    payload_to_internal_numbering,
    validate_structure,
    validate_sst,
    write_shifted_sst,
)
from moire_design_core.structure_worker import finalize


class Stage2Z2ModeTests(unittest.TestCase):

    def test_centered_z2_overlap_survives_all_structure_stages(self):
        cases = (
            (192, 64, 128, [112, 112], [7, 7]),
            (400, 160, 80, [64, 64], [4, 4]),
        )
        for first, spacing, second, expected_overlap, expected_columns in cases:
            with self.subTest(values=(first, spacing, second)), \
                    tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                sst = write_shifted_sst(
                    str(root / "sst.json"), first, spacing, second,
                    layers_design_sequence_identical=False)
                scaffold = root / "scaffold.json"
                final = root / "final.json"
                generate_scaffold_review(str(scaffold), str(sst))
                finalize(str(scaffold), str(final))

                for path in (sst, scaffold, final):
                    payload = payload_to_internal_numbering(json.loads(
                        path.read_text(encoding="utf-8")))
                    layout = payload["moire_structure_metadata"][
                        "variable_length_layout"]
                    self.assertEqual(
                        layout["seed_sst_overlap_bp"], expected_overlap)
                    self.assertEqual(
                        layout["capture_columns_by_layer"], expected_columns)
                final_layout = payload["moire_structure_metadata"][
                    "variable_length_layout"]
                self.assertEqual(
                    [final_layout["seed_z1_actual_bp"],
                     final_layout["seed_z3_actual_bp"]], expected_overlap)

    def test_scaffold_recovers_nominal_z2_range_from_layer_ranges(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sst = write_shifted_sst(
                str(root / "sst.json"), 96, 40, 96,
                layers_design_sequence_identical=False,
                mean_indel_per_helix=6.0)
            payload = json.loads(sst.read_text(encoding="utf-8"))
            layout = payload["moire_structure_metadata"][
                "variable_length_layout"]
            layout.pop("spacing_range", None)
            layout.pop("seed_z2_range", None)
            sst.write_text(json.dumps(payload), encoding="utf-8")

            generate_scaffold_review(
                str(root / "scaffold.json"), str(sst))
            result = payload_to_internal_numbering(json.loads(
                (root / "scaffold.json").read_text(encoding="utf-8")))
            result_layout = result["moire_structure_metadata"][
                "variable_length_layout"]
            first, second = result_layout["layer_ranges"]
            expected = [first[1] + 1, second[0] - 1]
            self.assertEqual(result_layout["spacing_range"], expected)
            self.assertEqual(result_layout["seed_z2_range"], expected)
            self.assertEqual(expected[1] - expected[0] + 1, 40)

    def test_mixed_lattice_preserves_indel_and_seed_overlap_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sst = write_shifted_sst(
                str(root / "sst.json"), 96, 40, 96,
                lattice_type="square_kagome",
                layers_design_sequence_identical=False,
                mean_indel_per_helix=6.0)
            scaffold_path = root / "scaffold.json"
            final_path = root / "final.json"
            generate_scaffold_review(str(scaffold_path), str(sst))
            finalize(str(scaffold_path), str(final_path))

            for path in (scaffold_path, final_path):
                payload = payload_to_internal_numbering(json.loads(
                    path.read_text(encoding="utf-8")))
                layout = payload["moire_structure_metadata"][
                    "variable_length_layout"]
                self.assertEqual(layout["seed_z1_actual_bp"], 96)
                self.assertEqual(layout["seed_z3_actual_bp"], 96)
                self.assertEqual(layout["seed_sst_overlap_bp"], [96, 96])
                self.assertEqual(layout["seed_z2_range"], [168, 207])
                self.assertAlmostEqual(
                    layout["mean_indel_per_helix_requested"], 6.0)
                self.assertAlmostEqual(
                    layout["mean_indel_per_helix_actual"], 6.0)
                self.assertAlmostEqual(
                    layout["actual_z2_spacing_bp"], 46.0)

            report = validate_structure(
                str(final_path), require_staples=True)
            self.assertEqual(report["seed_z1_overlap_bp"], 96.0)
            self.assertEqual(report["seed_z3_overlap_bp"], 96.0)
            self.assertAlmostEqual(report["seed_z2_actual_bp"], 46.0)
            self.assertEqual(
                report["seed_scaffold_lengths"],
                sorted(payload["moire_structure_metadata"][
                    "variable_length_layout"][
                        "seed_scaffold_lengths_after_indel"]))

    def test_deletions_are_balanced_across_nominal_8bp_domains(self):
        for spacing, deletion in ((32, -12.0), (64, -24.0)):
            with self.subTest(spacing=spacing), \
                    tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                sst = write_shifted_sst(
                    str(root / "sst.json"), 128, spacing, 128,
                    layers_design_sequence_identical=False,
                    mean_indel_per_helix=deletion)
                generate_scaffold_review(
                    str(root / "scaffold.json"), str(sst))
                payload = payload_to_internal_numbering(json.loads(
                    (root / "scaffold.json").read_text(encoding="utf-8")))
                layout = payload["moire_structure_metadata"][
                    "variable_length_layout"]
                placements = layout["seed_z2_indel_placements"]
                per_helix = Counter(
                    item["helix"] for item in placements)
                per_domain = Counter(
                    (item["helix"], item["z2_domain"])
                    for item in placements)
                self.assertEqual(len(placements), int(-deletion * 48))
                self.assertEqual(set(per_helix.values()), {int(-deletion)})
                self.assertLessEqual(max(per_domain.values()), 3)
                self.assertEqual(layout["seed_z2_domain_bp"], 8)
                self.assertEqual(
                    layout["maximum_seed_deletions_per_domain"], 3)
                # Every deletion is assigned to a unique equal-width Z2 bin,
                # so a legal result cannot cluster at one side of a domain.
                for helix in range(48):
                    helix_placements = sorted(
                        (item for item in placements
                         if item["helix"] == helix),
                        key=lambda item: item["z2_bin"])
                    self.assertEqual(
                        [item["z2_bin"] for item in helix_placements],
                        list(range(1, int(-deletion) + 1)))
                    self.assertTrue(all(
                        item["z2_bin_count"] == int(-deletion)
                        for item in helix_placements))
                self.assertLessEqual(
                    layout["seed_z2_indel_distribution"][
                        "maximum_distance_from_bin_center"], 1.0)

    def test_insertions_are_centered_in_equal_z2_partitions(self):
        for spacing, insertion in ((80, 2.0), (80, 5.0), (32, 9.0)):
            with self.subTest(spacing=spacing, insertion=insertion), \
                    tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                sst = write_shifted_sst(
                    str(root / "sst.json"), 128, spacing, 128,
                    layers_design_sequence_identical=False,
                    mean_indel_per_helix=insertion)
                generate_scaffold_review(
                    str(root / "scaffold.json"), str(sst))
                payload = payload_to_internal_numbering(json.loads(
                    (root / "scaffold.json").read_text(encoding="utf-8")))
                layout = payload["moire_structure_metadata"][
                    "variable_length_layout"]
                placements = layout["seed_z2_indel_placements"]
                z2_low, z2_high = layout["seed_z2_indel_range"]
                width = (z2_high-z2_low+1) / insertion
                per_domain = Counter(
                    (item["helix"], item["z2_domain"])
                    for item in placements)
                self.assertLessEqual(max(per_domain.values()), 3)
                self.assertEqual(
                    layout["maximum_seed_indels_per_domain"], 3)
                for helix in range(48):
                    helix_placements = sorted(
                        (item for item in placements
                         if item["helix"] == helix),
                        key=lambda item: item["z2_bin"])
                    self.assertEqual(len(helix_placements), int(insertion))
                    for rank, item in enumerate(helix_placements):
                        self.assertEqual(item["z2_bin"], rank+1)
                        self.assertEqual(
                            item["z2_bin_count"], int(insertion))
                        # In unconstrained cases the legal base remains
                        # inside its intended equal partition.  The 32-bp,
                        # +9 stress case may use the documented nearest-safe
                        # fallback to retain the hard three-per-domain cap.
                        if spacing == 80:
                            self.assertGreaterEqual(
                                item["base"], z2_low + rank*width)
                            self.assertLess(
                                item["base"], z2_low + (rank+1)*width)
                distribution = layout["seed_z2_indel_distribution"]
                self.assertEqual(
                    set(distribution["per_helix_counts"].values()),
                    {int(insertion)})
                self.assertLessEqual(
                    distribution["maximum_distance_from_bin_center"],
                    2.0 if spacing == 80 else 8.0)

    def test_fractional_mean_insertions_balance_zero_and_nonzero_helices(self):
        """Totals below 48 must keep zero-selection DP costs well typed."""
        for total in (1, 3, 47, 48, 49):
            with self.subTest(total=total), \
                    tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                sst = write_shifted_sst(
                    str(root / "sst.json"), 128, 120, 128,
                    layers_design_sequence_identical=False,
                    mean_indel_per_helix=total / 48.0)
                generate_scaffold_review(
                    str(root / "scaffold.json"), str(sst))
                payload = payload_to_internal_numbering(json.loads(
                    (root / "scaffold.json").read_text(encoding="utf-8")))
                layout = payload["moire_structure_metadata"][
                    "variable_length_layout"]
                placements = layout["seed_z2_indel_placements"]
                distribution = layout["seed_z2_indel_distribution"]
                counts = [
                    int(distribution["per_helix_counts"][str(helix)])
                    for helix in range(48)]

                self.assertEqual(len(placements), total)
                self.assertEqual(sum(counts), total)
                self.assertLessEqual(max(counts)-min(counts), 1)
                self.assertEqual(
                    sorted(Counter(counts).items()),
                    sorted(Counter(
                        [total // 48 + int(index < total % 48)
                         for index in range(48)]).items()))
                self.assertTrue(all(item["value"] == 1
                                    for item in placements))
                self.assertTrue(all(
                    layout["seed_z2_indel_range"][0] <= item["base"] <=
                    layout["seed_z2_indel_range"][1]
                    for item in placements))
                self.assertIn(
                    "no forced stagger", distribution["method"])
                self.assertTrue(all(
                    length <= 7557 for length in
                    layout["seed_scaffold_lengths_after_indel"]))

    def test_scaffold_capacity_precedes_insertion_uniformity(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sst = write_shifted_sst(
                str(root / "sst.json"), 128, 80, 128,
                layers_design_sequence_identical=False,
                mean_indel_per_helix=9.0)
            generate_scaffold_review(
                str(root / "scaffold.json"), str(sst))
            payload = payload_to_internal_numbering(json.loads(
                (root / "scaffold.json").read_text(encoding="utf-8")))
            layout = payload["moire_structure_metadata"][
                "variable_length_layout"]
            distribution = layout["seed_z2_indel_distribution"]
            self.assertEqual(len(layout["seed_z2_indel_placements"]), 9*48)
            self.assertTrue(all(
                length <= 7557 for length in
                layout["seed_scaffold_lengths_after_indel"]))
            self.assertEqual(
                distribution["scaffold_lengths_nt"],
                layout["seed_scaffold_lengths_after_indel"])
            # This fixture reaches one scaffold's hard capacity and therefore
            # deliberately accepts a larger spatial displacement.
            self.assertIn(7557, distribution["scaffold_lengths_nt"])
            self.assertGreater(
                distribution["maximum_distance_from_bin_center"], 2.0)

    def test_structure_worker_rejects_deletion_below_dynamic_limit(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sst = write_shifted_sst(
                str(root / "sst.json"), 128, 32, 128,
                layers_design_sequence_identical=False,
                mean_indel_per_helix=-12.1)
            with self.assertRaisesRegex(
                    RuntimeError, "Each 8-bp domain permits at most 3"):
                generate_scaffold_review(
                    str(root / "scaffold.json"), str(sst))

    def test_identical_layers_require_linked_phase(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            good = write_shifted_sst(
                str(root / "good.json"), 112, 48, 112,
                layers_design_sequence_identical=True)
            self.assertTrue(validate_sst(str(good))["valid"])
            self.assertTrue(phase_is_compatible(112, 48, 112))
            with self.assertRaisesRegex(ValueError, "32-bp"):
                write_shifted_sst(
                    str(root / "bad.json"), 112, 40, 112,
                    layers_design_sequence_identical=True)

    def test_independent_layers_keep_separate_reviewed_fixture_phases(self):
        with tempfile.TemporaryDirectory() as directory:
            target = write_shifted_sst(
                str(Path(directory) / "independent.json"),
                96, 48, 128, layers_design_sequence_identical=False)
            validation = validate_sst(str(target))
            self.assertTrue(validation["valid"], validation["errors"])
            payload = payload_to_internal_numbering(json.loads(
                target.read_text(encoding="utf-8")))
            layout = payload["moire_structure_metadata"][
                "variable_length_layout"]
            self.assertFalse(layout["layers_design_sequence_identical"])
            self.assertEqual(layout["sst_scaffold_ranges"],
                             [[64, 159], [208, 335]])
            self.assertEqual(layout[
                "sst_complementary_chain_ranges"],
                [[56, 167], [200, 343]])
            self.assertEqual(layout["seed_partition_lengths_bp"],
                             [112, 48, 128])
            self.assertEqual(layout["seed_sst_overlap_bp"], [96, 128])
            self.assertEqual(len(layout["capture_positions_by_layer"]), 2)
            self.assertFalse(phase_is_compatible(96, 48, 128))
            snapshot = build_output_sst_snapshot_payload(
                payload, "output.json")
            self.assertEqual(
                snapshot["moire_structure_metadata"]["export_role"],
                "output")


if __name__ == "__main__":
    unittest.main()
