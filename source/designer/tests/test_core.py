import math
import tempfile
import unittest
from pathlib import Path

from moire_design_core import (
    MAX_SEED_DELETIONS_PER_DOMAIN,
    SquareBilayerSettings,
    angle_from_moire_period,
    calibrated_indel_for_twist_per_base,
    calibrated_twist_per_base,
    elastic_calibrated_twist_for_cross_section,
    compatible_growth_values,
    compatible_z2_values,
    load_project,
    minimum_seed_deletion_per_helix,
    moire_period_from_angle,
    preview_seed_partition,
    save_project,
    solve_square_bilayer,
    validate_structure,
    validate_sst,
    write_shifted_sst,
    finalize_structure,
)
from moire_design_core.structure import (
    CAPTURE_PAIR_COLORS,
    build_output_sst_snapshot_payload,
    generate_scaffold_review,
    payload_to_internal_numbering,
)
from moire_design_core.project import export_capture_map
from moire_design_core.template import export_reference_seed, reference_seed_path


class CoreTests(unittest.TestCase):
    def test_angle_period_round_trip(self):
        period = moire_period_from_angle(3.8, 2.8)
        self.assertAlmostEqual(angle_from_moire_period(period, 2.8), 3.8)

    def test_symmetry_uses_its_lattice_constant_and_period_availability(self):
        cases = (
            ("square_square_c4", 2.8, 2.8, True),
            ("kagome_kagome", 5.4, 5.4, True),
            ("square_kagome", 2.8, 5.4, False),
        )
        for symmetry, first_a, second_a, period_available in cases:
            project = solve_square_bilayer(SquareBilayerSettings(
                target_mode="angle",
                target_angle_deg=3.3,
                lattice_symmetry=symmetry,
                lattice_constant_nm=first_a,
                layer1_lattice_constant_nm=first_a,
                layer2_lattice_constant_nm=second_a,
            ))
            self.assertEqual(
                project.prediction["layer_lattice_constants_nm"],
                [first_a, second_a])
            self.assertEqual(
                project.prediction["period_available"], period_available)
            if period_available:
                self.assertAlmostEqual(
                    project.prediction["predicted_moire_period_nm"],
                    moire_period_from_angle(3.3, first_a))
            else:
                self.assertIsNone(
                    project.prediction["predicted_moire_period_nm"])

    def test_calibrated_preset(self):
        project = solve_square_bilayer(SquareBilayerSettings())
        self.assertEqual(project.settings.spacer_bp_z2, 32)
        self.assertAlmostEqual(
            project.prediction["predicted_local_surface_angle_deg"],
            3.2967555036483183)
        self.assertEqual(project.seed_plan["occupied_helices"], 48)

    def test_phase_compatible_series(self):
        self.assertEqual(compatible_z2_values(128)[:3], [0, 32, 64])
        self.assertEqual(compatible_z2_values(120)[:3], [8, 40, 72])
        self.assertEqual(compatible_z2_values(136)[:3], [24, 56, 88])
        self.assertEqual(compatible_z2_values(112)[:3], [16, 48, 80])
        self.assertIn(120, compatible_growth_values(40))

    def test_z2_does_not_change_target_angle(self):
        settings = SquareBilayerSettings(target_mode="angle")
        first = solve_square_bilayer(settings)
        settings.spacer_bp_z2 = 64
        second = solve_square_bilayer(settings)
        self.assertAlmostEqual(first.prediction["reported_angle_deg"],
                               second.prediction["reported_angle_deg"])
        self.assertNotAlmostEqual(first.settings.mean_indel_per_helix,
                                  second.settings.mean_indel_per_helix)

    def test_seed_indel_limits_are_spacing_dependent(self):
        too_much_insertion = solve_square_bilayer(
            SquareBilayerSettings(
                target_mode="angle", target_angle_deg=45.0,
                spacer_bp_z2=32))
        self.assertGreater(
            too_much_insertion.settings.mean_indel_per_helix, 10.0)
        self.assertTrue(too_much_insertion.prediction[
            "seed_insertion_limit_exceeded"])
        self.assertFalse(any(
            item["title"] == "Seed insertion超过上限"
            for item in too_much_insertion.validation))

        excessive_deletion = solve_square_bilayer(
            SquareBilayerSettings(
                target_mode="indel", mean_indel_per_helix=-12.1,
                spacer_bp_z2=32))
        self.assertTrue(excessive_deletion.prediction[
            "seed_deletion_limit_exceeded"])
        self.assertTrue(excessive_deletion.prediction[
            "seed_indel_limit_exceeded"])
        self.assertFalse(excessive_deletion.prediction[
            "seed_insertion_limit_exceeded"])
        self.assertEqual(minimum_seed_deletion_per_helix(32), -12.0)
        self.assertEqual(minimum_seed_deletion_per_helix(64), -24.0)
        self.assertEqual(MAX_SEED_DELETIONS_PER_DOMAIN, 3)

        boundary = solve_square_bilayer(
            SquareBilayerSettings(
                target_mode="indel", mean_indel_per_helix=-24.0,
                spacer_bp_z2=64))
        self.assertFalse(boundary.prediction[
            "seed_deletion_limit_exceeded"])
        self.assertFalse(boundary.prediction[
            "seed_indel_limit_exceeded"])

    def test_twist_prediction_normalizes_indel_by_selected_length(self):
        cells = tuple((row, column) for row in range(8)
                      for column in range(8)
                      if not (2 <= row <= 5 and 2 <= column <= 5))
        per_base_96 = elastic_calibrated_twist_for_cross_section(
            9.0, 96.0, cells)
        per_base_32_same_density = (
            elastic_calibrated_twist_for_cross_section(
                3.0, 32.0, cells))
        self.assertAlmostEqual(per_base_96, per_base_32_same_density,
                               places=12)
        self.assertAlmostEqual(per_base_96*96.0, 32.3931, places=3)

        # Reusing the total +9 edits inside only 32 bp is a threefold larger
        # edit density and must not reproduce the old ~9.4-degree lookup.
        compact = solve_square_bilayer(SquareBilayerSettings(
            target_mode="indel", mean_indel_per_helix=9.0,
            spacer_bp_z2=32))
        self.assertGreater(
            compact.prediction["predicted_local_surface_angle_deg"], 25.0)
        self.assertEqual(
            compact.prediction["twist_prediction_model"],
            "effective-pitch-elastic-then-SNUPI")
        self.assertAlmostEqual(
            compact.prediction["calibration_equivalent_indel_per_96bp"],
            27.0)
        self.assertTrue(compact.prediction[
            "seed_twist_calibration_domain_exceeded"])

    def test_large_twist_is_preserved_until_feasibility_validation(self):
        project = solve_square_bilayer(
            SquareBilayerSettings(
                target_mode="angle", target_angle_deg=90.0,
                spacer_bp_z2=32))
        self.assertAlmostEqual(
            project.prediction["reported_angle_deg"], 90.0, places=6)
        self.assertGreater(project.settings.mean_indel_per_helix, 10.0)
        self.assertTrue(project.prediction[
            "seed_insertion_limit_exceeded"])

    def test_fixed_seed_preview_partition_linked_and_independent(self):
        linked = preview_seed_partition(40, 120, 120, True)
        self.assertEqual(
            (linked["z1_bp"], linked["z2_bp"], linked["z3_bp"]),
            (120, 40, 128))
        self.assertTrue(linked["phase_compatible"])
        self.assertEqual(
            linked["z1_bp"]+linked["z2_bp"]+linked["z3_bp"], 288)

        independent = preview_seed_partition(40, 128, 96, False)
        self.assertEqual(
            (independent["z1_bp"], independent["z2_bp"],
             independent["z3_bp"]), (128, 40, 120))
        self.assertEqual(
            (independent["sst_overlap_z1_bp"],
             independent["sst_overlap_z3_bp"]), (128, 96))

    def test_project_exposes_fixed_seed_preview_partition(self):
        project = solve_square_bilayer(SquareBilayerSettings(
            sst_growth_bp_z1=120, spacer_bp_z2=40,
            sst_growth_bp_z3=120,
            layers_design_sequence_identical=True))
        partition = project.prediction["preview_seed_partition"]
        self.assertEqual(partition["total_bp"], 288)
        self.assertEqual(
            partition["z1_bp"]+partition["z2_bp"]+partition["z3_bp"],
            288)

    def test_legacy_seed_lengths_are_ignored(self):
        settings = SquareBilayerSettings(
            growth_bp_z1=120,
            spacer_bp_z2=32,
            growth_bp_z3=136,
            layers_design_sequence_identical=False,
        )
        project = solve_square_bilayer(settings)
        self.assertEqual(project.settings.growth_bp_z1, 128)
        self.assertEqual(project.settings.spacer_bp_z2, 32)
        self.assertEqual(project.settings.growth_bp_z3, 128)
        self.assertIn(project.validation[1]["level"], ("info", "pass"))

    def test_seed_support_is_fixed_reference(self):
        settings = SquareBilayerSettings(
            growth_bp_z1=121,
            layers_design_sequence_identical=False,
        )
        project = solve_square_bilayer(settings)
        self.assertEqual(project.settings.growth_bp_z1, 128)

    def test_calibration_inverse(self):
        for indel in (-10, -8, -6, -4, -2, 0, 2, 4, 6):
            twist = calibrated_twist_per_base(indel)
            self.assertAlmostEqual(
                calibrated_indel_for_twist_per_base(twist), indel)

    def test_extended_s8_r4x4c_calibration_metadata(self):
        project = solve_square_bilayer(SquareBilayerSettings())
        calibration = project.prediction["calibration"]
        self.assertEqual(len(calibration["points"]), 9)
        self.assertEqual(
            calibration["validated_indel_per_helix_range"], [-10.0, 6.0])
        self.assertEqual(
            calibration["failed_indel_per_helix_points"], [8.0, 10.0])
        self.assertTrue(calibration["failed_points_excluded_from_fit"])
        self.assertAlmostEqual(
            calibrated_twist_per_base(-10), -0.16367829784857435)
        self.assertAlmostEqual(
            calibrated_twist_per_base(6), 0.24969345505563809)

    def test_extended_calibration_marks_only_out_of_domain_values(self):
        boundary = solve_square_bilayer(SquareBilayerSettings(
            target_mode="indel", mean_indel_per_helix=2.0,
            spacer_bp_z2=32))
        outside = solve_square_bilayer(SquareBilayerSettings(
            target_mode="indel", mean_indel_per_helix=8.0/3.0,
            spacer_bp_z2=32))
        self.assertFalse(boundary.prediction[
            "seed_twist_calibration_domain_exceeded"])
        self.assertTrue(outside.prediction[
            "seed_twist_calibration_domain_exceeded"])
        self.assertTrue(any(
            item["title"] == "校准外推" for item in outside.validation))

    def test_save_round_trip(self):
        project = solve_square_bilayer(SquareBilayerSettings())
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"test.moire.json"
            save_project(project, str(path))
            restored = load_project(str(path))
        self.assertEqual(restored.settings.project_name,
                         project.settings.project_name)

    def test_reference_exports(self):
        if reference_seed_path() is None:
            self.skipTest("reference Seed-S JSON is not installed")
        project = solve_square_bilayer(SquareBilayerSettings())
        with tempfile.TemporaryDirectory() as folder:
            seed = export_reference_seed(project, str(Path(folder)/"seed.json"))
            capture = export_capture_map(project, str(Path(folder)/"capture.csv"))
            self.assertGreater(seed.stat().st_size, 100000)
            self.assertIn("capture-0 + capture-1",
                          capture.read_text(encoding="utf-8-sig"))

    def test_sst_template_is_shifted_to_base_48(self):
        with tempfile.TemporaryDirectory() as folder:
            target = write_shifted_sst(str(Path(folder)/"sst.json"))
            payload = __import__("json").loads(target.read_text())
        rows = payload["vstrands"]
        for row in rows:
            occupied = [index for index, record in enumerate(row["scaf"])
                        if record != [-1, -1, -1, -1]]
            self.assertEqual((min(occupied), max(occupied)), (48, 335))
            self.assertTrue(all(
                48 <= index <= 175 or 208 <= index <= 335
                for index in occupied))
        self.assertEqual(
            payload["moire_structure_metadata"]["base_shift_bp"], 32)

    def test_scaffold_stage_uses_capacity_safe_balanced_bands(self):
        with tempfile.TemporaryDirectory() as folder:
            sst = write_shifted_sst(str(Path(folder)/"sst.json"))
            self.assertTrue(validate_sst(str(sst))["valid"])
            target = Path(folder)/"scaffold.json"
            report = generate_scaffold_review(str(target), str(sst))
            validation = validate_structure(str(target))
            scaffold_payload = payload_to_internal_numbering(
                __import__("json").loads(target.read_text()))
            final = Path(folder)/"structure.json"
            finalize_structure(str(target), str(final))
            final_validation = validate_structure(
                str(final), require_staples=True)
        self.assertEqual(report["sst_first_base"], 48)
        self.assertEqual(len(report["seed_scaffold_lengths"]), 2)
        self.assertLessEqual(max(report["seed_scaffold_lengths"]), 7557)
        self.assertGreaterEqual(
            min(report["seed_scaffold_lengths"]) * 2,
            max(report["seed_scaffold_lengths"]))
        seed_rows = {int(row["num"]): row
                     for row in scaffold_payload["vstrands"]
                     if int(row["num"]) < 48}
        actual_ranges = {}
        for number, row in seed_rows.items():
            occupied = [index for index, record in enumerate(row["scaf"])
                        if record != [-1, -1, -1, -1]]
            actual_ranges[number] = (min(occupied), max(occupied))
        self.assertLessEqual(
            max(high for unused_low, high in actual_ranges.values()) -
            min(high for unused_low, high in actual_ranges.values()), 21)
        self.assertLessEqual(
            max(low for low, unused_high in actual_ranges.values()) -
            min(low for low, unused_high in actual_ranges.values()), 21)
        self.assertTrue(validation["valid"])
        self.assertEqual(validation["seed_scaffold_single_nick_count"], 2)
        self.assertTrue(all(
            component["helices"] == list(range(48))
            for component in validation["seed_scaffold_components"]))
        routing = scaffold_payload["moire_structure_metadata"][
            "seed_edge_routing"]
        self.assertIn("capacity-safe Path-view bands", routing["partition"])
        self.assertEqual(
            routing["scope"],
            "Moiré Designer only; cadnano AutoCS unchanged")
        self.assertEqual(len(routing["band_ranges"]), 2)
        for number in range(48):
            bands = [item[str(number)] for item in routing["band_ranges"]]
            self.assertTrue(all(
                left[1] + 1 == right[0]
                for left, right in zip(bands, bands[1:])))
        self.assertTrue(final_validation["valid"])
        self.assertEqual(
            final_validation["capture_bridge_component_count"], 64)
        # Capture colors are copied from the accepted 2L template rather
        # than recolored by a post-finalization optimizer.
        self.assertGreaterEqual(final_validation["capture_color_count"], 4)

    def test_output_snapshot_closes_gaps_without_mutating_capture_design(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            sst = write_shifted_sst(str(root/"sst.json"))
            scaffold = root/"scaffold.json"
            final = root/"structure.json"
            generate_scaffold_review(str(scaffold), str(sst))
            finalize_structure(str(scaffold), str(final))
            source = __import__("json").loads(final.read_text())
            original = __import__("json").dumps(source, sort_keys=True)
            output = payload_to_internal_numbering(
                build_output_sst_snapshot_payload(source, "output.json"))
        self.assertEqual(original,
                         __import__("json").dumps(source, sort_keys=True))
        crossing = []
        for row in output["vstrands"]:
            if int(row["num"]) >= 48:
                continue
            for index, record in enumerate(row["stap"]):
                for offset in (0, 2):
                    if int(record[offset]) >= 48:
                        crossing.append((row["num"], index))
        self.assertEqual(crossing, [])
        occupied = [
            index for row in output["vstrands"] if int(row["num"]) >= 48
            for index, record in enumerate(row["scaf"])
            if record != [-1, -1, -1, -1]]
        self.assertEqual((min(occupied), max(occupied)), (48, 335))
        self.assertEqual(
            output["moire_structure_metadata"]["export_role"], "output")


if __name__ == "__main__":
    unittest.main()
