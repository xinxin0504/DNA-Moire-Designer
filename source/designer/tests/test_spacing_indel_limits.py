import math
import unittest

from moire_design_core import (
    SquareBilayerSettings,
    maximum_seed_insertion_per_helix,
    minimum_seed_deletion_per_helix,
    seed_indel_limits_for_spacing,
    solve_square_bilayer,
)


class SpacingIndelLimitTests(unittest.TestCase):
    def test_exact_spacing_limit_table(self):
        expected = {
            0: (0.0, 0.0),
            8: (-3.0, 3.0),
            16: (-6.0, 6.0),
            24: (-9.0, 9.0),
            32: (-12.0, 10.0),
            40: (-15.0, 10.0),
            64: (-24.0, 10.0),
            160: (-60.0, 10.0),
        }
        for spacing, limits in expected.items():
            with self.subTest(spacing=spacing):
                actual = seed_indel_limits_for_spacing(spacing)
                self.assertEqual(tuple(float(value) for value in actual),
                                 limits)
                self.assertEqual(
                    minimum_seed_deletion_per_helix(spacing), limits[0])
                self.assertEqual(
                    maximum_seed_insertion_per_helix(spacing), limits[1])

    def test_invalid_spacing_is_rejected_by_canonical_limit_function(self):
        for spacing in (-8, 1, 7, 9, 31):
            with self.subTest(spacing=spacing):
                with self.assertRaises(ValueError):
                    seed_indel_limits_for_spacing(spacing)

    def test_solver_flags_both_sides_of_every_compact_boundary(self):
        for spacing, lower, upper in (
                (8, -3.0, 3.0),
                (16, -6.0, 6.0),
                (24, -9.0, 9.0),
                (32, -12.0, 10.0)):
            base = dict(
                target_mode="indel",
                spacer_bp_z2=spacing,
                layers_design_sequence_identical=False,
            )
            for value in (lower, upper):
                project = solve_square_bilayer(SquareBilayerSettings(
                    mean_indel_per_helix=value, **base))
                self.assertFalse(
                    project.prediction["seed_indel_limit_exceeded"],
                    (spacing, value))
            below = solve_square_bilayer(SquareBilayerSettings(
                mean_indel_per_helix=lower-0.1, **base))
            above = solve_square_bilayer(SquareBilayerSettings(
                mean_indel_per_helix=upper+0.1, **base))
            self.assertTrue(
                below.prediction["seed_deletion_limit_exceeded"])
            self.assertTrue(
                above.prediction["seed_insertion_limit_exceeded"])

    def test_zero_spacing_forces_zero_twist_period_and_indel(self):
        for target_mode in ("angle", "period", "indel"):
            project = solve_square_bilayer(SquareBilayerSettings(
                target_mode=target_mode,
                target_angle_deg=37.0,
                target_period_nm=12.0,
                mean_indel_per_helix=9.0,
                spacer_bp_z2=0,
                layers_design_sequence_identical=False,
            ))
            with self.subTest(target_mode=target_mode):
                self.assertEqual(project.settings.target_mode, "angle")
                self.assertEqual(project.settings.target_angle_deg, 0.0)
                self.assertEqual(project.settings.mean_indel_per_helix, 0.0)
                self.assertEqual(
                    project.prediction["reported_angle_deg"], 0.0)
                self.assertTrue(math.isinf(
                    project.prediction["predicted_moire_period_nm"]))
                self.assertEqual(
                    project.prediction[
                        "maximum_seed_insertion_per_helix"], 0.0)
                self.assertFalse(
                    project.prediction["seed_indel_limit_exceeded"])


if __name__ == "__main__":
    unittest.main()
