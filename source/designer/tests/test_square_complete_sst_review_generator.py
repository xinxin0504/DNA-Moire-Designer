"""Regression tests for the complete-SST Square review generator."""

import importlib.util
from pathlib import Path
import unittest


SCRIPT = (Path(__file__).resolve().parents[1] / "tools" /
          "generate_sst_absolute_phase_review.py")
SPEC = importlib.util.spec_from_file_location("square_review_generator", SCRIPT)
GENERATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GENERATOR)


class SquareCompleteSSTReviewGeneratorTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.sst = GENERATOR.load(GENERATOR.SST_REFERENCE)
        cls.seed = GENERATOR.load(GENERATOR.SEED_REFERENCE)
        cls.source_rows = {
            int(row["num"]): row for row in cls.sst["vstrands"]
            if int(row["num"]) in GENERATOR.SST_SOURCE_HELICES}

    def test_review_cases_are_complete_phase_correct_square_sst(self):
        for case in GENERATOR.CASES:
            with self.subTest(case=case):
                unused_filename, unused_payload, metadata = (
                    GENERATOR.build_case(*case, self.sst, self.seed))
                self.assertTrue(metadata["audit"]["passed"])
                self.assertFalse(metadata["capture_gaps_reserved"])
                self.assertEqual(
                    metadata["audit"]["non_32nt_component_count"], 0)
                self.assertEqual(
                    metadata["audit"]["unintended_duplex_inside_Z2_count"],
                    0)
                self.assertIn(
                    metadata["placement"]["global_32bp_canvas_shift"],
                    (0, 32))

    def test_parameter_grid_has_unique_complete_u_boundaries(self):
        lengths = (64, 80, 96, 104, 112, 120, 128, 136, 144, 152, 160)
        spacings = (0, 32, 40, 48, 56, 64, 72, 80, 96, 128, 160)
        for first in lengths:
            for spacing in spacings:
                for second in (first, 128):
                    with self.subTest(values=(first, spacing, second)):
                        duplex, placement = GENERATOR.desired_duplex_ranges(
                            first, spacing, second)
                        resolved = GENERATOR.resolve_complete_sst_ranges(
                            self.source_rows, duplex, placement)
                        self.assertEqual(len(resolved[0]), 2)

    def test_centre_changes_in_discrete_eight_bp_steps(self):
        duplex_32, unused = GENERATOR.desired_duplex_ranges(128, 32, 128)
        duplex_48, unused = GENERATOR.desired_duplex_ranges(112, 48, 112)
        self.assertEqual(duplex_32, ((48, 175), (208, 335)))
        self.assertEqual(duplex_48, ((56, 167), (216, 327)))
        self.assertEqual(
            (duplex_48[0][0]-duplex_32[0][0]) % 8, 0)

    def test_invalid_linked_phase_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "32-bp"):
            GENERATOR.build_case(
                "linked", 112, 40, 112, True, self.sst, self.seed)

    def test_canvas_shift_translates_seed_without_changing_topology(self):
        unused_name, payload, metadata = GENERATOR.build_case(
            "linked", 160, 96, 160, True, self.sst, self.seed)
        shift = metadata["placement"]["global_32bp_canvas_shift"]
        self.assertEqual(shift, 32)
        output = {int(row["num"]): row for row in payload["vstrands"]}
        for source in self.seed["vstrands"]:
            source_helix = int(source["num"])
            if source_helix not in GENERATOR.SEED_SOURCE_HELICES:
                continue
            target_helix = source_helix+GENERATOR.SEED_TARGET_SHIFT
            for source_base, record in enumerate(source.get("scaf", [])):
                expected = []
                for side in (0, 2):
                    partner, partner_base = record[side:side+2]
                    if partner < 0:
                        expected.extend((-1, -1))
                    else:
                        expected.extend((
                            partner+GENERATOR.SEED_TARGET_SHIFT,
                            partner_base+shift))
                self.assertEqual(
                    output[target_helix]["scaf"][source_base+shift],
                    expected)


if __name__ == "__main__":
    unittest.main()
