import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools" / "generate_kagome_centered_sst_review.py"
SPEC = importlib.util.spec_from_file_location("kagome_centered_review", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class KagomeCenteredSSTReviewTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.seed = MODULE.load(MODULE.SEED_REFERENCE)
        unused, cls.source_rows = MODULE.kagome._source_rows()

    def test_all_review_cases_are_complete_and_valid(self):
        for case in MODULE.CASES:
            filename, payload, metadata = MODULE.build_case(
                *case, self.seed, self.source_rows)
            self.assertEqual(payload["name"], filename)
            self.assertTrue(metadata["audit"]["passed"])
            self.assertFalse(metadata["capture_gaps_reserved"])
            self.assertGreaterEqual(metadata["sst_duplex_ranges"][0][0], 32)
            self.assertEqual(metadata["audit"]["capture_nick_count"], 0)
            self.assertEqual(
                metadata["hole_kagome_sst_helices_public"], [5, 7, 13, 15])
            self.assertTrue(set(map(int, metadata["audit"]
                                    ["scaffold_component_lengths"])) <= {32, 48})
            self.assertTrue(set(map(int, metadata["audit"]
                                    ["complement_component_lengths"])) <= {32, 48})

    def test_small_spacing_review_stops_before_known_overlap_cases(self):
        spacings = {
            int(case[2]) for case in MODULE.CASES
            if int(case[1]) == int(case[3]) == 128}
        self.assertIn(16, spacings)
        self.assertIn(24, spacings)
        self.assertNotIn(8, spacings)
        self.assertNotIn(0, spacings)


if __name__ == "__main__":
    unittest.main()
