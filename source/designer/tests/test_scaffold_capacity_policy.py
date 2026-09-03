"""Hard threshold tests for the Moiré Seed scaffold allocation policy."""

import unittest

from moire_design_core.structure import scaffold_capacity_plan


class ScaffoldCapacityPolicyTests(unittest.TestCase):

    def test_exact_capacity_boundaries(self):
        expected = {
            1: (1, 8064),
            8064: (1, 8064),
            8065: (2, 7557),
            15114: (2, 7557),
            15115: (3, 7557),
            22671: (3, 7557),
        }
        for total, (count, per_scaffold) in expected.items():
            with self.subTest(total=total):
                plan = scaffold_capacity_plan(total)
                self.assertEqual(plan["count"], count)
                self.assertEqual(
                    plan["per_scaffold_capacity_nt"], per_scaffold)

    def test_over_capacity_is_rejected(self):
        with self.assertRaises(ValueError):
            scaffold_capacity_plan(22672)


if __name__ == "__main__":
    unittest.main()
