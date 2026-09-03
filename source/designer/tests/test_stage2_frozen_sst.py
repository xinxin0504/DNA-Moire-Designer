"""Square SST lengths must preserve the reviewed complete-U topology."""

import unittest

from moire_design_core.structure import (
    build_shifted_sst_payload,
)


class FrozenSSTTests(unittest.TestCase):

    def test_square_lengths_are_available_through_shared_complete_u_router(self):
        for length in (64, 80, 96, 104, 112, 120, 128, 136, 152, 160):
            with self.subTest(length=length):
                payload = build_shifted_sst_payload(
                    "square_%d.json" % length,
                    z1_bp=length, z2_bp=32, z3_bp=length,
                    lattice_type="square")
                layout = payload["moire_structure_metadata"][
                    "variable_length_layout"]
                self.assertEqual(layout["z1_bp"], length)
                self.assertEqual(layout["z3_bp"], length)

    def test_non_domain_square_length_is_rejected(self):
        with self.assertRaises(ValueError):
            build_shifted_sst_payload(
                "unsupported.json", z1_bp=100, z2_bp=32, z3_bp=100,
                lattice_type="square")


if __name__ == "__main__":
    unittest.main()
