"""Regression tests for the Z2=0/8 auxiliary review convention."""

import importlib.util
from pathlib import Path
import unittest


SCRIPT = (Path(__file__).resolve().parents[1] / "tools" /
          "generate_zero_z2_auxiliary_review.py")
SPEC = importlib.util.spec_from_file_location("aux_review", SCRIPT)
GENERATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GENERATOR)


class ZeroZ2AuxiliaryReviewTests(unittest.TestCase):

    def test_all_eight_files_keep_connected_partial_detours(self):
        for lattice in ("square", "kagome"):
            for case in GENERATOR.CASES:
                for capture_ready in (False, True):
                    with self.subTest(
                            lattice=lattice, case=case,
                            capture_ready=capture_ready):
                        unused_name, payload = GENERATOR.build_case(
                            lattice, *case, capture_ready=capture_ready)
                        metadata = payload["moire_structure_metadata"]
                        audit = metadata["audit"]
                        self.assertTrue(audit["passed"])
                        self.assertEqual(audit["reciprocity_error_count"], 0)
                        self.assertGreater(
                            audit["detoured_component_count"], 0)
                        self.assertGreater(
                            audit["auxiliary_nonempty_helix_count"], 0)
                        self.assertEqual(audit["seed_staple_record_count"], 0)
                        boundary_link_count = 0
                        for item in metadata["auxiliary_review_policy"][
                                "detoured_components"]:
                            self.assertGreater(item["detoured_node_count"], 0)
                            boundary_link_count += len(
                                item["primary_auxiliary_boundary_links"])
                            # A component can be wholly inside the occupied
                            # interval at Z2=0.  Otherwise it must expose the
                            # direct primary/auxiliary return link.
                            if not item["primary_auxiliary_boundary_links"]:
                                self.assertEqual(
                                    item["detoured_node_count"],
                                    item["component_length_nt"])
                        self.assertGreater(boundary_link_count, 0)

    def test_manual_first_pair_examples_are_reproduced(self):
        examples = (
            ("kagome", (64, 0, 64), False,
             "kagome_SST64_Z2_0_SST64_seed_scaffold_aux16_complete_sst_review.json"),
            ("square", (64, 0, 64), True,
             "square_SST64_Z2_0_SST64_seed_scaffold_aux16_capture_display_review.json"),
            ("square", (88, 8, 88), False,
             "square_SST88_Z2_8_SST88_seed_scaffold_aux16_complete_sst_review.json"),
        )
        desktop = Path.home() / "Desktop"
        for lattice, case, capture_ready, filename in examples:
            reference = desktop / filename
            if not reference.is_file():
                continue
            unused_name, payload = GENERATOR.build_case(
                lattice, *case, capture_ready=capture_ready)
            generated = {int(row["num"]): row
                         for row in payload["vstrands"]}
            import json
            manual = {int(row["num"]): row for row in
                      json.loads(reference.read_text())["vstrands"]}
            for helix in (0, 1, 64, 65):
                self.assertEqual(generated[helix]["scaf"],
                                 manual[helix]["scaf"])
                self.assertEqual(generated[helix]["stap"],
                                 manual[helix]["stap"])


if __name__ == "__main__":
    unittest.main()
