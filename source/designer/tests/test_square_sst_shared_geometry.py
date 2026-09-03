"""Square JSON generation and 1B preview share one centred geometry."""

import json
from pathlib import Path
import tempfile
import unittest

from moire_design_core.calculations import preview_seed_partition
from moire_design_core.square_sst_geometry import (
    centered_square_sst_geometry,
    complete_square_polymer_ranges,
)
from moire_design_core.structure import validate_sst, write_shifted_sst


class SquareSSTSharedGeometryTests(unittest.TestCase):

    CASES = (
        (128, 32, 128, True),
        (120, 40, 120, True),
        (112, 48, 112, True),
        (128, 40, 96, False),
        (96, 56, 128, False),
        (136, 56, 136, True),
        (152, 72, 152, True),
        (160, 64, 160, True),
        (96, 64, 96, True),
        (80, 48, 80, True),
        (136, 32, 136, False),
        (80, 32, 80, False),
        (160, 96, 160, True),
    )

    def test_preview_and_json_use_identical_absolute_ranges(self):
        with tempfile.TemporaryDirectory() as directory:
            for first, spacing, second, linked in self.CASES:
                with self.subTest(values=(first, spacing, second, linked)):
                    geometry = centered_square_sst_geometry(
                        first, spacing, second)
                    preview = preview_seed_partition(
                        spacing, first, second, linked)
                    self.assertEqual(preview["sst_layer_ranges"],
                                     geometry["layer_ranges"])
                    self.assertEqual(preview["seed_partition_ranges"],
                                     geometry["seed_partition_ranges"])
                    target = write_shifted_sst(
                        str(Path(directory) / (
                            "%d_%d_%d.json" % (first, spacing, second))),
                        first, spacing, second,
                        layers_design_sequence_identical=linked)
                    payload = json.loads(target.read_text(encoding="utf-8"))
                    layout = payload["moire_structure_metadata"][
                        "variable_length_layout"]
                    self.assertEqual(layout["layer_ranges"],
                                     geometry["layer_ranges"])
                    self.assertEqual(layout["sst_scaffold_ranges"],
                                     geometry["scaffold_ranges"])
                    self.assertEqual(layout[
                        "sst_complementary_chain_ranges"],
                        geometry["complement_ranges"])
                    validation = validate_sst(str(target))
                    self.assertTrue(validation["valid"], validation["errors"])
                    self.assertTrue(validation["complete_32nt_units"])

    def test_asymmetric_long_layers_keep_z2_centered(self):
        geometry = centered_square_sst_geometry(192, 64, 128)
        self.assertEqual(
            geometry["seed_partition_lengths_bp"], [112, 64, 112])
        self.assertEqual(geometry["optimized_seed_overlap_bp"], [112, 112])
        self.assertEqual(geometry["maximin_overlap_bp"], 112)
        self.assertEqual(
            geometry["seed_capture_positions_by_layer"], [
                [120, 136, 152, 168, 184, 200, 216],
                [296, 312, 328, 344, 360, 376, 392],
            ])

    def test_symmetric_short_layers_use_center_as_the_tiebreak(self):
        geometry = centered_square_sst_geometry(80, 64, 80)
        self.assertEqual(
            geometry["seed_partition_lengths_bp"], [112, 64, 112])
        self.assertEqual(geometry["optimized_seed_overlap_bp"], [80, 80])
        self.assertEqual(geometry["envelope_center_offset_bp"], 0.0)

    def test_short_layer_reassigns_unused_support_to_long_layer(self):
        cases = (
            (80, 64, 400, [80, 64, 144], [80, 144]),
            (400, 64, 80, [144, 64, 80], [144, 80]),
            (192, 64, 64, [160, 64, 64], [160, 64]),
        )
        for first, spacing, second, partition, overlap in cases:
            with self.subTest(values=(first, spacing, second)):
                geometry = centered_square_sst_geometry(
                    first, spacing, second)
                self.assertEqual(
                    geometry["seed_partition_lengths_bp"], partition)
                self.assertEqual(
                    geometry["optimized_seed_overlap_bp"], overlap)

    def test_maximum_spacing_always_leaves_64_bp_per_seed_support(self):
        # The immutable Seed is 288 bp long.  At the maximum 160-bp Z2,
        # max-min placement must divide the remaining 128 bp equally.  Since
        # both SST layers have a 64-bp minimum, the actual overlap is always
        # exactly 64/64, regardless of linkage or unequal layer lengths.
        for first, second in (
                (64, 64), (80, 400), (400, 80), (192, 128), (400, 400)):
            with self.subTest(first=first, second=second):
                geometry = centered_square_sst_geometry(first, 160, second)
                self.assertEqual(
                    geometry["seed_partition_lengths_bp"], [64, 160, 64])
                self.assertEqual(
                    geometry["optimized_seed_overlap_bp"], [64, 64])
                self.assertEqual(geometry["maximin_overlap_bp"], 64)

    def test_maximin_translation_preserves_reviewed_template_phase(self):
        for first, spacing, second in (
                (192, 64, 128), (80, 160, 88), (400, 160, 80)):
            with self.subTest(values=(first, spacing, second)):
                geometry = centered_square_sst_geometry(
                    first, spacing, second)
                self.assertTrue(geometry["template_phase_preserved"])
                self.assertEqual(geometry["coordinate_shift_bp"] % 32, 0)
                origin = geometry["capture_phase_reference_origin"]
                for index, layer in enumerate(geometry["layer_ranges"]):
                    scaffold, complement = \
                        complete_square_polymer_ranges(layer)
                    self.assertEqual(
                        geometry["scaffold_ranges"][index], scaffold)
                    self.assertEqual(
                        geometry["complement_ranges"][index], complement)
                    self.assertTrue(all(
                        (position-origin) % 16 == 0
                        for position in geometry[
                            "seed_capture_positions_by_layer"][index]))

    def test_every_supported_geometry_is_overlap_ranked_then_z2_centered(self):
        for first in range(64, 401, 8):
            for spacing in range(0, 161, 8):
                for second in (64, 80, 96, 128, 160, 192, 256, 400):
                    geometry = centered_square_sst_geometry(
                        first, spacing, second)
                    remaining = 288-spacing
                    candidates = []
                    for z1 in range(0, remaining+1, 8):
                        z3 = remaining-z1
                        overlap = (min(z1, first), min(z3, second))
                        candidates.append((
                            min(overlap), sum(overlap), -abs(z1-z3)))
                    expected = max(candidates)
                    z1, unused_z2, z3 = geometry[
                        "seed_partition_lengths_bp"]
                    actual = (
                        geometry["maximin_overlap_bp"],
                        geometry["total_overlap_bp"], -abs(z1-z3))
                    self.assertEqual(
                        actual, expected,
                        (first, spacing, second,
                         geometry["seed_partition_lengths_bp"]))


if __name__ == "__main__":
    unittest.main()
