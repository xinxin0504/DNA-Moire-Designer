"""Regression tests for 1.2-preview colours in pure-cylinder exports."""

import unittest

from moire_design_core.sequence_workflow_worker import (
    CYLINDER_PREVIEW_COLORS,
    _bild_color_command,
    _cylindrical_model_bild,
    _cylindrical_preview_color,
    _project_auxiliary_geometry,
    _seed_preview_support_ranges,
)


class _VirtualHelix:
    def __init__(self, number, coord=None, part=None):
        self._number = number
        self._coord = coord
        self._part = part

    def number(self):
        return self._number

    def coord(self):
        return self._coord

    def part(self):
        return self._part


class _Strand:
    def __init__(self, number=0, coord=None, part=None):
        self._helix = _VirtualHelix(number, coord, part)

    def virtualHelix(self):
        return self._helix


def _records(last_index=7):
    strand = _Strand()
    return [{
        "strand": strand,
        "count": 1,
        "direction": 1,
        "idx": index,
        "sub": 0,
        "pos": (0.0, 0.0, float(index)),
        "a1": (1.0, 0.0, 0.0),
    } for index in range(last_index + 1)]


class CylindricalPreviewColorTests(unittest.TestCase):
    def test_seed_export_uses_current_preview_partition_not_reference_ranges(self):
        layout = {
            "seed_partition_ranges": [[80, 167], [168, 271], [272, 367]],
            # Frozen reference bookkeeping; this is deliberately different.
            "seed_layer_ranges": [[80, 207], [240, 367]],
        }
        self.assertEqual(
            _seed_preview_support_ranges(layout),
            [[80, 167], [272, 367]])

    def test_old_design_falls_back_to_legacy_seed_support_ranges(self):
        self.assertEqual(
            _seed_preview_support_ranges({
                "seed_layer_ranges": [[48, 175], [208, 335]],
            }),
            [[48, 175], [208, 335]])

    def test_auxiliary_projection_moves_dat_and_atomic_coordinates_only(self):
        class _Part:
            def radius(self):
                return 1.0

            def latticeCoordToPositionXY(self, row, col):
                return float(row), float(col)

        part = _Part()
        strand = _Strand(64, (10, 0), part)
        record = {
            "strand": strand,
            "pos": (20.0, 4.0, 1.0),
            "output_pos": (25.0, 7.0, 2.0),
        }
        payload = {
            "vstrands": [
                {"num": 0, "row": 0, "col": 0},
                {"num": 64, "row": 10, "col": 0},
            ],
            "moire_structure_metadata": {
                "helix_numbering": "sst_first",
                "auxiliary_sst_routing": {"enabled": True},
            },
        }
        count = _project_auxiliary_geometry([record], payload)
        self.assertEqual(count, 1)
        self.assertEqual(record["cadnano_helix"], 64)
        self.assertEqual(record["ideal_geometry_helix"], 0)
        # PDB/mmCIF ``output_pos`` and oxDNA/BILD ``pos`` receive the same
        # spatial translation; their original centring offset is preserved.
        self.assertAlmostEqual(
            record["output_pos"][0] - record["pos"][0], 5.0)
        self.assertAlmostEqual(
            record["output_pos"][1] - record["pos"][1], 3.0)
        self.assertEqual(record["output_pos"][2] - record["pos"][2], 1.0)

    def test_seed_regions_match_preview_palette(self):
        ranges = ((0, 2), (5, 7))
        self.assertEqual(CYLINDER_PREVIEW_COLORS["seed_z2"], "#d9dee3")
        self.assertEqual(
            _cylindrical_preview_color("seed", 1, ranges),
            CYLINDER_PREVIEW_COLORS["seed_z1"])
        self.assertEqual(
            _cylindrical_preview_color("seed", 3, ranges),
            CYLINDER_PREVIEW_COLORS["seed_z2"])
        self.assertEqual(
            _cylindrical_preview_color("seed", 6, ranges),
            CYLINDER_PREVIEW_COLORS["seed_z3"])

        bild = _cylindrical_model_bild(
            _records(), "seed", radial_subdivisions=16,
            seed_support_ranges=ranges)
        for region in ("seed_z1", "seed_z2", "seed_z3"):
            self.assertIn(
                _bild_color_command(CYLINDER_PREVIEW_COLORS[region]), bild)

    def test_each_sst_layer_uses_its_preview_colour(self):
        layer_1 = _cylindrical_model_bild(
            _records(2), "sst_layer_1", radial_subdivisions=16)
        layer_2 = _cylindrical_model_bild(
            _records(2), "sst_layer_2", radial_subdivisions=16)
        self.assertIn(
            _bild_color_command(CYLINDER_PREVIEW_COLORS["sst_layer_1"]),
            layer_1)
        self.assertNotIn(
            _bild_color_command(CYLINDER_PREVIEW_COLORS["sst_layer_2"]),
            layer_1)
        self.assertIn(
            _bild_color_command(CYLINDER_PREVIEW_COLORS["sst_layer_2"]),
            layer_2)
        self.assertNotIn(
            _bild_color_command(CYLINDER_PREVIEW_COLORS["sst_layer_1"]),
            layer_2)


if __name__ == "__main__":
    unittest.main()
