import unittest

from cadnano2.model.guideddesign import (base_count_for_width,
                                         boolean_runs,
                                         cross_section_coords,
                                         estimate_scaffold_length,
                                         profile_count_for_height)


class GuidedDesignTests(unittest.TestCase):
    def test_physical_dimensions(self):
        self.assertEqual(base_count_for_width(34), 100)
        self.assertEqual(profile_count_for_height(28), 10)

    def test_honeycomb_z_layers_alternate_without_reversing_profiles(self):
        coords = cross_section_coords(3, 4, 'honeycomb', 'z')
        by_layer = {}
        for row, col, profile, layer in coords:
            by_layer.setdefault(layer, []).append((row, col, profile))
        self.assertEqual([item[2] for item in by_layer[0]], [0, 1, 2])
        self.assertEqual([item[2] for item in by_layer[1]], [0, 1, 2])
        self.assertEqual([item[0] for item in by_layer[0]], [1, 1, 1])
        self.assertEqual([item[0] for item in by_layer[1]], [2, 2, 2])

    def test_boolean_scanline_runs_and_length(self):
        runs = boolean_runs([False, True, True, False, True, True], 12)
        self.assertEqual(runs, [(2, 5), (8, 11)])
        self.assertEqual(estimate_scaffold_length([runs], 3), 24)


if __name__ == '__main__':
    unittest.main()
