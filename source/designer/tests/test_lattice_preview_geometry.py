import math
import unittest

from moire_designer.lattice_preview_geometry import (
    lattice_graph, rotated_graph_outside_square)


class RotatedSeedFootprintTests(unittest.TestCase):
    def test_exclusion_footprint_rotates_with_layer(self):
        points, edges = lattice_graph("square", 2.8, 30.8)
        angle_deg = 45.0
        half_extent = 12.2
        visible, unused_edges = rotated_graph_outside_square(
            points, edges, angle_deg, half_extent)

        angle = math.radians(-angle_deg)
        cosine, sine = math.cos(angle), math.sin(angle)
        for x, y in visible:
            local_x = cosine*x-sine*y
            local_y = sine*x+cosine*y
            self.assertFalse(
                abs(local_x) <= half_extent and
                abs(local_y) <= half_extent)

    def test_large_twist_keeps_same_number_of_visible_nodes(self):
        points, edges = lattice_graph("square", 2.8, 30.8)
        counts = []
        for angle_deg in (0.0, 20.0, 35.0, 45.0):
            visible, unused_edges = rotated_graph_outside_square(
                points, edges, angle_deg, 12.2)
            counts.append(len(visible))
        self.assertEqual(len(set(counts)), 1)


if __name__ == "__main__":
    unittest.main()
