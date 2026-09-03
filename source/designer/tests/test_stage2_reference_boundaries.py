"""Hard boundaries between two-layer routing and capture-only references."""

import unittest
import json

from moire_design_core.structure import (
    KAGOME_CAPTURE_TOPOLOGY_REFERENCE,
    SEED_CAPTURE_REFERENCE,
    SEED_ROUTING_REFERENCE,
    SQUARE_CAPTURE_REFERENCE,
)
from moire_design_core.structure_rules import assert_supported_seed_preset
from moire_design_core.structure_worker import (
    _reference_seed_scaffold_payload,
    capacity,
)


class Stage2ReferenceBoundaryTests(unittest.TestCase):

    def test_only_two_layer_square_file_can_route_seed(self):
        self.assertEqual(SEED_CAPTURE_REFERENCE, SEED_ROUTING_REFERENCE)
        self.assertEqual(SQUARE_CAPTURE_REFERENCE, SEED_ROUTING_REFERENCE)
        self.assertEqual(SEED_ROUTING_REFERENCE.name,
                         "Square_Seed_2L_newtemplate.json")
        self.assertNotIn("3L", SEED_ROUTING_REFERENCE.name)

    def test_kagome_three_layer_file_is_capture_catalogue_only(self):
        self.assertEqual(KAGOME_CAPTURE_TOPOLOGY_REFERENCE.name,
                         "Kagome_Seed_Ka-seed-pore_3L.json")
        self.assertNotEqual(KAGOME_CAPTURE_TOPOLOGY_REFERENCE,
                            SEED_ROUTING_REFERENCE)

    def test_only_supported_seed_preset_is_accepted(self):
        assert_supported_seed_preset("s8_r4x4")
        assert_supported_seed_preset("8x8_minus_4x4_pore")
        with self.assertRaises(ValueError):
            assert_supported_seed_preset("unsupported_cross_section")

    def test_frozen_two_layer_reference_has_two_capacity_safe_scaffolds(self):
        payload = json.loads(SEED_ROUTING_REFERENCE.read_text(encoding="utf-8"))
        rows = {int(row["num"]): row for row in payload["vstrands"]
                if 0 <= int(row["num"]) < 48}
        nodes = set()
        graph = {}
        for helix, row in rows.items():
            for base, record in enumerate(row["scaf"]):
                if record == [-1, -1, -1, -1]:
                    continue
                node = (helix, base)
                nodes.add(node)
                graph.setdefault(node, set())
                for slot in (0, 2):
                    if 0 <= int(record[slot]) < 48:
                        graph[node].add((int(record[slot]),
                                         int(record[slot + 1])))
        lengths = []
        unseen = set(nodes)
        while unseen:
            stack = [next(iter(unseen))]
            component = set()
            while stack:
                node = stack.pop()
                if node in component:
                    continue
                component.add(node)
                stack.extend(graph.get(node, ()) - component)
            unseen -= component
            lengths.append(len(component))
        self.assertEqual(sorted(lengths), [7300, 7336])

    def test_legacy_seed_lengths_cannot_change_fixed_reference(self):
        payload, lengths = _reference_seed_scaffold_payload({
            "z1_bp": 120,
            "z2_bp": 32,
            "z3_bp": 120,
        })
        self.assertEqual(lengths, [7300, 7336])
        self.assertTrue(payload["moire_edge_metadata"][
            "seed_routing_is_frozen_reference"])
        self.assertEqual(payload["moire_edge_metadata"][
            "seed_z1_edge_growth_bp"], 0)
        self.assertEqual(payload["moire_edge_metadata"][
            "seed_z3_edge_growth_bp"], 0)

    def test_capacity_uses_reference_transform_not_retired_router(self):
        report = capacity()
        self.assertEqual(report["seed_scaffold_count"], 2)
        self.assertEqual(report["planned_balanced_lengths"], [7300, 7336])
        self.assertEqual(report["seed_routing_reference"],
                         "Square_Seed_2L_newtemplate.json")


if __name__ == "__main__":
    unittest.main()
