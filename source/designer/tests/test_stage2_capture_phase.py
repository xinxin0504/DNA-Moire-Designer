"""Capture columns retain the immutable 128-bp phase at every SST length."""

import json
from pathlib import Path
import tempfile
import unittest

from moire_design_core.structure import (
    build_capture_ready_sst_payload,
    build_shifted_sst_payload,
    finalize_structure,
    fixed_seed_overlap_layout,
    generate_scaffold_review,
    payload_to_internal_numbering,
    structure_layout,
    validate_structure,
)


def _sst_crossover_edges(payload):
    internal = payload_to_internal_numbering(payload)
    rows = {int(row["num"]): row for row in internal["vstrands"]}
    edges = set()
    for helix, row in rows.items():
        if helix < 48:
            continue
        for base, record in enumerate(row["stap"]):
            for offset in (0, 2):
                partner, partner_base = map(
                    int, record[offset:offset + 2])
                if partner >= 48 and partner != helix:
                    edges.add(tuple(sorted(
                        ((helix, base), (partner, partner_base)))))
    return edges


def _seed_sst_edges(payload):
    internal = payload_to_internal_numbering(payload)
    rows = {int(row["num"]): row for row in internal["vstrands"]}
    edges = set()
    for helix in range(48):
        for base, record in enumerate(rows[helix]["stap"]):
            for offset in (0, 2):
                partner, partner_base = map(
                    int, record[offset:offset + 2])
                if partner >= 48:
                    edges.add(((helix, base), (partner, partner_base)))
    return edges


class Stage2CapturePhaseTests(unittest.TestCase):

    def test_square_128bp_support_follows_each_current_duplex(self):
        canonical = fixed_seed_overlap_layout(
            [[48, 175], [208, 335]])
        shifted = fixed_seed_overlap_layout(
            [[96, 223], [256, 383]])
        # Both layers remain exactly 128 bp.  The current duplex establishes
        # support, while the immutable Seed template still decides whether a
        # physical contact exists at each supported coordinate.
        self.assertEqual(canonical["capture_columns_by_layer"], [8, 8])
        self.assertEqual(shifted["capture_columns_by_layer"], [8, 5])
        self.assertNotEqual(canonical["actual_capture_positions"],
                            shifted["actual_capture_positions"])
        self.assertIn("real duplex start-to-end support",
                      shifted["capture_count_semantics"])

    def test_square_capture_grid_is_never_rephased_by_length(self):
        for length in (96, 104, 112, 120, 128, 136, 160):
            with self.subTest(length=length):
                layout = structure_layout(length, 32, length)
                origin = int(layout["capture_grid_origin"])
                self.assertEqual(origin % 16, 8)
                for layer in layout["capture_positions_by_layer"]:
                    self.assertTrue(layer)
                    low, high = layout["layer_ranges"][
                        layout["capture_positions_by_layer"].index(layer)]
                    self.assertTrue(all(low <= int(position) <= high
                                        for position in layer))
                    self.assertTrue(all(
                        (int(position) - origin) % 16 == 0
                        for position in layer))

    def test_short_square_sst_still_keeps_two_capture_pairs(self):
        for length in (64, 72, 80):
            with self.subTest(length=length):
                layout = structure_layout(length, 32, length)
                self.assertTrue(all(count >= 2 for count in
                                    layout["pair_count_by_layer"]))

    def test_accepted_layout_keeps_two_pairs_per_layer(self):
        for length in (96, 128):
            layout = structure_layout(length, 32, length)
            self.assertTrue(all(count >= 2 for count in
                                layout["pair_count_by_layer"]))

    def test_all_domain_aligned_square_lengths_use_complete_u_router(self):
        for length in (64, 72, 80, 96, 104, 112, 120, 128, 136, 160):
            with self.subTest(length=length):
                payload = build_shifted_sst_payload(
                    "frozen.json", z1_bp=length, z2_bp=32,
                    z3_bp=length, lattice_type="square")
                self.assertEqual(
                    payload["moire_structure_metadata"]
                    ["variable_length_layout"]["z1_bp"], length)
        for length in (65, 100):
            with self.subTest(length=length), self.assertRaises(ValueError):
                build_shifted_sst_payload(
                    "disabled.json", z1_bp=length, z2_bp=32,
                    z3_bp=length, lattice_type="square")

    def test_capture_ready_opens_only_two_u_edges_per_physical_column(self):
        complete = build_shifted_sst_payload(
            "complete.json", False,
            z1_bp=248, z2_bp=160, z3_bp=64,
            seed_z1_bp=128, seed_z3_bp=128,
            lattice_type="square",
            layers_design_sequence_identical=False)
        ready = build_capture_ready_sst_payload(
            complete, "capture_ready.json")
        layout = ready["moire_structure_metadata"]["variable_length_layout"]
        removed = (_sst_crossover_edges(complete) -
                   _sst_crossover_edges(ready))
        actual_positions = {
            int(position)
            for layer in layout["capture_positions_by_layer"]
            for position in layer}
        self.assertEqual(len(removed), 2 * len(actual_positions))
        self.assertEqual(
            {node[1] for edge in removed for node in edge},
            actual_positions)
        self.assertEqual(
            ready["moire_structure_metadata"]
            ["capture_gap_endpoint_count"],
            4 * len(actual_positions))

    def test_final_square_has_no_unconnected_or_outside_overlap_sst_gaps(self):
        complete = build_shifted_sst_payload(
            "complete.json", False,
            z1_bp=248, z2_bp=160, z3_bp=64,
            seed_z1_bp=128, seed_z3_bp=128,
            lattice_type="square",
            layers_design_sequence_identical=False)
        ready = build_capture_ready_sst_payload(
            complete, "capture_ready.json")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ready_path = root / "capture_ready.json"
            scaffold_path = root / "scaffold.json"
            final_path = root / "final.json"
            ready_path.write_text(json.dumps(ready), encoding="utf-8")
            generate_scaffold_review(
                str(scaffold_path), str(ready_path))
            finalize_structure(str(scaffold_path), str(final_path))
            final = json.loads(final_path.read_text(encoding="utf-8"))
            report = validate_structure(
                str(final_path), require_staples=True)

        self.assertTrue(report["valid"], report["errors"])
        removed = (_sst_crossover_edges(complete) -
                   _sst_crossover_edges(final))
        seed_sst = _seed_sst_edges(final)
        connected_sst_nodes = {sst_node for unused_seed, sst_node in seed_sst}
        self.assertEqual(len(removed), len(seed_sst) // 2)
        self.assertTrue(all(
            node in connected_sst_nodes
            for edge in removed for node in edge))
        actual_positions = {
            int(position)
            for layer in final["moire_structure_metadata"]
            ["capture_positions_by_layer"]
            for position in layer}
        self.assertEqual(
            {node[1] for edge in removed for node in edge},
            actual_positions)


if __name__ == "__main__":
    unittest.main()
