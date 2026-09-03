"""Production-path regressions for low-spacing auxiliary SST routing."""

import json
from pathlib import Path
import tempfile
import unittest

from moire_design_core.structure import (
    build_capture_ready_sst_payload,
    build_shifted_sst_payload,
    capture_site_assignments,
    finalize_structure,
    generate_scaffold_review,
    payload_to_internal_numbering,
    validate_structure,
    validate_sst,
)
from moire_design_core.kagome_sst import (
    _line_intervals,
    validate_kagome_sst_payload,
)
from moire_design_core.sst_auxiliary_routing import actual_helix


LENGTHS = (64, 80, 96, 104, 112, 120, 128, 136, 144, 152, 160)
KAGOME_SPACINGS = (0, 8, 16, 24, 32, 40, 48, 56, 64, 96, 160)


def _assert_reciprocal(test, payload):
    rows = {int(row["num"]): row for row in payload["vstrands"]}
    for field in ("scaf", "stap"):
        for helix, row in rows.items():
            for base, record in enumerate(row[field]):
                if record == [-1, -1, -1, -1]:
                    continue
                for slot in (0, 2):
                    partner, partner_base = map(
                        int, record[slot:slot + 2])
                    if partner < 0:
                        continue
                    test.assertIn(partner, rows)
                    reverse = rows[partner][field][partner_base]
                    test.assertIn([helix, base],
                                  (reverse[0:2], reverse[2:4]))


def _assert_cadnano_color_anchors(test, payload):
    """Every legacy colour marker must resolve to an existing strand."""
    empty = [-1, -1, -1, -1]
    for row in payload["vstrands"]:
        for topology_name, colors_name in (
                ("stap", "stap_colors"), ("scaf", "scaf_colors")):
            topology = row.get(topology_name, [])
            for item in row.get(colors_name, []):
                test.assertGreaterEqual(len(item), 2)
                index = int(item[0])
                test.assertGreaterEqual(index, 0)
                test.assertLess(index, len(topology))
                test.assertNotEqual(topology[index], empty)


def _sst_staple_edges(payload):
    """Return every SST--SST staple edge, including linear connections."""
    internal = payload_to_internal_numbering(payload)
    rows = {int(row["num"]): row for row in internal["vstrands"]}
    edges = set()
    for helix, row in rows.items():
        if helix < 48:
            continue
        for base, record in enumerate(row["stap"]):
            for slot in (0, 2):
                partner, partner_base = map(
                    int, record[slot:slot + 2])
                if partner >= 48:
                    edges.add(tuple(sorted(
                        ((helix, base), (partner, partner_base)))))
    return edges


def _seed_sst_staple_edges(payload):
    internal = payload_to_internal_numbering(payload)
    rows = {int(row["num"]): row for row in internal["vstrands"]}
    edges = set()
    for helix in range(48):
        for base, record in enumerate(rows[helix]["stap"]):
            for slot in (0, 2):
                partner, partner_base = map(
                    int, record[slot:slot + 2])
                if partner >= 48:
                    edges.add(((helix, base),
                               (partner, partner_base)))
    return edges


class ProductionSSTAuxiliaryMatrixTests(unittest.TestCase):

    def test_all_lattices_keep_centered_seed_and_only_physical_gaps(self):
        """Square, Kagome and mixed designs share one placement contract."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for lattice in ("square", "kagome", "square_kagome"):
                with self.subTest(lattice=lattice):
                    complete = build_shifted_sst_payload(
                        "complete.json", False,
                        z1_bp=192, z2_bp=96, z3_bp=88,
                        seed_z1_bp=128, seed_z3_bp=128,
                        lattice_type=lattice,
                        layers_design_sequence_identical=False)
                    complete_layout = complete[
                        "moire_structure_metadata"][
                            "variable_length_layout"]
                    geometry = complete_layout["square_centered_geometry"]
                    self.assertEqual(
                        complete_layout["coordinate_shift_bp"],
                        geometry["coordinate_shift_bp"])
                    self.assertEqual(
                        complete_layout["seed_layer_ranges"],
                        geometry["seed_layer_ranges"])

                    ready = build_capture_ready_sst_payload(
                        complete, "capture_ready.json")
                    ready_path = root / (lattice + "_ready.json")
                    scaffold_path = root / (lattice + "_scaffold.json")
                    final_path = root / (lattice + "_final.json")
                    ready_path.write_text(
                        json.dumps(ready), encoding="utf-8")
                    generate_scaffold_review(
                        str(scaffold_path), str(ready_path))
                    finalize_structure(str(scaffold_path), str(final_path))
                    final = json.loads(final_path.read_text(encoding="utf-8"))
                    report = validate_structure(
                        str(final_path), require_staples=True)
                    self.assertTrue(report["valid"], report["errors"])
                    self.assertFalse(report["capture_anchor_missing"])

                    final_layout = final["moire_structure_metadata"][
                        "variable_length_layout"]
                    self.assertEqual(
                        final_layout["coordinate_shift_bp"],
                        geometry["coordinate_shift_bp"])
                    self.assertEqual(
                        final_layout["seed_layer_ranges"],
                        geometry["seed_layer_ranges"])
                    self.assertEqual(
                        final_layout["overlap_ranges"],
                        geometry["optimized_overlap_ranges"])

                    removed = (_sst_staple_edges(complete) -
                               _sst_staple_edges(final))
                    bridges = _seed_sst_staple_edges(final)
                    physical_sst_nodes = {
                        sst_node for unused_seed, sst_node in bridges}
                    # Square removes a complete U edge whose two endpoints
                    # both become captures.  Kagome can instead start from a
                    # nick, crossover or linear connection, so requiring at
                    # least one physical endpoint is the shared invariant.
                    self.assertTrue(all(
                        any(node in physical_sst_nodes for node in edge)
                        for edge in removed))
                    # A Kagome linear capture opens the edge from the capture
                    # base to its immediate predecessor; that neighbour is
                    # intentionally not another capture base.  Every removed
                    # edge must nevertheless touch the exact physical target
                    # and may span at most one base when it is longitudinal.
                    for edge in removed:
                        target_nodes = [
                            node for node in edge
                            if node in physical_sst_nodes]
                        self.assertTrue(target_nodes)
                        if edge[0][1] != edge[1][1]:
                            self.assertEqual(
                                abs(edge[0][1] - edge[1][1]), 1)

    def test_square_zero_and_eight_spacing_matrix(self):
        cases = [(first, spacing, second)
                 for first in LENGTHS for spacing in (0, 8)
                 for second in (first, 128)]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for index, values in enumerate(cases):
                with self.subTest(values=values):
                    payload = build_shifted_sst_payload(
                        "square.json", False,
                        z1_bp=values[0], z2_bp=values[1], z3_bp=values[2],
                        lattice_type="square",
                        layers_design_sequence_identical=False)
                    path = root / ("square_%d.json" % index)
                    path.write_text(json.dumps(payload), encoding="utf-8")
                    report = validate_sst(str(path))
                    self.assertTrue(report["valid"], report["errors"])
                    self.assertTrue(payload["moire_structure_metadata"]
                                    ["auxiliary_sst_routing"]["enabled"])
                    _assert_reciprocal(self, payload)

    def test_kagome_242_parameter_combinations(self):
        cases = [(first, spacing, second)
                 for first in LENGTHS for spacing in KAGOME_SPACINGS
                 for second in (first, 128)]
        template_phase_edge_cases = []
        for values in cases:
            with self.subTest(values=values):
                payload = build_shifted_sst_payload(
                    "kagome.json", False,
                    z1_bp=values[0], z2_bp=values[1], z3_bp=values[2],
                    lattice_type="kagome",
                    layers_design_sequence_identical=False)
                metadata = payload["moire_structure_metadata"]
                ranges = metadata["sst_duplex_ranges"]
                self.assertEqual(ranges[0][1]-ranges[0][0]+1, values[0])
                self.assertEqual(ranges[1][1]-ranges[1][0]+1, values[2])
                self.assertEqual(ranges[1][0]-ranges[0][1]-1, values[1])
                shifts = metadata["linear_nick_phase_shifts_by_layer"]
                self.assertFalse(any(
                    value for layer in shifts for value in layer.values()))
                internal = payload_to_internal_numbering(payload)
                report = validate_kagome_sst_payload(internal)
                self.assertTrue(report["valid"], report["errors"])
                if report["components"]["stap"]["length_counts"].get(16):
                    template_phase_edge_cases.append(values)
                self.assertNotIn(
                    16, report["components"]["scaf"]["length_counts"])
                _assert_reciprocal(self, payload)
        # A 16-nt complete-SST component is permitted only for the immutable
        # template's one 64-nt non-scaffold topology: 16+32+16 becomes 48+16.
        # Several spacing combinations can place that same topology on one
        # layer, so assert its structural condition rather than hard-coding a
        # short parameter list.  Nick phase remains zero in every case.
        self.assertIn((64, 16, 64), template_phase_edge_cases)
        self.assertIn((64, 48, 64), template_phase_edge_cases)
        self.assertTrue(all(first == 64 or second == 64
                            for first, unused_spacing, second in
                            template_phase_edge_cases))

    def test_kagome_identical_layers_keep_shared_32_bp_phase(self):
        payload = build_shifted_sst_payload(
            "kagome_linked.json", False,
            z1_bp=128, z2_bp=32, z3_bp=128,
            lattice_type="kagome",
            layers_design_sequence_identical=True)
        shifts = payload["moire_structure_metadata"][
            "linear_nick_phase_shifts_by_layer"]
        self.assertFalse(any(value for layer in shifts
                             for value in layer.values()))
        for spacing in (16, 48):
            with self.subTest(spacing=spacing):
                with self.assertRaisesRegex(ValueError, "32-bp"):
                    build_shifted_sst_payload(
                        "kagome_linked_bad.json", False,
                        z1_bp=64, z2_bp=spacing, z3_bp=64,
                        lattice_type="kagome",
                        layers_design_sequence_identical=True)

    def test_kagome_template_phased_48_plus_16_reaches_final_structure(self):
        cases = ((64, 16, 64), (64, 48, 64))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for index, values in enumerate(cases):
                with self.subTest(values=values):
                    complete = build_shifted_sst_payload(
                        "kagome.json", False,
                        z1_bp=values[0], z2_bp=values[1],
                        z3_bp=values[2], lattice_type="kagome",
                        layers_design_sequence_identical=False)
                    ready = build_capture_ready_sst_payload(
                        complete, "kagome_ready.json")
                    ready_path = root / ("kagome_ready_%d.json" % index)
                    scaffold_path = root / (
                        "kagome_scaffold_%d.json" % index)
                    final_path = root / ("kagome_final_%d.json" % index)
                    ready_path.write_text(
                        json.dumps(ready), encoding="utf-8")
                    generate_scaffold_review(
                        str(scaffold_path), str(ready_path))
                    finalize_structure(str(scaffold_path), str(final_path))
                    report = validate_structure(
                        str(final_path), require_staples=True)
                    self.assertTrue(report["valid"], report["errors"])
                    self.assertFalse(report["capture_mapping_missing"])
                    self.assertFalse(report["capture_mapping_unexpected"])
                    final = json.loads(final_path.read_text(encoding="utf-8"))
                    _assert_cadnano_color_anchors(self, final)

    def test_mixed_square_kagome_layers_and_capture_ready_aux(self):
        cases = ((64, 0, 64), (88, 8, 88),
                 (128, 32, 128), (96, 56, 128))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for case_index, values in enumerate(cases):
              with self.subTest(values=values):
                complete = build_shifted_sst_payload(
                    "mixed.json", False,
                    z1_bp=values[0], z2_bp=values[1], z3_bp=values[2],
                    lattice_type="square_kagome",
                    layers_design_sequence_identical=False)
                metadata = complete["moire_structure_metadata"]
                layout = metadata["variable_length_layout"]
                self.assertEqual(metadata["lattice_by_layer"],
                                 ["square", "kagome"])
                self.assertEqual(layout["actual_spacing_bp"], values[1])
                self.assertEqual(layout["layer_ranges"][0][1] -
                                 layout["layer_ranges"][0][0] + 1,
                                 values[0])
                self.assertEqual(layout["layer_ranges"][1][1] -
                                 layout["layer_ranges"][1][0] + 1,
                                 values[2])
                complete_path = root / ("mixed_%d.json" % case_index)
                complete_path.write_text(
                    json.dumps(complete), encoding="utf-8")
                sst_report = validate_sst(str(complete_path))
                self.assertTrue(sst_report["valid"], sst_report["errors"])
                ready = build_capture_ready_sst_payload(
                    complete, "mixed_ready.json")
                self.assertTrue(ready["moire_structure_metadata"]
                                ["capture_gaps_reserved"])
                _assert_reciprocal(self, ready)
                ready_path = root / ("mixed_ready_%d.json" % case_index)
                scaffold_path = root / ("mixed_scaffold_%d.json" %
                                        case_index)
                final_path = root / ("mixed_final_%d.json" % case_index)
                ready_path.write_text(json.dumps(ready), encoding="utf-8")
                generate_scaffold_review(
                    str(scaffold_path), str(ready_path))
                finalize_structure(str(scaffold_path), str(final_path))
                final_report = validate_structure(
                    str(final_path), require_staples=True)
                self.assertTrue(final_report["valid"],
                                final_report["errors"])
                final_payload = json.loads(final_path.read_text())
                final_metadata = final_payload[
                    "moire_structure_metadata"]
                self.assertEqual(final_metadata["lattice_type"],
                                 "square_kagome")
                self.assertEqual(final_metadata["lattice_by_layer"],
                                 ["square", "kagome"])

    def test_mixed_square_kagome_242_parameter_combinations(self):
        cases = [(first, spacing, second)
                 for first in LENGTHS for spacing in KAGOME_SPACINGS
                 for second in (first, 128)]
        for values in cases:
            with self.subTest(values=values):
                payload = build_shifted_sst_payload(
                    "mixed.json", False,
                    z1_bp=values[0], z2_bp=values[1], z3_bp=values[2],
                    lattice_type="square_kagome",
                    layers_design_sequence_identical=False)
                metadata = payload["moire_structure_metadata"]
                ranges = metadata["sst_duplex_ranges"]
                self.assertEqual(ranges[0][1] - ranges[0][0] + 1,
                                 values[0])
                self.assertEqual(ranges[1][1] - ranges[1][0] + 1,
                                 values[2])
                self.assertEqual(ranges[1][0] - ranges[0][1] - 1,
                                 values[1])
                self.assertEqual(metadata["lattice_by_layer"],
                                 ["square", "kagome"])
                _assert_reciprocal(self, payload)

    def test_mixed_kagome_layer_keeps_template_phased_64_and_72_edges(self):
        expected = {
            # At this absolute phase the 64-nt Kagome non-scaffold range is
            # already 32+32; it must not be changed merely because another
            # 64-nt phase has the exceptional 48+16 topology.
            (64, 16, 64): [32, 32],
            # This is the one 64-nt topology that clips as 16+32+16; only one
            # side is merged and the immutable result is 48+16.
            (64, 48, 64): [48, 16],
            # A 72-bp Kagome duplex uses an 80-nt non-scaffold range.  It has
            # ordinary 32+48 components and no 16/24-nt exception.
            (72, 16, 72): [32, 48],
            (72, 48, 72): [48, 32],
        }
        for values, component_lengths in expected.items():
            with self.subTest(values=values):
                payload = build_shifted_sst_payload(
                    "mixed-edge.json", False,
                    z1_bp=values[0], z2_bp=values[1], z3_bp=values[2],
                    lattice_type="square_kagome",
                    layers_design_sequence_identical=False)
                metadata = payload["moire_structure_metadata"]
                layout = metadata["variable_length_layout"]
                self.assertEqual(metadata["lattice_by_layer"],
                                 ["square", "kagome"])
                kagome_staple_range = layout["staple_ranges"][1]
                intervals = _line_intervals(
                    kagome_staple_range[0], kagome_staple_range[1], 8)
                self.assertEqual(
                    [high-low+1 for low, high in intervals],
                    component_lengths)
                self.assertNotIn(24, component_lengths)
                self.assertEqual(layout["actual_spacing_bp"], values[1])
                _assert_reciprocal(self, payload)
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "mixed-edge.json"
                    path.write_text(json.dumps(payload), encoding="utf-8")
                    report = validate_sst(str(path))
                self.assertTrue(report["valid"], report["errors"])

    def test_zero_and_eight_spacing_capture_edges_follow_auxiliary_paths(self):
        cases = (("square", (64, 0, 64)),
                 ("square", (88, 8, 88)),
                 ("kagome", (64, 0, 64)),
                 ("kagome", (88, 8, 88)))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for case_index, (lattice, values) in enumerate(cases):
                with self.subTest(lattice=lattice, values=values):
                    payload = build_shifted_sst_payload(
                        "ready.json", True,
                        z1_bp=values[0], z2_bp=values[1],
                        z3_bp=values[2], lattice_type=lattice,
                        layers_design_sequence_identical=False)
                    internal = payload_to_internal_numbering(payload)
                    layout = internal["moire_structure_metadata"][
                        "variable_length_layout"]
                    assignments = capture_site_assignments(layout)
                    routing = layout["auxiliary_sst_routing"]
                    self.assertTrue(routing["enabled"])
                    self.assertTrue(any(
                        component["primary_auxiliary_boundary_links"]
                        for component in routing["detoured_components"]))
                    # Capture may sit exactly after the detoured single-chain
                    # interval (and therefore remain on the primary helix),
                    # or directly on an auxiliary node.  In both cases the
                    # chosen physical endpoint must be resolved from the real
                    # routed graph rather than assumed to be h48--63.
                    for assignment in assignments:
                        for bridge in assignment["bridges"]:
                            self.assertEqual(
                                int(bridge["sst_helix"]),
                                actual_helix(
                                    layout, int(assignment["layer"]),
                                    "stap",
                                    int(bridge["logical_sst_helix"]),
                                    int(assignment["position"])))
                    sst_path = root / ("low_sst_%d.json" % case_index)
                    scaffold_path = root / (
                        "low_scaffold_%d.json" % case_index)
                    final_path = root / ("low_final_%d.json" % case_index)
                    sst_path.write_text(json.dumps(payload), encoding="utf-8")
                    generate_scaffold_review(
                        str(scaffold_path), str(sst_path))
                    finalize_structure(str(scaffold_path), str(final_path))
                    report = validate_structure(
                        str(final_path), require_staples=True)
                    self.assertTrue(report["valid"], report["errors"])
                    self.assertFalse(report["capture_mapping_missing"])
                    self.assertFalse(report["capture_mapping_unexpected"])
                    final = json.loads(final_path.read_text(encoding="utf-8"))
                    _assert_cadnano_color_anchors(self, final)


if __name__ == "__main__":
    unittest.main()
