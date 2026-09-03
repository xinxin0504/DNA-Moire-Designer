"""The fixed 2L Seed must remain identical to its accepted template."""

import json
import tempfile
import unittest
from pathlib import Path

from moire_design_core.structure import (
    CAPTURE_COLUMN_COLORS,
    CAPTURE_DIRECT_POSITIONS,
    CAPTURE_PHASE_MAPPINGS,
    CAPTURE_REFERENCE_COLUMN_BY_COLOR,
    SEED_ROUTING_REFERENCE,
    capture_column_color,
    payload_to_internal_numbering,
    validate_structure,
    write_shifted_sst,
)
from moire_design_core.structure_worker import finalize, scaffold
from moire_design_core import structure_worker as runtime
from cadnano2.model.enum import LatticeType
from cadnano2.model.io.legacydecoder import import_legacy_dict


EMPTY = [-1, -1, -1, -1]


class FrozenSeedTrimTests(unittest.TestCase):

    def test_fixed_seed_preserves_complete_template_topology(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            sst = write_shifted_sst(directory / "sst.json")
            scaffold_path = directory / "scaffold.json"
            final_path = directory / "final.json"
            scaffold(str(scaffold_path), str(sst))
            finalize(str(scaffold_path), str(final_path))

            final = payload_to_internal_numbering(json.loads(
                final_path.read_text(encoding="utf-8")))
            reference = json.loads(SEED_ROUTING_REFERENCE.read_text(
                encoding="utf-8"))
            final_rows = {int(row["num"]): row
                          for row in final["vstrands"]}
            reference_rows = {int(row["num"]): row
                              for row in reference["vstrands"]}

            # Scaffold is byte-for-byte frozen on every Seed helix.
            for number in range(48):
                reference_length = len(reference_rows[number]["scaf"])
                self.assertEqual(
                    final_rows[number]["scaf"][:reference_length],
                    reference_rows[number]["scaf"])
                self.assertTrue(all(
                    record == EMPTY for record in
                    final_rows[number]["scaf"][reference_length:]))

            # Ordinary Seed staple and capture-core topology is also frozen;
            # no trim/growth/length-repair branch may touch it.
            def seed_edges(rows):
                edges = set()
                for number in range(48):
                    for base in range(64, 320):
                        record = rows[number]["stap"][base]
                        for slot in (0, 2):
                            partner, partner_base = map(
                                int, record[slot:slot + 2])
                            if 0 <= partner < 48 and partner != number:
                                edges.add(tuple(sorted(((number, base),
                                                        (partner,
                                                         partner_base)))))
                return edges

            self.assertEqual(seed_edges(final_rows),
                             seed_edges(reference_rows))
            for number in range(48):
                reference_length = len(reference_rows[number]["stap"])
                self.assertEqual(
                    final_rows[number]["stap"][:reference_length],
                    reference_rows[number]["stap"])
                self.assertTrue(all(
                    record == EMPTY for record in
                    final_rows[number]["stap"][reference_length:]))

            # Generated design colours normalize ordinary support staples to
            # black and give all eighteen immutable Capture columns their
            # own colour, including candidates and the two Z2 reserves.
            ordinary_template_colors = {0x000000, 0x60C9F6, 0xF49AE5}
            reference_markers = {
                (number, int(base)): int(color)
                for number, row in reference_rows.items()
                for base, color in row.get("stap_colors", [])}
            final_markers = {
                (number, int(base)): int(color)
                for number, row in final_rows.items()
                for base, color in row.get("stap_colors", [])}
            components, unused_labels, unused_adjacency = \
                runtime._staple_components_from_rows(reference_rows)
            layout = final["moire_structure_metadata"][
                "variable_length_layout"]
            expected_counts = {"ordinary": 0, "capture": 0}
            colors_by_column = {}
            for component in components:
                if not any(node[0] < 48 for node in component):
                    continue
                source_colors = {
                    reference_markers[node] for node in component
                    if node in reference_markers}
                self.assertEqual(len(source_colors), 1)
                source_color = next(iter(source_colors))
                generated_colors = {
                    final_markers[node] for node in component
                    if node in final_markers}
                self.assertEqual(len(generated_colors), 1)
                generated_color = next(iter(generated_colors))
                if source_color in ordinary_template_colors:
                    self.assertEqual(generated_color, 0x000000)
                    expected_counts["ordinary"] += 1
                else:
                    candidate_columns = (
                        (184, 200) if source_color == 0x999999 else
                        (CAPTURE_REFERENCE_COLUMN_BY_COLOR[source_color],))
                    matched_columns = []
                    for column in candidate_columns:
                        unit_index = (
                            int(column) - CAPTURE_DIRECT_POSITIONS[0]) // 16
                        phase = "B" if unit_index % 2 == 0 else "A"
                        capture_helices = {
                            int(seed_helix)
                            for cycle in (phase + "0", phase + "1")
                            for unused_sst, seed_helix in
                            CAPTURE_PHASE_MAPPINGS[cycle]}
                        if any((helix, int(column)) in component
                               for helix in capture_helices):
                            matched_columns.append(int(column))
                    self.assertEqual(len(matched_columns), 1)
                    column = matched_columns[0]
                    self.assertEqual(
                        generated_color, capture_column_color(column, layout))
                    colors_by_column.setdefault(column, set()).add(
                        generated_color)
                    expected_counts["capture"] += 1
            self.assertEqual(expected_counts, {
                "ordinary": 280, "capture": 144})
            self.assertEqual(len(colors_by_column), 18)
            self.assertTrue(all(len(colors) == 1
                                for colors in colors_by_column.values()))
            self.assertEqual(
                {next(iter(colors)) for colors in colors_by_column.values()},
                set(CAPTURE_COLUMN_COLORS))
            for component in components:
                if any(node[0] < 48 for node in component):
                    continue
                generated_colors = {
                    final_markers[node] for node in component
                    if node in final_markers}
                self.assertEqual(generated_colors, {0x000000})

            self.assertEqual(layout["pair_count_by_layer"], [4, 4])
            self.assertFalse(layout["seed_length_adjustment_enabled"])
            self.assertEqual(layout["seed_geometry_policy"],
                             "immutable_2L_reference")
            self.assertFalse(final["moire_structure_metadata"]["autobreak"]
                             ["rerun"])
            self.assertNotIn("capture_length_policy",
                             final["moire_structure_metadata"])

            report = validate_structure(str(final_path),
                                        require_staples=True)
            self.assertTrue(report["valid"], report["errors"])
            self.assertEqual(report["seed_scaffold_lengths"], [7300, 7336])
            self.assertEqual(report["seed_staple_missing_base_count"], 0)
            self.assertEqual(report["invalid_short_staple_count"], 0)
            self.assertEqual(report["expected_capture_bridge_count"], 64)
            self.assertEqual(report["capture_bridge_component_count"], 64)

            # A removed edge staple must not leave a color marker on an
            # empty base.  caDNAno's legacy decoder dereferences every
            # ``stap_colors`` coordinate and cannot open such a payload.
            public_payload = json.loads(final_path.read_text(
                encoding="utf-8"))
            for row in public_payload["vstrands"]:
                for base, unused_color in row.get("stap_colors", []):
                    self.assertNotEqual(row["stap"][int(base)], EMPTY)

            document = runtime.Document()
            part = import_legacy_dict(
                document, public_payload, LatticeType.Square,
                forceLatticeType=True)
            self.assertEqual(len(part.getVirtualHelices()), 64)

    def test_long_sst_translates_the_complete_frozen_seed_by_32_bp(self):
        """A canvas shift moves topology; it never regenerates the Seed."""
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            sst = write_shifted_sst(
                directory / "sst_shifted.json", 160, 96, 160,
                layers_design_sequence_identical=True)
            scaffold_path = directory / "scaffold_shifted.json"
            final_path = directory / "final_shifted.json"
            scaffold(str(scaffold_path), str(sst))
            finalize(str(scaffold_path), str(final_path))

            final = payload_to_internal_numbering(json.loads(
                final_path.read_text(encoding="utf-8")))
            reference = json.loads(SEED_ROUTING_REFERENCE.read_text(
                encoding="utf-8"))
            final_rows = {int(row["num"]): row
                          for row in final["vstrands"]}
            reference_rows = {int(row["num"]): row
                              for row in reference["vstrands"]}
            shift = 32

            def translated_record(record):
                result = list(record)
                for offset in (0, 2):
                    if int(result[offset]) >= 0:
                        result[offset + 1] = int(result[offset + 1]) + shift
                return result

            for number in range(48):
                self.assertTrue(all(
                    record == EMPTY for record in
                    final_rows[number]["scaf"][:shift]))
                for base, record in enumerate(
                        reference_rows[number]["scaf"]):
                    self.assertEqual(
                        final_rows[number]["scaf"][base + shift],
                        translated_record(record))

            layout = final["moire_structure_metadata"][
                "variable_length_layout"]
            self.assertEqual(layout["coordinate_shift_bp"], shift)
            self.assertEqual(layout["seed_layer_ranges"],
                             [[80, 207], [240, 367]])
            self.assertEqual(layout["seed_staple_physical_range"],
                             [80, 367])
            self.assertEqual(layout["seed_staple_required_coverage_range"],
                             [96, 351])
            report = validate_structure(str(final_path),
                                        require_staples=True)
            self.assertTrue(report["valid"], report["errors"])
            self.assertEqual(report["seed_scaffold_lengths"], [7300, 7336])
            self.assertEqual(report["capture_bridge_component_count"], 48)


if __name__ == "__main__":
    unittest.main()
