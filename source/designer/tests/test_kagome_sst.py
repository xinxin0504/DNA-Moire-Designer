"""Regression tests for Kagome SST-only and SST-side capture anchors."""

import json
import tempfile
import unittest
from pathlib import Path

from moire_design_core import validate_sst, write_shifted_sst
from moire_design_core.kagome_sst import (
    ACTIVE_HELICES,
    HOLE_HELICES,
    _line_intervals,
    build_kagome_sst_payload,
    kagome_capture_anchor_candidates,
    kagome_actual_capture_anchor_candidates,
    prepare_kagome_capture_sites,
    required_global_shift,
    validate_kagome_sst_payload,
)
from moire_design_core.structure import build_shifted_sst_payload


REFERENCE_ROOT = Path.home() / "Desktop"
BUNDLED_REFERENCE_ROOT = Path(__file__).resolve().parents[1] / \
    "moire_design_core" / "resources"


def _normalized_record(record, base_shift):
    output = []
    for slot in (0, 2):
        partner, partner_base = map(int, record[slot:slot + 2])
        if partner < 0:
            output.extend((-1, -1))
        else:
            output.extend((partner - 32, partner_base - base_shift))
    return output


def _normalized_sst_first_record(record, base_shift):
    output = []
    for slot in (0, 2):
        partner, partner_base = map(int, record[slot:slot + 2])
        if partner < 0:
            output.extend((-1, -1))
        else:
            output.extend((partner - 48, partner_base - base_shift))
    return output


def _first_occupied_across(rows, field):
    return min(index for row in rows.values()
               for index, record in enumerate(row[field])
               if record != [-1, -1, -1, -1])


class KagomeSSTTests(unittest.TestCase):
    def test_linear_sst_preserves_template_nick_and_capture_phase(self):
        self.assertEqual(_line_intervals(48, 175, 16), [
            (48, 79), (80, 111), (112, 143), (144, 175)])
        self.assertEqual(_line_intervals(56, 183, 8), [
            (56, 103), (104, 135), (136, 183)])
        self.assertEqual(_line_intervals(48, 191, 16), [
            (48, 79), (80, 111), (112, 143), (144, 191)])
        self.assertEqual(_line_intervals(56, 199, 8), [
            (56, 103), (104, 135), (136, 167), (168, 199)])

    def test_five_validated_lengths_match_reference_routing(self):
        for length in (96, 104, 112, 120, 128):
            reference_path = REFERENCE_ROOT / (
                "kagome_resource_%dbp.json" % length)
            if not reference_path.is_file():
                reference_path = REFERENCE_ROOT / (
                    "kagome_resource_%d.json" % length)
            if not reference_path.is_file():
                reference_path = BUNDLED_REFERENCE_ROOT / (
                    "kagome_sst_%dbp_fixture.json" % length)
            self.assertTrue(reference_path.is_file(), reference_path)
            reference = json.loads(reference_path.read_text())
            source_rows = {int(row["num"]): row
                           for row in reference["vstrands"]}
            source_is_bundled_sst_first = 16 not in source_rows
            source_offset = 0 if not source_is_bundled_sst_first else -16
            generated = build_kagome_sst_payload(
                "generated.json", length, 32, length)
            target_rows = {int(row["num"]): row
                           for row in generated["vstrands"]}
            shift = required_global_shift(length, 32, length)
            bundled_shifts = {}
            if source_is_bundled_sst_first:
                bundled_shifts = {
                    field: (_first_occupied_across(target_rows, field) -
                            _first_occupied_across(source_rows, field))
                    for field in ("scaf", "stap")}
            for source_helix in range(16, 32):
                source = source_rows[source_helix + source_offset]
                target = target_rows[source_helix + 32]
                for field in ("scaf", "stap"):
                    if not source_is_bundled_sst_first:
                        actual = [
                            _normalized_record(
                                target[field][base + shift], shift)
                            for base in range(len(source[field]))]
                    if source_is_bundled_sst_first:
                        field_shift = bundled_shifts[field]
                        generated_local = [
                            _normalized_sst_first_record(
                                (target[field][base + field_shift]
                                 if 0 <= base + field_shift <
                                 len(target[field]) else
                                 [-1, -1, -1, -1]),
                                field_shift)
                            for base in range(len(source[field]))]
                        self.assertEqual(source[field], generated_local)
                    else:
                        self.assertEqual(
                            source[field], actual,
                            "%dbp h%d %s" %
                            (length, source_helix, field))
                actual_colors = sorted(
                    [int(base) - shift, int(color)]
                    for base, color in target.get("stap_colors", []))
                expected_colors = source.get("stap_colors", [])
                if source_is_bundled_sst_first:
                    color_shift = bundled_shifts["stap"]
                    actual_colors = sorted(
                        [int(base) - color_shift, int(color)]
                        for base, color in
                        target.get("stap_colors", []))
                self.assertEqual(expected_colors, actual_colors)
            self.assertTrue(validate_kagome_sst_payload(generated)["valid"])

    def test_capture_anchors_use_three_original_kagome_sst_topologies(self):
        complete = build_kagome_sst_payload(
            "complete.json", 128, 32, 128)
        anchors = kagome_capture_anchor_candidates(complete)
        counts = {}
        for anchor in anchors:
            counts[anchor["origin_type"]] = \
                counts.get(anchor["origin_type"], 0) + 1
            self.assertIsNone(anchor["seed_helix"])
            self.assertTrue(anchor["seed_mapping_pending"])
        self.assertEqual(counts, {
            "sst_crossover": 32,
            "preexisting_nick": 12,
            "linear_continuous": 4,
        })
        prepared = prepare_kagome_capture_sites(complete)
        metadata = prepared["moire_structure_metadata"]
        self.assertEqual(metadata["kagome_capture_extension_counts"], {
            "16": 36, "32": 12})
        self.assertTrue(validate_kagome_sst_payload(prepared)["valid"])
        # The complete design passed to the derivation remains untouched.
        self.assertFalse(complete["moire_structure_metadata"]
                         ["capture_gaps_reserved"])

    def test_capture_grid_comes_from_duplex_and_absolute_seed_phase(self):
        payload = build_kagome_sst_payload(
            "clipped.json", 136, 56, 136)
        anchors = kagome_capture_anchor_candidates(payload)
        by_layer = {}
        for item in anchors:
            by_layer.setdefault(item["layer"], {}).setdefault(
                item["position"], []).append(item)
        duplex = payload["moire_structure_metadata"]["sst_duplex_ranges"]
        for layer, bounds in enumerate(duplex, 1):
            expected = list(range(
                int(bounds[0]) + ((8-int(bounds[0])) % 16),
                int(bounds[1]) + 1, 16))
            self.assertEqual(sorted(by_layer[layer]), expected)
            for position in expected:
                items = by_layer[layer][position]
                if position % 32 == 24:
                    self.assertEqual(len(items), 4)
                    self.assertTrue(all(item["capture_family"] ==
                                        "u_shaped_16nt" for item in items))
                else:
                    self.assertEqual(len(items), 2)
                    self.assertTrue(all(item["capture_family"].startswith(
                                        "linear_") for item in items))

    def test_same_sst_length_can_start_with_either_family_by_seed_overlap(self):
        payload = build_kagome_sst_payload("same-sst.json", 80, 32, 80)
        duplex = payload["moire_structure_metadata"]["sst_duplex_ranges"]
        theoretical = kagome_capture_anchor_candidates(payload)
        positions = sorted({int(item["position"]) for item in theoretical
                            if int(item["layer"]) == 1})
        linear_position = next(value for value in positions
                               if value % 32 == 8)
        u_position = next(value for value in positions
                          if value % 32 == 24)
        seed_ranges = [duplex[0], duplex[1]]
        linear_first = kagome_actual_capture_anchor_candidates(
            payload, [[linear_position], []], seed_ranges)
        u_first = kagome_actual_capture_anchor_candidates(
            payload, [[u_position], []], seed_ranges)
        self.assertEqual({item["capture_family"] for item in linear_first},
                         {"linear_32nt_or_right_edge_16nt"})
        self.assertEqual({item["capture_family"] for item in u_first},
                         {"u_shaped_16nt"})
        self.assertEqual(len(linear_first), 2)
        self.assertEqual(len(u_first), 4)

    def test_holes_remain_empty_and_spacing_is_independent(self):
        payload = build_kagome_sst_payload(
            "independent.json", 120, 40, 104)
        rows = {int(row["num"]): row for row in payload["vstrands"]}
        for helix in HOLE_HELICES:
            for field in ("scaf", "stap"):
                self.assertTrue(all(record == [-1, -1, -1, -1]
                                    for record in rows[helix][field]))
        for field in ("scaf", "stap"):
            occupied = {
                helix for helix, row in rows.items()
                if any(record != [-1, -1, -1, -1]
                       for record in row[field])}
            self.assertEqual(occupied, set(ACTIVE_HELICES))
        duplex = payload["moire_structure_metadata"]["sst_duplex_ranges"]
        self.assertEqual(duplex[1][0] - duplex[0][1] - 1, 40)
        self.assertTrue(validate_kagome_sst_payload(payload)["valid"])

    def test_public_writer_and_validator_dispatch_to_kagome(self):
        with tempfile.TemporaryDirectory() as folder:
            path = write_shifted_sst(
                str(Path(folder) / "kagome.json"), 128, 32, 128,
                lattice_type="kagome")
            payload = json.loads(path.read_text())
            self.assertEqual([int(row["num"]) for row in payload["vstrands"]],
                             list(range(16)))
            report = validate_sst(str(path))
        self.assertTrue(report["valid"], report["errors"])
        self.assertEqual(report["lattice_type"], "kagome")

    def test_square_default_and_explicit_square_are_identical(self):
        default = build_shifted_sst_payload("same.json")
        explicit = build_shifted_sst_payload(
            "same.json", lattice_type="square")
        self.assertEqual(default, explicit)


if __name__ == "__main__":
    unittest.main()
