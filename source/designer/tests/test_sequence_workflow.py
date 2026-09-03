"""Regression tests for the staged scaffold/SST/final-export workflow."""

import json
from pathlib import Path
import tempfile
import unittest
from zipfile import ZipFile

from cadnano2.model.io.sequencexlsx import _write_workbook

from moire_design_core import (
    finalize_structure,
    generate_scaffold_review,
    write_shifted_sst,
)
from moire_design_core.sequence_workflow_worker import (
    BASES,
    CADNANO_SCAFFOLDS,
    _add_scaffold_metadata_columns,
    _apply_seed_template_capture_colors,
    _build_seed_sst_assembled_payload,
    _design_targets,
    _load,
    _mapped_complete_sst_input_records,
    _mapped_input_records,
    _scaffold_base_map,
    _scaffold_export_summary,
    _sequence_sheets,
    _sst_physical_layer,
    analyze,
    assign_standard_scaffold,
    build_sequenced,
    export_template,
    import_template,
    scaffold_catalog,
)
from moire_design_core.sequence_workflow import build_sequenced_design
from moire_design_core.sequence_export_worker import (
    _add_capture_staple_columns,
    _capture_map_start_end,
    _capture_manifest_rows,
)


class SequenceWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        sst = write_shifted_sst(
            str(cls.root / "sst.json"), 128, 32, 128, 128, 128)
        cls.sst = sst
        generate_scaffold_review(str(cls.root / "scaffold.json"), str(sst))
        cls.scaffold = cls.root / "scaffold.json"
        finalize_structure(
            str(cls.root / "scaffold.json"), str(cls.root / "final.json"))
        cls.report = analyze(str(cls.root / "final.json"))
        cls.assignments = []
        for group in cls.report["targets"].values():
            for target in group:
                sequence = ("ACGT" * ((target["length"] + 3) // 4))[
                    :target["length"]]
                cls.assignments.append({
                    "target_id": target["id"], "sequence": sequence})
        cls.sequenced = cls.root / "with_sequences.json"
        cls.build_report = build_sequenced(
            str(cls.root / "final.json"), str(cls.sequenced),
            cls.assignments)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_detects_scaffolds_and_both_sst_layers(self):
        self.assertGreaterEqual(
            self.report["summary"]["seed_scaffold"]["count"], 1)
        self.assertEqual(
            self.report["summary"]["sst_input_layer_1"]["count"], 64)
        self.assertEqual(
            self.report["summary"]["sst_input_layer_2"]["count"], 64)
        self.assertEqual(
            self.report["summary"]["sst_input_layer_1"]["lengths"],
            {"32": 64})

    def test_large_assignments_round_trip_through_response_file(self):
        serialized = json.dumps(self.assignments, ensure_ascii=False)
        self.assertGreater(len(serialized), 8000)
        output = self.root / "with_sequences_via_response_file.json"
        report = build_sequenced_design(
            str(self.root / "final.json"), str(output), self.assignments)
        self.assertEqual(Path(report["path"]).resolve(), output.resolve())
        self.assertTrue(output.is_file())
        self.assertEqual(
            _load(output)["scaffold_sequences"],
            _load(self.sequenced)["scaffold_sequences"])

    def test_physical_layer_uses_real_sst_ranges_not_layer_centres(self):
        """Asymmetric layers must not inherit a Seed/Z2-like split point."""
        layer_ranges = [[24, 423], [432, 495]]
        routing = {
            "enabled": True,
            "layer": 2,
            "auxiliary_internal_helices": list(range(64, 80)),
        }
        # These strands are far to the right of the long Layer-1 centre, but
        # still lie physically inside Layer 1.
        self.assertEqual(_sst_physical_layer(
            [(1, 352, 367), (2, 352, 367)],
            layer_ranges, routing), 1)
        self.assertEqual(_sst_physical_layer(
            [(0, 408, 423), (7, 408, 423)],
            layer_ranges, routing), 1)
        # The first ordinary strand wholly inside the real Layer-2 interval.
        self.assertEqual(_sst_physical_layer(
            [(1, 432, 447), (0, 432, 447)],
            layer_ranges, routing), 2)
        # A short-spacing auxiliary detour is explicitly Layer 2 even when
        # its physical base indices touch the Layer-1 boundary.
        for auxiliary_helix in range(64, 80):
            self.assertEqual(_sst_physical_layer(
                [(0, 432, 439),
                 (auxiliary_helix, 424, 431),
                 (15, 432, 439)],
                layer_ranges, routing), 2)

    def test_asymmetric_generated_design_uses_physical_sst_boundary(self):
        """Regression for the 400/8/64 case previously split at base 367."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sst = write_shifted_sst(
                str(root / "sst.json"), 400, 8, 64, 128, 128,
                lattice_type="square",
                layers_design_sequence_identical=False)
            generate_scaffold_review(
                str(root / "scaffold.json"), str(sst))
            finalize_structure(
                str(root / "scaffold.json"), str(root / "final.json"))
            report = analyze(str(root / "final.json"))
        self.assertEqual(
            report["metadata"]["variable_length_layout"]["layer_ranges"],
            [[24, 423], [432, 495]])
        layer_1_starts = {
            item["start"] for item in report["targets"][
                "sst_input_layer_1"]}
        layer_2_starts = {
            item["start"] for item in report["targets"][
                "sst_input_layer_2"]}
        self.assertIn("1[367]", layer_1_starts)
        self.assertNotIn("1[367]", layer_2_starts)
        self.assertIn("1[447]", layer_2_starts)

    def test_identical_layer_template_maps_layer_one_to_layer_two(self):
        blank = self.root / "blank.xlsx"
        exported = export_template(
            str(self.root / "final.json"), str(blank), True)
        self.assertEqual(exported["row_count"], 64)
        rows = []
        for target in self.report["targets"]["sst_input_layer_1"]:
            rows.append((
                target["start"], target["end"], "AGCT" * 8,
                target["length"], target["color"]))
        filled = self.root / "filled.xlsx"
        _write_workbook(str(filled), (("input", rows),),
                        use_row_colors=True)
        imported = import_template(
            str(self.root / "final.json"), str(filled), True)
        self.assertEqual(len(imported["assignments"]), 128)
        self.assertEqual(sum(
            bool(item.get("copied_from_layer_1"))
            for item in imported["assignments"]), 64)

    def test_final_sequences_have_no_unknown_bases_and_capture_is_sorted(self):
        self.assertEqual(self.build_report["unresolved_output_bases"], 0)
        inputs, outputs, manifest = _sequence_sheets(_load(self.sequenced))
        for rows in list(inputs.values()) + list(outputs.values()):
            self.assertTrue(all("?" not in row[2] for row in rows))
        sort_keys = [
            (int(item["column_base"]), int(item["capture_base"]),
             int(item["capture_seed_helix"]))
            for item in manifest]
        self.assertEqual(sort_keys, sorted(sort_keys))
        by_column = {}
        for item in manifest:
            by_column.setdefault(int(item["column_base"]), set()).add(
                str(item["color"]).lower())
        self.assertEqual(len(by_column), 18)
        self.assertTrue(all(len(colors) == 1
                            for colors in by_column.values()))
        self.assertEqual(
            len({next(iter(colors)) for colors in by_column.values()}),
            len(by_column))

    def test_sst_output_uses_complete_routing_and_one_black_color(self):
        """Capture gaps must not split purchasing strands into 16-mers."""
        unused_inputs, outputs, unused_manifest = _sequence_sheets(
            _load(self.sequenced))
        for layer in (1, 2):
            rows = outputs["sst_output_layer_%d" % layer]
            self.assertTrue(rows)
            self.assertNotIn(16, {int(row[3]) for row in rows})
            self.assertEqual(
                {str(row[4]).lower() for row in rows}, {"#000000"})
        # The canonical 128/32/128 Square resource is made entirely from
        # complete 32-nt U-shaped SST output strands.
        self.assertEqual(
            {int(row[3]) for row in outputs["sst_output_layer_1"]}, {32})
        self.assertEqual(
            {int(row[3]) for row in outputs["sst_output_layer_2"]}, {32})

    def test_sst_output_prefers_the_saved_standalone_sst_design(self):
        """The accepted SST process file is the authoritative topology."""
        source = _load(self.sequenced)
        standalone = _load(self.sst)
        standalone["scaffold_sequences"] = \
            _mapped_complete_sst_input_records(source, self.sst)
        unused_inputs, outputs, unused_manifest = _sequence_sheets(
            source, complete_sst_payload=standalone)
        for layer in (1, 2):
            rows = outputs["sst_output_layer_%d" % layer]
            self.assertTrue(rows)
            self.assertEqual({int(row[3]) for row in rows}, {32})
            self.assertEqual(
                {str(row[4]).lower() for row in rows}, {"#000000"})

    def test_zero_spacing_complete_sst_maps_by_physical_base(self):
        """Auxiliary h64-79 routing may repartition, but not lose, inputs."""
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            sst = write_shifted_sst(
                str(root / "sst.json"), 64, 0, 64, 128, 128,
                lattice_type="square_kagome")
            generate_scaffold_review(
                str(root / "scaffold.json"), str(sst))
            finalize_structure(
                str(root / "scaffold.json"), str(root / "final.json"))
            report = analyze(str(root / "final.json"))
            assignments = []
            for group in report["targets"].values():
                for target in group:
                    sequence = ("ACGT" * (
                        (target["length"] + 3) // 4))[:target["length"]]
                    assignments.append({
                        "target_id": target["id"], "sequence": sequence})
            sequenced = root / "with_sequences.json"
            build_sequenced(
                str(root / "final.json"), str(sequenced), assignments)
            source = _load(sequenced)
            records = _mapped_complete_sst_input_records(source, sst)
            self.assertEqual(len(records), 68)
            self.assertEqual(
                sorted(len(item["sequence"]) for item in records),
                [32] * 64 + [48] * 4)
            standalone = _load(sst)
            standalone["scaffold_sequences"] = records
            unused_inputs, outputs, unused_manifest = _sequence_sheets(
                source, complete_sst_payload=standalone)
            for layer in (1, 2):
                rows = outputs["sst_output_layer_%d" % layer]
                self.assertTrue(rows)
                self.assertTrue(all("?" not in str(row[2]) for row in rows))

    def test_zero_eight_sixteen_spacing_sequence_matrix(self):
        """Low-spacing auxiliary routing must preserve every sequence."""
        for lattice_type in ("square", "kagome", "square_kagome"):
            for spacing in (0, 8, 16):
                with self.subTest(lattice=lattice_type, spacing=spacing), \
                        tempfile.TemporaryDirectory() as folder:
                    root = Path(folder)
                    sst = write_shifted_sst(
                        str(root / "sst.json"), 64, spacing, 64, 128, 128,
                        lattice_type=lattice_type)
                    generate_scaffold_review(
                        str(root / "scaffold.json"), str(sst))
                    finalize_structure(
                        str(root / "scaffold.json"),
                        str(root / "final.json"))
                    report = analyze(str(root / "final.json"))
                    assignments = []
                    for group in report["targets"].values():
                        for target in group:
                            offset = sum(map(ord, target["id"])) % 4
                            sequence = "".join(
                                "ACGT"[(offset + index) % 4]
                                for index in range(target["length"]))
                            assignments.append({
                                "target_id": target["id"],
                                "sequence": sequence,
                            })
                    sequenced = root / "with_sequences.json"
                    build_report = build_sequenced(
                        str(root / "final.json"), str(sequenced),
                        assignments)
                    self.assertEqual(
                        build_report["unresolved_output_bases"], 0)
                    source = _load(sequenced)
                    standalone = _load(sst)
                    standalone["scaffold_sequences"] = \
                        _mapped_complete_sst_input_records(source, sst)

                    coordinates = {
                        (int(row["row"]), int(row["col"]))
                        for row in standalone["vstrands"]}
                    self.assertEqual(
                        _scaffold_base_map(source, coordinates),
                        _scaffold_base_map(standalone))

                    unused_inputs, outputs, manifest = _sequence_sheets(
                        source, complete_sst_payload=standalone)
                    for rows in outputs.values():
                        for row in rows:
                            self.assertEqual(len(row[2]), int(row[3]))
                            self.assertLessEqual(set(row[2]), BASES)
                    self.assertEqual(len(manifest), 144)
                    for item in manifest:
                        self.assertEqual(
                            len(item["sequence"]), int(item["length"]))
                        self.assertLessEqual(set(item["sequence"]), BASES)

                    routing = source["moire_structure_metadata"].get(
                        "auxiliary_sst_routing", {})
                    self.assertEqual(
                        bool(routing.get("enabled")), spacing in (0, 8))
                    auxiliary = set(map(int, routing.get(
                        "auxiliary_internal_helices", [])))
                    self.assertFalse(any(
                        int(item["start_vh"]) in auxiliary
                        for item in report["targets"][
                            "sst_input_layer_1"]))

    def test_eight_spacing_complete_sst_mapping_ignores_deleted_seed(self):
        """Seed deletions must not enter the SST-only coordinate transfer."""
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            sst = write_shifted_sst(
                str(root / "sst.json"), 88, 8, 88, 128, 128,
                lattice_type="square_kagome",
                mean_indel_per_helix=-2.0)
            generate_scaffold_review(
                str(root / "scaffold.json"), str(sst))
            finalize_structure(
                str(root / "scaffold.json"), str(root / "final.json"))
            report = analyze(str(root / "final.json"))
            assignments = []
            for group in report["targets"].values():
                for target in group:
                    sequence = ("ACGT" * (
                        (target["length"] + 3) // 4))[:target["length"]]
                    assignments.append({
                        "target_id": target["id"], "sequence": sequence})
            sequenced = root / "with_sequences.json"
            build_sequenced(
                str(root / "final.json"), str(sequenced), assignments)
            source = _load(sequenced)
            records = _mapped_complete_sst_input_records(source, sst)
            self.assertTrue(records)
            standalone = _load(sst)
            standalone["scaffold_sequences"] = records
            coordinates = {
                (int(row["row"]), int(row["col"]))
                for row in standalone["vstrands"]}
            self.assertEqual(
                _scaffold_base_map(source, coordinates),
                _scaffold_base_map(standalone))

    def test_zero_spacing_capture_extensions_follow_sst_inputs(self):
        """Capture extensions must remain complementary to routed inputs."""
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            sst = write_shifted_sst(
                str(root / "sst.json"), 64, 0, 64, 128, 128,
                lattice_type="square_kagome")
            generate_scaffold_review(
                str(root / "scaffold.json"), str(sst))
            finalize_structure(
                str(root / "scaffold.json"), str(root / "final.json"))
            report = analyze(str(root / "final.json"))
            manifests = []
            for label, sst_base in (("a", "A"), ("c", "C")):
                assignments = []
                for group_name, group in report["targets"].items():
                    for target in group:
                        base = ("G" if group_name == "seed_scaffold"
                                else sst_base)
                        assignments.append({
                            "target_id": target["id"],
                            "sequence": base * target["length"],
                        })
                sequenced = root / (label + ".json")
                build_sequenced(
                    str(root / "final.json"), str(sequenced), assignments)
                unused_inputs, unused_outputs, manifest = _sequence_sheets(
                    _load(sequenced))
                manifests.append(manifest)

            def extension_sequence(item):
                return "".join(
                    item["sequence"][int(run["start"]):int(run["end"])]
                    for run in item.get("sequence_color_runs", [])
                    if run.get("role") == "capture_extension")

            second_by_target = {
                (int(item["seed_helix"]), int(item["base"]),
                 item["translation"]): item
                for item in manifests[1]}
            compared = 0
            for first in manifests[0]:
                if first.get("connection_role") != \
                        "physical Seed-SST connection":
                    continue
                second = second_by_target[(
                    int(first["seed_helix"]), int(first["base"]),
                    first["translation"])]
                first_extension = extension_sequence(first)
                second_extension = extension_sequence(second)
                if not first_extension:
                    continue
                compared += 1
                self.assertEqual(set(first_extension), {"T"})
                self.assertEqual(set(second_extension), {"G"})
            self.assertGreater(compared, 0)

    def test_capture_staple_rows_include_map_coordinates_and_template_color(
            self):
        inputs, outputs, manifest = _sequence_sheets(_load(self.sequenced))
        rows = outputs["staple_capture"]
        map_rows = _capture_manifest_rows(manifest)
        self.assertEqual(len(rows), len(manifest))
        self.assertEqual(len(rows), len(map_rows))
        for row, item, map_row in zip(rows, manifest, map_rows):
            self.assertEqual(len(row), 7)
            self.assertEqual(row[4], item["color"])
            self.assertEqual(row[5], map_row[0])
            self.assertEqual(row[6], map_row[1])

        # Every immutable Capture column has one row/map/extension colour,
        # and adjacent columns never reuse a cooperative-pair colour.
        by_column = {}
        for item in manifest:
            column_base = int(item["column_base"])
            by_column.setdefault(column_base, set()).add(
                str(item["color"]).lower())
            self.assertEqual(
                str(item["color"]).lower(),
                str(item["capture_color"]).lower())
            if item.get("capture_extension_color"):
                self.assertEqual(
                    str(item["color"]).lower(),
                    str(item["capture_extension_color"]).lower())
        self.assertEqual(len(by_column), 18)
        self.assertTrue(all(len(colors) == 1
                            for colors in by_column.values()))
        self.assertEqual(
            len({next(iter(colors)) for colors in by_column.values()}), 18)

        # Capture Map uses the actual extending endpoint, not the whole
        # core oligo's 5' start (the historical Seed 52[63] defect).
        valid_capture_face_helices = set(range(16, 24)) | set(range(40, 48))
        z2_rows = [item for item in manifest if item.get("phase") == "Z2"]
        self.assertTrue(z2_rows)
        for item in z2_rows:
            self.assertIn(
                int(item["capture_seed_helix"]),
                valid_capture_face_helices)
            self.assertNotEqual(
                (int(item["capture_seed_helix"]),
                 int(item["capture_base"])), (52, 63))
            capture_start, unused_end = _capture_map_start_end(item)
            self.assertIn(
                "Seed %d[%d]" % (
                    int(item["capture_seed_helix"]),
                    int(item["capture_base"])),
                capture_start)

        workbook = self.root / "capture_staple_columns.xlsx"
        _write_workbook(
            str(workbook), (("staple_capture", rows),),
            use_row_colors=True)
        _add_capture_staple_columns(workbook, 1)
        with ZipFile(workbook) as archive:
            sheet_xml = archive.read(
                "xl/worksheets/sheet1.xml").decode("utf-8")
        self.assertIn('dimension ref="A1:G', sheet_xml)
        self.assertIn("Capture Start", sheet_xml)
        self.assertIn("Capture End", sheet_xml)
        self.assertIn(rows[0][4], sheet_xml)
        self.assertIn(rows[0][5], sheet_xml)
        self.assertIn(rows[0][6], sheet_xml)

    def test_capture_end_uses_s_for_both_square_phases(self):
        base_item = {
            "capture_seed_helix": 0,
            "capture_base": 72,
            "seed_helix": 0,
            "base": 72,
            "sst_helix": 48,
            "face": "face1",
            "translation": "origin",
            "export_only": False,
        }
        for phase in ("A", "B"):
            item = dict(base_item, phase=phase)
            unused_start, capture_end = _capture_map_start_end(item)
            self.assertIn(" / S / origin / physical", capture_end)
            self.assertNotIn(" / %s /" % phase, capture_end)
        unused_start, kagome_end = _capture_map_start_end(
            dict(base_item, phase="K"))
        self.assertIn(" / K / origin / physical", kagome_end)

    def test_all_144_nonblack_seed_cores_export_as_staple_capture(self):
        unused_inputs, outputs, manifest = _sequence_sheets(
            _load(self.sequenced))
        self.assertEqual(len(outputs["staple_capture"]), 144)
        self.assertEqual(len(manifest), 144)
        self.assertEqual(sum(
            item.get("connection_role", "").startswith(
                "physical Seed-SST") for item in manifest), 64)
        self.assertEqual(sum(
            item.get("connection_role", "").startswith(
                "export-only translated") for item in manifest), 64)
        self.assertEqual(sum(
            item.get("connection_role", "").startswith(
                "immutable potential Z2") for item in manifest), 16)

    def test_capture_export_uses_only_template_capture_components(self):
        unused_inputs, outputs, manifest = _sequence_sheets(
            _load(self.sequenced))
        self.assertEqual(len(outputs["staple_capture"]), 144)
        self.assertNotIn(
            "#888888", {str(item["color"]).lower() for item in manifest})
        for item in manifest:
            self.assertEqual(item["staple_core_color"].lower(), "#000000")
            if item.get("connection_role", "").startswith((
                    "physical Seed-SST", "export-only translated")):
                core_runs = [
                    run for run in item.get("sequence_color_runs", ())
                    if run.get("role") == "staple_core"]
                self.assertTrue(core_runs)
                self.assertTrue(all(
                    run["color"].lower() == "#000000"
                    for run in core_runs))
                extension_runs = [
                    run for run in item.get("sequence_color_runs", ())
                    if run.get("role") == "capture_extension"]
                self.assertTrue(extension_runs)
                self.assertTrue(all(
                    run["color"].lower() ==
                    item["capture_extension_color"].lower()
                    for run in extension_runs))

    def test_connected_potential_z2_core_is_a_valid_physical_capture(self):
        manifest = [{
            "seed_helix": 45,
            "base": 200,
            "capture_extension_color": "#12ab34",
            "sequence_color_runs": [
                {"role": "staple_core", "color": "#ffffff"},
                {"role": "capture_extension", "color": "#ffffff"},
            ],
        }]
        _apply_seed_template_capture_colors(
            manifest, {(45, 200): "#999999"})
        self.assertEqual(manifest[0]["color"], "#12ab34")
        self.assertEqual(manifest[0]["capture_color"], "#12ab34")
        self.assertEqual(manifest[0]["staple_core_color"], "#000000")
        self.assertEqual(
            manifest[0]["sequence_color_runs"][0]["color"], "#000000")
        self.assertEqual(
            manifest[0]["sequence_color_runs"][1]["color"], "#12ab34")

        with self.assertRaisesRegex(ValueError, "does not map"):
            _apply_seed_template_capture_colors(
                [{"seed_helix": 45, "base": 199}],
                {(45, 199): "#000000"})

    def test_short_spacing_connected_z2_core_exports_once(self):
        """Regression for the Square-Kagome 112/16/112 final export."""
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            sst = write_shifted_sst(
                str(root / "sst.json"), 112, 16, 112, 128, 128,
                lattice_type="square_kagome")
            generate_scaffold_review(
                str(root / "scaffold.json"), str(sst))
            finalize_structure(
                str(root / "scaffold.json"), str(root / "final.json"))
            report = analyze(str(root / "final.json"))
            assignments = []
            for group in report["targets"].values():
                for target in group:
                    sequence = ("ACGT" * (
                        (target["length"] + 3) // 4))[:target["length"]]
                    assignments.append({
                        "target_id": target["id"], "sequence": sequence})
            sequenced = root / "with_sequences.json"
            build_sequenced(
                str(root / "final.json"), str(sequenced), assignments)

            unused_inputs, outputs, manifest = _sequence_sheets(
                _load(sequenced))
            self.assertEqual(len(outputs["staple_capture"]), 144)
            self.assertEqual(len(manifest), 144)
            connected = [
                item for item in manifest
                if int(item.get("seed_helix", -1)) == 45 and
                int(item.get("base", -1)) == 200]
            self.assertEqual(len(connected), 1)
            self.assertEqual(
                connected[0]["connection_role"],
                "physical Seed-SST connection")
            self.assertEqual(connected[0]["color"], "#8b4513")
            self.assertEqual(
                connected[0]["capture_extension_color"], "#8b4513")
            self.assertEqual(sum(
                item.get("connection_role", "").startswith(
                    "immutable potential Z2") for item in manifest), 12)

    def test_sequences_remap_to_each_stage_helix_numbering(self):
        source = _load(self.sequenced)
        sst_records = _mapped_input_records(
            source, self.sst,
            ("sst_input_layer_1", "sst_input_layer_2"))
        scaffold_records = _mapped_input_records(
            source, self.scaffold,
            ("seed_scaffold", "sst_input_layer_1", "sst_input_layer_2"))
        self.assertEqual(len(sst_records), 128)
        self.assertEqual(len(scaffold_records), len(self.assignments))
        self.assertTrue(all(set(item["sequence"]) <= set("ACGT")
                            for item in sst_records + scaffold_records))

    def test_assembled_structure_covers_both_eight_helix_faces(self):
        payload = _build_seed_sst_assembled_payload(_load(self.sequenced))
        metadata = payload["moire_structure_metadata"]
        self.assertEqual(metadata["assembled_sst_units_per_face"], 2)
        self.assertEqual(len(metadata["sst_helix_numbers"]), 32)
        self.assertEqual(len(payload["vstrands"]), 80)
        unused_document, part, unused_targets = _design_targets(payload)
        self.assertGreater(len(part.oligos()), 0)

    def test_single_scaffold_catalog_filters_by_routing_length(self):
        names = {
            item["name"] for item in
            scaffold_catalog(7400, multiple=False)["scaffolds"]}
        self.assertNotIn("M13mp18", names)
        self.assertEqual(names, {"CS3L", "CS4", "P7560"})
        self.assertTrue(all(
            item["length"] >= 7400 for item in
            scaffold_catalog(7400, multiple=False)["scaffolds"]))

    def test_multiple_scaffold_catalog_limits_choices_and_prevents_reuse(self):
        names = {
            item["name"] for item in
            scaffold_catalog(7400, multiple=True)["scaffolds"]}
        self.assertEqual(names, {"CS3L", "CS4", "P7560"})
        remaining = {
            item["name"] for item in scaffold_catalog(
                7400, multiple=True, used_names=("CS3L",))["scaffolds"]}
        self.assertEqual(remaining, {"CS4", "P7560"})

    def test_standard_scaffold_assignment_is_exact_and_rejects_duplicate(self):
        target = self.report["targets"]["seed_scaffold"][0]
        assignment = assign_standard_scaffold(
            target, "CS3L", multiple=True)
        self.assertEqual(assignment["scaffold_name"], "CS3L")
        self.assertEqual(assignment["source"],
                         "caDNAno built-in: CS3L_7559")
        self.assertEqual(
            assignment["sequence"],
            CADNANO_SCAFFOLDS["CS3L_7559"][:target["length"]])
        self.assertEqual(len(assignment["sequence"]), target["length"])
        self.assertEqual(assignment["length"], target["length"])
        with self.assertRaises(ValueError):
            assign_standard_scaffold(
                target, "CS3L", multiple=True, used_names=("CS3L",))

    def test_corrected_scaffold_names_keep_exact_cadnano_sequences(self):
        target = dict(self.report["targets"]["seed_scaffold"][0])
        target["length"] = 7400
        cases = (
            ("CS3L", "CS3L_7559", 7559),
            ("CS4", "CS4_7557", 7557),
        )
        for name, key, source_length in cases:
            assignment = assign_standard_scaffold(
                target, name, multiple=True)
            self.assertEqual(assignment["scaffold_name"], name)
            self.assertEqual(assignment["scaffold_source_length"],
                             source_length)
            self.assertEqual(assignment["sequence"],
                             CADNANO_SCAFFOLDS[key][:7400])

        # Old project labels remain loadable but are normalized immediately.
        legacy = assign_standard_scaffold(target, "CS3", multiple=True)
        self.assertEqual(legacy["scaffold_name"], "CS3L")
        hidden = {item["name"] for item in scaffold_catalog(
            7400, multiple=True, used_names=("CS4-L",))["scaffolds"]}
        self.assertNotIn("CS4", hidden)

    def test_scaffold_name_total_length_and_used_length_are_exported(self):
        assignments = [dict(item) for item in self.assignments]
        scaffold_targets = {
            item["id"]: item
            for item in self.report["targets"]["seed_scaffold"]}
        expected = {}
        scaffold_index = 0
        for assignment in assignments:
            target = scaffold_targets.get(assignment["target_id"])
            if target is None:
                continue
            scaffold_index += 1
            name = "Test Scaffold %d" % scaffold_index
            total_length = int(target["length"]) + scaffold_index * 11
            assignment.update({
                "category": "seed_scaffold",
                "scaffold_name": name,
                "scaffold_source_length": total_length,
            })
            expected[(target["start_vh"], target["start_idx"])] = (
                name, total_length, int(target["length"]))

        output = self.root / "with_scaffold_metadata.json"
        build_sequenced(
            str(self.root / "final.json"), str(output), assignments)
        payload = _load(output)
        saved = {
            (item["start_vh"], item["start_idx"]): item
            for item in payload["scaffold_sequences"]
            if item.get("scaffold_name")}
        self.assertEqual(set(saved), set(expected))
        for key, values in expected.items():
            name, total_length, used_length = values
            self.assertEqual(saved[key]["scaffold_name"], name)
            self.assertEqual(
                saved[key]["scaffold_source_length"], total_length)
            self.assertEqual(
                saved[key]["scaffold_used_length"], used_length)

        summary = _scaffold_export_summary(payload)
        self.assertEqual(len(summary), len(expected))
        self.assertEqual(
            {item["scaffold_name"] for item in summary},
            {values[0] for values in expected.values()})
        self.assertTrue(all(
            item["total_scaffold_length_nt"] >=
            item["length_used_in_structure_nt"]
            for item in summary))

        inputs, unused_outputs, unused_manifest = _sequence_sheets(payload)
        workbook = self.root / "scaffold_metadata.xlsx"
        _write_workbook(
            str(workbook), (("scaffold", inputs["scaffold"]),),
            use_row_colors=True)
        _add_scaffold_metadata_columns(workbook, summary)
        with ZipFile(workbook) as archive:
            sheet_xml = archive.read(
                "xl/worksheets/sheet1.xml").decode("utf-8")
        self.assertIn('dimension ref="A1:H', sheet_xml)
        self.assertIn("Scaffold Name", sheet_xml)
        self.assertIn("Total Scaffold Length (nt)", sheet_xml)
        self.assertIn("Length Used in Structure (nt)", sheet_xml)
        for name, total_length, used_length in expected.values():
            self.assertIn(name, sheet_xml)
            self.assertIn("<v>%d</v>" % total_length, sheet_xml)
            self.assertIn("<v>%d</v>" % used_length, sheet_xml)


if __name__ == "__main__":
    unittest.main()
