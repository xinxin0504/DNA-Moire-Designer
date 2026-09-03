"""Windows packaging regressions for the English-only distribution."""

from __future__ import annotations

import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from moire_design_core import sequence_workflow, sequence_workflow_worker
from moire_design_core.models import MoireProject, SquareBilayerSettings
from moire_designer import i18n


ROOT = Path(__file__).parents[3]
DESIGNER_ROOT = ROOT / "source" / "designer"
COMPANION_ROOT = ROOT / "source" / "cadnano_companion"
AUDIT_PATH = ROOT / "build_windows" / "tools" / "audit_english_runtime.py"


def _load_audit_module():
    spec = importlib.util.spec_from_file_location(
        "windows_english_audit", AUDIT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class WindowsReleaseRegressionTests(unittest.TestCase):
    def tearDown(self):
        i18n.set_language("en")

    def test_windows_distribution_exposes_only_english(self):
        self.assertEqual(i18n.LANGUAGES, (("en", "English", "English"),))
        self.assertEqual(i18n.set_language("zh_CN"), "en")
        self.assertEqual(i18n.current_language(), "en")
        chooser = (DESIGNER_ROOT / "moire_designer" /
                   "project_session.py").read_text(encoding="utf-8")
        self.assertNotIn("language_selector", chooser)
        self.assertNotIn("QComboBox", chooser)

        project = MoireProject.from_dict(MoireProject(
            settings=SquareBilayerSettings(interface_language="zh_CN"),
            prediction={}, validation=[], capture_plan={}, seed_plan={}
        ).to_dict())
        self.assertEqual(project.settings.interface_language, "en")

    def test_public_runtime_has_no_unresolved_chinese_presentation_text(self):
        audit = _load_audit_module()
        unresolved = [
            item for item in audit.audit(DESIGNER_ROOT)
            if not item.get("resolved", False)
        ]
        unresolved.extend(audit.audit_companion(COMPANION_ROOT))
        self.assertEqual(unresolved, [])
        catalog_cjk = [
            (source, target)
            for source, target in i18n._catalogs.get("en", {}).items()
            if audit.CJK_RE.search(str(target))
        ]
        self.assertEqual(catalog_cjk, [])

    def test_large_sequence_payload_uses_short_response_file(self):
        assignments = [
            {"target_id": "target-%03d" % index, "sequence": "ACGT" * 128}
            for index in range(128)
        ]
        observed = {}

        def fake_worker_command(name, *arguments):
            observed["name"] = name
            observed["arguments"] = list(arguments)
            return ["worker", *arguments]

        def fake_run(command, **unused):
            observed["command"] = command
            response_path = Path(command[-1])
            self.assertTrue(response_path.is_file())
            values = json.loads(response_path.read_text(encoding="utf-8"))
            self.assertEqual(len(json.loads(values[2])), 128)
            return SimpleNamespace(
                returncode=0, stdout='{"path": "sequenced_design.json"}',
                stderr="")

        with patch.object(sequence_workflow, "worker_command",
                          side_effect=fake_worker_command), \
                patch.object(sequence_workflow.subprocess, "run",
                             side_effect=fake_run):
            report = sequence_workflow.build_sequenced_design(
                "design.json", "sequenced_design.json", assignments)

        self.assertEqual(report["path"], "sequenced_design.json")
        self.assertEqual(observed["name"], "sequence-workflow")
        self.assertEqual(observed["arguments"][0], "build-sequenced")
        self.assertEqual(observed["arguments"][1], "@arguments-file")

    def test_sequence_worker_expands_response_file_before_dispatch(self):
        values = [
            "design.json", "sequenced.json",
            json.dumps([{"target_id": "target-001", "sequence": "ACGT"}]),
        ]
        observed = {}
        with tempfile.TemporaryDirectory() as directory:
            response_file = Path(directory) / "args.json"
            response_file.write_text(json.dumps(values), encoding="utf-8")

            def fake_build(design, output, assignments):
                observed.update(
                    design=design, output=output, assignments=assignments)
                return {"path": output}

            argv = [
                "sequence_workflow_worker.py", "build-sequenced",
                "@arguments-file", str(response_file),
            ]
            with patch.object(sys, "argv", argv), \
                    patch.object(sequence_workflow_worker, "build_sequenced",
                                 side_effect=fake_build), \
                    patch("builtins.print"):
                sequence_workflow_worker.main()

        self.assertEqual(observed["design"], "design.json")
        self.assertEqual(observed["output"], "sequenced.json")
        self.assertEqual(
            observed["assignments"],
            [{"target_id": "target-001", "sequence": "ACGT"}])

    def test_final_export_report_is_safe_on_windows_cp1252_stdout(self):
        """Unicode paths must survive the frozen worker transport unchanged."""
        result = {
            "path": r"C:\\export\\PDB∕oxView files",
            "message": "Moiré export complete",
        }
        argv = [
            "sequence_workflow_worker.py", "final-export",
            "project.moire.json", "sequenced.json", "output",
            json.dumps({}),
        ]
        raw = io.BytesIO()
        stdout = io.TextIOWrapper(raw, encoding="cp1252", errors="strict")
        try:
            with patch.object(sys, "argv", argv), \
                    patch.object(sequence_workflow_worker, "final_export",
                                 return_value=result), \
                    patch.object(sys, "stdout", stdout):
                sequence_workflow_worker.main()
                stdout.flush()
            transported = raw.getvalue().decode("cp1252")
        finally:
            stdout.detach()

        self.assertEqual(json.loads(transported), result)
        self.assertIn(r"\u2215", transported)

    def test_every_json_worker_uses_ascii_safe_stdout(self):
        worker_names = (
            "image_analysis_worker.py",
            "scale_detection_worker.py",
            "sequence_export_worker.py",
            "sequence_workflow_worker.py",
            "structure_worker.py",
            "vector_export_worker.py",
        )
        core = DESIGNER_ROOT / "moire_design_core"
        for name in worker_names:
            with self.subTest(worker=name):
                text = (core / name).read_text(encoding="utf-8")
                self.assertNotIn("print(json.dumps(result, ensure_ascii=False))",
                                 text)
                self.assertNotIn("print(json.dumps(report, ensure_ascii=False))",
                                 text)
                self.assertNotIn(
                    'print(json.dumps({"svg": str(target)}, ensure_ascii=False))',
                    text)

    def test_companion_spec_bundles_dynamic_qtsvg_modules(self):
        spec_text = (ROOT / "build_windows" / "cadnano2_companion.spec").read_text(
            encoding="utf-8")
        self.assertIn('"PyQt6.QtSvg"', spec_text)
        self.assertIn('"PyQt6.QtSvgWidgets"', spec_text)

    def test_frozen_specs_bundle_language_catalog_and_cadnano_metadata(self):
        designer_spec = (ROOT / "build_windows" /
                         "DNA_Moire_Designer.spec").read_text(encoding="utf-8")
        companion_spec = (ROOT / "build_windows" /
                          "cadnano2_companion.spec").read_text(encoding="utf-8")
        build_script = (ROOT / "build_windows" /
                        "build.ps1").read_text(encoding="utf-8")
        self.assertIn('"translations.json"', designer_spec)
        self.assertIn('copy_metadata("cadnano2")', companion_spec)
        self.assertIn("cadnano2-2.4.13-py3-none-any.whl", build_script)

    def test_frozen_self_tests_import_real_modules_and_block_packaging(self):
        launcher = (DESIGNER_ROOT /
                    "run_moire_designer.py").read_text(encoding="utf-8")
        build_script = (ROOT / "build_windows" /
                        "build.ps1").read_text(encoding="utf-8")
        self.assertIn(
            "from moire_design_core import models, structure, template",
            launcher)
        self.assertNotIn(
            "from moire_design_core import model, structure, template",
            launcher)
        self.assertGreaterEqual(build_script.count("Start-Process"), 2)
        self.assertGreaterEqual(build_script.count("-Wait -PassThru"), 2)
        self.assertIn("$DesignerSelfTest.ExitCode", build_script)
        self.assertIn("$CadnanoSelfTest.ExitCode", build_script)


if __name__ == "__main__":
    unittest.main()
