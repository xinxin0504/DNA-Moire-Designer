#!/usr/bin/env python3
"""Static checks for the deliberately reduced Windows product surface."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MAIN = (ROOT / "source" / "designer" / "moire_designer" /
        "mainwindow.py").read_text(encoding="utf-8")
BULK = (ROOT / "source" / "designer" / "moire_designer" /
        "analysis_bulk.py").read_text(encoding="utf-8")
PROJECT_SESSION = (ROOT / "source" / "designer" / "moire_designer" /
                   "project_session.py").read_text(encoding="utf-8")
IMAGE_WORKER = (ROOT / "source" / "designer" / "moire_design_core" /
                "image_analysis_worker.py").read_text(encoding="utf-8")
CORE = ROOT / "source" / "designer" / "moire_design_core"

assert "analysis_particle_action" not in MAIN
assert "_build_particle_analysis_panel" not in BULK
assert "_build_gel_analysis_panel" not in BULK
assert "analysis_module_stack.addWidget(self._build_crystal_analysis_panel())" in BULK
assert "elif step == 3:" in MAIN
assert 'menu_bar.addMenu("语言")' not in MAIN
assert "language_selector" not in PROJECT_SESSION
assert "self.capture_results_splitter = QSplitter(" in MAIN
assert "self.capture_results_splitter.addWidget(\n            self.structure_preview_status)" in MAIN
assert "self.structure_preview_status.setVisible(bool(text))" in MAIN
assert 'self.capture_preview.set_path_report("")' in MAIN
assert "if float(value) > 0.0" in MAIN
assert '"structure_accepted_at"' in MAIN
assert '"sequence_scaffold_accepted_at"' in MAIN
for argument in (
        "--analysis-kind", "--output-dir", "--domain-edits",
        "--theoretical-a-nm", "--theoretical-symmetry",
        "--pixel-size-nm", "--scale-value-nm"):
    assert argument in IMAGE_WORKER, "missing image-worker option: " + argument

# Frozen GUI workers communicate with the main process through stdout.  Keep
# every JSON response ASCII-safe so Windows legacy code pages (for example
# cp1252) cannot crash on user-facing Unicode paths such as PDB∕oxView files.
for worker_name in (
        "image_analysis_worker.py",
        "scale_detection_worker.py",
        "sequence_export_worker.py",
        "sequence_workflow_worker.py",
        "structure_worker.py",
        "vector_export_worker.py"):
    worker_text = (CORE / worker_name).read_text(encoding="utf-8")
    assert "print(json.dumps(" in worker_text, (
        "missing JSON stdout response: " + worker_name)
    assert "ensure_ascii=True" in worker_text, (
        "Windows-unsafe JSON stdout response: " + worker_name)
    assert "print(json.dumps(result, ensure_ascii=False))" not in worker_text
    assert "print(json.dumps(report, ensure_ascii=False))" not in worker_text
print("release_surface=English-only Designer+Moiré Analysis; Particle/Gel UI excluded")
