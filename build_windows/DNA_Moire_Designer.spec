# -*- mode: python ; coding: utf-8 -*-
from pathlib import Path
import sys

SPEC_DIR = Path(SPECPATH)
ROOT = SPEC_DIR.parent
SOURCE = ROOT / "source" / "designer"
PRIVATE_ENGINE = SOURCE / "designer_vendor"
sys.path.insert(0, str(PRIVATE_ENGINE))
sys.path.insert(0, str(SOURCE))

hiddenimports = [
    "moire_design_core.image_analysis_worker",
    "moire_design_core.scale_detection_worker",
    "moire_design_core.structure_worker",
    "moire_design_core.sequence_workflow_worker",
    "moire_design_core.sequence_export_worker",
    "moire_design_core.vector_export_worker",
    "cadnano2.model.io.sequencexlsx",
    "cadnano2.model.io.oxdnaexport",
    "cadnano2.model.io.primer3analysis",
    "PyQt6.QtSvg",
    "PyQt6.QtPrintSupport",
    "primer3.thermoanalysis",
]

datas = [
    (str(SOURCE / "assets"), "assets"),
    (str(SOURCE / "moire_designer" / "assets"), "moire_designer/assets"),
    # i18n.py loads this file through Path(__file__).with_name().  Python
    # modules are discovered automatically by PyInstaller, but adjacent JSON
    # resources are not.  Without this entry the frozen application silently
    # falls back to the very small built-in bootstrap catalog and most of the
    # historical Chinese UI remains visible.
    (str(SOURCE / "moire_designer" / "translations.json"),
     "moire_designer"),
    (str(SOURCE / "moire_design_core" / "resources"),
     "moire_design_core/resources"),
]
# Preserve the validated private engine as an isolated runtime resource.
# Analysis.datas accepts (source, destination) pairs; passing Tree's internal
# three-field TOC entries here fails on PyInstaller 6.x on Windows.
datas.append((str(PRIVATE_ENGINE / "cadnano2"), "designer_vendor/cadnano2"))

a = Analysis(
    [str(SOURCE / "run_moire_designer.py")],
    pathex=[str(SOURCE), str(PRIVATE_ENGINE)],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        "moire_designer.particle_analysis_ui",
        "moire_designer.gel_analysis_ui",
        "moire_design_core.particle_analysis_worker",
        "moire_design_core.gel_analysis_core",
    ],
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="DNA_Moire_Designer",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    icon=str(SPEC_DIR / "assets" / "moire-designer.ico"),
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="DNA_Moire_Designer",
)
