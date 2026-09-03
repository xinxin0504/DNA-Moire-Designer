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
    "PyQt6.QtSvgWidgets",
    "PyQt6.QtPrintSupport",
    "primer3.thermoanalysis",
]

datas = [
    (str(SOURCE / "assets"), "assets"),
    (str(SOURCE / "moire_designer" / "assets"), "moire_designer/assets"),
    (str(SOURCE / "moire_designer" / "translations.json"),
     "moire_designer"),
    (str(SOURCE / "moire_design_core" / "resources"),
     "moire_design_core/resources"),
    (str(PRIVATE_ENGINE / "cadnano2"), "designer_vendor/cadnano2"),
]

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
    upx=False,
    console=False,
    target_arch="arm64",
    codesign_identity="-",
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="DNA_Moire_Designer",
)
app = BUNDLE(
    coll,
    name="DNA Moiré Designer 0.9.2.app",
    icon=str(SPEC_DIR / "assets" / "moire-design.icns"),
    bundle_identifier="org.dnamoire.designer.release092",
    info_plist={
        "CFBundleDisplayName": "DNA Moiré Designer 0.9.2",
        "CFBundleName": "DNA Moiré Designer",
        "CFBundleShortVersionString": "0.9.2",
        "CFBundleVersion": "0.9.2",
        "LSMinimumSystemVersion": "13.0",
        "LSMultipleInstancesProhibited": True,
        "NSHighResolutionCapable": True,
        "NSPrincipalClass": "NSApplication",
    },
)
