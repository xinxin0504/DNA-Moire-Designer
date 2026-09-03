# -*- mode: python ; coding: utf-8 -*-
from pathlib import Path
import sys

from PyInstaller.utils.hooks import (collect_data_files, collect_submodules,
                                     copy_metadata)

SPEC_DIR = Path(SPECPATH)
ROOT = SPEC_DIR.parent
SOURCE = ROOT / "source" / "cadnano_companion"
sys.path.insert(0, str(SOURCE))

hiddenimports = collect_submodules("cadnano2") + [
    "PyQt6.QtSvg",
    "PyQt6.QtSvgWidgets",
]
datas = collect_data_files("cadnano2", include_py_files=False)
datas += copy_metadata("cadnano2")

a = Analysis(
    [str(SOURCE / "run_cadnano_companion.py")],
    pathex=[str(SOURCE)],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="cadnano2",
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
    name="cadnano2",
)
app = BUNDLE(
    coll,
    name="caDNAno Companion.app",
    icon=str(SPEC_DIR / "assets" / "cadnano2.icns"),
    bundle_identifier="org.dnamoire.cadnano-companion.release092",
    info_plist={
        "CFBundleDisplayName": "caDNAno Companion",
        "CFBundleName": "caDNAno Companion",
        "CFBundleShortVersionString": "2.4.13",
        "CFBundleVersion": "2.4.13",
        "LSMinimumSystemVersion": "13.0",
        "LSMultipleInstancesProhibited": False,
        "NSHighResolutionCapable": True,
        "NSPrincipalClass": "NSApplication",
    },
)
