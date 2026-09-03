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
    # cadnano imports these through util.qtWrapImport(), which PyInstaller
    # cannot discover statically.  Omitting QtSvg makes SVGButton fall through
    # every Qt backend and terminate with util.py's final AssertionError.
    "PyQt6.QtSvg",
    "PyQt6.QtSvgWidgets",
]
datas = collect_data_files("cadnano2", include_py_files=False)
# ui_mainwindow.py displays the package version through
# importlib.metadata.version("cadnano2").  PyInstaller does not include
# distribution metadata unless requested explicitly; its absence aborts
# cadnano before the first window is created.
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
    upx=True,
    console=False,
    icon=str(SPEC_DIR / "assets" / "cadnano2.ico"),
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="cadnano2",
)
