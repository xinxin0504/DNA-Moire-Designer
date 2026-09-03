"""Cross-platform runtime paths for the standalone distributions.

The Windows package contains two deliberately isolated applications:

* DNA Moiré Designer uses its validated private design engine.
* The companion cadnano executable uses the official clean cadnano2 base.

Workers are ordinary Python scripts in a development checkout.  In a frozen
PyInstaller build they are dispatched through the Designer executable so no
system Python installation is required.
"""

from __future__ import annotations

import importlib
import os
from pathlib import Path
import shutil
import sys


WORKER_MODULES = {
    "image-analysis": "moire_design_core.image_analysis_worker",
    "scale-detection": "moire_design_core.scale_detection_worker",
    "structure": "moire_design_core.structure_worker",
    "sequence-workflow": "moire_design_core.sequence_workflow_worker",
    "sequence-export": "moire_design_core.sequence_export_worker",
    "vector-export": "moire_design_core.vector_export_worker",
}


def frozen() -> bool:
    return bool(getattr(sys, "frozen", False))


def source_root() -> Path:
    """Return the read-only source/resource root bundled by PyInstaller."""
    if frozen() and hasattr(sys, "_MEIPASS"):
        return Path(getattr(sys, "_MEIPASS")).resolve()
    return Path(__file__).resolve().parent


def application_root() -> Path:
    """Return the installed directory containing the application EXE."""
    if frozen():
        return Path(sys.executable).resolve().parent
    return source_root()


def application_resources_root() -> Path:
    """Return the platform-native directory used for bundled helper tools."""
    executable = Path(sys.executable).resolve()
    if frozen() and sys.platform == "darwin" and executable.parent.name == "MacOS":
        return executable.parents[1] / "Resources"
    return application_root()


def configure_designer_engine() -> Path:
    """Put the validated private cadnano engine first on Designer's path."""
    vendor_root = source_root() / "designer_vendor"
    if vendor_root.is_dir():
        value = str(vendor_root)
        if value not in sys.path:
            sys.path.insert(0, value)
    return vendor_root


def _development_worker_python() -> Path:
    override = os.environ.get("MOIRE_WORKER_PYTHON", "").strip()
    if override:
        return Path(override).expanduser().resolve()
    return Path(sys.executable).resolve()


def worker_command(name: str, *arguments: object) -> list[str]:
    """Build a command for a private worker without relying on system Python."""
    if name not in WORKER_MODULES:
        raise KeyError("Unknown Moiré worker: %s" % name)
    args = [str(item) for item in arguments]
    if frozen():
        return [str(Path(sys.executable).resolve()), "--moire-worker", name,
                *args]
    module_path = source_root().joinpath(*WORKER_MODULES[name].split("."))
    script = module_path.with_suffix(".py")
    if not script.is_file():
        raise FileNotFoundError("Worker not found: %s" % script)
    return [str(_development_worker_python()), str(script), *args]


def dispatch_worker(name: str, arguments: list[str]) -> int:
    """Execute one worker inside the frozen Designer process."""
    module_name = WORKER_MODULES.get(name)
    if module_name is None:
        raise SystemExit("Unknown Moiré worker: %s" % name)
    module = importlib.import_module(module_name)
    previous = list(sys.argv)
    try:
        sys.argv = [getattr(module, "__file__", module_name), *arguments]
        result = module.main()
        return int(result or 0)
    finally:
        sys.argv = previous


def cadnano_executable() -> Path:
    """Locate the isolated official-base companion cadnano executable."""
    override = os.environ.get("MOIRE_CADNANO_EXECUTABLE", "").strip()
    if override:
        return Path(override).expanduser().resolve()
    if frozen() and sys.platform == "darwin":
        return (application_resources_root() / "cadnano" /
                "caDNAno Companion.app" / "Contents" / "MacOS" /
                "cadnano2")
    if sys.platform == "win32" or frozen():
        return application_root() / "cadnano" / "cadnano2.exe"
    located = shutil.which("cadnano2")
    return Path(located).resolve() if located else Path("cadnano2")


def tool_executable(name: str) -> str | None:
    """Locate an optional bundled command-line tool such as Tesseract."""
    variable = "MOIRE_%s_EXECUTABLE" % name.upper().replace("-", "_")
    override = os.environ.get(variable, "").strip()
    if override:
        return str(Path(override).expanduser().resolve())
    executable = "%s.exe" % name if sys.platform == "win32" else name
    bundled = application_resources_root() / "tools" / name / executable
    if bundled.is_file():
        return str(bundled)
    return shutil.which(executable) or shutil.which(name)
