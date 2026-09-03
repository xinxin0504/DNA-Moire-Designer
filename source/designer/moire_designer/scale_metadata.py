"""Scale-value fallbacks for TEM files, independent from the Qt UI."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

from moire_runtime import tool_executable


_NUMBER = r"([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)"
_UNIT = r"(nm|nanometers?|nanometres?|µm|μm|um|micrometers?|micrometres?|m)"


def _to_nm(value, unit):
    unit = unit.lower().replace("μ", "µ")
    factor = 1.0
    if unit in ("µm", "um") or unit.startswith("micro"):
        factor = 1000.0
    elif unit == "m":
        factor = 1e9
    return float(value) * factor


def scale_nm_from_filename(path):
    """Read a leading scale label such as ``20 nm_sample.tif``."""
    stem = Path(path).stem
    match = re.match(
        r"^\s*" + _NUMBER + r"\s*" + _UNIT + r"(?=$|[\s_\-])",
        stem, flags=re.IGNORECASE)
    return _to_nm(match.group(1), match.group(2)) if match else None


def pixel_size_nm_from_metadata_text(text):
    """Extract common microscopy/OME physical-pixel metadata spellings."""
    text = str(text or "")
    # OME-TIFF stores value and unit in separate XML attributes.
    ome = re.search(
        r'PhysicalSizeX\s*=\s*["\']' + _NUMBER + r'["\'][^>]{0,240}?'
        r'PhysicalSizeXUnit\s*=\s*["\']' + _UNIT + r'["\']',
        text, flags=re.IGNORECASE | re.DOTALL)
    if ome:
        return _to_nm(ome.group(1), ome.group(2))
    patterns = (
        r"(?:pixel[_ ]?size(?:_x|_nm)?|pixelwidth|pixel width|physical pixel size)"
        r"\s*[:=]\s*" + _NUMBER + r"\s*" + _UNIT,
        r"(?:xscale|x scale)\s*[:=]\s*" + _NUMBER + r"\s*" + _UNIT,
    )
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return _to_nm(match.group(1), match.group(2))
    return None


def raw_pixel_size_nm(path):
    """Read physical pixel size from a TIFF's embedded textual metadata."""
    source = Path(path)
    if source.suffix.lower() not in (".tif", ".tiff"):
        return None
    texts = []
    try:
        tiffinfo = tool_executable("tiffinfo")
        if not tiffinfo:
            raise FileNotFoundError("tiffinfo is not installed")
        process = subprocess.run(
            [tiffinfo, str(source)],
            capture_output=True, text=True, timeout=15, check=False)
        texts.extend((process.stdout, process.stderr))
    except (OSError, subprocess.SubprocessError):
        pass
    # ImageDescription is often plain XML/JSON inside TIFF; this fallback also
    # works on machines without the tiffinfo utility.
    try:
        with source.open("rb") as handle:
            texts.append(handle.read(2_000_000).decode("latin-1", "ignore"))
    except OSError:
        pass
    for text in texts:
        value = pixel_size_nm_from_metadata_text(text)
        if value and value > 0:
            return value
    return None
