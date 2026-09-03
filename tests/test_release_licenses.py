#!/usr/bin/env python3
"""Verify the release's main and third-party licensing surface."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_COPYRIGHT = "2. Physics Institute, University of Stuttgart"

license_text = (ROOT / "LICENSE").read_text(encoding="utf-8")
copyright_text = (ROOT / "COPYRIGHT").read_text(encoding="utf-8")
readme_text = (ROOT / "README.md").read_text(encoding="utf-8")
notices_text = (ROOT / "THIRD_PARTY_NOTICES.md").read_text(encoding="utf-8")
readme_words = " ".join(readme_text.split())

assert "GNU GENERAL PUBLIC LICENSE" in license_text
assert "Version 3, 29 June 2007" in license_text
assert EXPECTED_COPYRIGHT in copyright_text
assert "GPL-3.0-only" in copyright_text
assert "Generated research outputs" in readme_text
assert "complete corresponding source" in readme_words
assert "SOURCE_CODE.md" not in readme_text
assert "local sequence-aware JSON input/output changes" in readme_words
assert "caDNAno 2.4.13" in notices_text
assert "MIT License" in notices_text

required_licenses = {
    "cadnano2-LICENSE.txt",
    "GPL-3.0.txt",
    "PyInstaller-license.txt",
    "PyQt6-GPL-3.0.txt",
    "PyQt6-sip-license.txt",
    "Python-PSF-license.txt",
    "Qt6-LGPL-3.0.txt",
    "Tesseract-Apache-2.0.txt",
    "primer3-license.txt",
}
assert required_licenses <= {
    path.name for path in (ROOT / "licenses").iterdir() if path.is_file()
}
assert not any(ROOT.rglob("SOURCE_CODE.md"))

windows_build = ROOT / "build_windows" / "build.ps1"
if windows_build.exists():
    build_text = windows_build.read_text(encoding="utf-8")
    installer_text = (ROOT / "build_windows" / "installer" /
                      "DNA_Moire_Designer.iss").read_text(encoding="utf-8")
    assert 'Join-Path $DesignerDist "legal"' in build_text
    assert '"python-packages"' in build_text
    assert "LicenseFile=..\\..\\LICENSE" in installer_text
    assert "AppPublisher=" + EXPECTED_COPYRIGHT in installer_text

macos_build = ROOT / "build_macos" / "build.sh"
if macos_build.exists():
    build_text = macos_build.read_text(encoding="utf-8")
    assert 'Contents/Resources/legal' in build_text
    assert 'PACKAGE_LICENSE_TARGET=' in build_text
    assert 'DMG_STAGE/Legal notices' in build_text

print("release_licenses=GPL-3.0-only main license and curated notices present")
