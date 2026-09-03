#!/bin/bash
set -euo pipefail

BUILD_DIR="$(cd "$(dirname "$0")" && pwd)"
RELEASE_ROOT="$(cd "$BUILD_DIR/.." && pwd)"
BUILD_PYTHON="${BUILD_PYTHON:-$RELEASE_ROOT/../build_venv/bin/python}"
DIST="$RELEASE_ROOT/dist"
BUILD="$RELEASE_ROOT/build"
ARTIFACTS="$RELEASE_ROOT/artifacts"
DESIGNER_APP="$DIST/DNA Moiré Designer 0.9.2.app"
COMPANION_APP="$DIST/caDNAno Companion.app"
export PYINSTALLER_CONFIG_DIR="$BUILD/pyinstaller-cache"
export MPLCONFIGDIR="$BUILD/matplotlib-cache"

if [[ "$(uname -m)" != "arm64" ]]; then
  echo "An Apple Silicon (arm64) build host is required." >&2
  exit 2
fi
if [[ ! -x "$BUILD_PYTHON" ]]; then
  echo "Build Python not found: $BUILD_PYTHON" >&2
  exit 2
fi

"$BUILD_PYTHON" -c 'import PyInstaller,PyQt6,numpy,scipy,matplotlib,PIL,primer3'
"$BUILD_PYTHON" "$BUILD_DIR/tools/make_macos_icons.py"
"$BUILD_PYTHON" "$RELEASE_ROOT/build_windows/tools/audit_english_runtime.py" \
  --root "$RELEASE_ROOT/source/designer" \
  --companion-root "$RELEASE_ROOT/source/cadnano_companion" \
  --unresolved-only
"$BUILD_PYTHON" "$RELEASE_ROOT/tests/test_release_surface.py"
"$BUILD_PYTHON" "$RELEASE_ROOT/tests/test_release_licenses.py"

rm -rf "$DIST" "$BUILD" "$ARTIFACTS"
mkdir -p "$DIST" "$BUILD" "$ARTIFACTS"

"$BUILD_PYTHON" -m PyInstaller --clean --noconfirm \
  --distpath "$DIST" --workpath "$BUILD/designer" \
  "$BUILD_DIR/DNA_Moire_Designer.spec"
"$BUILD_PYTHON" -m PyInstaller --clean --noconfirm \
  --distpath "$DIST" --workpath "$BUILD/cadnano" \
  "$BUILD_DIR/cadnano2_companion.spec"

if [[ ! -x "$DESIGNER_APP/Contents/MacOS/DNA_Moire_Designer" ]]; then
  echo "Designer application executable is missing." >&2
  exit 3
fi
if [[ ! -x "$COMPANION_APP/Contents/MacOS/cadnano2" ]]; then
  echo "caDNAno Companion executable is missing." >&2
  exit 3
fi

mkdir -p "$DESIGNER_APP/Contents/Resources/cadnano"
ditto "$COMPANION_APP" \
  "$DESIGNER_APP/Contents/Resources/cadnano/caDNAno Companion.app"

TESSERACT_SOURCE="${TESSERACT_SOURCE:-/opt/homebrew/bin/tesseract}"
TESSDATA_SOURCE="${TESSDATA_SOURCE:-/opt/homebrew/share/tessdata}"
TESSERACT_TARGET="$DESIGNER_APP/Contents/Resources/tools/tesseract"
if [[ ! -x "$TESSERACT_SOURCE" || ! -d "$TESSDATA_SOURCE" ]]; then
  echo "A local arm64 Tesseract installation is required for the OCR-enabled build." >&2
  exit 4
fi
"$BUILD_PYTHON" "$BUILD_DIR/tools/bundle_macho_dependencies.py" \
  "$TESSERACT_SOURCE" "$TESSERACT_TARGET"
ditto "$TESSDATA_SOURCE" "$TESSERACT_TARGET/tessdata"
cp /opt/homebrew/Cellar/tesseract/5.3.3/LICENSE \
  "$TESSERACT_TARGET/TESSERACT-LICENSE.txt"

LEGAL_TARGET="$DESIGNER_APP/Contents/Resources/legal"
mkdir -p "$LEGAL_TARGET/licenses"
cp "$RELEASE_ROOT/LICENSE" "$LEGAL_TARGET/"
cp "$RELEASE_ROOT/COPYRIGHT" "$LEGAL_TARGET/"
cp "$RELEASE_ROOT/README.md" "$LEGAL_TARGET/"
cp "$RELEASE_ROOT/THIRD_PARTY_NOTICES.md" "$LEGAL_TARGET/"
/bin/cp -R -X "$RELEASE_ROOT/licenses/." "$LEGAL_TARGET/licenses/"

# Preserve the exact license metadata installed for every frozen Python
# dependency. This complements the curated, human-readable license inventory.
SITE_PACKAGES="$($BUILD_PYTHON -c 'import site; print(site.getsitepackages()[0])')"
PACKAGE_LICENSE_TARGET="$LEGAL_TARGET/licenses/python-packages"
mkdir -p "$PACKAGE_LICENSE_TARGET"
for DIST_INFO in "$SITE_PACKAGES"/*.dist-info; do
  if find "$DIST_INFO" -maxdepth 3 -type f \
      \( -iname 'LICENSE*' -o -iname 'LICENCE*' -o -iname 'COPYING*' \
         -o -iname 'NOTICE*' -o -iname 'AUTHORS*' \) | grep -q .; then
    /bin/cp -R -X "$DIST_INFO" \
      "$PACKAGE_LICENSE_TARGET/$(basename "$DIST_INFO")"
  fi
done
PACKAGE_LICENSE_COUNT="$(find "$PACKAGE_LICENSE_TARGET" -type f \
    \( -iname 'LICENSE*' -o -iname 'LICENCE*' -o -iname 'COPYING*' \
       -o -iname 'NOTICE*' -o -iname 'AUTHORS*' \) | wc -l | tr -d ' ')"
if [[ "$PACKAGE_LICENSE_COUNT" -eq 0 ]]; then
  echo "No Python dependency license or notice files were packaged." >&2
  exit 6
fi
echo "Packaged Python dependency license/notice files: $PACKAGE_LICENSE_COUNT"

PACKAGE_STAGE="$(mktemp -d /private/tmp/dna-moire-macos-v13.XXXXXX)"
cleanup_package_stage() {
  rm -rf "$PACKAGE_STAGE"
}
trap cleanup_package_stage EXIT
STAGED_APP="$PACKAGE_STAGE/$(basename "$DESIGNER_APP")"
ditto --noqtn "$DESIGNER_APP" "$STAGED_APP"

chmod -R u+w "$STAGED_APP"
# Remove only metadata that codesign rejects.  Do not use `xattr -c` here:
# Homebrew binaries carry a protected com.apple.provenance attribute on newer
# macOS releases, and attempting to clear that unrelated attribute aborts the
# otherwise reproducible workspace-only build.
xattr -dr com.apple.FinderInfo "$STAGED_APP" 2>/dev/null || true
xattr -dr com.apple.ResourceFork "$STAGED_APP" 2>/dev/null || true
xattr -dr com.apple.quarantine "$STAGED_APP" 2>/dev/null || true
codesign --force --deep --sign - \
  "$STAGED_APP/Contents/Resources/cadnano/caDNAno Companion.app"
codesign --force --deep --sign - "$STAGED_APP"
codesign --verify --deep --strict --verbose=2 "$STAGED_APP"

QT_QPA_PLATFORM=offscreen \
  "$STAGED_APP/Contents/Resources/cadnano/caDNAno Companion.app/Contents/MacOS/cadnano2" \
  --self-test > "$ARTIFACTS/cadnano-self-test.json"
QT_QPA_PLATFORM=offscreen \
  "$STAGED_APP/Contents/MacOS/DNA_Moire_Designer" \
  --self-test > "$ARTIFACTS/designer-self-test.json"

PORTABLE="$ARTIFACTS/DNA-Moire-Designer-0.9.2-macOS-Apple-Silicon-Portable.zip"
ditto -c -k --norsrc --keepParent "$STAGED_APP" "$PORTABLE"

DMG_STAGE="$PACKAGE_STAGE/dmg-stage"
mkdir -p "$DMG_STAGE"
ditto --noqtn "$STAGED_APP" "$DMG_STAGE/$(basename "$STAGED_APP")"
ln -s /Applications "$DMG_STAGE/Applications"
mkdir -p "$DMG_STAGE/Legal notices"
cp "$RELEASE_ROOT/LICENSE" "$DMG_STAGE/Legal notices/"
cp "$RELEASE_ROOT/COPYRIGHT" "$DMG_STAGE/Legal notices/"
cp "$RELEASE_ROOT/THIRD_PARTY_NOTICES.md" "$DMG_STAGE/Legal notices/"
DMG="$ARTIFACTS/DNA-Moire-Designer-0.9.2-macOS-Apple-Silicon.dmg"
hdiutil create -volname "DNA Moiré Designer 0.9.2" \
  -srcfolder "$DMG_STAGE" -ov -format UDZO "$DMG"

cp "$RELEASE_ROOT/README_MACOS_APPLE_SILICON.md" "$ARTIFACTS/"
cp "$RELEASE_ROOT/LICENSE" "$ARTIFACTS/"
cp "$RELEASE_ROOT/COPYRIGHT" "$ARTIFACTS/"
cp "$RELEASE_ROOT/THIRD_PARTY_NOTICES.md" "$ARTIFACTS/"
(
  cd "$ARTIFACTS"
  shasum -a 256 \
    "$(basename "$DMG")" \
    "$(basename "$PORTABLE")" > SHA256SUMS.txt
)

echo "Apple Silicon application archive source: $STAGED_APP"
echo "Release artifacts: $ARTIFACTS"
