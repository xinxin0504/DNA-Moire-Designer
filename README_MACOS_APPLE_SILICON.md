# DNA Moiré Designer 0.9.2 for Apple Silicon

This package is the macOS Apple Silicon build of the validated cross-platform
0.9.2 release source. It contains the English-only DNA design workflow,
sequence workflow, final export, Moiré analysis, and the isolated official-base
caDNAno Companion.

## Install

Open the DMG and drag **DNA Moiré Designer 0.9.2.app** to Applications, or
extract the Portable ZIP and run the app from any writable folder.

This local validation build is ad-hoc signed. On another Mac, the first launch
may require Control-clicking the app, choosing **Open**, and confirming once.
Public distribution without that warning requires Apple Developer ID signing
and notarization.

## Compatibility

- Apple Silicon only (`arm64`)
- macOS 13 or later
- No system Python, Homebrew, Tesseract, or separate caDNAno installation is
  required at runtime
- Existing DNA Moiré Designer installations are not modified by this package

The embedded caDNAno Companion preserves Designer `scaffold_sequences` and
supported top-level metadata when opening and saving JSON files.

## Build isolation

The application was built in a workspace-only directory. The build process did
not write to `/Applications`, replace an installed application, or modify the
`cn24x` virtual environment.

## Gatekeeper note

The build is ad-hoc signed and locally verified, but is not notarized. Before
external release, rebuild or re-sign it with a Developer ID Application
certificate and submit the DMG for Apple notarization.
