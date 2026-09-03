# DNA Moiré Designer — Windows x64 release kit

This kit builds one Windows installer containing two isolated applications:

1. **DNA Moiré Designer** — the validated English-only design, sequence, final-export, and **Moiré Analysis** workflows.
2. **cadnano2 Companion** — official cadnano2 2.4.13 plus the smallest JSON-I/O extension needed to open, edit, save, and reopen Designer JSON files with applied scaffold sequences.

The Designer continues to use its own validated private cadnano design engine. The companion cadnano installation is not used for automatic routing or sequence design.

## Build host

- Windows 10 or 11, x64
- 64-bit CPython 3.10
- Inno Setup 6 or 7
- Tesseract OCR for Windows, including `tesseract.exe` and its `tessdata` directory
- Internet access for the first Python dependency installation

The final Windows executables must be frozen on Windows. PyInstaller does not cross-compile a Windows executable from macOS.

## One-command build

Open PowerShell in this directory and run:

```powershell
Set-ExecutionPolicy -Scope Process Bypass
.\build_windows\build.ps1 -TesseractPath "C:\Program Files\Tesseract-OCR"
```

If Tesseract is already installed at the standard location, omit `-TesseractPath`.
If Inno Setup is installed in a non-standard directory, supply its compiler
explicitly:

```powershell
.\build_windows\build.ps1 `
  -TesseractPath "C:\Program Files\Tesseract-OCR" `
  -InnoSetupPath "C:\path\to\Inno Setup\ISCC.exe"
```

The script:

1. validates 64-bit Python 3.10;
2. creates a private build virtual environment;
3. installs pinned build/runtime dependencies;
4. rejects any untranslated public runtime text with the English-only source audit;
5. creates both PyInstaller onedir applications;
6. nests the clean-base cadnano companion under `cadnano\`;
7. copies Tesseract under `tools\tesseract\`;
8. initializes the complete Designer and cadnano GUIs off-screen as frozen self-tests, waits for both GUI processes to finish, verifies the packaged English catalog, scans the Designer window for CJK text, verifies cadnano package metadata, and aborts immediately if either process returns a non-zero exit code;
9. creates a portable ZIP and an Inno Setup installer.

Expected outputs:

- `artifacts\DNA-Moire-Designer-0.9.2-Windows-x64-Setup.exe`
- `artifacts\DNA-Moire-Designer-0.9.2-Windows-x64-Portable.zip`
- `artifacts\SHA256SUMS.txt`

Use `-SkipInstaller` only to create and test the onedir applications. Use `-SkipTesseract` only for packaging diagnostics: Moiré OCR functionality will then be incomplete.

## Required Windows acceptance checks

After installation, run these in PowerShell:

```powershell
& "$env:LOCALAPPDATA\Programs\DNA Moire Designer\DNA_Moire_Designer.exe" --self-test
& "$env:LOCALAPPDATA\Programs\DNA Moire Designer\cadnano\cadnano2.exe" --self-test
```

Then perform one GUI acceptance cycle:

1. create/open a `.moire.json` project;
2. generate and accept the final DNA design;
3. complete scaffold and SST-input sequence assignment;
4. export the final package and inspect sequence tables;
5. open the sequence-bearing JSON in the bundled cadnano companion;
6. save it under a new name, close it, reopen it, and verify the same scaffold sequences remain applied;
7. switch to Analysis Mode, load a TEM image, verify scale detection, and run Moiré Analysis.

The installed product has no language selector. Older project files carrying
`interface_language: zh_CN` are opened in English, and generated reports,
dialogs, spreadsheets, CSV files, and SVG labels remain English.

## Companion cadnano scope

The clean-base companion is pinned to official cadnano2 2.4.13. Its extension is intentionally limited to:

- loading and saving top-level `scaffold_sequences` by exact 5′ anchor;
- preserving Designer top-level metadata and scaffold colours;
- accepting a positional JSON path from the Designer launcher.

No Designer routing, SST/capture generation, orthogonal-sequence design, or analysis code is added to the companion.

## Distribution note

Review `THIRD_PARTY_NOTICES.md` and the licenses installed by Python dependencies before public redistribution. PyQt6 distribution is subject to its GPL/commercial licensing terms. Tesseract binaries are supplied by the Windows build operator and are not included in this source kit.
