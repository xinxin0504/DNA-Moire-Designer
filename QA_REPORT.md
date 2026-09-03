# Windows release-kit QA report

Date: 2026-09-03

## Scope verified on macOS source host

- The release snapshot is independent of the active Designer and cn24x installation.
- No file in the active installation was overwritten.
- The Windows Designer exposes the validated design, sequence, final-export, and Moiré Analysis workflows.
- The Windows release is English-only. The language selector has been removed, and projects saved with a legacy Chinese-language preference are opened in English.
- The Designer's validated private cadnano engine and the official-base companion cadnano are stored in separate package roots.
- The companion package differs from the official cadnano2 package only in `legacydecoder.py` and `legacyencoder.py`; its separate launcher adds positional JSON opening and self-test support.
- The latest accepted-state invalidation chain is included: changes to accepted Moiré parameters, DNA design, scaffold assignments, or SST-input assignments revoke the affected acceptance and all dependent downstream state.
- The Step 2 design report now occupies a full-width, vertically resizable bottom row and uses the same regular-weight report typography as the Step 3 SST-input report.
- The normal-staple histogram retains every nonzero percentage bin, including sub-1% edge bins.
- Square–Kagome linked layer lengths retain independent SST-input assignments; shared input mapping remains restricted to equal-topology Square–Square and Kagome–Kagome layers.
- Tesseract is resolved through the packaged cross-platform runtime. The Windows build therefore uses the bundled `tools\tesseract` executable and its `tessdata` directory rather than a host-specific macOS path.
- Spacing, Twist, and Seed Z2 insertion/deletion settings are linked to the actual number of nominal 8-bp spacing domains. For `N = spacing / 8`, the deletion bound is `-3N` and the insertion bound is `min(3N, 10)`: 0 bp permits only 0; 8 bp permits -3 through +3; 16 bp permits -6 through +6; 24 bp permits -9 through +9; and 32 bp permits -12 through +10.
- At 0-bp spacing, Twist and mean insertion/deletion are fixed at 0. The interface, parameter solver, accepted-project validation, and structure worker all enforce the same rule.
- Fixed-scaffold insertion headroom is checked before structure generation. A nominal domain bound is therefore necessary but not sufficient when the selected scaffold routes contain fewer available nucleotides than the requested total insertion count.
- The validated auxiliary-helix routing for 0/8-bp spacing, physical SST-layer attribution, complete/SST-only input-output-capture sequence mapping, and moment-balanced insertion/deletion allocation were preserved.
- Before invoking Inno Setup, the assembled PyInstaller application is copied to a unique short temporary source root. This prevents the legacy Windows path-length boundary from blocking `Setup.exe` creation when the build kit was extracted below repeated or deeply nested directory names; the temporary staging copy is removed after compilation.
- The Windows and macOS builds both package dependency distribution metadata and fail if no Python-package license or notice files are present. This covers modern wheels that store notices in nested `.dist-info/licenses` directories as well as older root-level layouts.

## Automated results

### Release-facing source regressions

Windows release-license and release-surface checks passed.

The complete source suite ran `196 tests`: `181` passed and `15` PyQt UI tests were skipped on the macOS source host because the local test interpreter does not contain PyQt. There were no failures or errors. A focused `48-test` release and workflow set also passed before the complete suite.

Passing coverage includes core Moiré parameter/calibration logic, fixed-Seed partitioning, spacing-dependent insertion/deletion limits, Square/Kagome SST behavior, strict accepted-input sequence export, Moiré Analysis labels and FFT/scale behavior, cylindrical-preview colour assignment, workflow structure, English UI/report strings, orthogonal-sequence workbook formatting, and the previously validated 0/8/16-bp spacing cases.

The focused v14.3 spacing regression covered the 0/8/16/24/32-bp domain boundaries, both insertion and deletion directions, UI control state and suffixes, stale-project preflight rejection, fixed-scaffold headroom reporting, real Square and Kagome-family structure generation, and final-export encoding. The frozen Apple Silicon worker additionally generated an 8-bp/+3 structure containing 144 insertion placements across 48 Seed helices, with no 8-bp domain receiving more than three edits and with the established deterministic first-moment-balanced allocation method. The separate 242-case Kagome matrix was intentionally not run for this focused correction.

Scaffold identity correction:

- the Designer display name `CS3L` maps exactly to cadnano key `CS3L_7559` (7559 nt; SHA-256 `1cde42cdddcdbb232536045e25e8021cab458baffda51021ed3362dd368bbf71`)
- the Designer display name `CS4` maps exactly to cadnano key `CS4_7557` (7557 nt; SHA-256 `ddb6cefc88a2f587d64296e56d70379bb27a9486e52a32b771fe9193c0fa579a`)
- the former labels `CS3` and `CS4-L` are accepted only as legacy project aliases and normalize to `CS3L` and `CS4`
- the complete scaffold/SST sequence-workflow regression and four focused catalogue/assignment regressions pass with the corrected names

### English-only presentation audit

- Designer runtime strings containing CJK characters audited: 878
- unresolved presentation strings: 0
- English translation-catalog values containing CJK characters: 0
- companion cadnano presentation strings containing CJK characters: 0
- the same audit is run automatically before every Windows build and stops the build if an unresolved presentation string is found
- `translations.json` is explicitly included in the PyInstaller data bundle; the earlier frozen build omitted this adjacent JSON resource even though the source-host audit could read it
- the frozen Designer self-test verifies representative full-catalog translations and scans the constructed window for CJK text before the installer is created
- the Designer self-test imports the real `moire_design_core.models` module (plural), matching the production package
- both GUI-subsystem executables are launched with `Start-Process -Wait -PassThru`; a non-zero self-test exit now aborts the build instead of allowing a false-positive installer

The audit includes the visible pages, dynamically opened dialogs and reports, exported CSV/XLSX/SVG text, sequence tables, warnings, progress messages, and restored-project state.

The current terminology review also verifies sentence-case interface text,
`caDNAno`, `oxView`, and `oxDNA` capitalization, consistent `SST sublattice`
wording, separate lattice-constant values for both layers, and the final-export
folder names `Input parameters`, `caDNAno design files`, `Oligonucleotide
sequences`, and `PDB/oxView files` (represented with a filesystem-safe Unicode
slash in the physical directory name).

### Windows SST-input command-length regression

- a synthetic assignment payload larger than 40 KB completed through the real sequence worker
- 128-assignment mock regression passed with exact argument recovery
- oversized payloads are now transferred through a short temporary response-file path rather than the Windows command line
- accepted intermediate output uses the fixed short filename `sequenced_design.json`; final export names are unchanged

### Synthetic companion sequence round trip

- input records: 3
- output records: 3
- input nucleotides: 14,668
- output nucleotides: 14,668
- exact records: yes
- Designer metadata preserved: yes
- input/output SHA-256: `41e1faad4089c3f4d8eecef1c187da1367ee7368d67c30da33f9eff896d4bb3e`

### Actual Designer export round trip

Fixture used for local QA:
`sst_scaffold_staple_capture_with_sequence.json`

- input records: 106
- output records: 106
- input nucleotides: 18,012
- output nucleotides: 18,012
- exact 5′ anchors and sequence records: yes
- Designer metadata preserved: yes
- input/output SHA-256: `3a465e6ed8ea55d22b5d8e67b767acc48f41a361b42095c4ca08b4eaf4d292a0`

### Syntax/static build checks

- Python compileall: passed
- both PyInstaller spec files compile as Python: passed
- the cadnano PyInstaller specification explicitly collects all `cadnano2` submodules and includes `PyQt6.QtSvg` and `PyQt6.QtSvgWidgets`
- the official `cadnano2` wheel metadata is installed into the isolated build environment and copied into the frozen companion, satisfying `importlib.metadata.version("cadnano2")`
- official wheel SHA-256 verified: `d877282c8a782b079070248b2743016de2ecae523939be857df20d3672b0e039`

## Windows-only acceptance still required

PyInstaller does not cross-compile Windows executables on macOS. The provided PowerShell build performs a complete offscreen GUI initialization of both frozen executables before compiling the installer. This specifically exercises the dynamically imported cadnano QtSvg modules. A final Windows GUI acceptance pass must confirm:

1. startup and workflow navigation;
2. one end-to-end design, Accept Current SST Input with a full-size assignment set, and final export;
3. sequence-bearing JSON open/save/reopen in companion cadnano;
4. Tesseract-backed scale OCR and one Moiré Analysis run;
5. installer install/uninstall and both shortcuts.
