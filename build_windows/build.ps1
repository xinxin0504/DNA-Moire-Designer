[CmdletBinding()]
param(
    [string]$TesseractPath = "",
    [string]$InnoSetupPath = "",
    [switch]$SkipTesseract,
    [switch]$SkipInstaller
)

$ErrorActionPreference = "Stop"
$BuildDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$ReleaseRoot = Split-Path -Parent $BuildDir
$Venv = Join-Path $BuildDir ".venv-build"
$Dist = Join-Path $ReleaseRoot "dist"
$Artifacts = Join-Path $ReleaseRoot "artifacts"

if (-not [Environment]::Is64BitOperatingSystem) {
    throw "A 64-bit Windows build host is required."
}

$Python = Get-Command py -ErrorAction SilentlyContinue
if ($Python) {
    & py -3.10 -c "import struct,sys; assert struct.calcsize('P') == 8; print(sys.version)"
    if ($LASTEXITCODE -ne 0) { throw "Install 64-bit Python 3.10." }
    & py -3.10 -m venv $Venv
} else {
    & python -c "import sys; assert sys.version_info[:2] == (3,10)"
    if ($LASTEXITCODE -ne 0) { throw "Install 64-bit Python 3.10." }
    & python -m venv $Venv
}

$BuildPython = Join-Path $Venv "Scripts\python.exe"
& $BuildPython -m pip install --upgrade pip wheel
& $BuildPython -m pip install -r (Join-Path $BuildDir "requirements-build.txt")
# Install only the official wheel metadata into the isolated build
# environment.  The companion's patched source directory remains first on
# PyInstaller's pathex; this installation lets copy_metadata("cadnano2")
# bundle the official name/version consumed by ui_mainwindow.py.
& $BuildPython -m pip install --no-deps `
    (Join-Path $ReleaseRoot "vendor\cadnano2-2.4.13-py3-none-any.whl")
if ($LASTEXITCODE -ne 0) {
    throw "Could not install the pinned cadnano2 metadata for packaging."
}
& $BuildPython (Join-Path $BuildDir "tools\make_windows_icons.py")
& $BuildPython (Join-Path $BuildDir "tools\audit_english_runtime.py") `
    --root (Join-Path $ReleaseRoot "source\designer") `
    --companion-root (Join-Path $ReleaseRoot "source\cadnano_companion") `
    --unresolved-only
if ($LASTEXITCODE -ne 0) {
    throw "English-only runtime audit failed. Resolve every reported presentation string before building."
}
& $BuildPython (Join-Path $ReleaseRoot "tests\test_release_surface.py")
if ($LASTEXITCODE -ne 0) { throw "Release-surface validation failed." }
& $BuildPython (Join-Path $ReleaseRoot "tests\test_release_licenses.py")
if ($LASTEXITCODE -ne 0) { throw "Release-license validation failed." }

if (Test-Path $Dist) { Remove-Item $Dist -Recurse -Force }
New-Item -ItemType Directory -Force -Path $Dist, $Artifacts | Out-Null

Push-Location $ReleaseRoot
try {
    & $BuildPython -m PyInstaller --clean --noconfirm `
        --distpath $Dist --workpath (Join-Path $ReleaseRoot "build\designer") `
        (Join-Path $BuildDir "DNA_Moire_Designer.spec")
    if ($LASTEXITCODE -ne 0) {
        throw "DNA Moire Designer PyInstaller build failed with exit code $LASTEXITCODE."
    }
    $DesignerExe = Join-Path $Dist "DNA_Moire_Designer\DNA_Moire_Designer.exe"
    if (-not (Test-Path $DesignerExe)) {
        throw "DNA Moire Designer build reported success, but the executable is missing: $DesignerExe"
    }
    & $BuildPython -m PyInstaller --clean --noconfirm `
        --distpath $Dist --workpath (Join-Path $ReleaseRoot "build\cadnano") `
        (Join-Path $BuildDir "cadnano2_companion.spec")
    if ($LASTEXITCODE -ne 0) {
        throw "cadnano2 companion PyInstaller build failed with exit code $LASTEXITCODE."
    }
    $CadnanoExe = Join-Path $Dist "cadnano2\cadnano2.exe"
    if (-not (Test-Path $CadnanoExe)) {
        throw "cadnano2 build reported success, but the executable is missing: $CadnanoExe"
    }
} finally {
    Pop-Location
}

$DesignerDist = Join-Path $Dist "DNA_Moire_Designer"
$CadnanoTarget = Join-Path $DesignerDist "cadnano"
Copy-Item (Join-Path $Dist "cadnano2") $CadnanoTarget -Recurse -Force

if (-not $SkipTesseract) {
    if (-not $TesseractPath) {
        $Candidates = @(
            "$env:ProgramFiles\Tesseract-OCR",
            "${env:ProgramFiles(x86)}\Tesseract-OCR"
        )
        $TesseractPath = $Candidates | Where-Object {
            Test-Path (Join-Path $_ "tesseract.exe")
        } | Select-Object -First 1
    }
    if (-not $TesseractPath -or
        -not (Test-Path (Join-Path $TesseractPath "tesseract.exe"))) {
        throw "Tesseract OCR was not found. Supply -TesseractPath or use -SkipTesseract explicitly."
    }
    $TesseractTarget = Join-Path $DesignerDist "tools\tesseract"
    Copy-Item $TesseractPath $TesseractTarget -Recurse -Force
}

# Keep the main GPL terms, institutional copyright notice, and all curated
# third-party notices inside both the installed application and portable ZIP.
$LegalTarget = Join-Path $DesignerDist "legal"
$LicenseTarget = Join-Path $LegalTarget "licenses"
New-Item -ItemType Directory -Force -Path $LegalTarget, $LicenseTarget | Out-Null
Copy-Item (Join-Path $ReleaseRoot "LICENSE") $LegalTarget -Force
Copy-Item (Join-Path $ReleaseRoot "COPYRIGHT") $LegalTarget -Force
Copy-Item (Join-Path $ReleaseRoot "README.md") $LegalTarget -Force
Copy-Item (Join-Path $ReleaseRoot "THIRD_PARTY_NOTICES.md") $LegalTarget -Force
Copy-Item (Join-Path $ReleaseRoot "licenses\*") $LicenseTarget -Recurse -Force

# Preserve the exact license metadata installed for every frozen Python
# dependency. This complements the curated, human-readable license inventory.
$SitePackagesOutput = @(& $BuildPython -c "import sysconfig; print(sysconfig.get_path('purelib'))")
if ($LASTEXITCODE -ne 0 -or $SitePackagesOutput.Count -eq 0) {
    throw "Unable to resolve the build environment's Python package directory."
}
$SitePackages = ([string]$SitePackagesOutput[-1]).Trim()
if (-not (Test-Path -LiteralPath $SitePackages -PathType Container)) {
    throw "The resolved Python package directory does not exist: $SitePackages"
}
$PackageLicenseTarget = Join-Path $LicenseTarget "python-packages"
New-Item -ItemType Directory -Force -Path $PackageLicenseTarget | Out-Null
$DistInfoDirectories = @(
    Get-ChildItem -LiteralPath $SitePackages -Directory -ErrorAction Stop |
        Where-Object { $_.Name -like "*.dist-info" }
)
if ($DistInfoDirectories.Count -eq 0) {
    throw "No installed Python distribution metadata was found for license packaging."
}
foreach ($DistInfo in $DistInfoDirectories) {
    # Preserve the complete metadata directory. Modern wheels commonly keep
    # license files below a nested ``licenses`` folder, while older wheels put
    # them at the dist-info root. Copying the complete directory handles both
    # layouts and retains the package name/version/license declarations.
    Copy-Item $DistInfo.FullName `
        (Join-Path $PackageLicenseTarget $DistInfo.Name) -Recurse -Force
}
$PackagedLicenseFiles = @(Get-ChildItem $PackageLicenseTarget -Recurse -File |
    Where-Object {
        $_.Name -match '^(LICENSE|LICENCE|COPYING|NOTICE|AUTHORS)'
    })
if ($PackagedLicenseFiles.Count -eq 0) {
    throw "Python package metadata was copied, but no license or notice files were found."
}
Write-Host "Packaged Python dependency license/notice files:" `
    $PackagedLicenseFiles.Count

# PowerShell does not reliably wait for GUI-subsystem executables invoked with
# the call operator.  Start-Process -Wait returns the real process exit code,
# so a traceback can no longer be followed by a nominally successful package.
$CadnanoSelfTest = Start-Process `
    -FilePath (Join-Path $CadnanoTarget "cadnano2.exe") `
    -ArgumentList "--self-test" -Wait -PassThru
if ($CadnanoSelfTest.ExitCode -ne 0) {
    throw "cadnano companion frozen GUI self-test failed with exit code $($CadnanoSelfTest.ExitCode)."
}
$DesignerSelfTest = Start-Process `
    -FilePath (Join-Path $DesignerDist "DNA_Moire_Designer.exe") `
    -ArgumentList "--self-test" -Wait -PassThru
if ($DesignerSelfTest.ExitCode -ne 0) {
    throw "DNA Moiré Designer frozen GUI self-test failed with exit code $($DesignerSelfTest.ExitCode)."
}

$Portable = Join-Path $Artifacts "DNA-Moire-Designer-0.9.2-Windows-x64-Portable.zip"
$ArchiveCreated = $false
for ($Attempt = 1; $Attempt -le 4; $Attempt++) {
    if (Test-Path $Portable) { Remove-Item $Portable -Force }
    try {
        Compress-Archive -Path (Join-Path $DesignerDist "*") `
            -DestinationPath $Portable -CompressionLevel Optimal `
            -ErrorAction Stop
        $ArchiveCreated = $true
        break
    } catch {
        if ($Attempt -eq 4) { throw }
        Write-Warning "Portable archive attempt $Attempt failed because a packaged file is temporarily busy. Retrying."
        Start-Sleep -Seconds (2 * $Attempt)
    }
}
if (-not $ArchiveCreated -or -not (Test-Path $Portable)) {
    throw "Portable archive was not created."
}

if (-not $SkipInstaller) {
    if ($InnoSetupPath) {
        if (-not (Test-Path $InnoSetupPath -PathType Leaf)) {
            throw "The supplied Inno Setup compiler does not exist: $InnoSetupPath"
        }
        $Iscc = (Resolve-Path $InnoSetupPath).Path
    } else {
        $IsccCandidates = @(
            "$env:ProgramFiles\Inno Setup 7\ISCC.exe",
            "${env:ProgramFiles(x86)}\Inno Setup 7\ISCC.exe",
            "$env:LOCALAPPDATA\Programs\Inno Setup 7\ISCC.exe",
            "$env:ProgramFiles\Inno Setup 6\ISCC.exe",
            "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe",
            "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe"
        )
        $Iscc = $IsccCandidates | Where-Object { Test-Path $_ } |
            Select-Object -First 1
    }
    if (-not $Iscc) { throw "Inno Setup 6 or 7 was not found." }
    # Inno Setup still encounters the legacy Windows path-length boundary
    # while recursively enumerating a deeply extracted PyInstaller onedir.
    # Compile from a short temporary source root so the kit works regardless
    # of how many parent directories the user created during extraction.
    $InnoStage = Join-Path $env:TEMP `
        ("DNA-Moire-Inno-" + [Guid]::NewGuid().ToString("N").Substring(0, 8))
    $InnoReleaseRoot = Join-Path $InnoStage "app"
    try {
        New-Item -ItemType Directory -Force -Path $InnoReleaseRoot | Out-Null
        Copy-Item (Join-Path $DesignerDist "*") $InnoReleaseRoot `
            -Recurse -Force
        if (-not (Test-Path (Join-Path $InnoReleaseRoot `
                    "DNA_Moire_Designer.exe"))) {
            throw "The short Inno Setup staging copy is incomplete."
        }
        & $Iscc "/DReleaseRoot=$InnoReleaseRoot" `
            "/DOutputDir=$Artifacts" `
            (Join-Path $BuildDir "installer\DNA_Moire_Designer.iss")
        if ($LASTEXITCODE -ne 0) { throw "Inno Setup compilation failed." }
    } finally {
        if (Test-Path $InnoStage) {
            Remove-Item $InnoStage -Recurse -Force
        }
    }
}

$HashFiles = Get-ChildItem $Artifacts -File | Where-Object {
    $_.Extension -in ".exe", ".zip"
}
$HashLines = foreach ($File in $HashFiles) {
    $Hash = Get-FileHash $File.FullName -Algorithm SHA256
    "$($Hash.Hash.ToLower())  $($File.Name)"
}
Set-Content -Path (Join-Path $Artifacts "SHA256SUMS.txt") `
    -Value $HashLines -Encoding ascii

Write-Host "Windows release assembled at: $DesignerDist"
Write-Host "Installer output: $Artifacts"
