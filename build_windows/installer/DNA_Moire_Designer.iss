#ifndef ReleaseRoot
  #error ReleaseRoot must point to the assembled Designer onedir directory.
#endif
#ifndef OutputDir
  #define OutputDir "..\..\artifacts"
#endif

#define AppName "DNA Moiré Designer"
#define AppVersion "0.9.2"

[Setup]
AppId={{C4CB72A8-F3B4-4EF6-AEA6-DFE6C0C4B95B}
AppName={#AppName}
AppVersion={#AppVersion}
AppPublisher=2. Physics Institute, University of Stuttgart
DefaultDirName={autopf}\DNA Moire Designer
DefaultGroupName={#AppName}
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
PrivilegesRequired=lowest
OutputDir={#OutputDir}
OutputBaseFilename=DNA-Moire-Designer-{#AppVersion}-Windows-x64-Setup
SetupIconFile=..\assets\moire-designer.ico
UninstallDisplayIcon={app}\DNA_Moire_Designer.exe
LicenseFile=..\..\LICENSE

[Files]
Source: "{#ReleaseRoot}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{autoprograms}\DNA Moiré Designer"; Filename: "{app}\DNA_Moire_Designer.exe"
Name: "{autoprograms}\cadnano2 Companion"; Filename: "{app}\cadnano\cadnano2.exe"
Name: "{autodesktop}\DNA Moiré Designer"; Filename: "{app}\DNA_Moire_Designer.exe"; Tasks: desktopicon

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Additional shortcuts:"

[Run]
Filename: "{app}\DNA_Moire_Designer.exe"; Description: "Launch DNA Moiré Designer"; Flags: nowait postinstall skipifsilent
