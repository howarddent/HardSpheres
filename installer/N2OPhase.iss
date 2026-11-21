; This file is part of HardSpheres.
; Copyright (C) 2026 Howard Dent
; SPDX-License-Identifier: GPL-3.0-or-later
;
; Inno Setup script for the N2OPhase Windows installer. Build the program
; first (lazbuild N2O\N2OPhase.lpi), then:
;
;   iscc /DAppVersion=1.0 installer\N2OPhase.iss
;
; which writes installer\out\N2OPhase-<version>-windows-x64-setup.exe. The
; GitHub workflow .github/workflows/windows-installers.yml does both.
;
; N2OPhase is one executable: the OpenCL kernels (Sticky/sticky.cl) and the
; reference data are compiled in. OpenGL comes with Windows and OpenCL with
; the graphics driver, so nothing else is installed. It installs per user by
; default (no administrator rights needed); the first page offers an
; install for all users instead.

#ifndef AppVersion
  #define AppVersion "0.0.0"
#endif
#define AppName "N2OPhase"
#define AppExe "N2OPhase.exe"

[Setup]
AppId={{7BFCB18D-992F-412E-80BC-536E43F899B4}
AppName={#AppName}
AppVersion={#AppVersion}
AppVerName={#AppName} {#AppVersion}
AppPublisher=Howard Dent
AppPublisherURL=https://github.com/howarddent/HardSpheres
AppSupportURL=https://github.com/howarddent/HardSpheres/issues
DefaultDirName={autopf}\{#AppName}
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesAllowed=x64
ArchitecturesInstallIn64BitMode=x64
LicenseFile=..\LICENSE
OutputDir=out
OutputBaseFilename={#AppName}-{#AppVersion}-windows-x64-setup
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
UninstallDisplayIcon={app}\{#AppExe}

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "..\N2O\{#AppExe}"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\LICENSE"; DestDir: "{app}"; DestName: "LICENSE.txt"

[Icons]
Name: "{autoprograms}\{#AppName}"; Filename: "{app}\{#AppExe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\{#AppExe}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#AppExe}"; Description: "{cm:LaunchProgram,{#AppName}}"; Flags: nowait postinstall skipifsilent
