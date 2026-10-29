{
  This file is part of HardSpheres.
  Copyright (C) 2026 Howard Dent
  SPDX-License-Identifier: GPL-3.0-or-later

  HardSpheres is free software: you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your
  option) any later version. It is distributed WITHOUT ANY WARRANTY; see
  the file LICENSE at the top of the repository for details.
}

unit un2omain;

{******************************************************************************
  N2OPhase - the sticky-sphere liquid and vapour of ../Sticky beside the
  pressure-volume isotherm of nitrous oxide, as a teaching aid (anaesthesia:
  isotherms and the pressure of a saturated vapour).

  Left, the simulation as in StickySpheres (stepped Lennard-Jones
  spheres in a box with gravity, an attractive floor and a piston lid;
  see ../Sticky/ustickymain.pas for the controls). Right, the N2O
  pressure-volume isotherm at the current temperature (un2opv.pas), with
  a red dot for the simulation's state and a trail of where it has been;
  N2O's saturation data and vapour pressure come from un2ophase.pas
  (whose pressure-temperature diagram control is no longer shown: the
  user found it didn't help). The dot's molar volume is the lever rule
  between N2O's saturated volumes while liquid is present, otherwise the
  simulation's compressibility factor PV/NkT is kept.

  The state is mapped onto N2O by corresponding states:
    - temperature: T_N2O = (T/eps) x Tc(N2O)/Tc(model), so the model's
      critical temperature lands on N2O's 36.4 C;
    - pressure: while liquid and vapour coexist (3 % of the spheres or
      more are in droplets of 32 or more), the system is a
      saturated vapour and its pressure is N2O's vapour pressure at that
      temperature: the dot sits on the flat line. Otherwise (all vapour,
      or above the critical temperature) the measured wall pressure is
      scaled by the model's own saturated vapour pressure at that
      temperature (unsaturated vapour), or by the model's critical
      pressure above the critical temperature.
  The model's critical temperature and vapour-pressure curve were
  measured from the simulation (README).

  Command line: as StickySpheres (--n, --phi, --eps, --g, --floor,
  --temp, --nohold, --maxwell, --snapshot, --frames), and
  --compress=<f>: after 100 frames move the lid to f times its height;
  --speed=<piston speed>.
*******************************************************************************}

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Math, Forms, Controls, Graphics, Dialogs, ExtCtrls,
  StdCtrls, ComCtrls, Spin, LCLType, FPImage, FPWritePNG,
  GL, GLU, GLext, OpenGLContext, ustickysim, un2ophase, un2opv;

type

  { TfmMain }

  TfmMain = class(TForm)
    procedure FormCreate(Sender: TObject);
    procedure FormDestroy(Sender: TObject);
  private
    FGL: TOpenGLControl;
    FPanel: TPanel;
    FPhasePanel: TPanel;
    FPV: TN2OPVDiagram;
    memPhase: TMemo;
    cbN, cbColour, cbPot: TComboBox;
    sePhi, seLambda: TFloatSpinEdit;
    cbMaxwell, cbHold, cbShowBonds: TCheckBox;
    btnReset, btnPause: TButton;
    tbTemp, tbEps, tbGrav, tbFloor, tbPiston, tbSpeed: TTrackBar;
    lblTemp, lblEps, lblGrav, lblFloor, lblPiston, lblSpeed: TLabel;
    memInfo: TMemo;
    pbHist: TPaintBox;
    lblStatus: TLabel;
    FSim: TStickySim;
    FGLReady, FSimReady, FPaused, FUpdating: Boolean;
    FProg: GLuint;
    FVBO, FEBO: GLuint;
    FVBOCount, FBondLines: Integer;
    FLocScale, FLocMode: GLint;
    FYaw, FPitch, FDist: Double;
    FDragging: Boolean;
    FLastX, FLastY: Integer;
    FHist: array of Double;
    FHistMax: Double;
    FVel: array of Single;
    FNbr, FShells: array of LongInt;
    FPosCopy: array of Single;
    FOvershoot: Double;
    FLines: array of GLuint;
    FClusters: TClusterStats;
    FHaveClusters: Boolean;
    FFrames: Integer;
    FFpsTime: QWord;
    FFpsFrames: Integer;
    FFps: Double;
    FError: string;
    FSnapshot: string;
    FSnapFrames: Integer;
    procedure BuildControls;
    function AddLabel(const S: string; ATop: Integer): TLabel;
    function AddCheck(const S: string; ATop: Integer; AChecked: Boolean): TCheckBox;
    procedure InitGL;
    procedure ResetSim;
    procedure GLPaint(Sender: TObject);
    procedure GLMouseDown(Sender: TObject; Button: TMouseButton; Shift: TShiftState; X, Y: Integer);
    procedure GLMouseMove(Sender: TObject; Shift: TShiftState; X, Y: Integer);
    procedure GLMouseUp(Sender: TObject; Button: TMouseButton; Shift: TShiftState; X, Y: Integer);
    procedure GLMouseWheel(Sender: TObject; Shift: TShiftState; WheelDelta: Integer;
      MousePos: TPoint; var Handled: Boolean);
    procedure AppIdle(Sender: TObject; var Done: Boolean);
    procedure ResetClick(Sender: TObject);
    procedure PauseClick(Sender: TObject);
    procedure TempChange(Sender: TObject);
    procedure EpsChange(Sender: TObject);
    procedure GravChange(Sender: TObject);
    procedure FloorChange(Sender: TObject);
    procedure PistonChange(Sender: TObject);
    procedure SpeedChange(Sender: TObject);
    procedure HoldChange(Sender: TObject);
    procedure HistPaint(Sender: TObject);
    procedure DrawScene;
    procedure UpdateBonds(Clusters: Boolean);
    procedure UpdateHistogram;
    procedure UpdateInfo;
    function SliderTemp: Double;
    function SliderEps: Double;
    function SliderGrav: Double;
    function SliderFloor: Double;
    function SliderPiston: Double;
    function SliderSpeed: Double;
    function MinLz: Double;
    procedure SaveSnapshot;
    procedure UpdatePhase;
    function FloorPressure: Double;
    procedure PhasePanelResize(Sender: TObject);
  end;

var
  fmMain: TfmMain;

implementation

{$R *.lfm}

{$IFDEF LINUX}
function glXGetCurrentContext: Pointer; cdecl; external 'GL';
function glXGetCurrentDisplay: Pointer; cdecl; external 'GL';
{$ENDIF}

const
  // the stepped Lennard-Jones model's own liquid-vapour coexistence, in
  // reduced units (T/eps, P sigma^3/eps), measured from the simulation:
  // critical temperature and the vapour pressure ln P = A - B/T + C ln T (README)
  ModelTc = 1.035;
  ModelPsatA = 0.550;
  ModelPsatB = 3.193;
  ModelPsatC = 3.363;
  // ... and its critical density (spheres per sigma^3), for volumes
  ModelRhoC = 0.302;
  // liquid present: this fraction of the spheres in droplets (below)
  LiquidFraction = 0.03;
  // spheres in clusters of 2^LiquidBin or more count as liquid, the rest
  // as vapour (for the lever rule)
  LiquidBin = 5;
  PVH = 700;

function ModelVapourPressure(Tr: Double): Double;
begin
  Tr := Max(Tr, 0.01);
  Result := Exp(ModelPsatA - ModelPsatB / Tr + ModelPsatC * Ln(Tr));
end;

function ModelPc: Double;
begin
  Result := ModelVapourPressure(ModelTc);
end;

// the simulation's temperature (T/eps) as N2O's, kelvin
function AsN2OKelvin(Tr: Double): Double;
begin
  Result := Tr / ModelTc * N2OTc;
end;

// the kernel file: beside the program, else in ../Sticky (shared)
function KernelFile: string;
begin
  Result := ExtractFilePath(ParamStr(0)) + 'sticky.cl';
  if not FileExists(Result) then
    Result := ExtractFilePath(ParamStr(0)) + '..' + PathDelim + 'Sticky' + PathDelim + 'sticky.cl';
end;

const
  FovY = 35.0;
  HistBins = 50;
  NChoices: array[0..5] of Integer = (1000, 2000, 5000, 10000, 20000, 50000);

  // pos.w = bonds + min(speed / (2.2 v_rms), 0.999)
  VertexShader =
    '#version 120'#10 +
    'uniform float uScale;'#10 +
    'uniform float uMode;'#10 +       // 0: by bonds, 1: by speed
    'varying vec3 vColor;'#10 +
    'vec3 ramp(float t) {'#10 +
    '  t = clamp(t, 0.0, 1.0);'#10 +
    '  return clamp(vec3(1.5 - abs(4.0*t - 3.0), 1.5 - abs(4.0*t - 2.0), 1.5 - abs(4.0*t - 1.0)), 0.0, 1.0);'#10 +
    '}'#10 +
    'void main() {'#10 +
    '  vec4 eye = gl_ModelViewMatrix * vec4(gl_Vertex.xyz, 1.0);'#10 +
    '  gl_Position = gl_ProjectionMatrix * eye;'#10 +
    '  gl_PointSize = max(1.0, uScale / -eye.z);'#10 +
    '  float nb = floor(gl_Vertex.w);'#10 +
    '  if (uMode < 0.5)'#10 +
    '    vColor = nb < 0.5 ? vec3(0.42, 0.50, 0.68) : ramp(0.35 + 0.65 * min((nb - 1.0) / 7.0, 1.0));'#10 +
    '  else'#10 +
    '    vColor = ramp(fract(gl_Vertex.w) / 0.999);'#10 +
    '}'#10;

  FragmentShader =
    '#version 120'#10 +
    'varying vec3 vColor;'#10 +
    'void main() {'#10 +
    '  vec2 c = gl_PointCoord * 2.0 - 1.0;'#10 +
    '  float r2 = dot(c, c);'#10 +
    '  if (r2 > 1.0) discard;'#10 +
    '  vec3 n = vec3(c.x, -c.y, sqrt(1.0 - r2));'#10 +
    '  vec3 l = normalize(vec3(-0.4, 0.5, 0.75));'#10 +
    '  float d = max(dot(n, l), 0.0);'#10 +
    '  float s = pow(max(dot(reflect(-l, n), vec3(0.0, 0.0, 1.0)), 0.0), 24.0);'#10 +
    '  gl_FragColor = vec4(vColor * (0.25 + 0.75 * d) + vec3(0.5 * s), 1.0);'#10 +
    '}'#10;

{ TfmMain }

procedure TfmMain.FormCreate(Sender: TObject);
var
  k, n: Integer;
begin
  DefaultFormatSettings.DecimalSeparator := '.';
  FYaw := -35; FPitch := 22; FDist := 2.9;
  FSnapFrames := StrToIntDef(Application.GetOptionValue('frames'), 200);
  FSnapshot := Application.GetOptionValue('snapshot');
  BuildControls;
  n := StrToIntDef(Application.GetOptionValue('n'), 5000);
  cbN.ItemIndex := 2;
  for k := 0 to High(NChoices) do
    if NChoices[k] = n then cbN.ItemIndex := k;
  if Application.HasOption('phi') then
    sePhi.Value := StrToFloatDef(Application.GetOptionValue('phi'), sePhi.Value);
  if Application.HasOption('lambda') then
  begin
    seLambda.Value := StrToFloatDef(Application.GetOptionValue('lambda'), seLambda.Value);
    cbPot.ItemIndex := 1;
  end;
  if Application.HasOption('eps') then
    tbEps.Position := Round(50 * StrToFloatDef(Application.GetOptionValue('eps'), 1));
  if Application.HasOption('g') then
    tbGrav.Position := Round(1000 * StrToFloatDef(Application.GetOptionValue('g'), 0.05));
  if Application.HasOption('floor') then
    tbFloor.Position := Round(20 * StrToFloatDef(Application.GetOptionValue('floor'), 1.5));
  if Application.HasOption('temp') then
    tbTemp.Position := Round(100 + 100 * Log10(StrToFloatDef(Application.GetOptionValue('temp'), 1)));
  if Application.HasOption('maxwell') then cbMaxwell.Checked := True;
  if Application.HasOption('speed') then
    tbSpeed.Position := Round(100 + 66.67 * Log10(StrToFloatDef(Application.GetOptionValue('speed'), 0.3)));
  if Application.HasOption('nohold') then cbHold.Checked := False;
  FSim := TStickySim.Create;
  FSim.AllowIntegrated := Application.HasOption('allow-igpu');
  Application.OnIdle := @AppIdle;
end;

procedure TfmMain.FormDestroy(Sender: TObject);
begin
  Application.OnIdle := nil;
  FreeAndNil(FSim);
end;

function TfmMain.AddLabel(const S: string; ATop: Integer): TLabel;
begin
  Result := TLabel.Create(Self);
  Result.Parent := FPanel;
  Result.Left := 10;
  Result.Top := ATop;
  Result.Caption := S;
end;

function TfmMain.AddCheck(const S: string; ATop: Integer; AChecked: Boolean): TCheckBox;
begin
  Result := TCheckBox.Create(Self);
  Result.Parent := FPanel;
  Result.Left := 10;
  Result.Top := ATop;
  Result.Caption := S;
  Result.Checked := AChecked;
end;

procedure TfmMain.BuildControls;
var
  y, k: Integer;
  function Track(ATop, AMax, APos: Integer; AOnChange: TNotifyEvent): TTrackBar;
  begin
    Result := TTrackBar.Create(Self);
    Result.Parent := FPanel;
    Result.Left := 6; Result.Top := ATop; Result.Width := FPanel.Width - 16;
    Result.Height := 28;
    Result.Min := 0; Result.Max := AMax; Result.Position := APos;
    Result.TickStyle := tsNone;
    Result.OnChange := AOnChange;
  end;
  function Spin(ATop: Integer; const ACaption: string; ADec: Integer; AInc, AMin, AMax, AVal: Double): TFloatSpinEdit;
  begin
    AddLabel(ACaption, ATop + 4);
    Result := TFloatSpinEdit.Create(Self);
    Result.Parent := FPanel;
    Result.Left := 150; Result.Top := ATop; Result.Width := 80;
    Result.DecimalPlaces := ADec; Result.Increment := AInc;
    Result.MinValue := AMin; Result.MaxValue := AMax; Result.Value := AVal;
  end;
begin
  Caption := 'Nitrous oxide: liquid, vapour and the P-V isotherm (OpenCL + OpenGL)';
  Width := 1930; Height := 1180;
  FPhasePanel := TPanel.Create(Self);
  FPhasePanel.Parent := Self;
  FPhasePanel.Align := alRight;
  FPhasePanel.Width := 480;
  FPhasePanel.BevelOuter := bvNone;
  FPhasePanel.OnResize := @PhasePanelResize;
  FPV := TN2OPVDiagram.Create(Self);
  FPV.Parent := FPhasePanel;
  FPV.SetBounds(0, 0, 480, PVH);
  memPhase := TMemo.Create(Self);
  memPhase.Parent := FPhasePanel;
  memPhase.SetBounds(6, PVH + 4, 468, 300);
  memPhase.ReadOnly := True;
  memPhase.ScrollBars := ssAutoVertical;
  memPhase.WordWrap := True;
  memPhase.Font.Size := 9;
  FPanel := TPanel.Create(Self);
  FPanel.Parent := Self;
  FPanel.Align := alRight;
  FPanel.Width := 430;
  FPanel.BevelOuter := bvNone;

  FGL := TOpenGLControl.Create(Self);
  FGL.Parent := Self;
  FGL.Align := alClient;
  FGL.DepthBits := 24;
  FGL.OnPaint := @GLPaint;
  FGL.OnMouseDown := @GLMouseDown;
  FGL.OnMouseMove := @GLMouseMove;
  FGL.OnMouseUp := @GLMouseUp;
  FGL.OnMouseWheel := @GLMouseWheel;

  y := 6;
  AddLabel('spheres:', y + 4);
  cbN := TComboBox.Create(Self);
  cbN.Parent := FPanel;
  cbN.Left := 150; cbN.Top := y; cbN.Width := 100;
  cbN.Style := csDropDownList;
  for k := 0 to High(NChoices) do cbN.Items.Add(IntToStr(NChoices[k]));
  btnReset := TButton.Create(Self);
  btnReset.Parent := FPanel;
  btnReset.Left := 262; btnReset.Top := y; btnReset.Width := 75;
  btnReset.Caption := 'Reset';
  btnReset.OnClick := @ResetClick;
  btnPause := TButton.Create(Self);
  btnPause.Parent := FPanel;
  btnPause.Left := 343; btnPause.Top := y; btnPause.Width := 75;
  btnPause.Caption := 'Pause';
  btnPause.OnClick := @PauseClick;
  y := y + 32;
  sePhi := Spin(y, 'packing fraction:', 3, 0.001, 0.001, 0.45, 0.003);   // a dilute, unsaturated vapour
  y := y + 30;
  AddLabel('potential:', y + 4);
  cbPot := TComboBox.Create(Self);
  cbPot.Parent := FPanel;
  cbPot.Left := 80; cbPot.Top := y; cbPot.Width := 200;
  cbPot.Style := csDropDownList;
  cbPot.Items.Add('stepped Lennard-Jones');
  cbPot.Items.Add('square well, width:');
  cbPot.ItemIndex := 0;
  seLambda := TFloatSpinEdit.Create(Self);
  seLambda.Parent := FPanel;
  seLambda.Left := 290; seLambda.Top := y; seLambda.Width := 75;
  seLambda.DecimalPlaces := 2; seLambda.Increment := 0.05;
  seLambda.MinValue := 1.05; seLambda.MaxValue := 2.0; seLambda.Value := 1.5;
  y := y + 30;
  cbMaxwell := AddCheck('start Maxwell-Boltzmann (else equal speeds)', y, False);
  y := y + 30;

  lblTemp := AddLabel('', y);
  tbTemp := Track(y + 16, 200, 90, @TempChange);    // T 0.79: about -37 °C as N2O
  y := y + 46;
  cbHold := AddCheck('hold T (thermostat)', y, True);
  cbHold.OnChange := @HoldChange;
  y := y + 26;
  lblEps := AddLabel('', y);
  tbEps := Track(y + 16, 200, 50, @EpsChange);
  y := y + 46;
  lblGrav := AddLabel('', y);
  tbGrav := Track(y + 16, 200, 5, @GravChange);      // weak: the vapour fills the box
  y := y + 46;
  lblFloor := AddLabel('', y);
  tbFloor := Track(y + 16, 80, 30, @FloorChange);     // 1.5: the liquid wets the floor
  y := y + 46;
  lblPiston := AddLabel('', y);
  tbPiston := Track(y + 16, 1000, 500, @PistonChange);
  y := y + 46;
  lblSpeed := AddLabel('', y);
  tbSpeed := Track(y + 16, 200, 65, @SpeedChange);   // 0.3: slow enough to stay near equilibrium
  y := y + 48;
  AddLabel('colour by:', y + 4);
  cbColour := TComboBox.Create(Self);
  cbColour.Parent := FPanel;
  cbColour.Left := 80; cbColour.Top := y; cbColour.Width := 110;
  cbColour.Style := csDropDownList;
  cbColour.Items.Add('bonds');
  cbColour.Items.Add('speed');
  cbColour.ItemIndex := 0;
  cbShowBonds := TCheckBox.Create(Self);
  cbShowBonds.Parent := FPanel;
  cbShowBonds.Left := 210; cbShowBonds.Top := y + 2;
  cbShowBonds.Caption := 'draw bonds';
  cbShowBonds.Checked := True;
  y := y + 34;

  memInfo := TMemo.Create(Self);
  memInfo.Parent := FPanel;
  memInfo.Left := 6; memInfo.Top := y; memInfo.Width := FPanel.Width - 12; memInfo.Height := 278;
  memInfo.ReadOnly := True;
  memInfo.ScrollBars := ssNone;
  memInfo.Font.Name := 'Monospace';
  memInfo.Font.Size := 8;
  y := y + 282;

  pbHist := TPaintBox.Create(Self);
  pbHist.Parent := FPanel;
  pbHist.Left := 6; pbHist.Top := y; pbHist.Width := FPanel.Width - 12;
  pbHist.Height := 230;
  pbHist.OnPaint := @HistPaint;
  y := y + 232;

  lblStatus := AddLabel('', y);
  lblStatus.WordWrap := True;
  lblStatus.AutoSize := False;
  lblStatus.Width := FPanel.Width - 16;
  lblStatus.Height := 34;
  lblStatus.Font.Size := 8;
  TempChange(nil);
  EpsChange(nil);
  GravChange(nil);
  FloorChange(nil);
  SpeedChange(nil);
end;

// temperature 0.1 .. 10, log scale
function TfmMain.SliderTemp: Double;
begin
  Result := Power(10, (tbTemp.Position - 100) / 100);
end;

// binding energy 0 .. 4
function TfmMain.SliderEps: Double;
begin
  Result := tbEps.Position / 50;
end;

// gravity 0 .. 0.2 diameters per unit time^2, along -z
function TfmMain.SliderFloor: Double;
begin
  Result := tbFloor.Position / 20;
end;

// the floor's attraction, a well one diameter deep above it: from about
// 1.5 eps the liquid wets the floor and spreads over it as a film
procedure TfmMain.FloorChange(Sender: TObject);
begin
  if SliderFloor > 0 then
    lblFloor.Caption := Format('floor attraction eps_f = %.2f  (eps_f/eps %.2f; wets from ~1.5)',
      [SliderFloor, SliderFloor / Max(SliderEps, 1e-9)])
  else
    lblFloor.Caption := 'floor attraction eps_f = 0  (hard floor, not wetted)';
  if FSim <> nil then FSim.FloorEps := SliderFloor;
end;

function TfmMain.SliderGrav: Double;
begin
  Result := tbGrav.Position / 1000;
end;

function TfmMain.SliderSpeed: Double;
begin
  Result := Power(10, (tbSpeed.Position - 100) / 66.67);
end;

// the piston's (the top face's) travel: from packing fraction 0.5 to its farthest
function TfmMain.MinLz: Double;
begin
  Result := FSim.Count * Pi / 6 / 0.5 / ((FSim.Lx - 1) * (FSim.Ly - 1)) + 1;
end;

function TfmMain.SliderPiston: Double;
begin
  if not FSimReady then Exit(0);
  Result := MinLz + (FSim.LzMax - MinLz) * tbPiston.Position / tbPiston.Max;
end;

procedure TfmMain.TempChange(Sender: TObject);
begin
  if (cbPot <> nil) and (cbPot.ItemIndex = 0) and (SliderEps > 0) then
    lblTemp.Caption := Format('temperature T = %.3g  (as N2O: %.0f °C)',
      [SliderTemp, AsN2OKelvin(SliderTemp / SliderEps) - 273.15])
  else
    lblTemp.Caption := Format('temperature T = %.3g  (rms speed %.3g)', [SliderTemp, Sqrt(3 * SliderTemp)]);
  if lblGrav <> nil then GravChange(nil);
  if lblEps <> nil then EpsChange(nil);
  if FUpdating or not FSimReady then Exit;
  FSim.SetTemperature(SliderTemp);
end;

procedure TfmMain.EpsChange(Sender: TObject);
begin
  lblEps.Caption := Format('well depth eps = %.2f   (eps/T = %.2f)', [SliderEps, SliderEps / SliderTemp]);
  if FSim <> nil then FSim.Eps := SliderEps;
  if (lblFloor <> nil) and (tbFloor <> nil) then FloorChange(nil);
end;

// the scale height T/g of a dilute gas, against the box's height
procedure TfmMain.GravChange(Sender: TObject);
begin
  if SliderGrav > 0 then
    lblGrav.Caption := Format('gravity g = %.3f  (gas scale height T/g %.0f)', [SliderGrav, SliderTemp / SliderGrav])
  else
    lblGrav.Caption := 'gravity g = 0';
  if FSim <> nil then FSim.Gravity := SliderGrav;
end;

procedure TfmMain.HoldChange(Sender: TObject);
begin
  if FSimReady then
  begin
    FSim.Thermostat := cbHold.Checked;
    FSim.TargetT := SliderTemp;
  end;
end;

procedure TfmMain.PistonChange(Sender: TObject);
begin
  if not FSimReady then Exit;
  lblPiston.Caption := Format('piston (top) at Lz = %.1f  (now %.1f; floor %.0f x %.0f)',
    [SliderPiston, FSim.Lz, FSim.Lx, FSim.Ly]);
  if FUpdating or (Sender = nil) then Exit;
  FSim.LzTarget := SliderPiston;
end;

procedure TfmMain.SpeedChange(Sender: TObject);
begin
  lblSpeed.Caption := Format('piston speed %.3g', [SliderSpeed]);
  if FSimReady then FSim.PistonSpeed := SliderSpeed;
end;

procedure TfmMain.ResetClick(Sender: TObject);
begin
  if FGLReady then
  begin
    FGL.MakeCurrent;
    ResetSim;
  end;
end;

procedure TfmMain.PauseClick(Sender: TObject);
begin
  FPaused := not FPaused;
  if FPaused then btnPause.Caption := 'Run' else btnPause.Caption := 'Pause';
end;

function CompileShader(Kind: GLenum; const Src: string): GLuint;
var
  P: PGLchar;
  ok, n: GLint;
  Log: string;
begin
  Result := glCreateShader(Kind);
  P := PGLchar(Src);
  glShaderSource(Result, 1, @P, nil);
  glCompileShader(Result);
  glGetShaderiv(Result, GL_COMPILE_STATUS, @ok);
  if ok = 0 then
  begin
    SetLength(Log, 4096);
    n := 0;
    glGetShaderInfoLog(Result, 4096, @n, PGLchar(Log));
    SetLength(Log, n);
    raise Exception.Create('shader: ' + Log);
  end;
end;

procedure TfmMain.InitGL;
var
  ok: GLint;
begin
  if not (Load_GL_version_1_5 and Load_GL_version_2_0) then
    raise Exception.Create('OpenGL 2.0 is needed');
  FProg := glCreateProgram();
  glAttachShader(FProg, CompileShader(GL_VERTEX_SHADER, VertexShader));
  glAttachShader(FProg, CompileShader(GL_FRAGMENT_SHADER, FragmentShader));
  glLinkProgram(FProg);
  glGetProgramiv(FProg, GL_LINK_STATUS, @ok);
  if ok = 0 then raise Exception.Create('shader program did not link');
  FLocScale := glGetUniformLocation(FProg, 'uScale');
  FLocMode := glGetUniformLocation(FProg, 'uMode');
  glGenBuffers(1, @FVBO);
  glGenBuffers(1, @FEBO);
  {$IFDEF LINUX}
  FSim.InitCL(KernelFile, glXGetCurrentContext, glXGetCurrentDisplay);
  {$ELSE}
  FSim.InitCL(KernelFile, nil, nil);
  {$ENDIF}
  FGLReady := True;
end;

procedure TfmMain.ResetSim;
var
  n: Integer;
begin
  FSimReady := False;
  if FPV <> nil then FPV.ClearState;
  n := NChoices[cbN.ItemIndex];
  glFinish;
  glBindBuffer(GL_ARRAY_BUFFER, FVBO);
  glBufferData(GL_ARRAY_BUFFER, 16 * n, nil, GL_DYNAMIC_DRAW);
  FVBOCount := n;
  glFinish;
  FSim.Eps := SliderEps;
  FSim.Gravity := SliderGrav;
  FSim.FloorEps := SliderFloor;
  if cbPot.ItemIndex = 1 then
    FSim.Setup(n, sePhi.Value, SliderTemp, SquareWellPotential(seLambda.Value), not cbMaxwell.Checked, FVBO)
  else
    FSim.Setup(n, sePhi.Value, SliderTemp, LennardJonesPotential, not cbMaxwell.Checked, FVBO);
  if not FSim.Shared then
    glBufferSubData(GL_ARRAY_BUFFER, 0, 16 * n, @FSim.HostPos[0]);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  FSim.PistonSpeed := SliderSpeed;
  FSim.Thermostat := cbHold.Checked;
  FUpdating := True;
  tbPiston.Position := Round(tbPiston.Max * (FSim.Lz - MinLz) / (FSim.LzMax - MinLz));
  FUpdating := False;
  FSimReady := True;
  FSim.LzTarget := FSim.Lz;
  PistonChange(nil);
  SetLength(FVel, 4 * n);
  SetLength(FNbr, n * MaxNeighbours);
  SetLength(FShells, n * MaxNeighbours);
  SetLength(FPosCopy, 4 * n);
  SetLength(FLines, n * MaxNeighbours);
  FBondLines := 0;
  SetLength(FHist, 0);
  FHaveClusters := False;
  FFrames := 0;
  if FSim.Shared then
    lblStatus.Caption := 'OpenCL on ' + FSim.DeviceName + ', sharing the GL buffer'
  else
    lblStatus.Caption := 'OpenCL on ' + FSim.DeviceName + '; GL on ' + string(PChar(glGetString(GL_RENDERER)))
      + ', positions copied each frame (run-nvidia.sh shares them)';
end;

// the bond lists from the GPU: the line indices to draw and, when asked,
// the clusters
procedure TfmMain.UpdateBonds(Clusters: Boolean);
var
  i, k, j, m: Integer;
begin
  FSim.ReadLists(FNbr, FShells);
  m := 0;
  for i := 0 to FSim.Count - 1 do
    for k := 0 to MaxNeighbours - 1 do
    begin
      j := FNbr[i * MaxNeighbours + k];
      if (j > i) and (FShells[i * MaxNeighbours + k] < FSim.Potential.BondEdge) then
      begin
        FLines[m] := i; FLines[m + 1] := j;
        Inc(m, 2);
      end;
    end;
  FBondLines := m div 2;
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, FEBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, Max(1, m) * 4, @FLines[0], GL_STREAM_DRAW);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  if Clusters then
  begin
    FClusters := FSim.ClusterStats(FNbr, FShells);
    FSim.ReadPositions(FPosCopy);
    FOvershoot := FSim.MaxOvershoot(FNbr, FShells, FPosCopy);
    FHaveClusters := True;
  end;
end;

procedure TfmMain.AppIdle(Sender: TObject; var Done: Boolean);
begin
  FGL.Invalidate;
  Done := False;
end;

procedure TfmMain.GLPaint(Sender: TObject);
var
  t: QWord;
  due: Boolean;
begin
  if not FGL.MakeCurrent then Exit;
  if FError <> '' then
  begin
    glClearColor(0.3, 0.05, 0.05, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    FGL.SwapBuffers;
    Exit;
  end;
  try
    if not FGLReady then InitGL;
    if not FSimReady then ResetSim;
    if not FPaused then
    begin
      glFinish;
      FSim.Frame;
      if not FSim.Shared then
      begin
        glBindBuffer(GL_ARRAY_BUFFER, FVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, 16 * FSim.Count, @FSim.HostPos[0]);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
      end;
      Inc(FFrames);
      // --compress=<f>: after 100 frames the lid moves to f times its height
      if (FFrames = 100) and Application.HasOption('compress') then
        FSim.LzTarget := Max(MinLz, FSim.Lz * StrToFloatDef(Application.GetOptionValue('compress'), 0.5));
      due := FFrames mod 10 = 1;
      if cbShowBonds.Checked or due then UpdateBonds(due);
      if due then UpdateHistogram;
    end;
    DrawScene;
    FGL.SwapBuffers;
    Inc(FFpsFrames);
    t := GetTickCount64;
    if t - FFpsTime >= 500 then
    begin
      FFps := FFpsFrames * 1000 / (t - FFpsTime);
      FFpsFrames := 0;
      FFpsTime := t;
      UpdateInfo;
      PistonChange(nil);
    end;
    if (FSnapshot <> '') and (FFrames >= FSnapFrames) then
    begin
      UpdateBonds(True);
      UpdateHistogram;
      UpdateInfo;
      SaveSnapshot;
      FSnapshot := '';
      Application.Terminate;
    end;
  except
    on E: Exception do
    begin
      FError := E.Message;
      lblStatus.Caption := 'ERROR: ' + FError;
      WriteLn('ERROR: ', FError);
      if FSnapshot <> '' then Application.Terminate;
    end;
  end;
end;

procedure TfmMain.DrawScene;
var
  w, h: Integer;
  L, x0, y0, z0, x1, y1, lz: Double;
begin
  w := FGL.Width; h := Max(1, FGL.Height);
  glViewport(0, 0, w, h);
  glClearColor(0.08, 0.08, 0.11, 1);
  glClear(GL_COLOR_BUFFER_BIT or GL_DEPTH_BUFFER_BIT);
  if not FSimReady then Exit;
  // framed on the floor (Lx = Ly = the starting height), not the current
  // height, so the floor and walls stay still on screen as the top moves
  L := Max(FSim.Lx, FSim.Ly);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity;
  gluPerspective(FovY, w / h, 0.05 * L, 20 * L);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity;
  glTranslatef(0, 0, -FDist * L);
  glRotatef(FPitch, 1, 0, 0);
  glRotatef(FYaw, 0, 1, 0);
  glRotatef(-90, 1, 0, 0);
  glTranslatef(-FSim.Lx / 2, -FSim.Ly / 2, -L / 2);
  glEnable(GL_DEPTH_TEST);

  glBindBuffer(GL_ARRAY_BUFFER, FVBO);
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 16, nil);
  // bonds, centre to centre (mostly hidden inside the spheres, so they
  // show the links between spheres of a cluster)
  if cbShowBonds.Checked and (FBondLines > 0) then
  begin
    glLineWidth(1.5);
    glColor3f(0.95, 0.85, 0.55);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, FEBO);
    glDrawElements(GL_LINES, 2 * FBondLines, GL_UNSIGNED_INT, nil);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  end;
  // spheres
  glUseProgram(FProg);
  glUniform1f(FLocScale, h / Tan(DegToRad(FovY / 2)) * 0.5);
  glUniform1f(FLocMode, cbColour.ItemIndex);
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
  glEnable(GL_POINT_SPRITE);
  glVertexPointer(4, GL_FLOAT, 16, nil);
  glDrawArrays(GL_POINTS, 0, FVBOCount);
  glDisableClientState(GL_VERTEX_ARRAY);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glUseProgram(0);
  glDisable(GL_POINT_SPRITE);
  glDisable(GL_VERTEX_PROGRAM_POINT_SIZE);

  // the floor, tinted blue as deep as its attraction
  if FSim.FloorEps > 0 then
  begin
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);
    glColor4f(0.3, 0.55, 1.0, Min(0.45, 0.12 * FSim.FloorEps));
    glBegin(GL_QUADS);
    glVertex3f(0, 0, 0); glVertex3f(FSim.Lx, 0, 0); glVertex3f(FSim.Lx, FSim.Ly, 0); glVertex3f(0, FSim.Ly, 0);
    glEnd;
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
  end;
  // the box: floor and side walls fixed, the top face (z = Lz) the piston
  x0 := 0; y0 := 0; z0 := 0; x1 := FSim.Lx; y1 := FSim.Ly;
  lz := FSim.Lz;
  glLineWidth(1.5);
  glColor3f(0.85, 0.85, 0.9);
  glBegin(GL_LINES);
  glVertex3f(x0, y0, z0); glVertex3f(x1, y0, z0);
  glVertex3f(x1, y0, z0); glVertex3f(x1, y1, z0);
  glVertex3f(x1, y1, z0); glVertex3f(x0, y1, z0);
  glVertex3f(x0, y1, z0); glVertex3f(x0, y0, z0);
  glVertex3f(x0, y0, z0); glVertex3f(x0, y0, lz);
  glVertex3f(x1, y0, z0); glVertex3f(x1, y0, lz);
  glVertex3f(x1, y1, z0); glVertex3f(x1, y1, lz);
  glVertex3f(x0, y1, z0); glVertex3f(x0, y1, lz);
  glEnd;
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDepthMask(GL_FALSE);
  glColor4f(0.9, 0.6, 0.2, 0.25);
  glBegin(GL_QUADS);
  glVertex3f(x0, y0, lz); glVertex3f(x1, y0, lz); glVertex3f(x1, y1, lz); glVertex3f(x0, y1, lz);
  glEnd;
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
  glColor3f(1.0, 0.7, 0.3);
  glLineWidth(2.5);
  glBegin(GL_LINE_LOOP);
  glVertex3f(x0, y0, lz); glVertex3f(x1, y0, lz); glVertex3f(x1, y1, lz); glVertex3f(x0, y1, lz);
  glEnd;
  glLineWidth(1);
end;

procedure TfmMain.UpdateHistogram;
var
  i, b, n: Integer;
  v, dv: Double;
begin
  if not FSimReady then Exit;
  FSim.ReadVelocities(FVel);
  n := FSim.Count;
  FHistMax := 4.0 * Sqrt(Max(FSim.Temperature, 1e-6));
  SetLength(FHist, HistBins);
  for b := 0 to HistBins - 1 do FHist[b] := 0;
  dv := FHistMax / HistBins;
  for i := 0 to n - 1 do
  begin
    v := Sqrt(Sqr(FVel[4 * i]) + Sqr(FVel[4 * i + 1]) + Sqr(FVel[4 * i + 2]));
    b := Trunc(v / dv);
    if b < HistBins then FHist[b] := FHist[b] + 1;
  end;
  for b := 0 to HistBins - 1 do FHist[b] := FHist[b] / (n * dv);
  pbHist.Invalidate;
end;

// top: the fraction of spheres in clusters of each size (powers of two);
// bottom: the speed distribution against Maxwell-Boltzmann
procedure TfmMain.HistPaint(Sender: TObject);
var
  C: TCanvas;
  w, h, b, x0, x1, yb, yt, k, nb, bw: Integer;
  T, fmax, v, f, sx: Double;
  s: string;
  function MB(v: Double): Double;
  begin
    Result := Sqrt(2 / Pi) * v * v / Power(T, 1.5) * Exp(-v * v / (2 * T));
  end;
begin
  C := pbHist.Canvas;
  w := pbHist.Width; h := pbHist.Height;
  C.Brush.Color := clWhite;
  C.FillRect(0, 0, w, h);
  C.Font.Size := 8;
  C.Font.Color := clBlack;
  if not FSimReady then Exit;

  // clusters
  C.TextOut(4, 2, 'fraction of spheres in clusters of size s');
  yt := 18; yb := 100;
  if FHaveClusters then
  begin
    nb := Max(1, FClusters.NBins);
    bw := (w - 10) div nb;
    C.Pen.Color := RGBToColor(150, 90, 40);
    for b := 0 to nb - 1 do
    begin
      x0 := 5 + b * bw;
      if b = 0 then C.Brush.Color := RGBToColor(110, 128, 175)
      else C.Brush.Color := RGBToColor(235, Max(60, 200 - 14 * b), 60);
      C.Rectangle(x0 + 1, yb - Round(FClusters.Bins[b] * (yb - yt)), x0 + bw - 1, yb + 1);
      if bw >= 14 then
      begin
        C.Brush.Style := bsClear;
        if (1 shl b) >= 1000 then s := Format('%dk', [(1 shl b) div 1000]) else s := IntToStr(1 shl b);
        C.TextOut(x0 + (bw - C.TextWidth(s)) div 2, yb + 2, s);
        C.Brush.Style := bsSolid;
      end;
    end;
    C.Pen.Color := clBlack;
    C.Line(5, yb, w - 5, yb);
  end;

  // speeds
  C.Brush.Color := clWhite;
  C.TextOut(4, 122, 'speed distribution; red: Maxwell-Boltzmann at T');
  if Length(FHist) = 0 then Exit;
  T := Max(FSim.Temperature, 1e-6);
  fmax := MB(Sqrt(2 * T)) * 1.25;
  for b := 0 to HistBins - 1 do fmax := Max(fmax, FHist[b] * 1.05);
  yt := 138; yb := h - 16;
  sx := (w - 10) / FHistMax;
  C.Brush.Color := RGBToColor(120, 160, 220);
  C.Pen.Color := RGBToColor(70, 100, 160);
  for b := 0 to HistBins - 1 do
  begin
    x0 := 5 + Round(b * FHistMax / HistBins * sx);
    x1 := 5 + Round((b + 1) * FHistMax / HistBins * sx);
    C.Rectangle(x0, yb - Round(FHist[b] / fmax * (yb - yt)), x1, yb + 1);
  end;
  C.Pen.Color := clRed;
  C.Pen.Width := 2;
  for k := 0 to 200 do
  begin
    v := FHistMax * k / 200;
    f := MB(v);
    if k = 0 then C.MoveTo(5 + Round(v * sx), yb - Round(f / fmax * (yb - yt)))
    else C.LineTo(5 + Round(v * sx), yb - Round(f / fmax * (yb - yt)));
  end;
  C.Pen.Width := 1;
  C.Pen.Color := clBlack;
  C.Line(5, yb, w - 5, yb);
  C.Brush.Style := bsClear;
  for k := 0 to 4 do
  begin
    v := Sqrt(3 * T) * k * 0.5;
    if v > FHistMax then Break;
    C.Line(5 + Round(v * sx), yb, 5 + Round(v * sx), yb + 4);
    C.TextOut(Max(0, 5 + Round(v * sx) - 8), yb + 3, Format('%.1f', [k * 0.5]));
  end;
  C.TextOut(w - 75, yt, 'speed / v_rms');
  C.Brush.Style := bsSolid;
end;

procedure TfmMain.UpdateInfo;
var
  L: TStrings;
  T, V, P, rho, n: Double;
begin
  if not FSimReady then Exit;
  n := FSim.Count;
  T := FSim.Temperature; V := FSim.Volume; P := FSim.Pressure;
  rho := n / V;
  L := memInfo.Lines;
  L.BeginUpdate;
  try
    L.Clear;
    L.Add(Format('N %d   t %.1f   %.0f frames/s, %d substeps', [FSim.Count, FSim.Time, FFps, FSim.Substeps]));
    L.Add(Format('T %.4f (target %.3g%s)   phi %.4f', [T, FSim.TargetT, BoolToStr(FSim.Thermostat, ', held', ''),
      FSim.PackingFraction]));
    L.Add(Format('%s, eps %.2f   eps/T %.2f', [FSim.Potential.Name, FSim.Eps, FSim.Eps / Max(T, 1e-9)]));
    L.Add(Format('event sweeps per substep %.1f%s', [FSim.MeanSweeps,
      BoolToStr(FSim.CappedSweeps > 0, Format(' (%d capped)', [FSim.CappedSweeps]), '')]));
    L.Add('');
    L.Add(Format('bonds (r < %.2f) %d, %.3f per sphere', [FSim.Potential.R[FSim.Potential.BondEdge], FSim.Bonds, 2 * FSim.Bonds / n]));
    L.Add(Format('energy per sphere: KE %.4f, gravity %.4f, pairs %.4f',
      [FSim.KineticEnergy / n, FSim.PotentialEnergy / n, FSim.PairEnergy / n]));
    if FSim.FloorEps > 0 then
      L.Add(Format('  floor %.4f (%.0f spheres in its well)', [FSim.FloorEnergy / n, -FSim.FloorEnergy / FSim.FloorEps]));
    L.Add(Format('  total %.4f,', [FSim.Energy / n]));
    L.Add(Format('  thermostat added %.4f', [FSim.ThermostatWork / n]));
    L.Add('per sphere per unit time:');
    L.Add(Format('  core collisions %.4g', [2 * FSim.CoreRate]));
    L.Add(Format('  bonds formed %.4g', [2 * FSim.FormRate]));
    L.Add(Format('  bonds broken %.4g', [2 * FSim.BreakRate]));
    if FHaveClusters then
    begin
      L.Add(Format('clusters (2+) %d, free spheres %.1f %%', [FClusters.Clusters, 100 * FClusters.Monomers / n]));
      L.Add(Format('  largest %d, mean seen by a sphere %.1f', [FClusters.Largest, FClusters.MeanSize]));
      L.Add(Format('  farthest a pair is outside its shell %.3f%s', [FOvershoot,
        BoolToStr(FOvershoot > 0.2, '  EVENTS MISSED', '')]));
    end;
    L.Add('');
    L.Add(Format('pressure, walls  %.5f   piston %.5f', [P, FSim.PistonPressure]));
    if FSim.Gravity > 0 then
      L.Add('  (in gravity P varies with height)')
    else
      L.Add(Format('PV/NT %.4f   1 + B2 rho %.4f', [P / (rho * T), 1 + SecondVirial(FSim.Potential, FSim.Eps, T) * rho]));
    if FSim.Overflows > 0 then
      L.Add(Format('WARNING: %d cell overflows (too dense)', [FSim.Overflows]));
  finally
    L.EndUpdate;
  end;
  UpdatePhase;
end;

procedure TfmMain.PhasePanelResize(Sender: TObject);
begin
  FPV.SetBounds(0, 0, FPhasePanel.ClientWidth, PVH);
  memPhase.SetBounds(6, PVH + 4, FPhasePanel.ClientWidth - 12, FPhasePanel.ClientHeight - PVH - 10);
end;

// the vapour's pressure at the floor, from the walls' mean. In gravity an
// isothermal vapour thins with height, P = P0 exp(-z/H), H = T/g, so the
// walls' mean is P0 [A (1 + e^-a) + p H (1 - e^-a)] / (2A + p h), with h
// the height open to the centres, a = h/H, A the floor's area and p its
// perimeter. The floor is where it first saturates. (The lid alone gives
// the same, but few spheres reach it: too noisy for a dilute vapour.)
function TfmMain.FloorPressure: Double;
var
  area, perim, h, H0, ah, f: Double;
begin
  Result := FSim.Pressure;
  if (FSim.Gravity <= 0) or (FSim.Temperature <= 0) then Exit;
  area := (FSim.Lx - 1) * (FSim.Ly - 1);
  perim := 2 * ((FSim.Lx - 1) + (FSim.Ly - 1));
  h := FSim.Lz - 1;
  H0 := FSim.Temperature / FSim.Gravity;
  ah := h / H0;
  f := (area * (1 + Exp(-ah)) + perim * H0 * (1 - Exp(-ah))) / (2 * area + perim * h);
  Result := FSim.Pressure / f;
end;

// the red dot: the simulation's state as N2O (see the header), and what
// it means, in the memo under the diagram
procedure TfmMain.UpdatePhase;
var
  L: TStrings;
  Tr, TK, TC, P, Pm, Psat, ratio, Z, V, Vl, Vv, xl: Double;
  liquid: Boolean;
  what: string;
  k: Integer;
begin
  if not FSimReady then Exit;
  L := memPhase.Lines;
  L.BeginUpdate;
  try
    L.Clear;
    if (cbPot.ItemIndex <> 0) or (FSim.Eps <= 0) then
    begin
      FPV.ClearState;
      L.Add('The mapping onto nitrous oxide is calibrated for the stepped '
        + 'Lennard-Jones potential with a well depth eps above 0.');
      Exit;
    end;
    Tr := FSim.Temperature / FSim.Eps;
    TK := AsN2OKelvin(Tr);
    TC := TK - 273.15;
    // the vapour's pressure: in gravity it varies with height, and where
    // it is highest, at the floor, is where it first saturates. The floor
    // carries the lid's push plus the weight of everything (exact in
    // equilibrium); without gravity, the walls' mean
    Pm := FloorPressure / FSim.Eps;
    if Pm <= 0 then Exit;     // not measured yet
    // liquid: droplets of 2^LiquidBin spheres or more holding LiquidFraction
    // of them. Waiting for one big cluster left the dot well above the
    // curve while many droplets were already growing
    xl := 0;
    if FHaveClusters then
      for k := LiquidBin to FClusters.NBins - 1 do xl := xl + FClusters.Bins[k];
    xl := EnsureRange(xl, 0, 1);
    liquid := xl >= LiquidFraction;
    Psat := N2OCoexistencePressure(TK);
    if Tr >= ModelTc then
    begin
      P := Pm / ModelPc * N2OPc;
      what := N2OPhaseName(TK, P) + ', above the critical temperature';
    end
    else if liquid then
    begin
      P := Psat;
      if TK < N2OTt then what := 'solid with its vapour (sublimation curve)'
      else what := 'saturated vapour over liquid';
    end
    else
    begin
      ratio := Pm / ModelVapourPressure(Tr);
      P := Psat * ratio;
      if ratio > 1 then what := 'supersaturated vapour (no liquid yet)'
      else what := 'unsaturated vapour';
    end;

    // the molar volume, for the P-V diagram. With liquid present, the
    // lever rule between N2O's saturated liquid and vapour volumes, by the
    // fraction of spheres in the liquid (clusters of 2^LiquidBin or more):
    // the dot slides along the flat part of the isotherm as vapour
    // condenses. Otherwise the simulation's compressibility factor
    // Z = PV/NkT kept, so an ideal vapour sits on the ideal-gas isotherm
    // and any departure is the model's own
    if liquid and (Tr < ModelTc) then
    begin
      if TK >= N2OTt then
        N2OSaturatedVolumes(TK, Vl, Vv)
      else
      begin
        N2OSaturatedVolumes(N2OTt, Vl, Vv);
        Vv := RBarLitre * TK / P;
      end;
      V := xl * Vl + (1 - xl) * Vv;
    end
    else
    begin
      Z := FSim.Pressure * FSim.Volume / (FSim.Count * Max(FSim.Temperature, 1e-9));
      V := Z * RBarLitre * TK / Max(P, 1e-9);
    end;
    FPV.SetState(TK, V, P, '');

    L.Add(Format('As nitrous oxide: %.1f °C, %.3g bar, %.3g L/mol - %s.', [TC, P, V, what]));
    if liquid and (Tr < ModelTc) then
      L.Add(Format('Liquid %.0f %% of the molecules, vapour %.0f %%.', [100 * xl, 100 * (1 - xl)]))
    else
      L.Add(Format('PV/RT = %.3f (1 for an ideal gas: the dashed line on the P-V diagram).', [Z]));
    L.Add('');
    if Tr >= ModelTc then
    begin
      L.Add('Above the critical temperature (36.4 °C for N2O) there is no '
        + 'flat part to the isotherm (P-V diagram): '
        + 'liquid surface and no vapour pressure: however hard the piston '
        + 'squeezes, the gas only becomes denser. There is no distinct liquid '
        + 'to form, so N2O above 36.4 °C is a gas, not a vapour.');
    end
    else if liquid then
    begin
      if TK < N2OTt then
        L.Add(Format('Below the triple point (-90.8 °C) the condensed phase of N2O '
          + 'would be solid; its vapour pressure is %.3g bar.', [Psat]))
      else
      begin
        L.Add(Format('Liquid and vapour coexist, so the vapour is saturated: its '
          + 'pressure is the vapour pressure at this temperature, %.3g bar, '
          + 'whatever the volume. Lower the piston and vapour condenses; '
          + 'raise it and liquid evaporates - the pressure stays the same '
          + 'while any liquid is left.', [Psat]));
        L.Add('');
        L.Add('This is why an N2O cylinder reads 50.5 bar at 20 °C until it is '
          + 'about 90 % empty (only vapour left): the gauge shows the vapour '
          + 'pressure, not the contents. The cylinder is weighed instead.');
      end;
      L.Add('');
      L.Add('');
      L.Add('On the P-V diagram the dot runs along the flat line: the overall '
        + 'volume shrinks as vapour condenses, at constant pressure. At the '
        + 'right-hand end the vapour has just become saturated; at the left '
        + 'it is all liquid, and the isotherm turns almost vertical because '
        + 'a liquid barely compresses.');
      L.Add('');
      L.Add('Raise the temperature and the flat line rises and shortens; '
        + 'at the critical point (the top of the dome) the liquid and vapour '
        + 'become one.');
    end
    else
    begin
      L.Add(Format('No liquid: the vapour is unsaturated and behaves nearly as '
        + 'an ideal gas. Lower the piston (constant temperature) and the '
        + 'pressure rises, PV ~ constant, until it reaches the vapour '
        + 'pressure, %.3g bar at this temperature. Then liquid forms and the '
        + 'pressure stops rising: the flat part of the isotherm.', [Psat]));
      L.Add('');
      L.Add('On the P-V diagram the dot follows the dashed ideal-gas curve '
        + '(Boyle''s law), dropping a little below it as the vapour nears '
        + 'saturation and the molecules'' attraction starts to tell.');
      if ratio > 1 then
      begin
        L.Add('');
        L.Add('Just now the vapour is above its vapour pressure but has not yet '
          + 'condensed: droplets need time to nucleate (a supersaturated, '
          + 'metastable vapour). On the P-V diagram it has run on past the '
          + 'dome; it drops onto the flat line once droplets form.');
      end;
    end;
    L.Add('');
    L.Add(Format('Simulation: T/eps %.3f, pressure at the floor %.4f (walls'' mean %.4f; reduced units); '
      + 'model critical temperature %.3f. Temperatures map by T/Tc, and '
      + 'while liquid is present the pressure is N2O''s vapour pressure at '
      + 'that temperature.', [Tr, Pm, FSim.Pressure / FSim.Eps, ModelTc]));
  finally
    L.EndUpdate;
  end;
end;

procedure TfmMain.GLMouseDown(Sender: TObject; Button: TMouseButton; Shift: TShiftState; X, Y: Integer);
begin
  FDragging := True;
  FLastX := X; FLastY := Y;
end;

procedure TfmMain.GLMouseMove(Sender: TObject; Shift: TShiftState; X, Y: Integer);
begin
  if not FDragging then Exit;
  FYaw := FYaw + (X - FLastX) * 0.4;
  FPitch := EnsureRange(FPitch + (Y - FLastY) * 0.4, -89, 89);
  FLastX := X; FLastY := Y;
  if FPaused then FGL.Invalidate;
end;

procedure TfmMain.GLMouseUp(Sender: TObject; Button: TMouseButton; Shift: TShiftState; X, Y: Integer);
begin
  FDragging := False;
end;

procedure TfmMain.GLMouseWheel(Sender: TObject; Shift: TShiftState; WheelDelta: Integer;
  MousePos: TPoint; var Handled: Boolean);
begin
  FDist := EnsureRange(FDist * Power(0.9, WheelDelta / 120), 0.3, 8);
  Handled := True;
end;

procedure TfmMain.SaveSnapshot;
var
  w, h, x, y: Integer;
  Buf: array of Byte;
  Img: TFPMemoryImage;
  Png: TFPWriterPNG;
  c: TFPColor;
begin
  FGL.MakeCurrent;
  DrawScene;
  w := FGL.Width; h := FGL.Height;
  SetLength(Buf, w * h * 4);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, @Buf[0]);
  Img := TFPMemoryImage.Create(w, h);
  Png := TFPWriterPNG.Create;
  try
    for y := 0 to h - 1 do
      for x := 0 to w - 1 do
      begin
        c.red := Buf[4 * ((h - 1 - y) * w + x)] * 257;
        c.green := Buf[4 * ((h - 1 - y) * w + x) + 1] * 257;
        c.blue := Buf[4 * ((h - 1 - y) * w + x) + 2] * 257;
        c.alpha := $FFFF;
        Img.Colors[x, y] := c;
      end;
    Img.SaveToFile(FSnapshot, Png);
  finally
    Png.Free;
    Img.Free;
  end;
  with GetFormImage do
  try
    SaveToFile(ChangeFileExt(FSnapshot, '_window.bmp'));
  finally
    Free;
  end;
  WriteLn(lblStatus.Caption);
  WriteLn(memInfo.Lines.Text);
  WriteLn(memPhase.Lines.Text);
end;

end.
