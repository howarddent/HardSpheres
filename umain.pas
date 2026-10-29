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

unit umain;

{******************************************************************************
  HardSpheres - an ideal (and, at higher density, a non-ideal) gas as hard
  spheres in a box, simulated on the GPU in OpenCL (ugassim.pas, gas.cl)
  and drawn with OpenGL from the same buffer.

  View: the spheres as shaded point sprites (one vertex each, a sphere
  painted in the fragment shader), coloured by speed against the rms
  speed; the box as lines and the piston (the top face, z = Lz) as a
  translucent plate; only the height changes with the volume. Drag to orbit, wheel to zoom.

  Controls:
    spheres, starting packing fraction, equal starting speeds or
      Maxwell-Boltzmann, Reset;
    temperature (sets it at once by rescaling every velocity) and
      "hold T" (a thermostat: compressions are then isothermal);
    piston position and piston speed: a slow piston compresses or expands
      the gas adiabatically, PV^(5/3) constant for an ideal gas;
    Pause.
  Readouts: T, V, packing fraction, pressure on the walls and on the
  piston, NT/V and the Carnahan-Starling pressure, PV/NT, PV^(5/3), the
  collision rate, and below them the speed histogram against the
  Maxwell-Boltzmann distribution at the current temperature.
  Far right, the P-V diagram (ugaspv.pas): the state as a red dot with
  its trail, against the ideal-gas and hard-sphere isotherms at T and,
  after a move without the thermostat, the adiabat; under it a short
  explanation of what the dot is doing.

  OpenCL-OpenGL sharing needs both on one GPU. On a hybrid-graphics laptop
  start the program with NVIDIA PRIME offload (run-nvidia.sh); otherwise
  it runs with the positions copied through the host each frame, and
  says so.

  Command line: --n=<spheres>, --phi=<packing>, --maxwell, --hold,
  --speed=<piston speed>, --compress=<f> (after 100 frames, move the
  piston to f times the box's height), --snapshot=<file.png> (runs
  --frames=<n> frames, default 200, saves the view and the window and
  prints the readouts).
*******************************************************************************}

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Math, Forms, Controls, Graphics, Dialogs, ExtCtrls,
  StdCtrls, ComCtrls, Spin, LCLType, FPImage, FPWritePNG,
  GL, GLU, GLext, OpenGLContext, ugassim, ugaspv;

type

  { TfmMain }

  TfmMain = class(TForm)
    procedure FormCreate(Sender: TObject);
    procedure FormDestroy(Sender: TObject);
  private
    FGL: TOpenGLControl;
    FPanel: TPanel;
    cbN: TComboBox;
    sePhi: TFloatSpinEdit;
    cbMaxwell: TCheckBox;
    btnReset, btnPause: TButton;
    tbTemp, tbPiston, tbSpeed: TTrackBar;
    lblTemp, lblPiston, lblSpeed: TLabel;
    cbHold: TCheckBox;
    memInfo: TMemo;
    pbHist: TPaintBox;
    lblStatus: TLabel;
    FPVPanel: TPanel;
    FPV: TGasPVDiagram;
    memPV: TMemo;
    FSim: TGasSim;
    FGLReady, FSimReady, FPaused, FUpdating: Boolean;
    FProg: GLuint;
    FVBO: GLuint;
    FVBOCount: Integer;
    FLocScale, FLocVScale: GLint;
    FYaw, FPitch, FDist: Double;
    FDragging: Boolean;
    FLastX, FLastY: Integer;
    FHist: array of Double;          // speed histogram, normalised
    FHistMax: Double;                // its speed range
    FVel: array of Single;
    FFrames: Integer;
    FFpsTime: QWord;
    FFpsFrames: Integer;
    FFps: Double;
    FStartPV, FStartTV: Double;
    FStartT, FStartV: Double;       // T and V per sphere before the move
    FError: string;
    FSnapshot: string;
    FSnapFrames: Integer;
    procedure BuildControls;
    function AddLabel(const S: string; ATop: Integer): TLabel;
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
    procedure PistonChange(Sender: TObject);
    procedure SpeedChange(Sender: TObject);
    procedure HoldChange(Sender: TObject);
    procedure HistPaint(Sender: TObject);
    procedure DrawScene;
    procedure UpdateHistogram;
    procedure UpdateInfo;
    procedure UpdatePV;
    procedure PVPanelResize(Sender: TObject);
    function SliderTemp: Double;
    function SliderPiston: Double;
    function SliderSpeed: Double;
    procedure SaveSnapshot;
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
  FovY = 35.0;
  HistBins = 60;
  NChoices: array[0..5] of Integer = (1000, 10000, 50000, 100000, 250000, 1000000);

  VertexShader =
    '#version 120'#10 +
    'uniform float uScale;'#10 +      // viewport height / tan(fov/2), times the radius 0.5
    'uniform float uVScale;'#10 +     // speed that maps to the top of the colour scale
    'varying vec3 vColor;'#10 +
    'vec3 ramp(float t) {'#10 +       // blue - cyan - green - yellow - red
    '  t = clamp(t, 0.0, 1.0);'#10 +
    '  return clamp(vec3(1.5 - abs(4.0*t - 3.0), 1.5 - abs(4.0*t - 2.0), 1.5 - abs(4.0*t - 1.0)), 0.0, 1.0);'#10 +
    '}'#10 +
    'void main() {'#10 +
    '  vec4 eye = gl_ModelViewMatrix * vec4(gl_Vertex.xyz, 1.0);'#10 +
    '  gl_Position = gl_ProjectionMatrix * eye;'#10 +
    '  gl_PointSize = max(1.0, uScale / -eye.z);'#10 +
    '  vColor = ramp(gl_Vertex.w / uVScale);'#10 +
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

const
  PVH = 520;     // the P-V diagram's height

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
  n := StrToIntDef(Application.GetOptionValue('n'), 100000);
  cbN.ItemIndex := 3;
  for k := 0 to High(NChoices) do
    if NChoices[k] = n then cbN.ItemIndex := k;
  if Application.HasOption('phi') then
    sePhi.Value := StrToFloatDef(Application.GetOptionValue('phi'), sePhi.Value);
  if Application.HasOption('maxwell') then cbMaxwell.Checked := True;
  if Application.HasOption('hold') then cbHold.Checked := True;
  if Application.HasOption('speed') then
    tbSpeed.Position := Round(100 + 66.67 * Log10(StrToFloatDef(Application.GetOptionValue('speed'), 1)));
  FSim := TGasSim.Create;
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

procedure TfmMain.BuildControls;
var
  y, k: Integer;
  function Track(ATop, AMax, APos: Integer; AOnChange: TNotifyEvent): TTrackBar;
  begin
    Result := TTrackBar.Create(Self);
    Result.Parent := FPanel;
    Result.Left := 6; Result.Top := ATop; Result.Width := FPanel.Width - 16;
    Result.Height := 30;
    Result.Min := 0; Result.Max := AMax; Result.Position := APos;
    Result.TickStyle := tsNone;
    Result.OnChange := AOnChange;
  end;
begin
  Caption := 'Hard-sphere gas (OpenCL + OpenGL)';
  Width := 1880; Height := 860;
  // the P-V diagram at the far right, created first so it sits outside
  // the controls
  FPVPanel := TPanel.Create(Self);
  FPVPanel.Parent := Self;
  FPVPanel.Align := alRight;
  FPVPanel.Width := 480;
  FPVPanel.BevelOuter := bvNone;
  FPVPanel.OnResize := @PVPanelResize;
  FPV := TGasPVDiagram.Create(Self);
  FPV.Parent := FPVPanel;
  FPV.SetBounds(0, 0, 480, PVH);
  memPV := TMemo.Create(Self);
  memPV.Parent := FPVPanel;
  memPV.SetBounds(6, PVH + 4, 468, 300);
  memPV.ReadOnly := True;
  memPV.ScrollBars := ssAutoVertical;
  memPV.WordWrap := True;
  memPV.Font.Size := 9;

  FPanel := TPanel.Create(Self);
  FPanel.Parent := Self;
  FPanel.Align := alRight;
  FPanel.Width := 400;
  FPanel.Left := FPVPanel.Left - 400;
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

  y := 8;
  AddLabel('spheres:', y + 4);
  cbN := TComboBox.Create(Self);
  cbN.Parent := FPanel;
  cbN.Left := 120; cbN.Top := y; cbN.Width := 120;
  cbN.Style := csDropDownList;
  for k := 0 to High(NChoices) do cbN.Items.Add(IntToStr(NChoices[k]));
  y := y + 32;
  AddLabel('packing fraction:', y + 4);
  sePhi := TFloatSpinEdit.Create(Self);
  sePhi.Parent := FPanel;
  sePhi.Left := 120; sePhi.Top := y; sePhi.Width := 80;
  sePhi.DecimalPlaces := 3; sePhi.Increment := 0.01;
  sePhi.MinValue := 0.001; sePhi.MaxValue := 0.45; sePhi.Value := 0.05;
  y := y + 32;
  cbMaxwell := TCheckBox.Create(Self);
  cbMaxwell.Parent := FPanel;
  cbMaxwell.Left := 10; cbMaxwell.Top := y;
  cbMaxwell.Caption := 'start Maxwell-Boltzmann (else equal speeds)';
  y := y + 28;
  btnReset := TButton.Create(Self);
  btnReset.Parent := FPanel;
  btnReset.Left := 10; btnReset.Top := y; btnReset.Width := 90;
  btnReset.Caption := 'Reset';
  btnReset.OnClick := @ResetClick;
  btnPause := TButton.Create(Self);
  btnPause.Parent := FPanel;
  btnPause.Left := 110; btnPause.Top := y; btnPause.Width := 90;
  btnPause.Caption := 'Pause';
  btnPause.OnClick := @PauseClick;
  y := y + 40;

  lblTemp := AddLabel('', y);
  tbTemp := Track(y + 16, 200, 100, @TempChange);
  y := y + 50;
  cbHold := TCheckBox.Create(Self);
  cbHold.Parent := FPanel;
  cbHold.Left := 10; cbHold.Top := y;
  cbHold.Caption := 'hold T (thermostat: isothermal)';
  cbHold.OnChange := @HoldChange;
  y := y + 28;
  lblPiston := AddLabel('', y);
  tbPiston := Track(y + 16, 1000, 500, @PistonChange);
  y := y + 50;
  lblSpeed := AddLabel('', y);
  tbSpeed := Track(y + 16, 200, 100, @SpeedChange);
  y := y + 54;

  memInfo := TMemo.Create(Self);
  memInfo.Parent := FPanel;
  memInfo.Left := 6; memInfo.Top := y; memInfo.Width := FPanel.Width - 12; memInfo.Height := 250;
  memInfo.ReadOnly := True;
  memInfo.ScrollBars := ssNone;
  memInfo.Font.Name := 'Monospace';
  memInfo.Font.Size := 9;
  y := y + 256;

  pbHist := TPaintBox.Create(Self);
  pbHist.Parent := FPanel;
  pbHist.Left := 6; pbHist.Top := y; pbHist.Width := FPanel.Width - 12;
  pbHist.Height := 200;
  pbHist.OnPaint := @HistPaint;
  y := y + 204;

  lblStatus := AddLabel('', y);
  lblStatus.WordWrap := True;
  lblStatus.AutoSize := False;
  lblStatus.Width := FPanel.Width - 16;
  lblStatus.Height := 60;
  TempChange(nil);
  SpeedChange(nil);
end;

// temperature 0.1 .. 10, log scale, 1 in the middle
function TfmMain.SliderTemp: Double;
begin
  Result := Power(10, (tbTemp.Position - 100) / 100);
end;

// piston speed 0.01 .. 10 diameters per unit time, log scale
function TfmMain.SliderSpeed: Double;
begin
  Result := Power(10, (tbSpeed.Position - 100) / 66.67);
end;

// piston (top face) position: from where the packing fraction would reach
// 0.5 to the farthest the piston may go (twice the starting height)
function TfmMain.SliderPiston: Double;
var
  lmin: Double;
begin
  if not FSimReady then Exit(0);
  lmin := FSim.Count * Pi / 6 / 0.5 / ((FSim.Lx - 1) * (FSim.Ly - 1)) + 1;
  Result := lmin + (FSim.LzMax - lmin) * tbPiston.Position / tbPiston.Max;
end;

procedure TfmMain.TempChange(Sender: TObject);
begin
  lblTemp.Caption := Format('temperature T = %.3g  (rms speed %.3g)', [SliderTemp, Sqrt(3 * SliderTemp)]);
  if FUpdating or not FSimReady then Exit;
  FSim.SetTemperature(SliderTemp);
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
  // a move from rest: PV^(5/3) is compared with its value now
  if Abs(FSim.Lz - FSim.LzTarget) < 1e-6 then
  begin
    FStartPV := FSim.Pressure * Power(FSim.Volume, 5 / 3);
    FStartTV := FSim.Temperature * Power(FSim.Volume, 2 / 3);
    FStartT := FSim.Temperature;
    FStartV := FSim.Volume / FSim.Count;
  end;
  FSim.LzTarget := SliderPiston;
end;

procedure TfmMain.SpeedChange(Sender: TObject);
begin
  // quasi-static (reversible adiabatic) needs the piston well below the
  // speed of sound, sqrt(5T/3) for a dilute gas, and the move slower
  // than a sound wave's crossing of the box
  if FSimReady then
    lblSpeed.Caption := Format('piston speed %.3g  (sound %.2f)', [SliderSpeed, Sqrt(5 / 3 * Max(FSim.Temperature, 1e-6))])
  else
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
  FLocVScale := glGetUniformLocation(FProg, 'uVScale');
  glGenBuffers(1, @FVBO);
  // the OpenCL context, sharing this GL context if the driver can
  {$IFDEF LINUX}
  FSim.InitCL(ExtractFilePath(ParamStr(0)) + 'gas.cl', glXGetCurrentContext, glXGetCurrentDisplay);
  {$ELSE}
  FSim.InitCL(ExtractFilePath(ParamStr(0)) + 'gas.cl', nil, nil);
  {$ENDIF}
  FGLReady := True;
end;

procedure TfmMain.ResetSim;
var
  n: Integer;
begin
  FSimReady := False;
  n := NChoices[cbN.ItemIndex];
  glFinish;
  // the GL buffer first, at its new size and with the starting
  // positions, then the simulation shares it
  glBindBuffer(GL_ARRAY_BUFFER, FVBO);
  glBufferData(GL_ARRAY_BUFFER, 16 * n, nil, GL_DYNAMIC_DRAW);
  FVBOCount := n;
  glFinish;
  FSim.Setup(n, sePhi.Value, SliderTemp, not cbMaxwell.Checked, FVBO);
  if not FSim.Shared then
    glBufferSubData(GL_ARRAY_BUFFER, 0, 16 * n, @FSim.HostPos[0]);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  FSim.PistonSpeed := SliderSpeed;
  FSim.Thermostat := cbHold.Checked;
  FUpdating := True;
  tbPiston.Position := Round(tbPiston.Max * (FSim.Lz - (FSim.Count * Pi / 6 / 0.5 / ((FSim.Lx - 1) * (FSim.Ly - 1)) + 1))
    / (FSim.LzMax - (FSim.Count * Pi / 6 / 0.5 / ((FSim.Lx - 1) * (FSim.Ly - 1)) + 1)));
  FUpdating := False;
  FSimReady := True;
  FSim.LzTarget := FSim.Lz;
  PistonChange(nil);
  SetLength(FVel, 4 * n);
  SetLength(FHist, 0);
  FStartPV := 0;
  FStartT := 0;
  FStartTV := 0;
  FFrames := 0;
  // the volume axis: to the farthest the piston can go, plus a margin
  FPV.VMax := (FSim.Lx - 1) * (FSim.Ly - 1) * (FSim.LzMax - 1) / FSim.Count * 1.1;
  FPV.ClearState;
  if FSim.Shared then
    lblStatus.Caption := 'OpenCL on ' + FSim.DeviceName + ', sharing the GL buffer'
  else
    lblStatus.Caption := 'OpenCL on ' + FSim.DeviceName + '; GL is on ' + string(PChar(glGetString(GL_RENDERER)))
      + ', so positions are copied each frame (run-nvidia.sh shares them)';
end;

procedure TfmMain.AppIdle(Sender: TObject; var Done: Boolean);
begin
  FGL.Invalidate;
  Done := False;
end;

procedure TfmMain.GLPaint(Sender: TObject);
var
  t: QWord;
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
      if (FFrames = 100) and Application.HasOption('compress') then
      begin
        FStartPV := FSim.Pressure * Power(FSim.Volume, 5 / 3);
        FStartTV := FSim.Temperature * Power(FSim.Volume, 2 / 3);
        FStartT := FSim.Temperature;
        FStartV := FSim.Volume / FSim.Count;
        FSim.LzTarget := FSim.Lz * StrToFloatDef(Application.GetOptionValue('compress'), 0.5);
      end;
      if FFrames mod 10 = 1 then
      begin
        UpdateHistogram;
        UpdatePV;
      end;
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
      UpdateHistogram;
      UpdateInfo;
      UpdatePV;
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
  L, cx, cy, cz, x0, y0, z0, x1, y1, lz: Double;
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
  glRotatef(-90, 1, 0, 0);          // z up
  cx := FSim.Lx / 2; cy := FSim.Ly / 2; cz := L / 2;
  glTranslatef(-cx, -cy, -cz);
  glEnable(GL_DEPTH_TEST);

  // spheres
  glUseProgram(FProg);
  glUniform1f(FLocScale, h / Tan(DegToRad(FovY / 2)) * 0.5);
  glUniform1f(FLocVScale, 2.2 * Sqrt(3 * Max(FSim.Temperature, 1e-6)));
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
  glEnable(GL_POINT_SPRITE);
  glBindBuffer(GL_ARRAY_BUFFER, FVBO);
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(4, GL_FLOAT, 16, nil);
  glDrawArrays(GL_POINTS, 0, FVBOCount);
  glDisableClientState(GL_VERTEX_ARRAY);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glUseProgram(0);
  glDisable(GL_POINT_SPRITE);
  glDisable(GL_VERTEX_PROGRAM_POINT_SIZE);

  // the box up to the piston (the top face)
  x0 := 0; y0 := 0; z0 := 0; x1 := FSim.Lx; y1 := FSim.Ly;
  lz := FSim.Lz;
  glLineWidth(1.5);
  glColor3f(0.85, 0.85, 0.9);
  glBegin(GL_LINES);
  // the floor
  glVertex3f(x0, y0, z0); glVertex3f(x1, y0, z0);
  glVertex3f(x1, y0, z0); glVertex3f(x1, y1, z0);
  glVertex3f(x1, y1, z0); glVertex3f(x0, y1, z0);
  glVertex3f(x0, y1, z0); glVertex3f(x0, y0, z0);
  // vertical edges up to the piston
  glVertex3f(x0, y0, z0); glVertex3f(x0, y0, lz);
  glVertex3f(x1, y0, z0); glVertex3f(x1, y0, lz);
  glVertex3f(x1, y1, z0); glVertex3f(x1, y1, lz);
  glVertex3f(x0, y1, z0); glVertex3f(x0, y1, lz);
  glEnd;
  // the piston: a translucent plate with a bright rim
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
  FHistMax := 4.0 * Sqrt(Max(FSim.Temperature, 1e-6));    // 2.3 v_rms
  SetLength(FHist, HistBins);
  for b := 0 to HistBins - 1 do FHist[b] := 0;
  dv := FHistMax / HistBins;
  for i := 0 to n - 1 do
  begin
    v := Sqrt(Sqr(FVel[4 * i]) + Sqr(FVel[4 * i + 1]) + Sqr(FVel[4 * i + 2]));
    b := Trunc(v / dv);
    if b < HistBins then FHist[b] := FHist[b] + 1;
  end;
  for b := 0 to HistBins - 1 do FHist[b] := FHist[b] / (n * dv);   // a density
  pbHist.Invalidate;
end;

procedure TfmMain.HistPaint(Sender: TObject);
var
  C: TCanvas;
  w, h, b, x0, x1, yb, k: Integer;
  T, fmax, v, f, sx: Double;
  function MB(v: Double): Double;   // Maxwell-Boltzmann speed density
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
  C.TextOut(4, 2, 'speed distribution; red: Maxwell-Boltzmann at T');
  if (Length(FHist) = 0) or not FSimReady then Exit;
  T := Max(FSim.Temperature, 1e-6);
  fmax := MB(Sqrt(2 * T)) * 1.25;
  for b := 0 to HistBins - 1 do fmax := Max(fmax, FHist[b] * 1.05);
  yb := h - 18;
  sx := (w - 10) / FHistMax;
  C.Brush.Color := RGBToColor(120, 160, 220);
  C.Pen.Color := RGBToColor(70, 100, 160);
  for b := 0 to HistBins - 1 do
  begin
    x0 := 5 + Round(b * FHistMax / HistBins * sx);
    x1 := 5 + Round((b + 1) * FHistMax / HistBins * sx);
    C.Rectangle(x0, yb - Round(FHist[b] / fmax * (yb - 18)), x1, yb + 1);
  end;
  C.Pen.Color := clRed;
  C.Pen.Width := 2;
  for k := 0 to 200 do
  begin
    v := FHistMax * k / 200;
    f := MB(v);
    if k = 0 then C.MoveTo(5 + Round(v * sx), yb - Round(f / fmax * (yb - 18)))
    else C.LineTo(5 + Round(v * sx), yb - Round(f / fmax * (yb - 18)));
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
    C.TextOut(Max(0, 5 + Round(v * sx) - 8), yb + 4, Format('%.1f', [k * 0.5]));
  end;
  C.TextOut(w - 70, 16, 'speed / v_rms');
  C.Brush.Style := bsSolid;
end;

procedure TfmMain.UpdateInfo;
var
  L: TStrings;
  T, V, P, phi, Pid, Pcs: Double;
begin
  if not FSimReady then Exit;
  T := FSim.Temperature; V := FSim.Volume; P := FSim.Pressure;
  phi := FSim.PackingFraction;
  Pid := FSim.Count * T / V;
  Pcs := Pid * CarnahanStarling(phi);
  L := memInfo.Lines;
  L.BeginUpdate;
  try
    L.Clear;
    L.Add(Format('N %d spheres, diameter 1, mass 1', [FSim.Count]));
    L.Add(Format('t %.1f   %.0f frames/s, %d substeps/frame', [FSim.Time, FFps, FSim.Substeps]));
    L.Add(Format('T  %.4f   (target %.3g%s)', [T, FSim.TargetT, BoolToStr(FSim.Thermostat, ', held', '')]));
    L.Add(Format('V  %.4g  (box %.0f x %.0f x %.1f)', [V, FSim.Lx, FSim.Ly, FSim.Lz]));
    L.Add(Format('packing fraction  %.4f', [phi]));
    L.Add('');
    L.Add(Format('pressure, walls   %.5f', [P]));
    L.Add(Format('pressure, piston  %.5f', [FSim.PistonPressure]));
    L.Add(Format('ideal gas  NT/V   %.5f', [Pid]));
    L.Add(Format('Carnahan-Starling %.5f', [Pcs]));
    L.Add(Format('PV/NT %.4f  (CS %.4f)', [P / Pid, CarnahanStarling(phi)]));
    if FStartPV > 0 then
    begin
      // both constant in a slow compression of an ideal gas
      L.Add(Format('PV^(5/3) %.4g  (%.3f x before the move)', [P * Power(V, 5 / 3), P * Power(V, 5 / 3) / FStartPV]));
      L.Add(Format('TV^(2/3) %.4g  (%.3f x before the move)', [T * Power(V, 2 / 3), T * Power(V, 2 / 3) / FStartTV]));
    end
    else
    begin
      L.Add(Format('PV^(5/3) %.4g', [P * Power(V, 5 / 3)]));
      L.Add(Format('TV^(2/3) %.4g', [T * Power(V, 2 / 3)]));
    end;
    L.Add(Format('collisions %.4g per sphere per unit time', [2 * FSim.CollisionRate]));
    L.Add(Format('  mean free path %.3g diameters', [Sqrt(8 * T / Pi) / Max(2 * FSim.CollisionRate, 1e-9)]));
    if FSim.Overflows > 0 then
      L.Add(Format('WARNING: %d cell overflows (too dense)', [FSim.Overflows]));
  finally
    L.EndUpdate;
  end;
end;

procedure TfmMain.PVPanelResize(Sender: TObject);
begin
  FPV.SetBounds(0, 0, FPVPanel.ClientWidth, PVH);
  memPV.SetBounds(6, PVH + 4, FPVPanel.ClientWidth - 12, FPVPanel.ClientHeight - PVH - 10);
end;

// the hard-sphere (Carnahan-Starling) excess free energy over NkT
function HSExcess(phi: Double): Double;
begin
  Result := (4 * phi - 3 * phi * phi) / Sqr(1 - phi);
end;

// the state on the P-V diagram, per sphere, and what it is doing
procedure TfmMain.UpdatePV;
var
  L: TStrings;
  T, v, P, Z, Zcs, T0: Double;
  moved: Boolean;
begin
  // the pressure is smoothed over ~20 frames: wait for it to settle
  if (not FSimReady) or (FFrames < 3 * FSim.PressureFrames) then Exit;
  T := FSim.Temperature;
  v := FSim.Volume / FSim.Count;
  P := FSim.Pressure;
  if (T <= 0) or (P <= 0) then Exit;
  moved := FStartPV > 0;
  if moved and not FSim.Thermostat then T0 := FStartT else T0 := 0;
  FPV.SetState(T, v, P, T0, FStartV);
  Z := P * v / T;
  Zcs := CarnahanStarling(FSim.PackingFraction);
  L := memPV.Lines;
  L.BeginUpdate;
  try
    L.Clear;
    L.Add(Format('Volume per sphere %.3g, pressure %.4g, T %.3g: PV/NkT = %.3f '
      + '(1 for an ideal gas; %.3f for hard spheres at this density).', [v, P, T, Z, Zcs]));
    L.Add('');
    if Zcs < 1.1 then
      L.Add('The spheres take up little of the box, so the gas is nearly ideal: '
        + 'the dot sits on the dashed ideal-gas isotherm, P = NkT/V.')
    else
      L.Add(Format('The spheres now fill %.0f %% of the box. They have no attraction, '
        + 'but their own volume is no longer negligible, so the pressure is above '
        + 'the ideal gas''s (the solid hard-sphere isotherm): less room to move in, '
        + 'more collisions with the walls.', [100 * FSim.PackingFraction]));
    L.Add('');
    if FSim.Thermostat then
      L.Add('Temperature held (isothermal): move the piston and the dot runs along '
        + 'the isotherm. Halve the volume and the pressure doubles (Boyle''s law).')
    else
      L.Add('Temperature not held: a slow piston compresses the gas adiabatically. '
        + 'The work done heats it, so the dot climbs across the isotherms along the '
        + 'adiabat, steeper than Boyle''s hyperbola. For an ideal gas that is '
        + 'PV^(5/3) constant (dotted); hard spheres follow the steeper solid one, '
        + 'since their own volume makes each step of compression cost more work. '
        + 'Tick "hold T" for an isothermal compression.');
    if moved and not FSim.Thermostat then
      L.Add(Format('Since the move began T has risen %.3f times: hard spheres %.3f, '
        + 'an ideal gas %.3f.',
        [T / FStartT, Power(FStartV / v, 2 / 3) * Exp(2 / 3 * (HSExcess(FSim.PackingFraction)
          - HSExcess(Pi / 6 / FStartV))), Power(FStartV / v, 2 / 3)]));
    L.Add('');
    L.Add('However far the piston goes there is no flat part: with no attraction '
      + 'between the spheres the gas never condenses, unlike a vapour '
      + '(compare the N2O program).');
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

// the GL view read back (before the swap, as TVMPlot2D does) beside the
// panel's readouts in the memo printed to stdout
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
  WriteLn(memPV.Lines.Text);
end;

end.
