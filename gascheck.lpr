program gascheck;

{******************************************************************************
  Headless checks of the hard-sphere gas (ugassim.pas) against theory:
    1. energy: T stays put with no thermostat and a still piston;
    2. thermalisation: equal speeds relax to Maxwell-Boltzmann, whose
       <v^4>/<v^2>^2 is 5/3 (equal speeds: 1);
    3. pressure: PV/NT against Carnahan-Starling, at packing 0.05 and 0.3;
    4. collision rate: pairs per sphere per unit time against Enskog,
       2 n sqrt(pi T) g(1), g(1) = (1 - phi/2)/(1 - phi)^3;
    5. adiabatic compression of a dilute gas to half its height: T rises
       as V^(-2/3).
  Build: fpc -O2 gascheck.lpr      Run: ./gascheck [N] [--allow-igpu]
*******************************************************************************}

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX} cthreads, {$ENDIF}
  SysUtils, Math, ugassim;

var
  N: Integer;
  AllowIGPU: Boolean;
  i: Integer;

function Kernel: string;
begin
  Result := ExtractFilePath(ParamStr(0)) + 'gas.cl';
end;

function SpeedRatio(S: TGasSim): Double;
var
  V: array of Single;
  i: Integer;
  v2, s2, s4: Double;
begin
  SetLength(V, 4 * S.Count);
  S.ReadVelocities(V);
  s2 := 0; s4 := 0;
  for i := 0 to S.Count - 1 do
  begin
    v2 := Sqr(V[4 * i]) + Sqr(V[4 * i + 1]) + Sqr(V[4 * i + 2]);
    s2 := s2 + v2; s4 := s4 + v2 * v2;
  end;
  s2 := s2 / S.Count; s4 := s4 / S.Count;
  Result := s4 / (s2 * s2);
end;

// spheres whose centres lie more than Margin inside the box, and that
// region's volume: the bulk, away from the layering at the walls
procedure Bulk(S: TGasSim; Margin: Double; out Count: Integer; out Vol: Double);
var
  i: Integer;
  x, y, z: Single;
begin
  Count := 0;
  for i := 0 to S.Count - 1 do
  begin
    x := S.HostPos[4 * i]; y := S.HostPos[4 * i + 1]; z := S.HostPos[4 * i + 2];
    if (x > Margin) and (x < S.Lx - Margin) and (y > Margin) and (y < S.Ly - Margin)
       and (z > Margin) and (z < S.Lz - Margin) then Inc(Count);
  end;
  Vol := (S.Lx - 2 * Margin) * (S.Ly - 2 * Margin) * (S.Lz - 2 * Margin);
end;

procedure Equilibrium(Phi: Double; Frames: Integer);
var
  S: TGasSim;
  f: Integer;
  P, T0, Tm, Z, Zcs, g1, nd, rate, k0, sumP, sumR, sumN, vb, nb, phib: Double;
  cnt, nin: Integer;
  w0: QWord;
begin
  S := TGasSim.Create;
  try
    S.AllowIntegrated := AllowIGPU;
    S.InitCL(Kernel, nil, nil);
    S.Setup(N, Phi, 1.0, True, 0);
    if Phi = 0.05 then WriteLn('device: ', S.DeviceName);
    k0 := SpeedRatio(S);
    S.Frame;
    T0 := S.Temperature;
    sumP := 0; sumR := 0; sumN := 0; cnt := 0;
    w0 := GetTickCount64;
    for f := 2 to Frames do
    begin
      S.Frame;
      if f > Frames div 2 then
      begin
        sumP := sumP + S.Pressure; sumR := sumR + S.CollisionRate; Inc(cnt);
        Bulk(S, 4, nin, vb);
        sumN := sumN + nin;
      end;
    end;
    Tm := S.Temperature;
    P := sumP / cnt;
    rate := sumR / cnt;
    // against the bulk density (4 diameters in from the walls), which the
    // wall layering leaves below N/V in a dense gas
    nb := sumN / cnt / vb;
    phib := nb * Pi / 6;
    Z := P / (nb * Tm);
    Zcs := CarnahanStarling(phib);
    nd := nb;
    g1 := (1 - phib / 2) / Power(1 - phib, 3);
    WriteLn(Format('phi %.3f (bulk %.4f), N %d, box %.1f, cells %.2f holding %d, %d frames (t = %.1f, %d substeps/frame), %.2f s',
      [S.PackingFraction, phib, S.Count, S.Lx, S.CellSize, S.CellCapacity, Frames, S.Time, S.Substeps, (GetTickCount64 - w0) / 1000]));
    WriteLn(Format('  T %.5f -> %.5f (drift %.2e)', [T0, Tm, Tm / T0 - 1]));
    WriteLn(Format('  <v^4>/<v^2>^2 %.3f at start -> %.3f (Maxwell 5/3 = 1.667)', [k0, SpeedRatio(S)]));
    WriteLn(Format('  Z = P/(n_bulk T) %.4f, Carnahan-Starling %.4f (%.1f %%)', [Z, Zcs, 100 * (Z / Zcs - 1)]));
    WriteLn(Format('  collisions %.4f per sphere per unit time, Enskog %.4f (%.1f %%)',
      [rate, 2 * nd * Sqrt(Pi * Tm) * g1, 100 * (rate / (2 * nd * Sqrt(Pi * Tm) * g1) - 1)]));
    WriteLn(Format('  cell overflows %d', [S.Overflows]));
  finally
    S.Free;
  end;
end;

procedure Adiabatic;
var
  S: TGasSim;
  V0, T0, V1, T1: Double;
  f: Integer;
begin
  S := TGasSim.Create;
  try
    S.AllowIntegrated := AllowIGPU;
    S.InitCL(Kernel, nil, nil);
    S.Setup(N, 0.01, 1.0, False, 0);
    for f := 1 to 50 do S.Frame;
    V0 := S.Volume; T0 := S.Temperature;
    S.PistonSpeed := 0.05;               // slow against v_rms = 1.7
    S.LzTarget := S.Lz / 2;
    f := 0;
    while (S.Lz > S.LzTarget + 1e-6) and (f < 100000) do
    begin
      S.Frame;
      Inc(f);
    end;
    for f := 1 to 20 do S.Frame;
    V1 := S.Volume; T1 := S.Temperature;
    WriteLn(Format('adiabatic, phi 0.01 -> %.3f, piston at 0.05 (v_rms %.2f), t = %.0f:',
      [S.PackingFraction, Sqrt(3 * T0), S.Time]));
    WriteLn(Format('  V %.0f -> %.0f, T %.4f -> %.4f: ratio %.4f, ideal (V0/V1)^(2/3) %.4f',
      [V0, V1, T0, T1, T1 / T0, Power(V0 / V1, 2 / 3)]));
  finally
    S.Free;
  end;
end;

begin
  DefaultFormatSettings.DecimalSeparator := '.';
  RandSeed := 20261020;
  N := StrToIntDef(ParamStr(1), 20000);
  AllowIGPU := False;
  for i := 1 to ParamCount do
    if ParamStr(i) = '--allow-igpu' then AllowIGPU := True;
  Equilibrium(0.05, 400);
  Equilibrium(0.30, 400);
  Adiabatic;
end.
