program stickycheck;

{******************************************************************************
  Headless checks of the sticky-sphere simulation (ustickysim.pas), on the
  OpenCL GPU without a window. SW: square well, width 1.5; LJ: stepped
  Lennard-Jones, cutoff 2.

  1. Eps = 0: plain hard spheres; PV/NT against Carnahan-Starling.
  2. Capture off, no thermostat: KE + PE_gravity + pair energy constant
     (SW, LJ, LJ in gravity, and with the floor's attraction).
  3. Capture off, thermostat, low density, equilibrated: pairs per sphere
     within the cutoff and within 1.5 (the bonds) against the Boltzmann
     counts rho sum_s e^(-Eps U_s/T) (4 pi/3)(R_s+1^3 - R_s^3) (a little
     high: spheres near the walls have fewer neighbours), and PV/NT
     against 1 + B2 rho for the stepped potential.
  4. Gravity, Eps = 0: the mean height against the barometric profile,
     ideal and with the hard-sphere (Carnahan-Starling) pressure.
  5. Capture on, thermostat: bonding, capture and breaking, the clusters,
     and the largest distance any listed pair is outside its shell (a
     health check: events missed).

  Build: fpc -O2 -Fu.. -FUlib stickycheck.lpr      Run: ./stickycheck [N [max sweeps]] [--allow-igpu]
*******************************************************************************}

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX} cthreads, {$ENDIF}
  SysUtils, Math, clgl, ustickysim;

var
  S: TStickySim;
  N, Sweeps: Integer;
  Nbr, Shells: array of LongInt;
  SW, LJ: TStepPotential;

procedure Run(Frames: Integer);
var
  f: Integer;
begin
  for f := 1 to Frames do S.Frame;
end;

procedure Configure(Eps, G: Double; Capture, Hold: Boolean);
begin
  S.Eps := Eps;
  S.Gravity := G;
  S.Capture := Capture;
  S.Thermostat := Hold;
end;

procedure HardSphereCheck(Phi: Double);
var
  T, P, rho: Double;
  f: Integer;
begin
  Configure(0, 0, False, False);
  S.Setup(N, Phi, 1.0, SW, False, 0);
  S.PressureFrames := 1;
  Run(200);
  T := 0; P := 0;
  for f := 1 to 400 do
  begin
    S.Frame;
    T := T + S.Temperature / 400; P := P + S.Pressure / 400;
  end;
  rho := N / S.Volume;
  WriteLn(Format('  phi %.3f  T %.4f  PV/NT %.4f  (Carnahan-Starling %.4f; walls lower it a little)',
    [Phi, T, P / (rho * T), (1 + Phi + Phi * Phi - Phi * Phi * Phi) / Power(1 - Phi, 3)]));
end;

procedure EnergyCheck(const Pot: TStepPotential; Phi, T0, G: Double; EF: Double = 0);
var
  E0, Emin, Emax: Double;
  f: Integer;
begin
  Configure(1, G, False, False);
  S.FloorEps := EF;
  S.Setup(N, Phi, T0, Pot, False, 0);
  S.Frame;
  E0 := S.Energy;
  Emin := E0; Emax := E0;
  for f := 1 to 1000 do
  begin
    S.Frame;
    Emin := Min(Emin, S.Energy); Emax := Max(Emax, S.Energy);
  end;
  S.FloorEps := 0;
  WriteLn(Format('  %s, phi %.2f, g %.2f, floor %.1f: T %.3f -> %.3f, %d bonds, pair energy/N %.3f;  E/N %.5f -> %.5f, range %.2g',
    [Pot.Name, Phi, G, EF, T0, S.Temperature, S.Bonds, S.PairEnergy / N, E0 / N, S.Energy / N, (Emax - Emin) / N]));
end;

procedure EquilibriumCheck(const Pot: TStepPotential; Phi, T: Double);
var
  P, rho, nAll, nBond: Double;
  f, i, k, m, mb: Integer;
begin
  Configure(1, 0, False, True);
  S.Setup(N, Phi, T, Pot, False, 0);
  S.PressureFrames := 1;
  Run(2000);
  P := 0; nAll := 0; nBond := 0;
  for f := 1 to 40 do
  begin
    for i := 1 to 50 do
    begin
      S.Frame;
      P := P + S.Pressure / 2000;
    end;
    S.ReadLists(Nbr, Shells);
    m := 0; mb := 0;
    for i := 0 to N * MaxNeighbours - 1 do
      if Nbr[i] >= 0 then
      begin
        Inc(m);
        if Shells[i] < Pot.BondEdge then Inc(mb);
      end;
    nAll := nAll + m / N / 40; nBond := nBond + mb / N / 40;
  end;
  rho := N / S.Volume;
  WriteLn(Format('  %s, phi %.3f, T %.2f:', [Pot.Name, Phi, T]));
  WriteLn(Format('     neighbours within the cutoff %.4f (low density %.4f), within %.2f %.4f (%.4f);  PV/NT %.4f (1 + B2 rho %.4f)',
    [nAll, LowDensityNeighbours(Pot, Pot.K, 1, T, rho), Pot.R[Pot.BondEdge], nBond,
     LowDensityNeighbours(Pot, Pot.BondEdge, 1, T, rho), P / (rho * T), 1 + SecondVirial(Pot, 1, T) * rho]));
end;

// the mean height of Nb hard spheres in gravity from hydrostatic
// equilibrium with the Carnahan-Starling chemical potential,
// ln n + mu_ex/T = c - g z/T, mu_ex/T = (8 phi - 9 phi^2 + 3 phi^3)/(1 - phi)^3
function HydrostaticMeanHeight(Nb: Integer; Wd, Lz, T, G: Double): Double;
const
  M = 2000;
var
  A, w, c, lo, hi, tot, sz, z, n: Double;
  i, it: Integer;
  function Dens(rhs: Double): Double;
  var
    a0, b0, mm, ph, f: Double;
    k: Integer;
  begin
    a0 := -60; b0 := 2;
    for k := 1 to 80 do
    begin
      mm := (a0 + b0) / 2;
      ph := Exp(mm) * Pi / 6;
      f := mm + (8 * ph - 9 * ph * ph + 3 * ph * ph * ph) / Power(1 - ph, 3) - rhs;
      if f > 0 then b0 := mm else a0 := mm;
    end;
    Result := Exp((a0 + b0) / 2);
  end;
begin
  A := Sqr(Wd - 1); w := Lz - 1;
  lo := -40; hi := 5;
  for it := 1 to 50 do
  begin
    c := (lo + hi) / 2;
    tot := 0;
    for i := 0 to M - 1 do tot := tot + Dens(c - G * w * (i + 0.5) / M / T);
    if tot * w / M * A > Nb then hi := c else lo := c;
  end;
  tot := 0; sz := 0;
  for i := 0 to M - 1 do
  begin
    z := w * (i + 0.5) / M;
    n := Dens((lo + hi) / 2 - G * z / T);
    tot := tot + n; sz := sz + z * n;
  end;
  Result := sz / tot;
end;

procedure BarometricCheck(Nb: Integer; Phi, T, G: Double);
var
  P: array of Single;
  f, i: Integer;
  zm, H, w, th: Double;
begin
  Configure(0, G, False, True);
  S.Setup(Nb, Phi, T, SW, False, 0);
  SetLength(P, 4 * Nb);
  Run(4000);
  zm := 0;
  for f := 1 to 400 do
  begin
    Run(20);
    S.ReadPositions(P);
    for i := 0 to Nb - 1 do zm := zm + P[4 * i + 2];
  end;
  zm := zm / (400 * Nb) - 0.5;
  H := T / G;
  w := S.Lz - 1;
  th := H - w * Exp(-w / H) / (1 - Exp(-w / H));
  WriteLn(Format('  N %d phi %.3f T %.2f g %.3f (scale height %.1f, box %.1f): mean height %.3f;  ideal barometric %.3f, with hard-sphere pressure %.3f',
    [Nb, Phi, T, G, H, S.Lz, zm, th, HydrostaticMeanHeight(Nb, S.Ly, S.Lz, T, G)]));
end;

procedure CaptureCheck(const Pot: TStepPotential; Phi, T0: Double; Frames: Integer);
var
  C: TClusterStats;
  P: array of Single;
  f, k: Integer;
  line: string;
begin
  Configure(1, 0, True, True);
  S.Setup(N, Phi, T0, Pot, False, 0);
  S.PressureFrames := 20;
  SetLength(P, 4 * N);
  for f := 1 to Frames do
  begin
    S.Frame;
    if (f mod (Frames div 4) = 0) then
    begin
      S.ReadLists(Nbr, Shells);
      S.ReadPositions(P);
      C := S.ClusterStats(Nbr, Shells);
      line := '';
      for k := 0 to Min(C.NBins - 1, 14) do line := line + Format(' %.2f', [C.Bins[k]]);
      WriteLn(Format('  t %3.0f  bonds/sphere %.2f  free %.3f  clusters %d  largest %d  per sphere per time: formed %.3g captured %.3g broken %.3g;  max outside its shell %.3f;  sweeps %.1f (%d capped), substeps %d',
        [S.Time, 2 * S.Bonds / N, C.Monomers / N, C.Clusters, C.Largest,
         2 * S.FormRate, 2 * S.CaptureRate, 2 * S.BreakRate, S.MaxOvershoot(Nbr, Shells, P),
         S.MeanSweeps, S.CappedSweeps, S.Substeps]));
      WriteLn('        fraction in clusters of 1, 2-3, 4-7, ...:', line);
    end;
  end;
end;

const
  Temps: array[0..2] of Double = (2.0, 1.0, 0.6);
var
  Msg: string;
  k, i: Integer;
  AllowIGPU: Boolean;
  Args: array of string;
begin
  // --allow-igpu anywhere; the rest are N and max sweeps, in order
  AllowIGPU := False;
  Args := nil;
  for i := 1 to ParamCount do
    if ParamStr(i) = '--allow-igpu' then AllowIGPU := True
    else Insert(ParamStr(i), Args, Length(Args));
  if Length(Args) >= 1 then N := StrToIntDef(Args[0], 20000) else N := 20000;
  if Length(Args) >= 2 then Sweeps := StrToInt(Args[1]) else Sweeps := 32;
  if not LoadOpenCL(Msg) then begin WriteLn(Msg); Halt(1); end;
  Randomize;
  SW := SquareWellPotential(1.5);
  LJ := LennardJonesPotential;
  S := TStickySim.Create;
  try
    S.MaxSweeps := Sweeps;
    S.AllowIntegrated := AllowIGPU;
    S.InitCL(ExtractFilePath(ParamStr(0)) + 'sticky.cl', nil, nil);
    WriteLn('OpenCL on ', S.DeviceName, ', N = ', N, ', at most ', Sweeps, ' sweeps per substep');
    Write('LJ shells:');
    for k := 0 to LJ.K - 1 do Write(Format(' [%.2f,%.2f) %.3f', [LJ.R[k], LJ.R[k + 1], LJ.U[k]]));
    WriteLn;
    SetLength(Nbr, N * MaxNeighbours);
    SetLength(Shells, N * MaxNeighbours);
    WriteLn('1. Eps = 0 (hard spheres)');
    HardSphereCheck(0.05);
    WriteLn('2. capture off, no thermostat: energy');
    EnergyCheck(SW, 0.05, 0.6, 0);
    EnergyCheck(LJ, 0.05, 0.6, 0);
    EnergyCheck(LJ, 0.05, 1.0, 0.05);
    EnergyCheck(LJ, 0.05, 1.0, 0.05, 1.5);
    WriteLn('3. capture off, thermostat: equilibrium at low density');
    EquilibriumCheck(SW, 0.01, 2.0);
    EquilibriumCheck(LJ, 0.01, 2.0);
    EquilibriumCheck(LJ, 0.01, 1.2);
    WriteLn('4. gravity, Eps = 0');
    BarometricCheck(2000, 0.002, 1.0, 0.03);
    BarometricCheck(2000, 0.002, 1.0, 0.1);
    WriteLn('5. capture on, thermostat, phi 0.05');
    for k := 0 to High(Temps) do
    begin
      WriteLn(Format(' %s, T = %.2f', [LJ.Name, Temps[k]]));
      CaptureCheck(LJ, 0.05, Temps[k], 800);
    end;
  finally
    S.Free;
  end;
end.
