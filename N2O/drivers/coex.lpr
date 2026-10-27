program coex;
{ Liquid-vapour coexistence and the critical temperature: 200k spheres
  at rho 0.31 (near rho_c), long runs. The wall pressure below Tc is the
  vapour pressure. The box interior is cut into cubes of two sizes;
  Binder's cumulant U4 = 1 - <m^4>/(3 <m^2>^2), m the cube's density less
  the mean, is ~0 for one Gaussian peak (above Tc) and ~2/3 for two (liquid
  and vapour, below Tc). The density histogram of the smaller cubes has its
  peaks at the coexisting densities.
  Usage: ./coex T* [T* ...]   Appends to coex.csv, histograms to
  histall_<T>.csv. }
{$mode objfpc}{$H+}
uses {$IFDEF UNIX} cthreads, {$ENDIF} Classes, SysUtils, Math, DateUtils, clgl, ustickysim;
const
  N = 200000;
  Rho = 0.31;
  TEquil = 300; TSample = 300; SampleEvery = 5;
  NGrids = 2;
  Grids: array[0..NGrids - 1] of Integer = (6, 14);   // cubes per side, the outer layer dropped
  HistBins = 60; HistMax = 0.9;
var
  S: TStickySim;
  Pos: array of Single;
  Cnt: array of Integer;
  T, P, side, mean, m, d: Double;
  M2, M4, U4: array[0..NGrids - 1] of Double;
  NS: array[0..NGrids - 1] of Int64;
  NH: Int64;
  Hist: array[0..HistBins - 1] of Double;
  a, g, f, fe, i, ix, iy, iz, nb, k: Integer;
  t0: TDateTime;
  OutF: TextFile;
begin
  DefaultFormatSettings.DecimalSeparator := '.';
  Randomize;
  S := TStickySim.Create;
  S.InitCL(ExtractFilePath(ParamStr(0)) + '../../Sticky/sticky.cl', nil, nil);
  S.Eps := 1; S.Gravity := 0; S.FloorEps := 0; S.Capture := False; S.Thermostat := True;
  S.SimTimePerFrame := 0.1;
  SetLength(Pos, 4 * N);
  for a := 1 to ParamCount do
  begin
    T := StrToFloat(ParamStr(a));
    t0 := Now;
    S.Setup(N, Pi / 6 * Rho, T, LennardJonesPotential, False, 0, 1.0);
    S.PressureFrames := 1;
    fe := Round(TEquil / S.SimTimePerFrame);
    for f := 1 to fe do S.Frame;
    for g := 0 to NGrids - 1 do begin M2[g] := 0; M4[g] := 0; NS[g] := 0; end;
    FillChar(Hist, SizeOf(Hist), 0); NH := 0;
    P := 0;
    fe := Round(TSample / S.SimTimePerFrame);
    for f := 1 to fe do
    begin
      S.Frame;
      P := P + S.Pressure / fe;
      if f mod Round(SampleEvery / S.SimTimePerFrame) <> 0 then Continue;
      S.ReadPositions(Pos);
      for g := 0 to NGrids - 1 do
      begin
        nb := Grids[g];
        side := (S.Lx - 1) / nb;              // centres lie in [0.5, L - 0.5]
        SetLength(Cnt, nb * nb * nb);
        FillChar(Cnt[0], Length(Cnt) * SizeOf(Integer), 0);
        for i := 0 to N - 1 do
        begin
          ix := Min(nb - 1, Max(0, Floor((Pos[4 * i] - 0.5) / side)));
          iy := Min(nb - 1, Max(0, Floor((Pos[4 * i + 1] - 0.5) / side)));
          iz := Min(nb - 1, Max(0, Floor((Pos[4 * i + 2] - 0.5) / side)));
          Inc(Cnt[(iz * nb + iy) * nb + ix]);
        end;
        // the mean over the interior cubes, then the moments about it
        mean := 0; k := 0;
        for iz := 1 to nb - 2 do for iy := 1 to nb - 2 do for ix := 1 to nb - 2 do
        begin
          mean := mean + Cnt[(iz * nb + iy) * nb + ix]; Inc(k);
        end;
        mean := mean / k;
        for iz := 1 to nb - 2 do for iy := 1 to nb - 2 do for ix := 1 to nb - 2 do
        begin
          m := (Cnt[(iz * nb + iy) * nb + ix] - mean) / Power(side, 3);
          M2[g] := M2[g] + m * m; M4[g] := M4[g] + m * m * m * m; Inc(NS[g]);
        end;
        if g = NGrids - 1 then
          for k := 0 to High(Cnt) do
          begin
            d := Cnt[k] / Power(side, 3);
            Hist[Min(HistBins - 1, Floor(d / HistMax * HistBins))] += 1;
            Inc(NH);
          end;
      end;
    end;
    for g := 0 to NGrids - 1 do
    begin
      M2[g] := M2[g] / NS[g]; M4[g] := M4[g] / NS[g];
      U4[g] := 1 - M4[g] / (3 * M2[g] * M2[g]);
    end;
    AssignFile(OutF, 'coex.csv');
    if FileExists('coex.csv') then Append(OutF) else
    begin
      Rewrite(OutF);
      WriteLn(OutF, '# Tstar,P,U4_cube', Grids[0], ',U4_cube', Grids[1], ',sd_cube', Grids[0], ',sd_cube', Grids[1], ',seconds');
    end;
    WriteLn(OutF, Format('%.4f,%.6f,%.4f,%.4f,%.5f,%.5f,%.0f',
      [T, P, U4[0], U4[1], Sqrt(M2[0]), Sqrt(M2[1]), SecondSpan(Now, t0)]));
    CloseFile(OutF);
    AssignFile(OutF, Format('histall_%.4f.csv', [T]));
    Rewrite(OutF);
    WriteLn(OutF, '# density (cube side ', (S.Lx - 1) / Grids[NGrids - 1]:0:2, '),fraction');
    for k := 0 to HistBins - 1 do
      WriteLn(OutF, Format('%.4f,%.5f', [(k + 0.5) * HistMax / HistBins, Hist[k] / NH]));
    CloseFile(OutF);
    WriteLn(Format('T* %.4f: P %.5f, U4 %.3f (cube %.1f) %.3f (cube %.1f), %.0f s',
      [T, P, U4[0], (S.Lx - 1) / Grids[0], U4[1], (S.Lx - 1) / Grids[1], SecondSpan(Now, t0)]));
    Flush(Output);
  end;
  S.Free;
end.
