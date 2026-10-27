program iso3;
{ The liquid and near-critical isotherm points again, with longer
  equilibration (60 units, for the box's breathing mode to die out), and
  the density measured in the box interior, away from the walls. Appends
  to iso3.csv: as iso2.csv plus rho_interior. Points already there are
  skipped. Usage: ./iso3 (run from the folder that holds iso3.csv). }
{$mode objfpc}{$H+}
uses {$IFDEF UNIX} cthreads, {$ENDIF} Classes, SysUtils, Math, DateUtils, clgl, ustickysim;
const
  N = 1000000;
  OutFile = 'iso3.csv';
  Margin = 10;      // interior: this far from every wall
  SnapEvery = 2;    // sim time units between interior counts
type
  TPoint = record T, Rho, Te, Tm: Double; end;
var
  Pts: array of TPoint;
  Done: TStringList;

procedure Add(T: Double; const Rhos: array of Double; Te, Tm: Double);
var r: Double;
begin
  for r in Rhos do
  begin
    SetLength(Pts, Length(Pts) + 1);
    Pts[High(Pts)].T := T; Pts[High(Pts)].Rho := r;
    Pts[High(Pts)].Te := Te; Pts[High(Pts)].Tm := Tm;
  end;
end;

function Key(T, Rho: Double): string;
begin
  Result := Format('%.4f,%.4f', [T, Rho]);
end;

var
  S: TStickySim;
  i, b, f, fe, fm: Integer;
  P, Blk: array[0..3] of Double;
  Pm, Tsum, RhoIn: Double;
  Pos: array of Single;
  ns, nin, j: Integer;
  t0: TDateTime;
  OutF: TextFile;
  SL: TStringList;
begin
  DefaultFormatSettings.DecimalSeparator := '.';
  Randomize;
  // liquid and near-critical points again: 60 units for the box's breathing
  // mode to die out, then the interior density measured
  Add(0.9134, [0.60, 0.615, 0.63, 0.645], 60, 28);
  Add(0.9803, [0.51, 0.525, 0.54, 0.555], 60, 28);
  Add(1.0350, [0.1, 0.13, 0.16, 0.2, 0.24, 0.29, 0.34, 0.39], 60, 40);
  Add(1.0350, [0.44, 0.49, 0.54], 60, 28);
  Add(1.2478, [0.13, 0.17, 0.22, 0.28], 60, 28);
  Add(1.2478, [0.005, 0.01, 0.02, 0.035, 0.06, 0.09], 20, 28);
  Add(0.9134, [0.06], 30, 28);
  Add(0.9803, [0.10], 30, 28);

  Done := TStringList.Create;
  if FileExists(OutFile) then
  begin
    SL := TStringList.Create;
    SL.LoadFromFile(OutFile);
    for i := 0 to SL.Count - 1 do
      if (SL[i] <> '') and (SL[i][1] <> '#') then
        Done.Add(Copy(SL[i], 1, Length(Key(1, 0.1))));
    SL.Free;
  end;
  AssignFile(OutF, OutFile);
  if FileExists(OutFile) then Append(OutF) else
  begin
    Rewrite(OutF);
    WriteLn(OutF, '# Tstar,rho,P,P1,P2,P3,P4,Tmeas,overflows,capped,seconds,rho_interior');
  end;

  S := TStickySim.Create;
  S.InitCL(ExtractFilePath(ParamStr(0)) + '../../Sticky/sticky.cl', nil, nil);
  WriteLn('OpenCL on ', S.DeviceName, ', ', Length(Pts), ' points, ', Done.Count, ' done');
  S.Eps := 1; S.Gravity := 0; S.FloorEps := 0; S.Thermostat := True;
  S.SimTimePerFrame := 0.1;
  for i := 0 to High(Pts) do
  begin
    if Done.IndexOf(Key(Pts[i].T, Pts[i].Rho)) >= 0 then Continue;
    t0 := Now;
    S.Setup(N, Pi / 6 * Pts[i].Rho, Pts[i].T, LennardJonesPotential, False, 0, 1.0);
    S.PressureFrames := 1;
    fe := Round(Pts[i].Te / S.SimTimePerFrame);
    fm := Round(Pts[i].Tm / S.SimTimePerFrame / 4);
    for f := 1 to fe do S.Frame;
    Tsum := 0; RhoIn := 0; ns := 0;
    SetLength(Pos, 4 * N);
    for b := 0 to 3 do
    begin
      Blk[b] := 0;
      for f := 1 to fm do
      begin
        S.Frame;
        Blk[b] := Blk[b] + S.Pressure / fm;
        Tsum := Tsum + S.Temperature / (4 * fm);
        if f mod Round(SnapEvery / S.SimTimePerFrame) = 0 then
        begin
          S.ReadPositions(Pos);
          nin := 0;
          for j := 0 to N - 1 do
            if (Pos[4 * j] > Margin) and (Pos[4 * j] < S.Lx - Margin)
               and (Pos[4 * j + 1] > Margin) and (Pos[4 * j + 1] < S.Ly - Margin)
               and (Pos[4 * j + 2] > Margin) and (Pos[4 * j + 2] < S.Lz - Margin) then Inc(nin);
          RhoIn := RhoIn + nin / ((S.Lx - 2 * Margin) * (S.Ly - 2 * Margin) * (S.Lz - 2 * Margin));
          Inc(ns);
        end;
      end;
    end;
    Pm := (Blk[0] + Blk[1] + Blk[2] + Blk[3]) / 4;
    RhoIn := RhoIn / Max(1, ns);
    WriteLn(OutF, Format('%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f,%d,%d,%.0f,%.5f',
      [Key(Pts[i].T, Pts[i].Rho), Pm, Blk[0], Blk[1], Blk[2], Blk[3], Tsum, S.Overflows, S.CappedSweeps,
       SecondSpan(Now, t0), RhoIn]));
    Flush(OutF);
    WriteLn(Format('T* %.4f rho %.3f: P %.5f (blocks %.5f %.5f %.5f %.5f), %.0f s',
      [Pts[i].T, Pts[i].Rho, Pm, Blk[0], Blk[1], Blk[2], Blk[3], SecondSpan(Now, t0)]));
  end;
  CloseFile(OutF);
  S.Free;
end.
