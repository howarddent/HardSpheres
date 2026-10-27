program iso2;
{ N2O isotherms from the stepped-LJ model, 1M spheres, fixed volume per
  point, at T* 0.9134, 0.9803, 1.0350 (Tc) and 1.2478: 0, 20, 36 and 100 C
  as N2O. Appends to iso2.csv: Tstar, rho, P mean, 4 block means, T,
  overflows, capped substeps, wall seconds. Points already there are skipped.
  Usage: ./iso2 (run from the folder that holds iso2.csv). }
{$mode objfpc}{$H+}
uses {$IFDEF UNIX} cthreads, {$ENDIF} Classes, SysUtils, Math, DateUtils, clgl, ustickysim;
const
  N = 1000000;
  OutFile = 'iso2.csv';
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
  Pm, Tsum: Double;
  t0: TDateTime;
  OutF: TextFile;
  SL: TStringList;
begin
  DefaultFormatSettings.DecimalSeparator := '.';
  Randomize;
  // Tc* = 1.035: 0 C, 20 C, critical, 100 C. gas | near-critical | liquid
  Add(0.9134, [0.005, 0.01, 0.02, 0.03, 0.04, 0.05], 15, 20);
  Add(0.9134, [0.63, 0.64, 0.65, 0.66, 0.68], 20, 12);
  Add(0.9803, [0.005, 0.01, 0.02, 0.03, 0.045, 0.06, 0.075, 0.088], 15, 20);
  Add(0.9803, [0.54, 0.555, 0.57, 0.59, 0.61, 0.64], 20, 12);
  Add(1.0350, [0.005, 0.01, 0.02, 0.03, 0.045, 0.06, 0.08, 0.1, 0.13], 15, 20);
  Add(1.0350, [0.16, 0.2, 0.24, 0.29, 0.34, 0.39, 0.44], 60, 40);
  Add(1.0350, [0.49, 0.54, 0.59, 0.64], 20, 12);
  Add(1.2478, [0.005, 0.01, 0.02, 0.03, 0.045, 0.06, 0.08, 0.1, 0.13, 0.17, 0.22], 15, 20);
  Add(1.2478, [0.28, 0.35, 0.45], 15, 12);

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
    WriteLn(OutF, '# Tstar,rho,P,P1,P2,P3,P4,Tmeas,overflows,capped,seconds');
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
    Tsum := 0;
    for b := 0 to 3 do
    begin
      Blk[b] := 0;
      for f := 1 to fm do
      begin
        S.Frame;
        Blk[b] := Blk[b] + S.Pressure / fm;
        Tsum := Tsum + S.Temperature / (4 * fm);
      end;
    end;
    Pm := (Blk[0] + Blk[1] + Blk[2] + Blk[3]) / 4;
    WriteLn(OutF, Format('%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f,%d,%d,%.0f',
      [Key(Pts[i].T, Pts[i].Rho), Pm, Blk[0], Blk[1], Blk[2], Blk[3], Tsum, S.Overflows, S.CappedSweeps,
       SecondSpan(Now, t0)]));
    Flush(OutF);
    WriteLn(Format('T* %.4f rho %.3f: P %.5f (blocks %.5f %.5f %.5f %.5f), %.0f s',
      [Pts[i].T, Pts[i].Rho, Pm, Blk[0], Blk[1], Blk[2], Blk[3], SecondSpan(Now, t0)]));
  end;
  CloseFile(OutF);
  S.Free;
end.
