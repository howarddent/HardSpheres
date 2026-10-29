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

program bench;
{ Speed check: 10 time units of 1M spheres at T* 0.953 and rho 0.02, 0.3
  and 0.7, with the substep and sweep counts. Usage: ./bench }
{$mode objfpc}{$H+}
uses {$IFDEF UNIX} cthreads, {$ENDIF} SysUtils, Math, DateUtils, clgl, ustickysim;
var
  S: TStickySim; i: Integer; t0: TDateTime; rho: Double;
  Rhos: array[0..2] of Double = (0.02, 0.3, 0.7);
  k: Integer;
begin
  Randomize;
  S := TStickySim.Create;
  S.InitCL(ExtractFilePath(ParamStr(0)) + '../../Sticky/sticky.cl', nil, nil);
  WriteLn(S.DeviceName);
  S.Eps := 1; S.Gravity := 0; S.FloorEps := 0; S.Thermostat := True;
  S.SimTimePerFrame := 0.1;
  for k := 0 to 2 do
  begin
    rho := Rhos[k];
    S.Setup(1000000, Pi / 6 * rho, 0.953, LennardJonesPotential, False, 0, 1.0);
    S.PressureFrames := 1;
    t0 := Now;
    for i := 1 to 100 do S.Frame;
    WriteLn(Format('rho %.3f: %.1f s per 10 time units, substeps %d, sweeps %.2f capped %d, overflow %d, P %.4f T %.3f',
      [rho, MilliSecondsBetween(Now, t0) / 1000, S.Substeps, S.MeanSweeps, S.CappedSweeps, S.Overflows, S.Pressure, S.Temperature]));
  end;
  S.Free;
end.
