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

program N2OPhase;

{******************************************************************************
  N2OPhase: the sticky-sphere liquid and vapour of ../Sticky (stepped
  Lennard-Jones spheres, gravity, an attractive floor, a piston lid)
  beside the phase diagram of nitrous oxide, with a red dot for the
  simulation's state mapped onto N2O by corresponding states - a teaching
  aid for the phase diagram and isotherms of a saturated vapour. See
  un2omain.pas (the window), un2ophase.pas (the diagram), and
  ../Sticky/ustickysim.pas and sticky.cl (the simulation, shared).
*******************************************************************************}

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}
  cthreads,
  {$ENDIF}
  Interfaces, Forms, un2omain;

begin
  RequireDerivedFormResource := True;
  Application.Scaled := True;
  Application.Initialize;
  Application.CreateForm(TfmMain, fmMain);
  Application.Run;
end.
