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

program Isotherms;

{******************************************************************************
  Isotherms: spheres in a box with the intermolecular attraction on (a
  real gas, which condenses) or off (an ideal gas), beside the
  pressure-volume isotherm of nitrous oxide, carbon dioxide or water, with
  a red dot for the simulation's state mapped onto the substance by
  corresponding states. See uisomain.pas (the window), uisopv.pas (the
  diagram), usubstance.pas (the substances' data), and
  ../Sticky/ustickysim.pas and sticky.cl (the simulation, shared).
*******************************************************************************}

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}
  cthreads,
  {$ENDIF}
  Interfaces, Forms, uisomain;

begin
  RequireDerivedFormResource := True;
  Application.Scaled := True;
  Application.Initialize;
  Application.CreateForm(TfmMain, fmMain);
  Application.Run;
end.
