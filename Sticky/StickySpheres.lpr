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

program StickySpheres;

{******************************************************************************
  StickySpheres: hard spheres with a square-well attraction, which bind
  into clusters and break apart again, simulated in OpenCL and drawn in
  OpenGL from the same GPU buffer. See ustickymain.pas (the window),
  ustickysim.pas (the simulation), sticky.cl (the kernels) and
  stickycheck.lpr (headless checks). Shares ../clgl.pas with HardSpheres.
*******************************************************************************}

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}
  cthreads,
  {$ENDIF}
  Interfaces, Forms, ustickymain;

begin
  RequireDerivedFormResource := True;
  Application.Scaled := True;
  Application.Initialize;
  Application.CreateForm(TfmMain, fmMain);
  Application.Run;
end.
