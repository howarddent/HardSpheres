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
