program HardSpheres;

{******************************************************************************
  HardSpheres: a hard-sphere gas in a box with a piston, simulated in
  OpenCL and drawn in OpenGL from the same GPU buffer. See umain.pas
  (the window), ugassim.pas (the simulation), gas.cl (the kernels) and
  gascheck.lpr (headless checks against kinetic theory).
*******************************************************************************}

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}
  cthreads,
  {$ENDIF}
  Interfaces, Forms, umain;


begin
  RequireDerivedFormResource := True;
  Application.Scaled := True;
  Application.Initialize;
  Application.CreateForm(TfmMain, fmMain);
  Application.Run;
end.
