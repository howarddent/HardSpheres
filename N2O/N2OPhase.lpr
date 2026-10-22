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
