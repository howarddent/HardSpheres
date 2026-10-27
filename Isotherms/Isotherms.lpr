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
