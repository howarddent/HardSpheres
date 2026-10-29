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

unit un2ophase;

{******************************************************************************
  The pressure-temperature phase diagram of nitrous oxide, drawn on a
  canvas, with a red dot for the current state of a simulation.

  Data:
    - liquid-vapour (vapour pressure) curve, triple point 182.33 K /
      0.878 bar to critical point 309.52 K / 72.45 bar: the NIST
      Chemistry WebBook's reference equation of state (n2odata.inc),
      interpolated as ln P linear in 1/T;
    - solid-vapour (sublimation) curve: Clausius-Clapeyron from the
      triple point with the enthalpy of sublimation, the enthalpy of
      vaporisation there (16.59 kJ/mol, NIST) plus that of fusion
      (6.54 kJ/mol);
    - solid-liquid (melting) curve: a straight Clapeyron line from the
      triple point, dP/dT = dH_fus/(T dv) ~ 44 bar/K with the solid's
      density taken as ~1600 kg/m3 (approximate: drawn for shape).

  Axes: temperature in °C, pressure in bar on a log scale, so the
  triple point (under 1 bar) and the critical point (72 bar) both show.
*******************************************************************************}

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Math, Controls, Graphics, Types;

{$I n2odata.inc}

const
  N2OTc = 309.52;            // K
  N2OPc = 72.447;            // bar
  N2OTt = 182.33;            // K
  N2OPt = 0.8784;            // bar
  N2OTb = 184.67;            // K, normal boiling point (1.01325 bar)
  N2OSubSlope = 2782;        // K: dH_sub/R, (16.59 + 6.54) kJ/mol / R
  N2OMeltSlope = 44;         // bar/K, approximate

// coexistence pressures, bar
function N2OVapourPressure(TK: Double): Double;      // liquid-vapour, Tt..Tc
function N2OSublimationPressure(TK: Double): Double; // solid-vapour, below Tt
function N2OCoexistencePressure(TK: Double): Double; // whichever applies (Tt..Tc or below)
function N2OMeltingTemperature(Pbar: Double): Double;
// the phase at (TK, Pbar), as a word
function N2OPhaseName(TK, Pbar: Double): string;

type
  TStatePoint = record
    T, P: Double;              // K, bar
  end;

  TN2OPhaseDiagram = class(TCustomControl)
  private
    FHasState: Boolean;
    FStateT, FStateP: Double;       // K, bar
    FStateText: string;
    FTrail: array of TStatePoint;        // (T K, P bar), oldest first
    FTrailLen: Integer;
    function XOf(TC: Double): Integer;
    function YOf(P: Double): Integer;
    procedure CurvePoints(var Pts: array of TPoint; var Count: Integer; FromT, ToT: Double;
      Kind: Integer);
  protected
    procedure Paint; override;
  public
    TMinC, TMaxC: Double;            // axis range, deg C
    PMin, PMax: Double;              // bar
    constructor Create(AOwner: TComponent); override;
    // the state to mark (K, bar) and a line describing it; the dot leaves a
    // short trail, so a slow change draws its path
    procedure SetState(TK, Pbar: Double; const ALabel: string);
    procedure ClearState;
  end;

implementation

const
  PlotL = 52; PlotR = 12; PlotT = 24; PlotB = 62;
  TrailMax = 400;

function N2OVapourPressure(TK: Double): Double;
var
  i: Integer;
  f: Double;
begin
  if TK <= N2OSat[0, 0] then Exit(N2OSat[0, 1]);
  if TK >= N2OSat[N2OSatCount - 1, 0] then Exit(N2OSat[N2OSatCount - 1, 1]);
  i := 0;
  while N2OSat[i + 1, 0] < TK do Inc(i);
  f := (1 / TK - 1 / N2OSat[i, 0]) / (1 / N2OSat[i + 1, 0] - 1 / N2OSat[i, 0]);
  Result := Exp(Ln(N2OSat[i, 1]) + f * (Ln(N2OSat[i + 1, 1]) - Ln(N2OSat[i, 1])));
end;

function N2OSublimationPressure(TK: Double): Double;
begin
  Result := N2OPt * Exp(-N2OSubSlope * (1 / TK - 1 / N2OTt));
end;

function N2OCoexistencePressure(TK: Double): Double;
begin
  if TK < N2OTt then Result := N2OSublimationPressure(TK)
  else Result := N2OVapourPressure(TK);
end;

function N2OMeltingTemperature(Pbar: Double): Double;
begin
  Result := N2OTt + (Pbar - N2OPt) / N2OMeltSlope;
end;

function N2OPhaseName(TK, Pbar: Double): string;
begin
  if (TK < N2OTt) and (Pbar > N2OSublimationPressure(TK)) then Result := 'solid'
  else if (TK >= N2OTt) and (TK < N2OMeltingTemperature(Pbar)) and (Pbar > N2OPt) then Result := 'solid'
  else if TK >= N2OTc then
  begin
    if Pbar >= N2OPc then Result := 'supercritical fluid' else Result := 'gas';
  end
  else if Pbar > N2OCoexistencePressure(TK) then Result := 'liquid'
  else Result := 'vapour';
end;

{ TN2OPhaseDiagram }

constructor TN2OPhaseDiagram.Create(AOwner: TComponent);
begin
  inherited Create(AOwner);
  TMinC := -120; TMaxC := 60;
  PMin := 0.1; PMax := 300;
  SetLength(FTrail, TrailMax);
  DoubleBuffered := True;
  Color := clWhite;
end;

function TN2OPhaseDiagram.XOf(TC: Double): Integer;
begin
  Result := PlotL + Round((TC - TMinC) / (TMaxC - TMinC) * (ClientWidth - PlotL - PlotR));
end;

function TN2OPhaseDiagram.YOf(P: Double): Integer;
var
  f: Double;
begin
  P := EnsureRange(P, PMin * 0.5, PMax * 2);
  f := (Log10(P) - Log10(PMin)) / (Log10(PMax) - Log10(PMin));
  Result := ClientHeight - PlotB - Round(f * (ClientHeight - PlotT - PlotB));
end;

// points along a curve between two temperatures (K). Kind 0: vapour
// pressure, 1: sublimation
procedure TN2OPhaseDiagram.CurvePoints(var Pts: array of TPoint; var Count: Integer;
  FromT, ToT: Double; Kind: Integer);
const
  Steps = 120;
var
  k: Integer;
  T, P: Double;
begin
  for k := 0 to Steps do
  begin
    T := FromT + (ToT - FromT) * k / Steps;
    if Kind = 0 then P := N2OVapourPressure(T) else P := N2OSublimationPressure(T);
    Pts[Count] := Point(XOf(T - 273.15), YOf(P));
    Inc(Count);
  end;
end;

procedure TN2OPhaseDiagram.SetState(TK, Pbar: Double; const ALabel: string);
begin
  if FHasState and (Abs(TK - FStateT) < 1e-9) and (Abs(Pbar - FStateP) < 1e-9) and (ALabel = FStateText) then Exit;
  FHasState := True;
  FStateT := TK; FStateP := Pbar; FStateText := ALabel;
  if FTrailLen = TrailMax then
  begin
    Move(FTrail[1], FTrail[0], (TrailMax - 1) * SizeOf(TStatePoint));
    Dec(FTrailLen);
  end;
  FTrail[FTrailLen].T := TK;
  FTrail[FTrailLen].P := Pbar;
  Inc(FTrailLen);
  Invalidate;
end;

procedure TN2OPhaseDiagram.ClearState;
begin
  FHasState := False;
  FTrailLen := 0;
  Invalidate;
end;

procedure TN2OPhaseDiagram.Paint;
var
  C: TCanvas;
  Pts: array of TPoint;
  n, k, x, y, x0, x1, y0, y1, xt, yt, xc, yc: Integer;
  TK, TminK, TmaxK, Tsub0, P, v: Double;
  s: string;
  dec: Integer;

  procedure Region(AColor: TColor);
  begin
    C.Brush.Color := AColor;
    C.Pen.Style := psClear;
    C.Polygon(Pts[0..n - 1]);
    C.Pen.Style := psSolid;
  end;

  procedure Add(px, py: Integer);
  begin
    Pts[n] := Point(px, py);
    Inc(n);
  end;

  procedure Lbl(px, py: Integer; const Txt: string; AColor: TColor; Bold: Boolean = False; Centre: Boolean = True);
  begin
    C.Brush.Style := bsClear;
    C.Font.Color := AColor;
    if Bold then C.Font.Style := [fsBold] else C.Font.Style := [];
    if Centre then
      C.TextOut(px - C.TextWidth(Txt) div 2, py - C.TextHeight(Txt) div 2, Txt)
    else
      C.TextOut(px, py, Txt);
    C.Font.Style := [];
    C.Brush.Style := bsSolid;
  end;

begin
  C := Canvas;
  C.Brush.Color := clWhite;
  C.FillRect(ClientRect);
  C.Font.Size := 8;
  x0 := XOf(TMinC); x1 := XOf(TMaxC); y0 := YOf(PMax); y1 := YOf(PMin);
  TminK := TMinC + 273.15; TmaxK := TMaxC + 273.15;
  SetLength(Pts, 1000);
  xt := XOf(N2OTt - 273.15); yt := YOf(N2OPt);
  xc := XOf(N2OTc - 273.15); yc := YOf(N2OPc);
  C.ClipRect := Rect(x0, y0, x1 + 1, y1 + 1);
  C.Clipping := True;

  // regions: vapour (whole plot first), then the others over it
  C.Brush.Color := RGBToColor(255, 248, 214);
  C.FillRect(Rect(x0, y0, x1 + 1, y1 + 1));
  // gas: T > Tc, P < Pc
  n := 0;
  Add(xc, yc); Add(x1 + 1, yc); Add(x1 + 1, y1 + 1); Add(xc, y1 + 1);
  Region(RGBToColor(255, 228, 200));
  // supercritical fluid: T > Tc, P > Pc
  n := 0;
  Add(xc, y0); Add(x1 + 1, y0); Add(x1 + 1, yc); Add(xc, yc);
  Region(RGBToColor(232, 220, 250));
  // liquid: above the vapour-pressure curve, right of the melting line, left of Tc
  n := 0;
  Add(XOf(N2OMeltingTemperature(PMax * 2) - 273.15), y0);
  Add(xc, y0);
  CurvePoints(Pts, n, N2OTc, N2OTt, 0);
  Region(RGBToColor(205, 228, 255));
  // solid: left of the sublimation curve and the melting line
  n := 0;
  Add(x0, y0);
  Add(XOf(N2OMeltingTemperature(PMax * 2) - 273.15), y0);
  Add(xt, yt);
  CurvePoints(Pts, n, N2OTt, TminK, 1);
  Add(x0, y1 + 1);
  Region(RGBToColor(215, 215, 222));

  // grid
  C.Pen.Color := RGBToColor(200, 200, 200);
  C.Pen.Width := 1;
  C.Pen.Style := psDot;
  k := Ceil(TMinC / 20) * 20;
  while k <= TMaxC do
  begin
    C.Line(XOf(k), y0, XOf(k), y1);
    Inc(k, 20);
  end;
  for dec := -1 to 2 do
  begin
    C.Line(x0, YOf(Power(10, dec)), x1, YOf(Power(10, dec)));
  end;
  C.Pen.Style := psSolid;

  // reference lines: Tc and Pc beyond the critical point, 1 atm, 20 deg C
  C.Pen.Color := RGBToColor(120, 90, 160);
  C.Pen.Style := psDash;
  C.Line(xc, y0, xc, y1);
  C.Line(xc, yc, x1, yc);
  C.Pen.Color := RGBToColor(110, 110, 110);
  C.Line(x0, YOf(1.01325), x1, YOf(1.01325));
  C.Line(XOf(20), y0, XOf(20), y1);
  C.Pen.Style := psSolid;

  // the curves
  C.Pen.Width := 2;
  C.Pen.Color := RGBToColor(20, 60, 160);
  n := 0;
  CurvePoints(Pts, n, N2OTt, N2OTc, 0);
  C.Polyline(Pts[0..n - 1]);
  C.Pen.Color := RGBToColor(60, 60, 60);
  n := 0;
  CurvePoints(Pts, n, TminK, N2OTt, 1);
  C.Polyline(Pts[0..n - 1]);
  C.Line(xt, yt, XOf(N2OMeltingTemperature(PMax * 2) - 273.15), YOf(PMax * 2));
  C.Pen.Width := 1;

  // region labels
  Lbl(XOf(-108), YOf(40), 'SOLID', RGBToColor(80, 80, 90), True);
  Lbl(XOf(-25), YOf(120), 'LIQUID', RGBToColor(20, 60, 160), True);
  Lbl(XOf(-20), YOf(0.4), 'VAPOUR', RGBToColor(150, 110, 0), True);
  Lbl(XOf(48), YOf(8), 'GAS', RGBToColor(170, 90, 20), True);
  Lbl(XOf(48), YOf(180), 'SUPER-', RGBToColor(100, 60, 150), True);
  Lbl(XOf(48), YOf(130), 'CRITICAL', RGBToColor(100, 60, 150), True);

  // special points
  C.Brush.Color := clBlack;
  C.Pen.Color := clBlack;
  C.Ellipse(xt - 4, yt - 4, xt + 5, yt + 5);
  C.Ellipse(xc - 4, yc - 4, xc + 5, yc + 5);
  x := XOf(20); y := YOf(N2OVapourPressure(293.15));
  C.Brush.Color := RGBToColor(0, 120, 60);
  C.Pen.Color := RGBToColor(0, 120, 60);
  C.Rectangle(x - 4, y - 4, x + 5, y + 5);
  x := XOf(N2OTb - 273.15); y := YOf(1.01325);
  C.Brush.Color := RGBToColor(110, 110, 110);
  C.Pen.Color := RGBToColor(110, 110, 110);
  C.Ellipse(x - 3, y - 3, x + 4, y + 4);
  Lbl(xt + 7, yt + 2, 'triple point −90.8 °C, 0.88 bar', clBlack, False, False);
  s := 'critical point 36.4 °C, 72.4 bar';
  Lbl(xc - 8 - C.TextWidth(s), yc - 20, s, clBlack, False, False);
  s := 'cylinder at 20 °C: 50.5 bar';
  Lbl(XOf(20) - 8 - C.TextWidth(s), YOf(N2OVapourPressure(293.15)) - 6, s, RGBToColor(0, 110, 50), False, False);
  s := '1 atm (boils at −88.5 °C)';
  Lbl(x1 - 4 - C.TextWidth(s), YOf(1.01325) - C.TextHeight(s) - 2, s, RGBToColor(90, 90, 90), False, False);

  // the state: a trail and a red dot
  if FHasState then
  begin
    C.Pen.Color := RGBToColor(240, 140, 140);
    C.Pen.Width := 2;
    for k := 1 to FTrailLen - 1 do
      C.Line(XOf(FTrail[k - 1].T - 273.15), YOf(FTrail[k - 1].P), XOf(FTrail[k].T - 273.15), YOf(FTrail[k].P));
    C.Pen.Width := 1;
    x := XOf(FStateT - 273.15); y := YOf(FStateP);
    C.Brush.Color := clRed;
    C.Pen.Color := RGBToColor(120, 0, 0);
    C.Ellipse(x - 6, y - 6, x + 7, y + 7);
  end;
  C.Clipping := False;

  // frame and axes
  C.Brush.Style := bsClear;
  C.Pen.Color := clBlack;
  C.Rectangle(x0, y0, x1 + 1, y1 + 1);
  C.Brush.Style := bsSolid;
  C.Font.Color := clBlack;
  k := Ceil(TMinC / 20) * 20;
  while k <= TMaxC do
  begin
    s := IntToStr(k);
    C.Line(XOf(k), y1, XOf(k), y1 + 4);
    Lbl(XOf(k), y1 + 12, s, clBlack);
    Inc(k, 20);
  end;
  for dec := -1 to 2 do
  begin
    v := Power(10, dec);
    if v < 1 then s := FormatFloat('0.#', v) else s := IntToStr(Round(v));
    y := YOf(v);
    C.Line(x0 - 4, y, x0, y);
    Lbl(x0 - 8 - C.TextWidth(s) div 2, y, s, clBlack);
    for k := 2 to 9 do
    begin
      P := v * k;
      if (P > PMin) and (P < PMax) then C.Line(x0 - 2, YOf(P), x0, YOf(P));
    end;
  end;
  Lbl((x0 + x1) div 2, y1 + 26, 'temperature (°C)', clBlack);
  C.Font.Orientation := 900;
  s := 'pressure (bar, log scale)';
  C.Brush.Style := bsClear;
  C.TextOut(4, (y0 + y1) div 2 + C.TextWidth(s) div 2, s);
  C.Font.Orientation := 0;
  Lbl((x0 + x1) div 2, 10, 'Nitrous oxide: phase diagram', clBlack, True);
  if FHasState then
    Lbl(x0, y1 + 38, FStateText, RGBToColor(170, 0, 0), True, False);
end;

end.
