unit ugaspv;

{******************************************************************************
  The pressure-volume diagram of the hard-sphere gas, drawn on a canvas,
  with a red dot for the simulation's state and a trail of where it has
  been (as the piston moves, the path the gas actually takes).

  Volume is per sphere, v = V/N (diameters cubed), and pressure in the
  simulation's units (kT = T, mass 1, diameter 1); both axes linear, so an
  ideal-gas isotherm is Boyle's hyperbola. Drawn at the current
  temperature T:
    - the ideal-gas isotherm P = T/v, dashed;
    - the hard-sphere isotherm P = (T/v) Z(phi), Z the Carnahan-Starling
      compressibility factor and phi = pi/(6v): it parts from the ideal
      one only when the spheres' own volume starts to matter;
    - once the piston has moved without the thermostat, the ideal-gas
      adiabat through the state before the move, P v^(5/3) constant,
      dotted: a fast enough (but still slow) compression follows it and
      climbs across the isotherms as the gas heats.
  There is no flat part: with no attraction the gas never condenses (the
  contrast with ../N2O).
*******************************************************************************}

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Math, Controls, Graphics, ugassim;

type
  TGasPVPoint = record
    V, P: Double;
  end;

  TGasPVDiagram = class(TCustomControl)
  private
    FHasState: Boolean;
    FT, FV, FP: Double;
    FAdiabat: Double;             // P v^(5/3) of the adiabat, 0 for none
    FTrail: array of TGasPVPoint;
    FTrailLen: Integer;
    FPMax: Double;
    function XOf(V: Double): Integer;
    function YOf(P: Double): Integer;
  protected
    procedure Paint; override;
  public
    VMax: Double;                 // the volume axis's extent, per sphere
    constructor Create(AOwner: TComponent); override;
    procedure SetState(T, V, P, Adiabat: Double);
    procedure ClearState;
  end;

implementation

const
  PlotL = 58; PlotR = 12; PlotT = 24; PlotB = 48;
  TrailMax = 800;
  VClose = Pi / 6 / 0.5;          // per-sphere volume at packing 0.5

constructor TGasPVDiagram.Create(AOwner: TComponent);
begin
  inherited Create(AOwner);
  VMax := 25;
  FPMax := 1;
  SetLength(FTrail, TrailMax);
  DoubleBuffered := True;
  Color := clWhite;
end;

function TGasPVDiagram.XOf(V: Double): Integer;
begin
  V := EnsureRange(V, -0.5 * VMax, 2 * VMax);
  Result := PlotL + Round(V / VMax * (ClientWidth - PlotL - PlotR));
end;

function TGasPVDiagram.YOf(P: Double): Integer;
begin
  P := EnsureRange(P, -0.5 * FPMax, 2 * FPMax);
  Result := ClientHeight - PlotB - Round(P / FPMax * (ClientHeight - PlotT - PlotB));
end;

procedure TGasPVDiagram.SetState(T, V, P, Adiabat: Double);
begin
  FHasState := True;
  FT := T; FV := V; FP := P; FAdiabat := Adiabat;
  if FTrailLen = TrailMax then
  begin
    Move(FTrail[1], FTrail[0], (TrailMax - 1) * SizeOf(TGasPVPoint));
    Dec(FTrailLen);
  end;
  FTrail[FTrailLen].V := V;
  FTrail[FTrailLen].P := P;
  Inc(FTrailLen);
  Invalidate;
end;

procedure TGasPVDiagram.ClearState;
begin
  FHasState := False;
  FTrailLen := 0;
  Invalidate;
end;

// a round axis maximum at least x
function NiceMax(x: Double): Double;
const
  Steps: array[0..5] of Double = (1, 2, 2.5, 5, 7.5, 10);
var
  e, m: Double;
  k: Integer;
begin
  e := Power(10, Floor(Log10(x)));
  m := x / e;
  for k := 0 to High(Steps) do
    if Steps[k] >= m - 1e-9 then Exit(Steps[k] * e);
  Result := 10 * e;
end;

procedure TGasPVDiagram.Paint;
const
  NPts = 300;
var
  C: TCanvas;
  Pts: array of TPoint;
  k, x, y, x0, x1, y0, y1: Integer;
  T, V, P, pm, step, tick, vstep: Double;
  s: string;

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

  // the curve P(v) from the packing-0.5 volume to the axis end
  procedure Curve(Kind: Integer);
  var
    j: Integer;
    vv, pp: Double;
  begin
    for j := 0 to NPts do
    begin
      vv := VClose + (VMax - VClose) * j / NPts;
      case Kind of
        0: pp := T / vv;
        1: pp := T / vv * CarnahanStarling(Pi / 6 / vv);
      else
        pp := FAdiabat / Power(vv, 5 / 3);
      end;
      Pts[j] := Point(XOf(vv), YOf(pp));
    end;
    C.Polyline(Pts);
  end;

begin
  C := Canvas;
  C.Brush.Color := clWhite;
  C.FillRect(ClientRect);
  C.Font.Size := 8;
  if FHasState then T := FT else T := 1;
  // the pressure range: room above everything the gas has been through,
  // and at least the isotherm at a third of the axis
  pm := 3 * T / (VMax / 3);
  if FHasState then
  begin
    pm := Max(pm, 1.3 * FP);
    for k := 0 to FTrailLen - 1 do pm := Max(pm, 1.3 * FTrail[k].P);
  end;
  FPMax := NiceMax(pm);
  x0 := XOf(0); x1 := XOf(VMax); y0 := YOf(FPMax); y1 := YOf(0);

  // grid
  C.Pen.Color := RGBToColor(205, 205, 205);
  C.Pen.Style := psDot;
  vstep := NiceMax(VMax / 6);
  tick := vstep;
  while tick < VMax - 1e-9 do
  begin
    C.Line(XOf(tick), y0, XOf(tick), y1);
    tick := tick + vstep;
  end;
  step := FPMax / 5;
  tick := step;
  while tick < FPMax - 1e-9 do
  begin
    C.Line(x0, YOf(tick), x1, YOf(tick));
    tick := tick + step;
  end;
  C.Pen.Style := psSolid;

  C.ClipRect := Rect(x0, y0, x1 + 1, y1 + 1);
  C.Clipping := True;
  SetLength(Pts, NPts + 1);
  // spheres can't be packed beyond about 0.5 as a fluid: shade that
  C.Brush.Color := RGBToColor(235, 235, 235);
  C.Pen.Style := psClear;
  C.Rectangle(x0, y0, XOf(VClose), y1 + 1);
  C.Pen.Style := psSolid;

  // the hard-sphere (Carnahan-Starling) and ideal-gas isotherms at T
  C.Pen.Color := RGBToColor(20, 60, 160);
  C.Pen.Width := 2;
  Curve(1);
  C.Pen.Color := RGBToColor(200, 120, 0);
  C.Pen.Width := 1;
  C.Pen.Style := psDash;
  Curve(0);
  C.Pen.Style := psSolid;
  // the adiabat through the state before the piston moved
  if FAdiabat > 0 then
  begin
    C.Pen.Color := RGBToColor(150, 0, 150);
    C.Pen.Style := psDot;
    Curve(2);
    C.Pen.Style := psSolid;
  end;

  // the state: trail and red dot
  if FHasState then
  begin
    C.Pen.Color := RGBToColor(240, 140, 140);
    C.Pen.Width := 2;
    for k := 1 to FTrailLen - 1 do
      C.Line(XOf(FTrail[k - 1].V), YOf(FTrail[k - 1].P), XOf(FTrail[k].V), YOf(FTrail[k].P));
    C.Pen.Width := 1;
    x := XOf(FV); y := YOf(FP);
    C.Brush.Color := clRed;
    C.Pen.Color := RGBToColor(120, 0, 0);
    C.Ellipse(x - 6, y - 6, x + 7, y + 7);
  end;
  C.Clipping := False;

  // legend, top right
  x := x1 - 200; y := y0 + 8;
  C.Brush.Color := clWhite;
  C.Pen.Color := RGBToColor(180, 180, 180);
  if FAdiabat > 0 then k := 5 else k := 4;
  C.Rectangle(x - 6, y - 4, x1 - 6, y + 16 * k + 2);
  C.Pen.Width := 2; C.Pen.Color := RGBToColor(20, 60, 160);
  C.Line(x, y + 7, x + 24, y + 7);
  Lbl(x + 30, y, 'hard spheres (isotherm at T)', RGBToColor(20, 60, 160), False, False);
  Inc(y, 16);
  C.Pen.Width := 1; C.Pen.Color := RGBToColor(200, 120, 0); C.Pen.Style := psDash;
  C.Line(x, y + 7, x + 24, y + 7);
  Lbl(x + 30, y, 'ideal gas, PV = NkT', RGBToColor(190, 110, 0), False, False);
  Inc(y, 16);
  if FAdiabat > 0 then
  begin
    C.Pen.Color := RGBToColor(150, 0, 150); C.Pen.Style := psDot;
    C.Line(x, y + 7, x + 24, y + 7);
    Lbl(x + 30, y, 'adiabat, PV^(5/3) const', RGBToColor(150, 0, 150), False, False);
    Inc(y, 16);
  end;
  C.Pen.Style := psSolid;
  C.Pen.Width := 2; C.Pen.Color := RGBToColor(240, 140, 140);
  C.Line(x, y + 7, x + 24, y + 7);
  C.Pen.Width := 1;
  C.Brush.Color := clRed; C.Pen.Color := RGBToColor(120, 0, 0);
  C.Ellipse(x + 8, y + 3, x + 17, y + 12);
  Lbl(x + 30, y, 'the simulation, and its path', RGBToColor(160, 0, 0), False, False);
  Inc(y, 16);
  C.Brush.Color := RGBToColor(235, 235, 235); C.Pen.Color := RGBToColor(200, 200, 200);
  C.Rectangle(x + 4, y + 2, x + 20, y + 12);
  Lbl(x + 30, y, 'packed past φ = 0.5', RGBToColor(110, 110, 110), False, False);

  // frame and axes
  C.Brush.Style := bsClear;
  C.Pen.Color := clBlack;
  C.Rectangle(x0, y0, x1 + 1, y1 + 1);
  C.Brush.Style := bsSolid;
  tick := 0;
  while tick <= VMax + 1e-9 do
  begin
    x := XOf(tick);
    C.Line(x, y1, x, y1 + 4);
    Lbl(x, y1 + 12, FormatFloat('0.##', tick), clBlack);
    tick := tick + vstep;
  end;
  tick := 0;
  while tick <= FPMax + 1e-12 do
  begin
    y := YOf(tick);
    C.Line(x0 - 4, y, x0, y);
    s := FormatFloat('0.####', tick);
    Lbl(x0 - 8 - C.TextWidth(s) div 2, y, s, clBlack);
    tick := tick + step;
  end;
  Lbl((x0 + x1) div 2, y1 + 26, 'volume per sphere V/N (diameters³)', clBlack);
  C.Font.Orientation := 900;
  s := 'pressure (simulation units)';
  C.Brush.Style := bsClear;
  C.TextOut(4, (y0 + y1) div 2 + C.TextWidth(s) div 2, s);
  C.Font.Orientation := 0;
  Lbl((x0 + x1) div 2, 10, Format('Hard-sphere gas: P against V at T = %.3g', [T]), clBlack, True);
end;

end.
