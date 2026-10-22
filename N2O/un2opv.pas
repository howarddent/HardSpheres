unit un2opv;

{******************************************************************************
  The pressure-volume diagram of nitrous oxide at the current temperature,
  drawn on a canvas, with a red dot for the simulation's state and a trail
  of where it has been (as the piston moves, the simulated isotherm).

  Drawn:
    - the liquid-vapour dome: the saturated liquid and vapour molar volumes
      from the triple point to the critical point (NIST, n2odata.inc via
      un2ophase), meeting at the critical point (0.0974 L/mol, 72.4 bar);
    - at the current temperature T, below the critical temperature: the
      isotherm's flat part, the tie line at the vapour pressure from the
      saturated liquid to the saturated vapour volume, and a steep liquid
      branch above it (schematic: a liquid is barely compressible);
    - the ideal-gas isotherm P = RT/V at T, dashed, for comparison.
  Volume is the molar volume in L/mol on a log scale, so the liquid
  (~0.04 L/mol) and the vapour (litres) both show; pressure is linear in
  bar, its range set from the vapour pressure at T.
*******************************************************************************}

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Math, Controls, Graphics, un2ophase;

const
  N2OMolarMass = 44.013;      // g/mol
  RBarLitre = 0.0831446;      // L bar / (mol K)

// saturated liquid and vapour molar volumes at TK (L/mol), Tt..Tc
procedure N2OSaturatedVolumes(TK: Double; out Vl, Vv: Double);
function N2OCriticalVolume: Double;

type
  TPVPoint = record
    V, P: Double;
  end;

  TN2OPVDiagram = class(TCustomControl)
  private
    FHasState: Boolean;
    FT, FV, FP: Double;           // K, L/mol, bar
    FLabel: string;
    FTrail: array of TPVPoint;
    FTrailLen: Integer;
    FPMax: Double;
    function XOf(V: Double): Integer;
    function YOf(P: Double): Integer;
  protected
    procedure Paint; override;
  public
    VMin, VMax: Double;           // L/mol
    constructor Create(AOwner: TComponent); override;
    procedure SetState(TK, V, P: Double; const ALabel: string);
    procedure ClearState;
  end;

implementation

const
  PlotL = 52; PlotR = 12; PlotT = 24; PlotB = 48;
  TrailMax = 600;

procedure N2OSaturatedVolumes(TK: Double; out Vl, Vv: Double);
var
  i: Integer;
  f, rl, rv: Double;
begin
  if TK <= N2OSat[0, 0] then i := 0
  else if TK >= N2OSat[N2OSatCount - 1, 0] then i := N2OSatCount - 2
  else
  begin
    i := 0;
    while N2OSat[i + 1, 0] < TK do Inc(i);
  end;
  f := EnsureRange((TK - N2OSat[i, 0]) / (N2OSat[i + 1, 0] - N2OSat[i, 0]), 0, 1);
  rl := N2OSat[i, 2] + f * (N2OSat[i + 1, 2] - N2OSat[i, 2]);
  // the vapour density grows roughly exponentially: interpolate its log
  rv := Exp(Ln(N2OSat[i, 3]) + f * (Ln(N2OSat[i + 1, 3]) - Ln(N2OSat[i, 3])));
  Vl := N2OMolarMass / rl;     // g/mol / (kg/m3) = L/mol
  Vv := N2OMolarMass / rv;
end;

function N2OCriticalVolume: Double;
begin
  Result := N2OMolarMass / N2OSat[N2OSatCount - 1, 2];
end;

constructor TN2OPVDiagram.Create(AOwner: TComponent);
begin
  inherited Create(AOwner);
  VMin := 0.01; VMax := 50;
  FPMax := 100;
  SetLength(FTrail, TrailMax);
  DoubleBuffered := True;
  Color := clWhite;
end;

function TN2OPVDiagram.XOf(V: Double): Integer;
var
  f: Double;
begin
  V := EnsureRange(V, VMin * 0.5, VMax * 2);
  f := (Ln(V) - Ln(VMin)) / (Ln(VMax) - Ln(VMin));
  Result := PlotL + Round(f * (ClientWidth - PlotL - PlotR));
end;

function TN2OPVDiagram.YOf(P: Double): Integer;
begin
  P := EnsureRange(P, -0.5 * FPMax, 2 * FPMax);
  Result := ClientHeight - PlotB - Round(P / FPMax * (ClientHeight - PlotT - PlotB));
end;

procedure TN2OPVDiagram.SetState(TK, V, P: Double; const ALabel: string);
begin
  if FHasState and (Abs(TK - FT) < 1e-9) and (Abs(V - FV) < 1e-12) and (Abs(P - FP) < 1e-9) and (ALabel = FLabel) then Exit;
  FHasState := True;
  FT := TK; FV := V; FP := P; FLabel := ALabel;
  if FTrailLen = TrailMax then
  begin
    Move(FTrail[1], FTrail[0], (TrailMax - 1) * SizeOf(TPVPoint));
    Dec(FTrailLen);
  end;
  FTrail[FTrailLen].V := V;
  FTrail[FTrailLen].P := P;
  Inc(FTrailLen);
  Invalidate;
end;

procedure TN2OPVDiagram.ClearState;
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

procedure TN2OPVDiagram.Paint;
var
  C: TCanvas;
  Pts: array of TPoint;
  n, k, x, y, x0, x1, y0, y1, xl, xv, yp: Integer;
  T, Vl, Vv, Psat, V, P, step, tick: Double;
  s: string;
  below: Boolean;

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
  if FHasState then T := FT else T := 230;
  below := T < N2OTc;
  // the pressure range: room above the vapour pressure (or the state)
  if below then Psat := N2OVapourPressure(Max(T, N2OTt)) else Psat := N2OPc;
  P := 1.6 * Psat;
  if FHasState then P := Max(P, 1.15 * FP);
  FPMax := NiceMax(P);
  x0 := XOf(VMin); x1 := XOf(VMax); y0 := YOf(FPMax); y1 := YOf(0);
  C.ClipRect := Rect(x0, y0, x1 + 1, y1 + 1);
  C.Clipping := True;
  SetLength(Pts, 2 * N2OSatCount + 2);

  // the dome: liquid side up to the critical point, then the vapour side down
  n := 0;
  for k := 0 to N2OSatCount - 1 do
  begin
    Pts[n] := Point(XOf(N2OMolarMass / N2OSat[k, 2]), YOf(N2OSat[k, 1]));
    Inc(n);
  end;
  for k := N2OSatCount - 2 downto 0 do
  begin
    Pts[n] := Point(XOf(N2OMolarMass / N2OSat[k, 3]), YOf(N2OSat[k, 1]));
    Inc(n);
  end;
  C.Brush.Color := RGBToColor(226, 238, 252);
  C.Pen.Style := psClear;
  C.Polygon(Pts[0..n - 1]);
  C.Pen.Style := psSolid;
  C.Pen.Color := RGBToColor(20, 60, 160);
  C.Pen.Width := 2;
  C.Polyline(Pts[0..n - 1]);
  C.Pen.Width := 1;

  // grid
  C.Pen.Color := RGBToColor(205, 205, 205);
  C.Pen.Style := psDot;
  V := 0.01;
  while V <= VMax * 1.001 do
  begin
    if V >= VMin then C.Line(XOf(V), y0, XOf(V), y1);
    V := V * 10;
  end;
  step := FPMax / 5;
  tick := step;
  while tick < FPMax - 1e-9 do
  begin
    C.Line(x0, YOf(tick), x1, YOf(tick));
    tick := tick + step;
  end;
  C.Pen.Style := psSolid;

  // the ideal-gas isotherm at T
  C.Pen.Color := RGBToColor(200, 120, 0);
  C.Pen.Style := psDash;
  C.Pen.Width := 1;
  n := 0;
  SetLength(Pts, 201);
  for k := 0 to 200 do
  begin
    V := Exp(Ln(VMin) + (Ln(VMax) - Ln(VMin)) * k / 200);
    Pts[n] := Point(XOf(V), YOf(RBarLitre * T / V));
    Inc(n);
  end;
  C.Polyline(Pts[0..n - 1]);
  C.Pen.Style := psSolid;

  // the real isotherm's flat part and liquid branch at T
  if below and (T >= N2OTt) then
  begin
    N2OSaturatedVolumes(T, Vl, Vv);
    xl := XOf(Vl); xv := XOf(Vv); yp := YOf(Psat);
    C.Pen.Color := RGBToColor(0, 120, 60);
    C.Pen.Width := 3;
    C.Line(xl, yp, xv, yp);
    // the liquid barely compresses: a near-vertical branch (schematic)
    C.Line(xl, yp, XOf(Vl * 0.97), y0);
    C.Pen.Width := 1;
    C.Brush.Color := RGBToColor(0, 120, 60);
    C.Pen.Color := RGBToColor(0, 120, 60);
    C.Ellipse(xl - 3, yp - 3, xl + 4, yp + 4);
    C.Ellipse(xv - 3, yp - 3, xv + 4, yp + 4);
    s := Format('saturated: %.3g bar', [Psat]);
    Lbl((xl + xv) div 2 - C.TextWidth(s) div 2, yp - C.TextHeight(s) - 3, s, RGBToColor(0, 110, 50), True, False);
  end;

  // critical point
  x := XOf(N2OCriticalVolume); y := YOf(N2OPc);
  C.Brush.Color := clBlack; C.Pen.Color := clBlack;
  C.Ellipse(x - 4, y - 4, x + 5, y + 5);

  // labels
  Lbl(XOf(0.6), YOf(0.3 * Psat) - 6, 'liquid + vapour', RGBToColor(20, 60, 160), True);
  Lbl(XOf(0.019), YOf(0.85 * FPMax), 'liquid', RGBToColor(20, 60, 160), True);
  if below then s := 'vapour' else s := 'gas';
  Lbl(XOf(15), YOf(0.85 * FPMax), s, RGBToColor(150, 110, 0), True);
  s := 'ideal gas, PV = RT';
  V := RBarLitre * T / (0.6 * FPMax);
  Lbl(XOf(V) + 6, YOf(0.6 * FPMax) - 4, s, RGBToColor(190, 110, 0), False, False);

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

  // frame and axes
  C.Brush.Style := bsClear;
  C.Pen.Color := clBlack;
  C.Rectangle(x0, y0, x1 + 1, y1 + 1);
  C.Brush.Style := bsSolid;
  V := 0.01;
  while V <= VMax * 1.001 do
  begin
    if V >= VMin then
    begin
      C.Line(XOf(V), y1, XOf(V), y1 + 4);
      if V < 1 then s := FormatFloat('0.###', V) else s := IntToStr(Round(V));
      Lbl(XOf(V), y1 + 12, s, clBlack);
    end;
    for k := 2 to 9 do
      if (V * k > VMin) and (V * k < VMax) then C.Line(XOf(V * k), y1, XOf(V * k), y1 + 2);
    V := V * 10;
  end;
  tick := 0;
  while tick <= FPMax + 1e-9 do
  begin
    y := YOf(tick);
    C.Line(x0 - 4, y, x0, y);
    s := FormatFloat('0.##', tick);
    Lbl(x0 - 8 - C.TextWidth(s) div 2, y, s, clBlack);
    tick := tick + step;
  end;
  Lbl((x0 + x1) div 2, y1 + 26, 'molar volume (L/mol, log scale)', clBlack);
  C.Font.Orientation := 900;
  s := 'pressure (bar)';
  C.Brush.Style := bsClear;
  C.TextOut(4, (y0 + y1) div 2 + C.TextWidth(s) div 2, s);
  C.Font.Orientation := 0;
  if FHasState then
    s := Format('Nitrous oxide: isotherm at %.1f °C', [FT - 273.15])
  else
    s := 'Nitrous oxide: pressure and volume';
  Lbl((x0 + x1) div 2, 10, s, clBlack, True);
end;

end.
