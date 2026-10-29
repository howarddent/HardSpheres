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

unit usubstance;

{******************************************************************************
  The substances the simulation can be mapped onto - nitrous oxide, carbon
  dioxide and water - with their critical and triple points, saturation
  curves and a few lines of teaching notes.

  Data:
    - liquid-vapour curve, triple point to critical point, with the
      saturated liquid and vapour densities: the NIST Chemistry WebBook's
      reference equations of state (n2odata.inc, co2data.inc,
      h2odata.inc), the pressure interpolated as ln P linear in 1/T, the
      liquid density linearly and the vapour density's log linearly in T;
    - solid-vapour (sublimation) curve: Clausius-Clapeyron from the triple
      point, with dH_sub/R chosen to pass through a known point (N2O: from
      dH_vap + dH_fus; CO2: dry ice at -78.5 C, 1 atm; water: 51.06 kJ/mol);
    - solid-liquid (melting) curve: a straight Clapeyron line from the
      triple point (water's slopes backwards: ice melts under pressure).
*******************************************************************************}

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Math;

const
  RBarLitre = 0.0831446;      // L bar / (mol K)

type
  TSatRow = array[0..3] of Double;   // T K, P bar, liquid and vapour density kg/m3

  TSubstance = class
  public
    Name: string;              // 'nitrous oxide'
    Formula: string;           // 'N2O'
    MolarMass: Double;         // g/mol
    Tc, Pc: Double;            // K, bar
    Tt, Pt: Double;            // K, bar
    SubSlope: Double;          // K: dH_sub/R
    MeltSlope: Double;         // bar/K
    Sat: array of TSatRow;
    // shown while liquid and vapour coexist: what the vapour pressure means
    // for this substance in practice
    SaturatedNote: string;
    // coexistence pressures, bar
    function VapourPressure(TK: Double): Double;      // liquid-vapour, Tt..Tc
    function SublimationPressure(TK: Double): Double; // solid-vapour, below Tt
    function CoexistencePressure(TK: Double): Double; // whichever applies
    function MeltingTemperature(Pbar: Double): Double;
    // the phase at (TK, Pbar), as a word
    function PhaseName(TK, Pbar: Double): string;
    // saturated liquid and vapour molar volumes at TK (L/mol), Tt..Tc
    procedure SaturatedVolumes(TK: Double; out Vl, Vv: Double);
    function CriticalVolume: Double;
  end;

const
  SubstanceCount = 3;

var
  Substances: array[0..SubstanceCount - 1] of TSubstance;

// kelvin as °C
function Celsius(TK: Double): Double;
// a pressure for the text: bar, and kPa as well when it is under 1 bar
function PressureText(Pbar: Double): string;

implementation

{$I n2odata.inc}
{$I co2data.inc}
{$I h2odata.inc}

function Celsius(TK: Double): Double;
begin
  Result := TK - 273.15;
end;

function PressureText(Pbar: Double): string;
begin
  if Pbar < 1 then Result := Format('%.3g bar (%.3g kPa)', [Pbar, 100 * Pbar])
  else Result := Format('%.3g bar', [Pbar]);
end;

function TSubstance.VapourPressure(TK: Double): Double;
var
  i, n: Integer;
  f: Double;
begin
  n := Length(Sat);
  if TK <= Sat[0, 0] then Exit(Sat[0, 1]);
  if TK >= Sat[n - 1, 0] then Exit(Sat[n - 1, 1]);
  i := 0;
  while Sat[i + 1, 0] < TK do Inc(i);
  f := (1 / TK - 1 / Sat[i, 0]) / (1 / Sat[i + 1, 0] - 1 / Sat[i, 0]);
  Result := Exp(Ln(Sat[i, 1]) + f * (Ln(Sat[i + 1, 1]) - Ln(Sat[i, 1])));
end;

function TSubstance.SublimationPressure(TK: Double): Double;
begin
  Result := Pt * Exp(-SubSlope * (1 / TK - 1 / Tt));
end;

function TSubstance.CoexistencePressure(TK: Double): Double;
begin
  if TK < Tt then Result := SublimationPressure(TK)
  else Result := VapourPressure(TK);
end;

function TSubstance.MeltingTemperature(Pbar: Double): Double;
begin
  Result := Tt + (Pbar - Pt) / MeltSlope;
end;

function TSubstance.PhaseName(TK, Pbar: Double): string;
var
  solid: Boolean;
begin
  // solid: above the sublimation curve below Tt, or on the cold side of
  // the melting line above Pt (for water that side lies below Tt)
  solid := (Pbar > Pt) and (TK < MeltingTemperature(Pbar));
  if TK < Tt then solid := solid or ((Pbar > SublimationPressure(TK)) and (Pbar <= Pt));
  if solid then Result := 'solid'
  else if TK >= Tc then
  begin
    if Pbar >= Pc then Result := 'supercritical fluid' else Result := 'gas';
  end
  else if Pbar > CoexistencePressure(TK) then Result := 'liquid'
  else Result := 'vapour';
end;

procedure TSubstance.SaturatedVolumes(TK: Double; out Vl, Vv: Double);
var
  i, n: Integer;
  f, rl, rv: Double;
begin
  n := Length(Sat);
  if TK <= Sat[0, 0] then i := 0
  else if TK >= Sat[n - 1, 0] then i := n - 2
  else
  begin
    i := 0;
    while Sat[i + 1, 0] < TK do Inc(i);
  end;
  f := EnsureRange((TK - Sat[i, 0]) / (Sat[i + 1, 0] - Sat[i, 0]), 0, 1);
  rl := Sat[i, 2] + f * (Sat[i + 1, 2] - Sat[i, 2]);
  // the vapour density grows roughly exponentially: interpolate its log
  rv := Exp(Ln(Sat[i, 3]) + f * (Ln(Sat[i + 1, 3]) - Ln(Sat[i, 3])));
  Vl := MolarMass / rl;     // g/mol / (kg/m3) = L/mol
  Vv := MolarMass / rv;
end;

function TSubstance.CriticalVolume: Double;
begin
  Result := MolarMass / Sat[High(Sat), 2];
end;

function Make(const AName, AFormula: string; AM, ATc, APc, ATt, APt, ASub, AMelt: Double;
  Table: Pointer; ACount: Integer; const ANote: string): TSubstance;
type
  PSatRow = ^TSatRow;
var
  k: Integer;
begin
  Result := TSubstance.Create;
  with Result do
  begin
    Name := AName; Formula := AFormula; MolarMass := AM;
    Tc := ATc; Pc := APc; Tt := ATt; Pt := APt;
    SubSlope := ASub; MeltSlope := AMelt;
    SetLength(Sat, ACount);
    for k := 0 to ACount - 1 do Sat[k] := PSatRow(Table)[k];
    SaturatedNote := ANote;
  end;
end;

procedure CreateSubstances;
begin
  Substances[0] := Make('nitrous oxide', 'N2O', 44.013, 309.52, 72.447, 182.33, 0.8784,
    2782,     // (16.59 + 6.54) kJ/mol / R
    44,       // approximate, solid density ~1600 kg/m3
    @N2OSat, N2OSatCount,
    'This is why an N2O cylinder reads 50.5 bar at 20 °C until it is '
    + 'about 90 % empty (only vapour left): the gauge shows the vapour '
    + 'pressure, not the contents. The cylinder is weighed instead.');
  Substances[1] := Make('carbon dioxide', 'CO2', 44.010, 304.128, 73.773, 216.592, 5.1796,
    3143,     // through dry ice's sublimation at 194.69 K, 1.01325 bar
    46.8,     // Span-Wagner's melting curve at the triple point
    @CO2Sat, CO2SatCount,
    'A CO2 cylinder reads 57.3 bar at 20 °C while any liquid is left, '
    + 'so, as with N2O, it is weighed to tell how full it is. CO2''s triple '
    + 'point is at 5.2 bar, above 1 atm: at atmospheric pressure there is '
    + 'no liquid, and dry ice turns straight to vapour at −78.5 °C.');
  Substances[2] := Make('water', 'H2O', 18.015, 647.096, 220.64, 273.16, 0.0061165,
    6141,     // 51.06 kJ/mol / R
    -135,     // ice is less dense than water: melts under pressure
    @H2OSat, H2OSatCount,
    'Water''s vapour pressure is 2.3 kPa at 20 °C and 6.3 kPa at 37 °C: '
    + 'air in the lungs is saturated, so 6.3 kPa of it is water vapour. '
    + 'A liquid boils when its vapour pressure reaches the pressure around '
    + 'it, which for water at 1 atm (101.3 kPa) is at 100 °C.');
end;

procedure FreeSubstances;
var
  k: Integer;
begin
  for k := 0 to SubstanceCount - 1 do FreeAndNil(Substances[k]);
end;

initialization
  CreateSubstances;
finalization
  FreeSubstances;
end.
