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

unit ugassim;

{******************************************************************************
  Hard-sphere gas on the GPU (OpenCL), time-stepped.

  Units: sphere diameter 1, mass 1, Boltzmann's constant 1, so the
  temperature is the mean kinetic energy per degree of freedom,
  T = 2 KE/(3N), and the rms speed is sqrt(3T).

  The box is [0, Lx] x [0, Ly] x [0, Lz]; the top face, z = Lz, is a
  piston, so only the height changes with the volume (the floor and side
  walls stay put). It moves towards LzTarget at PistonSpeed, so a slow move is a
  reversible adiabatic compression or expansion (spheres bouncing off a
  receding wall lose energy, off an advancing one gain it), and a fast
  one is not. A thermostat, when on, rescales the velocities each frame
  to hold T (isothermal).

  Each frame is SimTimePerFrame of simulated time, in substeps short
  enough that a typical sphere moves StepFrac diameters per substep (see
  gas.cl for the kernels). A hot gas therefore takes more substeps per
  frame and is seen to move faster.

  Measured:
    - Temperature: from the kinetic energy, averaged over the frame;
    - Pressure: momentum given to all six walls per unit time and area,
      PistonPressure the same for the piston alone, both smoothed over
      about PressureFrames frames;
    - the collision rate (pairs per sphere per unit time).
  The areas and volume are those open to the sphere centres, (L - 1)
  in each direction, so an ideal gas gives exactly P = N T / V.

  Positions live in a buffer that is, when the GUI can share it, the GL
  vertex buffer it draws (cl_khr_gl_sharing): nothing crosses the bus per
  frame. Otherwise the positions are read back to HostPos for the GUI to
  upload.
*******************************************************************************}

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Math, clgl;

type
  TGasSim = class
  private
    FCtx: cl_context;
    FQueue: cl_command_queue;
    FDevice: cl_device_id;
    FProg: cl_program;
    KMove, KClear, KBin, KPairs, KApply, KScale: cl_kernel;
    FPos, FVelA, FVelB, FCorr, FPartner, FCounts, FCells, FOverflow, FNColl, FPartial: cl_mem;
    FShared: Boolean;             // FPos is the GL buffer
    FDeviceName: string;
    FN, FGroups, FLocal: Integer;
    FGrid: array[0..3] of LongInt;
    FNCells: Integer;
    FCellSize: Single;
    FMaxP: LongInt;
    FLx, FLy, FLz, FLzMax: Single;
    FT0: Double;
    FTime: Double;
    FTemp, FPress, FPistonPress, FCollRate: Double;
    FHaveStats: Boolean;
    FOverflowCount: LongInt;
    FSubsteps: Integer;
    FHostPartial: array of Single;
    procedure BuildProgram(const Src: string);
    procedure ReleaseBuffers;
    procedure Enqueue(K: cl_kernel; Count: Integer; Local: Integer = 0);
    procedure Arg(K: cl_kernel; Index: Integer; Size: PtrUInt; Value: Pointer);
    procedure ScaleVelocities(s: Single);
  public
    HostPos: array of Single;     // x, y, z, |v| per sphere, when not shared
    LzTarget: Double;
    PistonSpeed: Double;          // diameters per unit time
    Thermostat: Boolean;
    TargetT: Double;
    StepFrac: Double;             // rms displacement per substep, diameters
    // run on an integrated GPU (memory shared with the CPU). Off by
    // default: on a laptop that GPU drives the display, and a long kernel
    // on it can hang the desktop past a failed GPU reset
    AllowIntegrated: Boolean;
    SimTimePerFrame: Double;
    PressureFrames: Integer;
    constructor Create;
    destructor Destroy; override;
    // an OpenCL context: sharing the current GL context when GLContext and
    // GLDisplay are given (GLX) and the driver allows, otherwise on the
    // GPU with the most compute units. Raises on failure.
    procedure InitCL(const KernelFile: string; GLContext, GLDisplay: Pointer);
    // N spheres at packing fraction Phi in a cube, temperature T; equal
    // speeds in random directions, or Maxwell-Boltzmann. GLBuffer: a GL
    // vertex buffer of N float4 to draw from (0: none); if the context
    // shares it, the positions are written into it here, otherwise the
    // caller uploads HostPos. LzMaxFactor: how far the piston may withdraw, in units of
    // the starting length.
    procedure Setup(N: Integer; Phi, T: Double; EqualSpeeds: Boolean;
      GLBuffer: LongWord; LzMaxFactor: Double = 2.0);
    // one frame; the caller must have finished with the GL buffer (glFinish)
    procedure Frame;
    // sets the temperature at once, by rescaling every velocity
    procedure SetTemperature(T: Double);
    // velocities, x, y, z, 0 per sphere
    procedure ReadVelocities(var V: array of Single);
    function Volume: Double;      // open to the centres
    function PackingFraction: Double;
    property Count: Integer read FN;
    property Lx: Single read FLx;
    property Ly: Single read FLy;
    property Lz: Single read FLz;
    property LzMax: Single read FLzMax;
    property Time: Double read FTime;
    property Temperature: Double read FTemp;
    property Pressure: Double read FPress;
    property PistonPressure: Double read FPistonPress;
    property CollisionRate: Double read FCollRate;
    property Overflows: LongInt read FOverflowCount;
    property CellSize: Single read FCellSize;
    property CellCapacity: LongInt read FMaxP;
    property Substeps: Integer read FSubsteps;
    property Shared: Boolean read FShared;
    property DeviceName: string read FDeviceName;
    property Context: cl_context read FCtx;
  end;

// Carnahan-Starling compressibility factor PV/NkT of hard spheres
function CarnahanStarling(Phi: Double): Double;

implementation

const
  LocalSize = 256;
  CellBudget = 8000000;

function CarnahanStarling(Phi: Double): Double;
begin
  Result := (1 + Phi + Phi * Phi - Phi * Phi * Phi) / Power(1 - Phi, 3);
end;

constructor TGasSim.Create;
begin
  StepFrac := 0.02;
  SimTimePerFrame := 0.05;
  PistonSpeed := 0.5;
  TargetT := 1;
  PressureFrames := 20;
end;

destructor TGasSim.Destroy;
begin
  ReleaseBuffers;
  if KMove <> nil then clReleaseKernel(KMove);
  if KClear <> nil then clReleaseKernel(KClear);
  if KBin <> nil then clReleaseKernel(KBin);
  if KPairs <> nil then clReleaseKernel(KPairs);
  if KApply <> nil then clReleaseKernel(KApply);
  if KScale <> nil then clReleaseKernel(KScale);
  if FProg <> nil then clReleaseProgram(FProg);
  if FQueue <> nil then clReleaseCommandQueue(FQueue);
  if FCtx <> nil then clReleaseContext(FCtx);
  inherited Destroy;
end;

procedure TGasSim.ReleaseBuffers;
  procedure R(var M: cl_mem);
  begin
    if M <> nil then clReleaseMemObject(M);
    M := nil;
  end;
begin
  R(FPos); R(FVelA); R(FVelB); R(FCorr); R(FPartner); R(FCounts); R(FCells);
  R(FOverflow); R(FNColl); R(FPartial);
end;

procedure TGasSim.InitCL(const KernelFile: string; GLContext, GLDisplay: Pointer);
var
  Msg, Src: string;
  np, nd, cu, bestCU: cl_uint;
  Plats: array of cl_platform_id;
  Devs: array of cl_device_id;
  p, d: Integer;
  Props: array[0..6] of cl_context_properties;
  GetGLInfo: TclGetGLContextInfoKHR;
  Dev: cl_device_id;
  Err: cl_int;
  bestP: cl_platform_id;
  SL: TStringList;
  Skipped: string;

  // an integrated GPU, refused unless AllowIntegrated; its name is noted
  function Refused(D: cl_device_id): Boolean;
  var
    unified: cl_uint;
    s: string;
  begin
    unified := 0;
    clGetDeviceInfo(D, CL_DEVICE_HOST_UNIFIED_MEMORY, SizeOf(unified), @unified, nil);
    Result := (unified <> 0) and not AllowIntegrated;
    if Result then
    begin
      s := DeviceString(D, CL_DEVICE_NAME);
      if Pos(s, Skipped) = 0 then
      begin
        if Skipped <> '' then Skipped := Skipped + ', ';
        Skipped := Skipped + s;
      end;
    end;
  end;

begin
  Skipped := '';
  if not LoadOpenCL(Msg) then raise Exception.Create(Msg);
  np := 0;
  if (clGetPlatformIDs(0, nil, @np) <> CL_SUCCESS) or (np = 0) then
    raise Exception.Create('no OpenCL platform');
  SetLength(Plats, np);
  clGetPlatformIDs(np, @Plats[0], nil);

  // first choice: the device driving the current GL context, shared
  FShared := False;
  if (GLContext <> nil) and (GLDisplay <> nil) and Assigned(clGetExtensionFunctionAddressForPlatform)
     and Assigned(clCreateFromGLBuffer) then
    for p := 0 to np - 1 do
    begin
      if Pos('cl_khr_gl_sharing', PlatformString(Plats[p], CL_PLATFORM_EXTENSIONS)) = 0 then Continue;
      GetGLInfo := TclGetGLContextInfoKHR(clGetExtensionFunctionAddressForPlatform(Plats[p], 'clGetGLContextInfoKHR'));
      if not Assigned(GetGLInfo) then Continue;
      Props[0] := CL_GL_CONTEXT_KHR;   Props[1] := cl_context_properties(GLContext);
      Props[2] := CL_GLX_DISPLAY_KHR;  Props[3] := cl_context_properties(GLDisplay);
      Props[4] := CL_CONTEXT_PLATFORM; Props[5] := cl_context_properties(Plats[p]);
      Props[6] := 0;
      Dev := nil;
      if (GetGLInfo(@Props[0], CL_CURRENT_DEVICE_FOR_GL_CONTEXT_KHR, SizeOf(Dev), @Dev, nil) <> CL_SUCCESS)
         or (Dev = nil) or Refused(Dev) then Continue;
      FCtx := clCreateContext(@Props[0], 1, @Dev, nil, nil, @Err);
      if (Err = CL_SUCCESS) and (FCtx <> nil) then
      begin
        FDevice := Dev;
        FShared := True;
        Break;
      end;
      FCtx := nil;
    end;

  // otherwise the GPU with the most compute units, positions copied
  if FCtx = nil then
  begin
    bestCU := 0;
    bestP := nil;
    for p := 0 to np - 1 do
    begin
      nd := 0;
      if (clGetDeviceIDs(Plats[p], CL_DEVICE_TYPE_GPU, 0, nil, @nd) <> CL_SUCCESS) or (nd = 0) then Continue;
      SetLength(Devs, nd);
      clGetDeviceIDs(Plats[p], CL_DEVICE_TYPE_GPU, nd, @Devs[0], nil);
      for d := 0 to nd - 1 do
      begin
        if Refused(Devs[d]) then Continue;
        cu := 0;
        clGetDeviceInfo(Devs[d], CL_DEVICE_MAX_COMPUTE_UNITS, SizeOf(cu), @cu, nil);
        if cu > bestCU then
        begin
          bestCU := cu;
          FDevice := Devs[d];
          bestP := Plats[p];
        end;
      end;
    end;
    if (FDevice = nil) and (Skipped <> '') then
      raise Exception.Create('only an integrated GPU (' + Skipped + ') is available to OpenCL, '
        + 'and it may be driving the display: is the NVIDIA driver loaded (nvidia-smi)? '
        + 'Pass --allow-igpu to run on it anyway');
    if FDevice = nil then raise Exception.Create('no OpenCL GPU device');
    Props[0] := CL_CONTEXT_PLATFORM; Props[1] := cl_context_properties(bestP); Props[2] := 0;
    FCtx := clCreateContext(@Props[0], 1, @FDevice, nil, nil, @Err);
    CheckCL(Err, 'clCreateContext');
  end;
  FDeviceName := DeviceString(FDevice, CL_DEVICE_NAME);
  FQueue := clCreateCommandQueue(FCtx, FDevice, 0, @Err);
  CheckCL(Err, 'clCreateCommandQueue');

  SL := TStringList.Create;
  try
    SL.LoadFromFile(KernelFile);
    Src := SL.Text;
  finally
    SL.Free;
  end;
  BuildProgram(Src);
end;

procedure TGasSim.BuildProgram(const Src: string);
var
  Err: cl_int;
  P: PAnsiChar;
  Log: string;
  n: PtrUInt;
  function K(const Name: string): cl_kernel;
  var
    E: cl_int;
  begin
    Result := clCreateKernel(FProg, PAnsiChar(Name), @E);
    CheckCL(E, 'clCreateKernel ' + Name);
  end;
begin
  P := PAnsiChar(Src);
  FProg := clCreateProgramWithSource(FCtx, 1, @P, nil, @Err);
  CheckCL(Err, 'clCreateProgramWithSource');
  Err := clBuildProgram(FProg, 1, @FDevice, '-cl-fast-relaxed-math', nil, nil);
  if Err <> CL_SUCCESS then
  begin
    n := 0;
    clGetProgramBuildInfo(FProg, FDevice, CL_PROGRAM_BUILD_LOG, 0, nil, @n);
    SetLength(Log, n);
    if n > 0 then clGetProgramBuildInfo(FProg, FDevice, CL_PROGRAM_BUILD_LOG, n, @Log[1], nil);
    raise Exception.Create('gas.cl did not build:'#10 + Log);
  end;
  KMove := K('move');
  KClear := K('clear');
  KBin := K('bin');
  KPairs := K('pairs');
  KApply := K('apply');
  KScale := K('scale');
end;

procedure TGasSim.Arg(K: cl_kernel; Index: Integer; Size: PtrUInt; Value: Pointer);
begin
  CheckCL(clSetKernelArg(K, Index, Size, Value), Format('clSetKernelArg %d', [Index]));
end;

procedure TGasSim.Enqueue(K: cl_kernel; Count: Integer; Local: Integer);
var
  g, l: PtrUInt;
begin
  if Local > 0 then
  begin
    l := Local;
    g := ((Count + Local - 1) div Local) * Local;
    CheckCL(clEnqueueNDRangeKernel(FQueue, K, 1, nil, @g, @l, 0, nil, nil), 'clEnqueueNDRangeKernel');
  end
  else
  begin
    g := Count;
    CheckCL(clEnqueueNDRangeKernel(FQueue, K, 1, nil, @g, nil, 0, nil, nil), 'clEnqueueNDRangeKernel');
  end;
end;

// a normal deviate (Box-Muller)
function Gauss: Double;
var
  u1, u2: Double;
begin
  repeat
    u1 := Random;
  until u1 > 1e-300;
  u2 := Random;
  Result := Sqrt(-2 * Ln(u1)) * Cos(2 * Pi * u2);
end;

procedure TGasSim.Setup(N: Integer; Phi, T: Double; EqualSpeeds: Boolean;
  GLBuffer: LongWord; LzMaxFactor: Double);
var
  L, a, jit, s, ke, vm: Double;
  ns, i, ix, iy, iz, k, slots: Integer;
  V: array of Single;
  cm: array[0..2] of Double;
  Err: cl_int;
  Zero: LongInt;
  perm: array of Integer;
begin
  ReleaseBuffers;
  FN := N;
  // the box: a cube of side L with N pi/6 = Phi (L - 1)^3, centres in L - 1
  L := Power(N * Pi / 6 / Phi, 1 / 3) + 1;
  FLx := L; FLy := L; FLz := L;
  FLzMax := L * LzMaxFactor;
  LzTarget := L;
  FT0 := T;
  TargetT := T;
  FTime := 0;
  FHaveStats := False;
  FOverflowCount := 0;

  // a simple cubic lattice of at least N sites inside [0.5, L - 0.5],
  // jittered without overlap, the sites taken in random order
  ns := Ceil(Power(N, 1 / 3) - 1e-9);
  a := (L - 1) / ns;
  if a < 1.0 then raise Exception.CreateFmt('packing fraction %.3g too high to start from a lattice', [Phi]);
  jit := 0.45 * (a - 1);
  slots := ns * ns * ns;
  SetLength(perm, slots);
  for i := 0 to slots - 1 do perm[i] := i;
  for i := slots - 1 downto 1 do
  begin
    k := Random(i + 1);
    ix := perm[i]; perm[i] := perm[k]; perm[k] := ix;
  end;
  SetLength(HostPos, 4 * N);
  SetLength(V, 4 * N);
  cm[0] := 0; cm[1] := 0; cm[2] := 0;
  for i := 0 to N - 1 do
  begin
    k := perm[i];
    ix := k mod ns; iy := (k div ns) mod ns; iz := k div (ns * ns);
    HostPos[4 * i] := 0.5 + (ix + 0.5) * a + jit * (2 * Random - 1);
    HostPos[4 * i + 1] := 0.5 + (iy + 0.5) * a + jit * (2 * Random - 1);
    HostPos[4 * i + 2] := 0.5 + (iz + 0.5) * a + jit * (2 * Random - 1);
    if EqualSpeeds then
    begin
      // a random direction, speed sqrt(3T)
      repeat
        V[4 * i] := 2 * Random - 1; V[4 * i + 1] := 2 * Random - 1; V[4 * i + 2] := 2 * Random - 1;
        s := Sqr(V[4 * i]) + Sqr(V[4 * i + 1]) + Sqr(V[4 * i + 2]);
      until (s > 1e-4) and (s <= 1);
      s := Sqrt(3 * T / s);
      V[4 * i] := V[4 * i] * s; V[4 * i + 1] := V[4 * i + 1] * s; V[4 * i + 2] := V[4 * i + 2] * s;
    end
    else
    begin
      s := Sqrt(T);
      V[4 * i] := Gauss * s; V[4 * i + 1] := Gauss * s; V[4 * i + 2] := Gauss * s;
    end;
    V[4 * i + 3] := 0;
    cm[0] := cm[0] + V[4 * i]; cm[1] := cm[1] + V[4 * i + 1]; cm[2] := cm[2] + V[4 * i + 2];
  end;
  // no net drift, and exactly temperature T
  ke := 0;
  for i := 0 to N - 1 do
  begin
    for k := 0 to 2 do V[4 * i + k] := V[4 * i + k] - cm[k] / N;
    ke := ke + 0.5 * (Sqr(V[4 * i]) + Sqr(V[4 * i + 1]) + Sqr(V[4 * i + 2]));
  end;
  s := Sqrt(T / (2 * ke / (3 * N)));
  for i := 0 to 4 * N - 1 do V[i] := V[i] * s;
  for i := 0 to N - 1 do
  begin
    vm := Sqrt(Sqr(V[4 * i]) + Sqr(V[4 * i + 1]) + Sqr(V[4 * i + 2]));
    HostPos[4 * i + 3] := vm;
  end;

  // buffers
  if FShared and (GLBuffer <> 0) then
  begin
    // the GL buffer (already sized by the caller, glFinish'ed) becomes
    // the position buffer, and the starting positions go into it
    FPos := clCreateFromGLBuffer(FCtx, CL_MEM_READ_WRITE, GLBuffer, @Err);
    CheckCL(Err, 'clCreateFromGLBuffer');
    CheckCL(clEnqueueAcquireGLObjects(FQueue, 1, @FPos, 0, nil, nil), 'acquire GL buffer');
    CheckCL(clEnqueueWriteBuffer(FQueue, FPos, CL_TRUE, 0, 16 * N, @HostPos[0], 0, nil, nil), 'write pos');
    CheckCL(clEnqueueReleaseGLObjects(FQueue, 1, @FPos, 0, nil, nil), 'release GL buffer');
    CheckCL(clFinish(FQueue), 'clFinish');
  end
  else
  begin
    FShared := FShared and (GLBuffer <> 0);
    FPos := clCreateBuffer(FCtx, CL_MEM_READ_WRITE, 16 * N, nil, @Err);
    CheckCL(Err, 'clCreateBuffer pos');
    CheckCL(clEnqueueWriteBuffer(FQueue, FPos, CL_TRUE, 0, 16 * N, @HostPos[0], 0, nil, nil), 'write pos');
  end;
  FVelA := clCreateBuffer(FCtx, CL_MEM_READ_WRITE, 16 * N, nil, @Err);
  CheckCL(Err, 'clCreateBuffer vel');
  FVelB := clCreateBuffer(FCtx, CL_MEM_READ_WRITE, 16 * N, nil, @Err);
  CheckCL(Err, 'clCreateBuffer vel');
  CheckCL(clEnqueueWriteBuffer(FQueue, FVelA, CL_TRUE, 0, 16 * N, @V[0], 0, nil, nil), 'write vel');
  FCorr := clCreateBuffer(FCtx, CL_MEM_READ_WRITE, 16 * N, nil, @Err);
  CheckCL(Err, 'clCreateBuffer corr');
  FillChar(V[0], 16 * N, 0);
  CheckCL(clEnqueueWriteBuffer(FQueue, FCorr, CL_TRUE, 0, 16 * N, @V[0], 0, nil, nil), 'write corr');
  FPartner := clCreateBuffer(FCtx, CL_MEM_READ_WRITE, 4 * N, nil, @Err);
  CheckCL(Err, 'clCreateBuffer partner');

  // cells over the largest box the piston allows: unit cubes (one
  // diameter, the reach of a collision) unless that makes more than
  // CellBudget, when they grow, and with them how many spheres one may
  // hold - enough for random close packing, 1.3 per unit volume, and a
  // margin. The memory then follows the box's volume, not N.
  FCellSize := Max(1.0, Power(FLx * FLy * FLzMax / CellBudget, 1 / 3));
  FGrid[0] := Ceil(FLx / FCellSize); FGrid[1] := Ceil(FLy / FCellSize);
  FGrid[2] := Ceil(FLzMax / FCellSize); FGrid[3] := 1;
  FNCells := FGrid[0] * FGrid[1] * FGrid[2];
  FMaxP := Max(8, Ceil(1.3 * Power(FCellSize, 3)) + 4);
  FCounts := clCreateBuffer(FCtx, CL_MEM_READ_WRITE, 4 * FNCells, nil, @Err);
  CheckCL(Err, 'clCreateBuffer counts');
  FCells := clCreateBuffer(FCtx, CL_MEM_READ_WRITE, 4 * PtrUInt(FMaxP) * PtrUInt(FNCells), nil, @Err);
  CheckCL(Err, 'clCreateBuffer cells');
  Zero := 0;
  FOverflow := clCreateBuffer(FCtx, CL_MEM_READ_WRITE, 4, nil, @Err);
  CheckCL(Err, 'clCreateBuffer overflow');
  CheckCL(clEnqueueWriteBuffer(FQueue, FOverflow, CL_TRUE, 0, 4, @Zero, 0, nil, nil), 'write');
  FNColl := clCreateBuffer(FCtx, CL_MEM_READ_WRITE, 4, nil, @Err);
  CheckCL(Err, 'clCreateBuffer ncoll');
  FLocal := LocalSize;
  FGroups := (N + FLocal - 1) div FLocal;
  FPartial := clCreateBuffer(FCtx, CL_MEM_READ_WRITE, 16 * FGroups, nil, @Err);
  CheckCL(Err, 'clCreateBuffer partial');
  SetLength(FHostPartial, 4 * FGroups);
end;

procedure TGasSim.ScaleVelocities(s: Single);
var
  n: LongInt;
begin
  n := FN;
  Arg(KScale, 0, SizeOf(cl_mem), @FVelA);
  Arg(KScale, 1, 4, @n);
  Arg(KScale, 2, 4, @s);
  Enqueue(KScale, FN);
end;

procedure TGasSim.SetTemperature(T: Double);
begin
  if FTemp > 0 then ScaleVelocities(Sqrt(T / FTemp))
  else if FT0 > 0 then ScaleVelocities(Sqrt(T / FT0));
  TargetT := T;
  FTemp := T;
end;

procedure TGasSim.Frame;
var
  vrms, dt, W, imp, impP, ke, ax, ay, az, tf, alpha, x: Double;
  sub, i, n: Integer;
  Lv: array[0..3] of Single;
  dtS, WS, rS, ics: Single;
  Zero, nc: LongInt;
  tmp: cl_mem;
begin
  n := FN;
  tf := SimTimePerFrame;
  vrms := Sqrt(3 * Max(FTemp, 1e-6));
  if not FHaveStats then vrms := Sqrt(3 * FT0);
  FSubsteps := Max(1, Ceil(tf * vrms / StepFrac));
  dt := tf / FSubsteps;

  if FShared then
    CheckCL(clEnqueueAcquireGLObjects(FQueue, 1, @FPos, 0, nil, nil), 'acquire GL buffer');
  FillChar(FHostPartial[0], Length(FHostPartial) * 4, 0);
  CheckCL(clEnqueueWriteBuffer(FQueue, FPartial, CL_FALSE, 0, 16 * FGroups, @FHostPartial[0], 0, nil, nil), 'write');
  Zero := 0;
  CheckCL(clEnqueueWriteBuffer(FQueue, FNColl, CL_FALSE, 0, 4, @Zero, 0, nil, nil), 'write');

  rS := 0.5;
  Arg(KClear, 0, SizeOf(cl_mem), @FCounts);
  Arg(KClear, 1, 4, @FNCells);
  Arg(KBin, 0, SizeOf(cl_mem), @FPos);
  Arg(KBin, 1, 4, @n);
  Arg(KBin, 2, 16, @FGrid[0]);
  ics := 1 / FCellSize;
  Arg(KBin, 3, 4, @ics);
  Arg(KBin, 4, 4, @FMaxP);
  Arg(KBin, 5, SizeOf(cl_mem), @FCounts);
  Arg(KBin, 6, SizeOf(cl_mem), @FCells);
  Arg(KBin, 7, SizeOf(cl_mem), @FOverflow);
  for sub := 1 to FSubsteps do
  begin
    // the piston, towards its target
    W := 0;
    if Abs(LzTarget - FLz) > 1e-9 then
    begin
      W := Sign(LzTarget - FLz) * PistonSpeed;
      if Abs(W * dt) > Abs(LzTarget - FLz) then W := (LzTarget - FLz) / dt;
    end;
    FLz := FLz + W * dt;
    Lv[0] := FLx; Lv[1] := FLy; Lv[2] := FLz; Lv[3] := 0;
    dtS := dt; WS := W;
    Arg(KMove, 0, SizeOf(cl_mem), @FPos);
    Arg(KMove, 1, SizeOf(cl_mem), @FVelA);
    Arg(KMove, 2, SizeOf(cl_mem), @FCorr);
    Arg(KMove, 3, 4, @n);
    Arg(KMove, 4, 4, @dtS);
    Arg(KMove, 5, 16, @Lv[0]);
    Arg(KMove, 6, 4, @WS);
    Arg(KMove, 7, 4, @rS);
    Arg(KMove, 8, SizeOf(cl_mem), @FPartial);
    Arg(KMove, 9, 16 * FLocal, nil);
    Enqueue(KMove, FN, FLocal);
    Enqueue(KClear, FNCells);
    Enqueue(KBin, FN);
    Arg(KPairs, 0, SizeOf(cl_mem), @FPos);
    Arg(KPairs, 1, SizeOf(cl_mem), @FVelA);
    Arg(KPairs, 2, 4, @n);
    Arg(KPairs, 3, 16, @FGrid[0]);
    Arg(KPairs, 4, 4, @ics);
    Arg(KPairs, 5, 4, @FMaxP);
    Arg(KPairs, 6, SizeOf(cl_mem), @FCounts);
    Arg(KPairs, 7, SizeOf(cl_mem), @FCells);
    Arg(KPairs, 8, SizeOf(cl_mem), @FPartner);
    Enqueue(KPairs, FN);
    Arg(KApply, 0, SizeOf(cl_mem), @FPos);
    Arg(KApply, 1, SizeOf(cl_mem), @FVelA);
    Arg(KApply, 2, SizeOf(cl_mem), @FVelB);
    Arg(KApply, 3, SizeOf(cl_mem), @FCorr);
    Arg(KApply, 4, 4, @n);
    Arg(KApply, 5, SizeOf(cl_mem), @FPartner);
    Arg(KApply, 6, SizeOf(cl_mem), @FNColl);
    Enqueue(KApply, FN);
    tmp := FVelA; FVelA := FVelB; FVelB := tmp;
  end;
  if not FShared then
    CheckCL(clEnqueueReadBuffer(FQueue, FPos, CL_FALSE, 0, 16 * FN, @HostPos[0], 0, nil, nil), 'read pos');
  CheckCL(clEnqueueReadBuffer(FQueue, FPartial, CL_FALSE, 0, 16 * FGroups, @FHostPartial[0], 0, nil, nil), 'read');
  CheckCL(clEnqueueReadBuffer(FQueue, FNColl, CL_FALSE, 0, 4, @nc, 0, nil, nil), 'read');
  CheckCL(clEnqueueReadBuffer(FQueue, FOverflow, CL_FALSE, 0, 4, @FOverflowCount, 0, nil, nil), 'read');
  if FShared then
    CheckCL(clEnqueueReleaseGLObjects(FQueue, 1, @FPos, 0, nil, nil), 'release GL buffer');
  CheckCL(clFinish(FQueue), 'clFinish');
  FTime := FTime + tf;

  imp := 0; impP := 0; ke := 0;
  for i := 0 to FGroups - 1 do
  begin
    imp := imp + FHostPartial[4 * i];
    impP := impP + FHostPartial[4 * i + 1];
    ke := ke + FHostPartial[4 * i + 2];
  end;
  ke := ke / FSubsteps;
  FTemp := 2 * ke / (3 * FN);
  ax := FLx - 1; ay := FLy - 1; az := FLz - 1;
  // smoothed over ~PressureFrames frames (the first frame taken as is)
  if FHaveStats then alpha := 1 / Max(1, PressureFrames) else alpha := 1;
  x := imp / (2 * (ay * az + ax * az + ax * ay) * tf);
  FPress := FPress + alpha * (x - FPress);
  x := impP / (ax * ay * tf);
  FPistonPress := FPistonPress + alpha * (x - FPistonPress);
  x := nc / (FN * tf);
  FCollRate := FCollRate + alpha * (x - FCollRate);
  FHaveStats := True;

  if Thermostat and (FTemp > 0) and (Abs(FTemp - TargetT) > 1e-6 * TargetT) then
  begin
    ScaleVelocities(Sqrt(TargetT / FTemp));
    clFinish(FQueue);
    FTemp := TargetT;
  end;
end;

procedure TGasSim.ReadVelocities(var V: array of Single);
begin
  CheckCL(clEnqueueReadBuffer(FQueue, FVelA, CL_TRUE, 0, 16 * FN, @V[0], 0, nil, nil), 'read vel');
end;

function TGasSim.Volume: Double;
begin
  Result := (FLx - 1) * (FLy - 1) * (FLz - 1);
end;

function TGasSim.PackingFraction: Double;
begin
  Result := FN * Pi / 6 / Volume;
end;

end.
