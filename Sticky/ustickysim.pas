unit ustickysim;

{******************************************************************************
  Sticky spheres: hard spheres with a stepped attraction - a square well,
  or a stepped approximation to Lennard-Jones - on the GPU (OpenCL),
  time-stepped as in HardSpheres' ugassim.pas.

  Units: sphere diameter 1 (the hard core, and the Lennard-Jones sigma),
  mass 1, Boltzmann's constant 1; T = 2 KE/(3N). Energies of the potential
  are Eps times the table U.

  The potential (TStepPotential): the hard core at R[0] = 1, then K flat
  shells, shell s on R[s] <= r < R[s+1] at Eps U[s], zero beyond the
  cutoff R[K]. Every force is an impulse at a step (see sticky.cl).
    SquareWellPotential(lam): one shell, U = -1, out to lam;
    LennardJonesPotential: nine shells out to 2 sigma, each at the volume
      average of 4 (r^-12 - r^-6) over the shell - the well's floor about
      -0.98 Eps near r = 1.12, a repulsive rise to 1, and the tail falling
      to -0.09 Eps at the cutoff (the full LJ there is -0.06). The hard core
      replaces LJ's repulsion below sigma, where LJ crosses zero.
  A "bond" is a pair closer than R[BondEdge] (the well's edge for the
  square well, 1.5 sigma for Lennard-Jones - the Stillinger criterion):
  bonds are counted, coloured, drawn and make the clusters.

  Capture on: on a downhill step inward a pair with less than Eps of
  relative kinetic energy along the line of centres steps in without
  speeding up, and the step's energy is given off (inelastic; Radiated).
  Capture off: energy-conserving, E = KE + PE_gravity + PE_pairs constant.
  Either way a pair climbs out of the well only with the energy to do so.

  The box and thermostat are those of HardSpheres, but the piston is the
  top face (z = Lz): only Lz changes with the volume, the floor and side
  walls stay put. Gravity (Gravity, along -z) pulls condensate to the
  floor.
*******************************************************************************}

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Math, clgl;

const
  MaxNeighbours = 48;   // per sphere, within the cutoff: ~41 within 2 at close packing
  MaxShells = 15;

type
  TStepPotential = record
    Name: string;
    K: Integer;                                // shells
    R: array[0..MaxShells] of Single;          // R[0] = 1 (core) .. R[K] (cutoff)
    U: array[0..MaxShells - 1] of Single;      // shell energies, units of Eps
    BondEdge: Integer;                         // a bond: r < R[BondEdge]
  end;

  TClusterStats = record
    Clusters: Integer;        // of two or more spheres
    Monomers: Integer;
    Largest: Integer;
    MeanSize: Double;         // seen by a sphere: sum s^2 / N
    // Bins[k]: fraction of spheres in clusters of size 2^k .. 2^(k+1) - 1
    Bins: array[0..20] of Double;
    NBins: Integer;
  end;

  TStickySim = class
  private
    FCtx: cl_context;
    FQueue: cl_command_queue;
    FDevice: cl_device_id;
    FProg: cl_program;
    KMove, KClear, KBin, KPairs, KApply, KScale, KFill: cl_kernel;
    FPos, FVelA, FVelB, FCorr, FPartner, FPType, FNbr, FShl, FNN, FRBuf, FUBuf,
      FCounts, FCells, FOverflow, FCnt, FPartial, FPending, FActA, FActB: cl_mem;
    FShared: Boolean;
    FDeviceName: string;
    FN, FGroups, FLocal: Integer;
    FGrid: array[0..3] of LongInt;
    FNCells: Integer;
    FCellSize: Single;
    FMaxP: LongInt;
    FLx, FLy, FLz, FLzMax: Single;
    FPot: TStepPotential;
    FT0: Double;
    FTime: Double;
    FTemp, FKE, FPE, FPairPE, FPress, FPistonPress: Double;
    FRates: array[0..3] of Double;      // core, formed, broken, captured per sphere per unit time
    FTotals: array[0..3] of Int64;
    FBondCount: Integer;
    FRadiated, FThermostatWork: Double;
    FHaveStats: Boolean;
    FOverflowCount: LongInt;
    FSubsteps: Integer;
    FMeanSweeps: Double;          // per substep, last frame
    FCappedSweeps: Integer;       // substeps that hit MaxSweeps, last frame
    FHostPartial: array of Single;
    procedure BuildProgram(const Src: string);
    procedure ReleaseBuffers;
    procedure Enqueue(K: cl_kernel; Count: Integer; Local: Integer = 0);
    procedure Arg(K: cl_kernel; Index: Integer; Size: PtrUInt; Value: Pointer);
    procedure ScaleVelocities(s: Single);
    function GetRate(Index: Integer): Double;
    function GetTotal(Index: Integer): Int64;
  public
    HostPos: array of Single;     // x, y, z, w per sphere, when not shared
    LzTarget: Double;             // the piston (the top face) moves towards this
    PistonSpeed: Double;
    Thermostat: Boolean;
    TargetT: Double;
    Eps: Double;                  // energy scale: the square well's depth, LJ's epsilon
    Gravity: Double;              // acceleration along -z, diameters per time^2
    Capture: Boolean;             // slow pairs stick, giving off the step energy
    StepFrac: Double;
    MaxSweeps: Integer;           // pairs/apply passes per substep, at most
    PendingTrace: array of LongInt;   // pending count per sweep, last substep (diagnostics)
    SimTimePerFrame: Double;
    PressureFrames: Integer;
    constructor Create;
    destructor Destroy; override;
    procedure InitCL(const KernelFile: string; GLContext, GLDisplay: Pointer);
    // N spheres at packing fraction Phi in a cube, temperature T, the
    // potential Pot; no pair inside the cutoff to start with where the
    // spacing allows
    procedure Setup(N: Integer; Phi, T: Double; const Pot: TStepPotential; EqualSpeeds: Boolean;
      GLBuffer: LongWord; LzMaxFactor: Double = 2.0);
    procedure Frame;
    procedure SetTemperature(T: Double);
    procedure ReadVelocities(var V: array of Single);
    // neighbour lists, MaxNeighbours per sphere (-1: empty), and shells
    procedure ReadLists(var Nbr, Shells: array of LongInt);
    procedure ReadPositions(var P: array of Single);
    // clusters joined by bonds (pairs closer than R[BondEdge])
    function ClusterStats(const Nbr, Shells: array of LongInt): TClusterStats;
    // the farthest any listed pair is outside its shell: ~0 when the
    // events keep up (a health check)
    function MaxOvershoot(const Nbr, Shells: array of LongInt; const P: array of Single): Double;
    function Volume: Double;
    function PackingFraction: Double;
    // KE + PE (gravity) + pair energy: constant with capture off and no thermostat or piston
    function Energy: Double;
    property Potential: TStepPotential read FPot;
    property Count: Integer read FN;
    property Lx: Single read FLx;
    property Ly: Single read FLy;
    property Lz: Single read FLz;
    property LzMax: Single read FLzMax;
    property Time: Double read FTime;
    property Temperature: Double read FTemp;
    property KineticEnergy: Double read FKE;
    property PotentialEnergy: Double read FPE;    // in gravity
    property PairEnergy: Double read FPairPE;     // of the pair potential
    property Pressure: Double read FPress;
    property PistonPressure: Double read FPistonPress;
    property Bonds: Integer read FBondCount;
    property CoreRate: Double index 0 read GetRate;
    property FormRate: Double index 1 read GetRate;
    property BreakRate: Double index 2 read GetRate;
    property CaptureRate: Double index 3 read GetRate;
    property CoreTotal: Int64 index 0 read GetTotal;
    property FormTotal: Int64 index 1 read GetTotal;
    property BreakTotal: Int64 index 2 read GetTotal;
    property CaptureTotal: Int64 index 3 read GetTotal;
    property Radiated: Double read FRadiated;
    property ThermostatWork: Double read FThermostatWork;
    property Overflows: LongInt read FOverflowCount;
    property CellSize: Single read FCellSize;
    property Substeps: Integer read FSubsteps;
    property MeanSweeps: Double read FMeanSweeps;
    property CappedSweeps: Integer read FCappedSweeps;
    property Shared: Boolean read FShared;
    property DeviceName: string read FDeviceName;
  end;

function SquareWellPotential(Lambda: Double): TStepPotential;
function LennardJonesPotential: TStepPotential;
// second virial coefficient of the stepped potential,
// (2 pi/3)(1 + sum_s (1 - e^(-Eps U[s]/T))(R[s+1]^3 - R[s]^3))
function SecondVirial(const Pot: TStepPotential; Eps, T: Double): Double;
// neighbours per sphere at low density rho, sum over shells below Edge of
// rho e^(-Eps U[s]/T) (4 pi/3)(R[s+1]^3 - R[s]^3)
function LowDensityNeighbours(const Pot: TStepPotential; Edge: Integer; Eps, T, Rho: Double): Double;

implementation

const
  LocalSize = 256;
  CellBudget = 8000000;
  SweepBatch = 4;

function SquareWellPotential(Lambda: Double): TStepPotential;
begin
  Result := Default(TStepPotential);
  Result.Name := Format('square well, width %.2f', [Lambda]);
  Result.K := 1;
  Result.R[0] := 1;
  Result.R[1] := Max(1.01, Lambda);
  Result.U[0] := -1;
  Result.BondEdge := 1;
end;

// shells chosen finest round the minimum (2^(1/6) = 1.122) and the steep
// inner wall, coarser in the tail
function LennardJonesPotential: TStepPotential;
const
  Edges: array[0..9] of Single = (1.0, 1.06, 1.12, 1.2, 1.3, 1.4, 1.5, 1.65, 1.8, 2.0);
  function F(r: Double): Double;    // integral of 4 (r^-12 - r^-6) r^2 dr
  begin
    Result := 4 * (-Power(r, -9) / 9 + Power(r, -3) / 3);
  end;
var
  k: Integer;
begin
  Result := Default(TStepPotential);
  Result.Name := 'stepped Lennard-Jones, cutoff 2';
  Result.K := High(Edges);
  for k := 0 to Result.K do Result.R[k] := Edges[k];
  for k := 0 to Result.K - 1 do
    Result.U[k] := (F(Edges[k + 1]) - F(Edges[k])) / ((Power(Edges[k + 1], 3) - Power(Edges[k], 3)) / 3);
  Result.BondEdge := 6;               // r < 1.5
end;

function SecondVirial(const Pot: TStepPotential; Eps, T: Double): Double;
var
  s: Integer;
begin
  Result := 1;
  for s := 0 to Pot.K - 1 do
    Result := Result + (1 - Exp(-Eps * Pot.U[s] / T)) * (Power(Pot.R[s + 1], 3) - Power(Pot.R[s], 3));
  Result := 2 * Pi / 3 * Result;
end;

function LowDensityNeighbours(const Pot: TStepPotential; Edge: Integer; Eps, T, Rho: Double): Double;
var
  s: Integer;
begin
  Result := 0;
  for s := 0 to Min(Edge, Pot.K) - 1 do
    Result := Result + Exp(-Eps * Pot.U[s] / T) * (Power(Pot.R[s + 1], 3) - Power(Pot.R[s], 3));
  Result := Result * Rho * 4 * Pi / 3;
end;

constructor TStickySim.Create;
begin
  StepFrac := 0.02;
  MaxSweeps := 32;
  SimTimePerFrame := 0.05;
  PistonSpeed := 0.5;
  TargetT := 1;
  PressureFrames := 20;
  Eps := 1;
  Capture := True;
  FPot := LennardJonesPotential;
end;

destructor TStickySim.Destroy;
begin
  ReleaseBuffers;
  if KMove <> nil then clReleaseKernel(KMove);
  if KClear <> nil then clReleaseKernel(KClear);
  if KBin <> nil then clReleaseKernel(KBin);
  if KPairs <> nil then clReleaseKernel(KPairs);
  if KApply <> nil then clReleaseKernel(KApply);
  if KScale <> nil then clReleaseKernel(KScale);
  if KFill <> nil then clReleaseKernel(KFill);
  if FProg <> nil then clReleaseProgram(FProg);
  if FQueue <> nil then clReleaseCommandQueue(FQueue);
  if FCtx <> nil then clReleaseContext(FCtx);
  inherited Destroy;
end;

procedure TStickySim.ReleaseBuffers;
  procedure R(var M: cl_mem);
  begin
    if M <> nil then clReleaseMemObject(M);
    M := nil;
  end;
begin
  R(FPos); R(FVelA); R(FVelB); R(FCorr); R(FPartner); R(FPType); R(FNbr); R(FShl); R(FNN);
  R(FRBuf); R(FUBuf); R(FCounts); R(FCells); R(FOverflow); R(FCnt); R(FPartial); R(FPending);
  R(FActA); R(FActB);
end;

function TStickySim.GetRate(Index: Integer): Double;
begin
  Result := FRates[Index];
end;

function TStickySim.GetTotal(Index: Integer): Int64;
begin
  Result := FTotals[Index];
end;

procedure TStickySim.InitCL(const KernelFile: string; GLContext, GLDisplay: Pointer);
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
begin
  if not LoadOpenCL(Msg) then raise Exception.Create(Msg);
  np := 0;
  if (clGetPlatformIDs(0, nil, @np) <> CL_SUCCESS) or (np = 0) then
    raise Exception.Create('no OpenCL platform');
  SetLength(Plats, np);
  clGetPlatformIDs(np, @Plats[0], nil);

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
         or (Dev = nil) then Continue;
      FCtx := clCreateContext(@Props[0], 1, @Dev, nil, nil, @Err);
      if (Err = CL_SUCCESS) and (FCtx <> nil) then
      begin
        FDevice := Dev;
        FShared := True;
        Break;
      end;
      FCtx := nil;
    end;

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

procedure TStickySim.BuildProgram(const Src: string);
var
  Err: cl_int;
  P: PAnsiChar;
  Log, Opts: string;
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
  Opts := Format('-cl-fast-relaxed-math -DMAXB=%d', [MaxNeighbours]);
  Err := clBuildProgram(FProg, 1, @FDevice, PAnsiChar(Opts), nil, nil);
  if Err <> CL_SUCCESS then
  begin
    n := 0;
    clGetProgramBuildInfo(FProg, FDevice, CL_PROGRAM_BUILD_LOG, 0, nil, @n);
    SetLength(Log, n);
    if n > 0 then clGetProgramBuildInfo(FProg, FDevice, CL_PROGRAM_BUILD_LOG, n, @Log[1], nil);
    raise Exception.Create('sticky.cl did not build:'#10 + Log);
  end;
  KMove := K('move');
  KClear := K('clear');
  KBin := K('bin');
  KPairs := K('pairs');
  KApply := K('apply');
  KScale := K('scale');
  KFill := K('fill');
end;

procedure TStickySim.Arg(K: cl_kernel; Index: Integer; Size: PtrUInt; Value: Pointer);
begin
  CheckCL(clSetKernelArg(K, Index, Size, Value), Format('clSetKernelArg %d', [Index]));
end;

procedure TStickySim.Enqueue(K: cl_kernel; Count: Integer; Local: Integer);
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

procedure TStickySim.Setup(N: Integer; Phi, T: Double; const Pot: TStepPotential; EqualSpeeds: Boolean;
  GLBuffer: LongWord; LzMaxFactor: Double);
var
  L, a, jit, s, ke, rc: Double;
  ns, i, ix, iy, iz, k, slots: Integer;
  V: array of Single;
  Empty: array of LongInt;
  cm: array[0..2] of Double;
  Err: cl_int;
  Zero: array[0..3] of LongInt;
  perm: array of Integer;
  function Buf(Size: PtrUInt; const What: string): cl_mem;
  var
    E: cl_int;
  begin
    Result := clCreateBuffer(FCtx, CL_MEM_READ_WRITE, Size, nil, @E);
    CheckCL(E, 'clCreateBuffer ' + What);
  end;
begin
  ReleaseBuffers;
  FN := N;
  FPot := Pot;
  rc := FPot.R[FPot.K];
  L := Power(N * Pi / 6 / Phi, 1 / 3) + 1;
  FLx := L; FLy := L; FLz := L;
  FLzMax := L * LzMaxFactor;
  LzTarget := L;
  FT0 := T;
  FTemp := T;
  TargetT := T;
  FTime := 0;
  FHaveStats := False;
  FOverflowCount := 0;
  FBondCount := 0;
  FPE := 0;
  FPairPE := 0;
  FRadiated := 0;
  FThermostatWork := 0;
  for k := 0 to 3 do
  begin
    FRates[k] := 0;
    FTotals[k] := 0;
  end;

  // a jittered simple cubic lattice, kept outside one another's cutoff
  // where the spacing allows, so the gas starts with empty lists
  ns := Ceil(Power(N, 1 / 3) - 1e-9);
  a := (L - 1) / ns;
  if a < 1.0 then raise Exception.CreateFmt('packing fraction %.3g too high to start from a lattice', [Phi]);
  if a > rc then jit := 0.45 * (a - rc) else jit := 0;
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
    HostPos[4 * i + 3] := 0.5;
    if EqualSpeeds then
    begin
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
  ke := 0;
  for i := 0 to N - 1 do
  begin
    for k := 0 to 2 do V[4 * i + k] := V[4 * i + k] - cm[k] / N;
    ke := ke + 0.5 * (Sqr(V[4 * i]) + Sqr(V[4 * i + 1]) + Sqr(V[4 * i + 2]));
  end;
  s := Sqrt(T / (2 * ke / (3 * N)));
  for i := 0 to 4 * N - 1 do V[i] := V[i] * s;
  FKE := 1.5 * N * T;

  if FShared and (GLBuffer <> 0) then
  begin
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
    FPos := Buf(16 * N, 'pos');
    CheckCL(clEnqueueWriteBuffer(FQueue, FPos, CL_TRUE, 0, 16 * N, @HostPos[0], 0, nil, nil), 'write pos');
  end;
  FVelA := Buf(16 * N, 'vel');
  FVelB := Buf(16 * N, 'vel');
  CheckCL(clEnqueueWriteBuffer(FQueue, FVelA, CL_TRUE, 0, 16 * N, @V[0], 0, nil, nil), 'write vel');
  FCorr := Buf(16 * N, 'corr');
  FillChar(V[0], 16 * N, 0);
  CheckCL(clEnqueueWriteBuffer(FQueue, FCorr, CL_TRUE, 0, 16 * N, @V[0], 0, nil, nil), 'write corr');
  FPartner := Buf(4 * N, 'partner');
  FPType := Buf(4 * N, 'ptype');
  FNN := Buf(4 * N, 'nn');
  CheckCL(clEnqueueWriteBuffer(FQueue, FNN, CL_TRUE, 0, 4 * N, @V[0], 0, nil, nil), 'write nn');
  FNbr := Buf(4 * PtrUInt(N) * MaxNeighbours, 'nbr');
  FShl := Buf(4 * PtrUInt(N) * MaxNeighbours, 'shl');
  SetLength(Empty, N * MaxNeighbours);
  for i := 0 to High(Empty) do Empty[i] := -1;
  CheckCL(clEnqueueWriteBuffer(FQueue, FNbr, CL_TRUE, 0, 4 * N * MaxNeighbours, @Empty[0], 0, nil, nil), 'write nbr');
  for i := 0 to High(Empty) do Empty[i] := 0;
  CheckCL(clEnqueueWriteBuffer(FQueue, FShl, CL_TRUE, 0, 4 * N * MaxNeighbours, @Empty[0], 0, nil, nil), 'write shl');
  FRBuf := Buf(4 * (MaxShells + 1), 'R');
  CheckCL(clEnqueueWriteBuffer(FQueue, FRBuf, CL_TRUE, 0, 4 * (MaxShells + 1), @FPot.R[0], 0, nil, nil), 'write R');
  FUBuf := Buf(4 * MaxShells, 'U');
  CheckCL(clEnqueueWriteBuffer(FQueue, FUBuf, CL_TRUE, 0, 4 * MaxShells, @FPot.U[0], 0, nil, nil), 'write U');

  // cells at least one cutoff wide, and room for a liquid
  FCellSize := Max(rc, Power(FLx * FLy * FLzMax / CellBudget, 1 / 3));
  FGrid[0] := Ceil(FLx / FCellSize); FGrid[1] := Ceil(FLy / FCellSize);
  FGrid[2] := Ceil(FLzMax / FCellSize); FGrid[3] := 1;
  FNCells := FGrid[0] * FGrid[1] * FGrid[2];
  FMaxP := Max(12, Ceil(1.5 * Power(FCellSize, 3)) + 6);
  FCounts := Buf(4 * FNCells, 'counts');
  FCells := Buf(4 * PtrUInt(FMaxP) * PtrUInt(FNCells), 'cells');
  FillChar(Zero, SizeOf(Zero), 0);
  FOverflow := Buf(4, 'overflow');
  CheckCL(clEnqueueWriteBuffer(FQueue, FOverflow, CL_TRUE, 0, 4, @Zero, 0, nil, nil), 'write');
  FCnt := Buf(4 * (4 + MaxShells + 1), 'counters');
  FPending := Buf(4, 'pending');
  FActA := Buf(4 * N, 'active');
  FActB := Buf(4 * N, 'active');
  FLocal := LocalSize;
  FGroups := (N + FLocal - 1) div FLocal;
  FPartial := Buf(32 * FGroups, 'partial');
  SetLength(FHostPartial, 8 * FGroups);
end;

procedure TStickySim.ScaleVelocities(s: Single);
var
  n: LongInt;
begin
  n := FN;
  Arg(KScale, 0, SizeOf(cl_mem), @FVelA);
  Arg(KScale, 1, 4, @n);
  Arg(KScale, 2, 4, @s);
  Enqueue(KScale, FN);
end;

procedure TStickySim.SetTemperature(T: Double);
begin
  if FTemp > 0 then
  begin
    ScaleVelocities(Sqrt(T / FTemp));
    clFinish(FQueue);
    FThermostatWork := FThermostatWork + 1.5 * FN * (T - FTemp);
  end;
  TargetT := T;
  FTemp := T;
  FKE := 1.5 * FN * T;
end;

procedure TStickySim.Frame;
var
  vchar, dt, W, imp, impP, ke, ends, sz, pe, ax, ay, az, tf, alpha, x, umin: Double;
  sub, sw, i, n, b, nsw: Integer;
  pend, one, zero1: LongInt;
  Lv: array[0..3] of Single;
  dtS, WS, rS, ics, epsS, vsS, tmS, gS: Single;
  capS, KS, beS: LongInt;
  Zero: array[0..4 + MaxShells] of LongInt;
  nc: array[0..4 + MaxShells] of LongInt;
  tmp: cl_mem;
  function Level(s: Integer): Double;     // shell energy, 0 outside
  begin
    if s < FPot.K then Result := FPot.U[s] else Result := 0;
  end;
begin
  n := FN;
  tf := SimTimePerFrame;
  // a sphere falling to the well's floor speeds up by up to 2 sqrt(depth),
  // and one falling half the box's height by sqrt(g Lz)
  umin := 0;
  for i := 0 to FPot.K - 1 do umin := Min(umin, FPot.U[i]);
  vchar := Sqrt(3 * Max(FTemp, 1e-6) + 2 * Max(Eps, 0) * Abs(umin) + Max(Gravity, 0) * FLz);
  FSubsteps := Max(1, Ceil(tf * vchar / StepFrac));
  dt := tf / FSubsteps;

  if FShared then
    CheckCL(clEnqueueAcquireGLObjects(FQueue, 1, @FPos, 0, nil, nil), 'acquire GL buffer');
  FillChar(FHostPartial[0], Length(FHostPartial) * 4, 0);
  CheckCL(clEnqueueWriteBuffer(FQueue, FPartial, CL_FALSE, 0, 32 * FGroups, @FHostPartial[0], 0, nil, nil), 'write');
  FillChar(Zero, SizeOf(Zero), 0);
  CheckCL(clEnqueueWriteBuffer(FQueue, FCnt, CL_FALSE, 0, SizeOf(Zero), @Zero, 0, nil, nil), 'write');

  rS := 0.5;
  ics := 1 / FCellSize;
  epsS := Max(Eps, 0);
  if Capture then capS := 1 else capS := 0;
  vsS := 2.2 * Sqrt(3 * Max(FTemp, 1e-6));
  tmS := 2 * dt;
  gS := Max(Gravity, 0);
  KS := FPot.K;
  beS := FPot.BondEdge;
  nsw := 0;
  FCappedSweeps := 0;
  Arg(KClear, 0, SizeOf(cl_mem), @FCounts);
  Arg(KClear, 1, 4, @FNCells);
  Arg(KBin, 0, SizeOf(cl_mem), @FPos);
  Arg(KBin, 1, 4, @n);
  Arg(KBin, 2, 16, @FGrid[0]);
  Arg(KBin, 3, 4, @ics);
  Arg(KBin, 4, 4, @FMaxP);
  Arg(KBin, 5, SizeOf(cl_mem), @FCounts);
  Arg(KBin, 6, SizeOf(cl_mem), @FCells);
  Arg(KBin, 7, SizeOf(cl_mem), @FOverflow);
  for sub := 1 to FSubsteps do
  begin
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
    Arg(KMove, 8, SizeOf(cl_mem), @FNN);
    Arg(KMove, 9, SizeOf(cl_mem), @FShl);
    Arg(KMove, 10, SizeOf(cl_mem), @FUBuf);
    Arg(KMove, 11, 4, @beS);
    Arg(KMove, 12, 4, @epsS);
    Arg(KMove, 13, 4, @vsS);
    Arg(KMove, 14, 4, @gS);
    Arg(KMove, 15, SizeOf(cl_mem), @FPartial);
    Arg(KMove, 16, 32 * FLocal, nil);
    Enqueue(KMove, FN, FLocal);
    Enqueue(KClear, FNCells);
    Enqueue(KBin, FN);
    SetLength(PendingTrace, 0);
    // every sphere is searched in the first sweep
    one := 1; zero1 := 0;
    Arg(KFill, 0, SizeOf(cl_mem), @FActA);
    Arg(KFill, 1, 4, @n);
    Arg(KFill, 2, 4, @one);
    Enqueue(KFill, FN);
    for sw := 1 to Max(1, MaxSweeps) do
    begin
      Arg(KFill, 0, SizeOf(cl_mem), @FActB);
      Arg(KFill, 1, 4, @n);
      Arg(KFill, 2, 4, @zero1);
      Enqueue(KFill, FN);
      pend := 0;
      CheckCL(clEnqueueWriteBuffer(FQueue, FPending, CL_FALSE, 0, 4, @pend, 0, nil, nil), 'write');
      Arg(KPairs, 0, SizeOf(cl_mem), @FPos);
      Arg(KPairs, 1, SizeOf(cl_mem), @FVelA);
      Arg(KPairs, 2, 4, @n);
      Arg(KPairs, 3, 16, @FGrid[0]);
      Arg(KPairs, 4, 4, @ics);
      Arg(KPairs, 5, 4, @FMaxP);
      Arg(KPairs, 6, SizeOf(cl_mem), @FCounts);
      Arg(KPairs, 7, SizeOf(cl_mem), @FCells);
      Arg(KPairs, 8, SizeOf(cl_mem), @FNbr);
      Arg(KPairs, 9, SizeOf(cl_mem), @FShl);
      Arg(KPairs, 10, SizeOf(cl_mem), @FNN);
      Arg(KPairs, 11, SizeOf(cl_mem), @FRBuf);
      Arg(KPairs, 12, 4, @KS);
      Arg(KPairs, 13, SizeOf(cl_mem), @FPartner);
      Arg(KPairs, 14, SizeOf(cl_mem), @FPType);
      Arg(KPairs, 15, SizeOf(cl_mem), @FPending);
      Arg(KPairs, 16, SizeOf(cl_mem), @FActA);
      Arg(KPairs, 17, SizeOf(cl_mem), @FActB);
      Enqueue(KPairs, FN);
      tmp := FActA; FActA := FActB; FActB := tmp;
      // nothing pending: this substep is done. Read every SweepBatch
      // sweeps only: each read waits for the GPU, and a spare sweep over
      // an empty backlog costs next to nothing (only active spheres search)
      Inc(nsw);
      if (sw mod SweepBatch = 0) or (sw = MaxSweeps) then
      begin
        CheckCL(clEnqueueReadBuffer(FQueue, FPending, CL_TRUE, 0, 4, @pend, 0, nil, nil), 'read pending');
        SetLength(PendingTrace, Length(PendingTrace) + 1);
        PendingTrace[High(PendingTrace)] := pend;
        if pend = 0 then Break;
        if sw = MaxSweeps then Inc(FCappedSweeps);
      end;
      Arg(KApply, 0, SizeOf(cl_mem), @FPos);
      Arg(KApply, 1, SizeOf(cl_mem), @FVelA);
      Arg(KApply, 2, SizeOf(cl_mem), @FVelB);
      Arg(KApply, 3, SizeOf(cl_mem), @FCorr);
      Arg(KApply, 4, 4, @n);
      Arg(KApply, 5, SizeOf(cl_mem), @FPartner);
      Arg(KApply, 6, SizeOf(cl_mem), @FPType);
      Arg(KApply, 7, SizeOf(cl_mem), @FNbr);
      Arg(KApply, 8, SizeOf(cl_mem), @FShl);
      Arg(KApply, 9, SizeOf(cl_mem), @FNN);
      Arg(KApply, 10, SizeOf(cl_mem), @FRBuf);
      Arg(KApply, 11, SizeOf(cl_mem), @FUBuf);
      Arg(KApply, 12, 4, @KS);
      Arg(KApply, 13, 4, @beS);
      Arg(KApply, 14, 4, @epsS);
      Arg(KApply, 15, 4, @capS);
      Arg(KApply, 16, 4, @tmS);
      Arg(KApply, 17, SizeOf(cl_mem), @FCnt);
      Enqueue(KApply, FN);
      tmp := FVelA; FVelA := FVelB; FVelB := tmp;
    end;
  end;
  if not FShared then
    CheckCL(clEnqueueReadBuffer(FQueue, FPos, CL_FALSE, 0, 16 * FN, @HostPos[0], 0, nil, nil), 'read pos');
  CheckCL(clEnqueueReadBuffer(FQueue, FPartial, CL_FALSE, 0, 32 * FGroups, @FHostPartial[0], 0, nil, nil), 'read');
  CheckCL(clEnqueueReadBuffer(FQueue, FCnt, CL_FALSE, 0, SizeOf(nc), @nc[0], 0, nil, nil), 'read');
  CheckCL(clEnqueueReadBuffer(FQueue, FOverflow, CL_FALSE, 0, 4, @FOverflowCount, 0, nil, nil), 'read');
  if FShared then
    CheckCL(clEnqueueReleaseGLObjects(FQueue, 1, @FPos, 0, nil, nil), 'release GL buffer');
  CheckCL(clFinish(FQueue), 'clFinish');
  FTime := FTime + tf;
  FMeanSweeps := nsw / FSubsteps;

  imp := 0; impP := 0; ke := 0; ends := 0; sz := 0; pe := 0;
  for i := 0 to FGroups - 1 do
  begin
    imp := imp + FHostPartial[8 * i];
    impP := impP + FHostPartial[8 * i + 1];
    ke := ke + FHostPartial[8 * i + 2];
    ends := ends + FHostPartial[8 * i + 3];
    sz := sz + FHostPartial[8 * i + 4];
    pe := pe + FHostPartial[8 * i + 5];
  end;
  // all averaged over the frame's substeps (as of each substep's move);
  // the height from the lowest a centre can be, z = 1/2
  ke := ke / FSubsteps;
  FKE := ke;
  FTemp := 2 * ke / (3 * FN);
  FBondCount := Round(ends / FSubsteps / 2);
  FPE := Max(Gravity, 0) * (sz / FSubsteps - 0.5 * FN);
  FPairPE := pe / FSubsteps;
  for i := 0 to 3 do FTotals[i] := FTotals[i] + nc[i];
  // a capture at boundary b gives off the step there
  for b := 1 to FPot.K do
    FRadiated := FRadiated + nc[4 + b] * Max(Eps, 0) * (Level(b) - Level(b - 1));
  ax := FLx - 1; ay := FLy - 1; az := FLz - 1;
  if FHaveStats then alpha := 1 / Max(1, PressureFrames) else alpha := 1;
  x := imp / (2 * (ay * az + ax * az + ax * ay) * tf);
  FPress := FPress + alpha * (x - FPress);
  x := impP / (ax * ay * tf);
  FPistonPress := FPistonPress + alpha * (x - FPistonPress);
  for i := 0 to 3 do
  begin
    x := nc[i] / (FN * tf);
    FRates[i] := FRates[i] + alpha * (x - FRates[i]);
  end;
  FHaveStats := True;

  if Thermostat and (FTemp > 0) and (Abs(FTemp - TargetT) > 1e-6 * TargetT) then
  begin
    ScaleVelocities(Sqrt(TargetT / FTemp));
    clFinish(FQueue);
    FThermostatWork := FThermostatWork + 1.5 * FN * (TargetT - FTemp);
    FTemp := TargetT;
    FKE := 1.5 * FN * TargetT;
  end;
end;

function TStickySim.Energy: Double;
begin
  Result := FKE + FPE + FPairPE;
end;

procedure TStickySim.ReadVelocities(var V: array of Single);
begin
  CheckCL(clEnqueueReadBuffer(FQueue, FVelA, CL_TRUE, 0, 16 * FN, @V[0], 0, nil, nil), 'read vel');
end;

procedure TStickySim.ReadLists(var Nbr, Shells: array of LongInt);
begin
  CheckCL(clEnqueueReadBuffer(FQueue, FNbr, CL_TRUE, 0, 4 * FN * MaxNeighbours, @Nbr[0], 0, nil, nil), 'read nbr');
  CheckCL(clEnqueueReadBuffer(FQueue, FShl, CL_TRUE, 0, 4 * FN * MaxNeighbours, @Shells[0], 0, nil, nil), 'read shl');
end;

procedure TStickySim.ReadPositions(var P: array of Single);
begin
  if FShared then
    CheckCL(clEnqueueAcquireGLObjects(FQueue, 1, @FPos, 0, nil, nil), 'acquire GL buffer');
  CheckCL(clEnqueueReadBuffer(FQueue, FPos, CL_TRUE, 0, 16 * FN, @P[0], 0, nil, nil), 'read pos');
  if FShared then
    CheckCL(clEnqueueReleaseGLObjects(FQueue, 1, @FPos, 0, nil, nil), 'release GL buffer');
  CheckCL(clFinish(FQueue), 'clFinish');
end;

// clusters: connected components of the bond graph (union-find)
function TStickySim.ClusterStats(const Nbr, Shells: array of LongInt): TClusterStats;
var
  Parent, Size: array of Integer;
  i, k, j, a, b2, s, bin: Integer;
  sum2: Double;
  function Find(x: Integer): Integer;
  begin
    while Parent[x] <> x do
    begin
      Parent[x] := Parent[Parent[x]];
      x := Parent[x];
    end;
    Result := x;
  end;
begin
  FillChar(Result, SizeOf(Result), 0);
  SetLength(Parent, FN);
  SetLength(Size, FN);
  for i := 0 to FN - 1 do
  begin
    Parent[i] := i;
    Size[i] := 1;
  end;
  for i := 0 to FN - 1 do
    for k := 0 to MaxNeighbours - 1 do
    begin
      j := Nbr[i * MaxNeighbours + k];
      if (j <= i) or (j >= FN) or (Shells[i * MaxNeighbours + k] >= FPot.BondEdge) then Continue;
      a := Find(i); b2 := Find(j);
      if a = b2 then Continue;
      if Size[a] < Size[b2] then begin s := a; a := b2; b2 := s; end;
      Parent[b2] := a;
      Size[a] := Size[a] + Size[b2];
    end;
  sum2 := 0;
  Result.NBins := 1;
  while (1 shl Result.NBins) <= FN do Inc(Result.NBins);
  Result.NBins := Min(Result.NBins, Length(Result.Bins));
  for i := 0 to FN - 1 do
    if Find(i) = i then
    begin
      s := Size[i];
      if s = 1 then Inc(Result.Monomers) else Inc(Result.Clusters);
      Result.Largest := Max(Result.Largest, s);
      sum2 := sum2 + Sqr(Double(s));
      bin := Min(Trunc(Log2(s) + 1e-9), Result.NBins - 1);
      Result.Bins[bin] := Result.Bins[bin] + s / FN;
    end;
  Result.MeanSize := sum2 / FN;
end;

function TStickySim.MaxOvershoot(const Nbr, Shells: array of LongInt; const P: array of Single): Double;
var
  i, k, j, s: Integer;
  r: Double;
begin
  Result := 0;
  for i := 0 to FN - 1 do
    for k := 0 to MaxNeighbours - 1 do
    begin
      j := Nbr[i * MaxNeighbours + k];
      if (j <= i) or (j >= FN) then Continue;
      s := Shells[i * MaxNeighbours + k];
      r := Sqrt(Sqr(P[4 * i] - P[4 * j]) + Sqr(P[4 * i + 1] - P[4 * j + 1]) + Sqr(P[4 * i + 2] - P[4 * j + 2]));
      Result := Max(Result, Max(r - FPot.R[s + 1], FPot.R[s] - r));
    end;
end;

function TStickySim.Volume: Double;
begin
  Result := (FLx - 1) * (FLy - 1) * (FLz - 1);
end;

function TStickySim.PackingFraction: Double;
begin
  Result := FN * Pi / 6 / Volume;
end;

end.
