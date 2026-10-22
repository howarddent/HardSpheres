unit clgl;

{******************************************************************************
  OpenCL bindings for HardSpheres: the core calls the simulation uses plus
  cl_khr_gl_sharing (a CL buffer made from a GL vertex buffer), loaded at
  run time with LoadLibrary, so the program starts on a machine without
  OpenCL and says so instead of failing to link.

  Only what ugassim.pas calls is declared. The handle types are opaque
  pointers, as in the C headers.

  GL sharing (Linux, GLX): the CL context is created with the current GLX
  context and display as properties, on the device clGetGLContextInfoKHR
  names for that GL context (an extension function, looked up per
  platform). On a hybrid-graphics laptop GL runs on the integrated GPU
  unless the program is started with NVIDIA PRIME offload
  (__NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia), and
  the sharing device is then that GPU.
*******************************************************************************}

{$mode objfpc}{$H+}

interface

uses
  SysUtils, DynLibs;

type
  cl_int = LongInt;
  cl_uint = LongWord;
  cl_ulong = QWord;
  cl_bitfield = cl_ulong;
  cl_platform_id = Pointer;
  cl_device_id = Pointer;
  cl_context = Pointer;
  cl_command_queue = Pointer;
  cl_mem = Pointer;
  cl_program = Pointer;
  cl_kernel = Pointer;
  cl_event = Pointer;
  Pcl_platform_id = ^cl_platform_id;
  Pcl_device_id = ^cl_device_id;
  Pcl_mem = ^cl_mem;
  Pcl_event = ^cl_event;
  Pcl_int = ^cl_int;
  Pcl_uint = ^cl_uint;
  cl_context_properties = PtrInt;
  Pcl_context_properties = ^cl_context_properties;

const
  CL_SUCCESS = 0;
  CL_TRUE = 1;
  CL_FALSE = 0;
  CL_DEVICE_TYPE_GPU = 1 shl 2;
  CL_DEVICE_TYPE_ALL = $FFFFFFFF;
  CL_PLATFORM_NAME = $0902;
  CL_PLATFORM_EXTENSIONS = $0904;
  CL_DEVICE_NAME = $102B;
  CL_DEVICE_MAX_COMPUTE_UNITS = $1002;
  CL_DEVICE_EXTENSIONS = $1030;
  CL_PROGRAM_BUILD_LOG = $1183;
  CL_MEM_READ_WRITE = 1 shl 0;
  CL_MEM_READ_ONLY = 1 shl 2;
  CL_CONTEXT_PLATFORM = $1084;
  // cl_khr_gl_sharing
  CL_GL_CONTEXT_KHR = $2008;
  CL_GLX_DISPLAY_KHR = $200A;
  CL_CURRENT_DEVICE_FOR_GL_CONTEXT_KHR = $2006;

type
  TclGetPlatformIDs = function(num_entries: cl_uint; platforms: Pcl_platform_id; num_platforms: Pcl_uint): cl_int; cdecl;
  TclGetPlatformInfo = function(platform_: cl_platform_id; param_name: cl_uint; size: PtrUInt;
    value: Pointer; size_ret: PPtrUInt): cl_int; cdecl;
  TclGetDeviceIDs = function(platform_: cl_platform_id; device_type: cl_bitfield; num_entries: cl_uint;
    devices: Pcl_device_id; num_devices: Pcl_uint): cl_int; cdecl;
  TclGetDeviceInfo = function(device: cl_device_id; param_name: cl_uint; size: PtrUInt;
    value: Pointer; size_ret: PPtrUInt): cl_int; cdecl;
  TclCreateContext = function(properties: Pcl_context_properties; num_devices: cl_uint;
    devices: Pcl_device_id; notify: Pointer; user_data: Pointer; errcode_ret: Pcl_int): cl_context; cdecl;
  TclCreateCommandQueue = function(context: cl_context; device: cl_device_id;
    properties: cl_bitfield; errcode_ret: Pcl_int): cl_command_queue; cdecl;
  TclCreateBuffer = function(context: cl_context; flags: cl_bitfield; size: PtrUInt; host_ptr: Pointer;
    errcode_ret: Pcl_int): cl_mem; cdecl;
  TclEnqueueWriteBuffer = function(queue: cl_command_queue; buffer: cl_mem; blocking: cl_uint;
    offset, size: PtrUInt; ptr: Pointer; num_events: cl_uint; wait_list: Pcl_event; event: Pcl_event): cl_int; cdecl;
  TclEnqueueReadBuffer = function(queue: cl_command_queue; buffer: cl_mem; blocking: cl_uint;
    offset, size: PtrUInt; ptr: Pointer; num_events: cl_uint; wait_list: Pcl_event; event: Pcl_event): cl_int; cdecl;
  TclCreateProgramWithSource = function(context: cl_context; count: cl_uint; strings: PPAnsiChar;
    lengths: PPtrUInt; errcode_ret: Pcl_int): cl_program; cdecl;
  TclBuildProgram = function(program_: cl_program; num_devices: cl_uint; devices: Pcl_device_id;
    options: PAnsiChar; notify: Pointer; user_data: Pointer): cl_int; cdecl;
  TclGetProgramBuildInfo = function(program_: cl_program; device: cl_device_id; param_name: cl_uint;
    size: PtrUInt; value: Pointer; size_ret: PPtrUInt): cl_int; cdecl;
  TclCreateKernel = function(program_: cl_program; name: PAnsiChar; errcode_ret: Pcl_int): cl_kernel; cdecl;
  TclSetKernelArg = function(kernel: cl_kernel; index: cl_uint; size: PtrUInt; value: Pointer): cl_int; cdecl;
  TclEnqueueNDRangeKernel = function(queue: cl_command_queue; kernel: cl_kernel; work_dim: cl_uint;
    global_offset, global_size, local_size: PPtrUInt; num_events: cl_uint; wait_list: Pcl_event;
    event: Pcl_event): cl_int; cdecl;
  TclFinish = function(queue: cl_command_queue): cl_int; cdecl;
  TclRelease = function(obj: Pointer): cl_int; cdecl;
  TclGetExtensionFunctionAddressForPlatform = function(platform_: cl_platform_id; name: PAnsiChar): Pointer; cdecl;
  // cl_khr_gl_sharing
  TclCreateFromGLBuffer = function(context: cl_context; flags: cl_bitfield; bufobj: LongWord;
    errcode_ret: Pcl_int): cl_mem; cdecl;
  TclEnqueueAcquireGLObjects = function(queue: cl_command_queue; num_objects: cl_uint; objects: Pcl_mem;
    num_events: cl_uint; wait_list: Pcl_event; event: Pcl_event): cl_int; cdecl;
  TclGetGLContextInfoKHR = function(properties: Pcl_context_properties; param_name: cl_uint;
    size: PtrUInt; value: Pointer; size_ret: PPtrUInt): cl_int; cdecl;

var
  clGetPlatformIDs: TclGetPlatformIDs;
  clGetPlatformInfo: TclGetPlatformInfo;
  clGetDeviceIDs: TclGetDeviceIDs;
  clGetDeviceInfo: TclGetDeviceInfo;
  clCreateContext: TclCreateContext;
  clCreateCommandQueue: TclCreateCommandQueue;
  clCreateBuffer: TclCreateBuffer;
  clEnqueueWriteBuffer: TclEnqueueWriteBuffer;
  clEnqueueReadBuffer: TclEnqueueReadBuffer;
  clCreateProgramWithSource: TclCreateProgramWithSource;
  clBuildProgram: TclBuildProgram;
  clGetProgramBuildInfo: TclGetProgramBuildInfo;
  clCreateKernel: TclCreateKernel;
  clSetKernelArg: TclSetKernelArg;
  clEnqueueNDRangeKernel: TclEnqueueNDRangeKernel;
  clFinish: TclFinish;
  clReleaseMemObject, clReleaseKernel, clReleaseProgram,
  clReleaseCommandQueue, clReleaseContext: TclRelease;
  clGetExtensionFunctionAddressForPlatform: TclGetExtensionFunctionAddressForPlatform;
  clCreateFromGLBuffer: TclCreateFromGLBuffer;
  clEnqueueAcquireGLObjects, clEnqueueReleaseGLObjects: TclEnqueueAcquireGLObjects;

// loads the OpenCL ICD loader once; False (and Msg) if it isn't there
function LoadOpenCL(out Msg: string): Boolean;

// the text of a platform or device string property
function PlatformString(P: cl_platform_id; Param: cl_uint): string;
function DeviceString(D: cl_device_id; Param: cl_uint): string;

// an exception carrying the call's name and the CL error code
procedure CheckCL(Err: cl_int; const What: string);

implementation

var
  LibHandle: TLibHandle = NilHandle;

function LoadOpenCL(out Msg: string): Boolean;
  function P(const Name: string): Pointer;
  begin
    Result := GetProcedureAddress(LibHandle, Name);
  end;
begin
  Msg := '';
  if LibHandle = NilHandle then
  begin
    {$IFDEF WINDOWS}
    LibHandle := LoadLibrary('OpenCL.dll');
    {$ELSE}
    LibHandle := LoadLibrary('libOpenCL.so.1');
    if LibHandle = NilHandle then LibHandle := LoadLibrary('libOpenCL.so');
    {$ENDIF}
    if LibHandle = NilHandle then
    begin
      Msg := 'the OpenCL library (libOpenCL.so.1 / OpenCL.dll) was not found';
      Exit(False);
    end;
    Pointer(clGetPlatformIDs) := P('clGetPlatformIDs');
    Pointer(clGetPlatformInfo) := P('clGetPlatformInfo');
    Pointer(clGetDeviceIDs) := P('clGetDeviceIDs');
    Pointer(clGetDeviceInfo) := P('clGetDeviceInfo');
    Pointer(clCreateContext) := P('clCreateContext');
    Pointer(clCreateCommandQueue) := P('clCreateCommandQueue');
    Pointer(clCreateBuffer) := P('clCreateBuffer');
    Pointer(clEnqueueWriteBuffer) := P('clEnqueueWriteBuffer');
    Pointer(clEnqueueReadBuffer) := P('clEnqueueReadBuffer');
    Pointer(clCreateProgramWithSource) := P('clCreateProgramWithSource');
    Pointer(clBuildProgram) := P('clBuildProgram');
    Pointer(clGetProgramBuildInfo) := P('clGetProgramBuildInfo');
    Pointer(clCreateKernel) := P('clCreateKernel');
    Pointer(clSetKernelArg) := P('clSetKernelArg');
    Pointer(clEnqueueNDRangeKernel) := P('clEnqueueNDRangeKernel');
    Pointer(clFinish) := P('clFinish');
    Pointer(clReleaseMemObject) := P('clReleaseMemObject');
    Pointer(clReleaseKernel) := P('clReleaseKernel');
    Pointer(clReleaseProgram) := P('clReleaseProgram');
    Pointer(clReleaseCommandQueue) := P('clReleaseCommandQueue');
    Pointer(clReleaseContext) := P('clReleaseContext');
    Pointer(clGetExtensionFunctionAddressForPlatform) := P('clGetExtensionFunctionAddressForPlatform');
    // core since OpenCL 1.0 under cl_khr_gl_sharing, exported by the ICD loader
    Pointer(clCreateFromGLBuffer) := P('clCreateFromGLBuffer');
    Pointer(clEnqueueAcquireGLObjects) := P('clEnqueueAcquireGLObjects');
    Pointer(clEnqueueReleaseGLObjects) := P('clEnqueueReleaseGLObjects');
  end;
  Result := Assigned(clGetPlatformIDs) and Assigned(clCreateContext)
    and Assigned(clEnqueueNDRangeKernel);
  if not Result then Msg := 'the OpenCL library lacks core entry points';
end;

function PlatformString(P: cl_platform_id; Param: cl_uint): string;
var
  n: PtrUInt;
begin
  Result := '';
  n := 0;
  if (clGetPlatformInfo(P, Param, 0, nil, @n) <> CL_SUCCESS) or (n = 0) then Exit;
  SetLength(Result, n);
  clGetPlatformInfo(P, Param, n, @Result[1], nil);
  SetLength(Result, StrLen(PChar(Result)));
end;

function DeviceString(D: cl_device_id; Param: cl_uint): string;
var
  n: PtrUInt;
begin
  Result := '';
  n := 0;
  if (clGetDeviceInfo(D, Param, 0, nil, @n) <> CL_SUCCESS) or (n = 0) then Exit;
  SetLength(Result, n);
  clGetDeviceInfo(D, Param, n, @Result[1], nil);
  SetLength(Result, StrLen(PChar(Result)));
end;

procedure CheckCL(Err: cl_int; const What: string);
begin
  if Err <> CL_SUCCESS then
    raise Exception.CreateFmt('%s failed: OpenCL error %d', [What, Err]);
end;

end.
