#!/bin/sh
# Runs Isotherms with OpenGL on the NVIDIA GPU (PRIME render offload), so
# OpenCL and OpenGL are on the same device and share the position buffer.
cd "$(dirname "$0")"
exec env __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia ./Isotherms "$@"
