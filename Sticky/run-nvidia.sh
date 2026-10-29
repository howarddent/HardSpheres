#!/bin/sh
# This file is part of HardSpheres.
# Copyright (C) 2026 Howard Dent
# SPDX-License-Identifier: GPL-3.0-or-later
#
# HardSpheres is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your
# option) any later version. It is distributed WITHOUT ANY WARRANTY; see
# the file LICENSE at the top of the repository for details.

# Runs StickySpheres with OpenGL on the NVIDIA GPU (PRIME render offload), so
# OpenCL and OpenGL are on the same device and share the position buffer.
cd "$(dirname "$0")"
exec env __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia ./StickySpheres "$@"
