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

# Build the N2O simulation drivers with Free Pascal: ./build.sh [name ...]
# (all of them by default). They use Sticky/ustickysim.pas and clgl.pas, and
# load Sticky/sticky.cl relative to the executable.
cd "$(dirname "$0")"
mkdir -p lib
for p in ${@:-iso2 iso3 coex bench}; do
  fpc -O2 -Fu../../Sticky -Fu../.. -FUlib "$p.lpr" | grep -iE 'error|fatal' && exit 1
  echo "built $p"
done
