#!/bin/sh
# Build the N2O simulation drivers with Free Pascal: ./build.sh [name ...]
# (all of them by default). They use Sticky/ustickysim.pas and clgl.pas, and
# load Sticky/sticky.cl relative to the executable.
cd "$(dirname "$0")"
mkdir -p lib
for p in ${@:-iso2 iso3 coex bench}; do
  fpc -O2 -Fu../../Sticky -Fu../.. -FUlib "$p.lpr" | grep -iE 'error|fatal' && exit 1
  echo "built $p"
done
