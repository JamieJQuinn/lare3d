#!/usr/bin/env bash

set -e

for folder in "$@"; do
  cd $folder
  latest_snapshot=$(ls Data/*.sdf | sort | tail -n 1 | grep -Eo '[0-9]{4}')
  sed -i -e 's/\(restart_snapshot = \).*/\1'${latest_snapshot}'/' src/control.f90
  #snapshot=5
  #sed -i -e 's/\(restart_snapshot = \).*/\1'${snapshot}'/' src/control.f90
  sed -i -e 's/\(initial = \).*/\1IC_RESTART/' src/control.f90
  sed -i -e 's/\(t_end = \).*/\1'600.0_num'/' src/control.f90
  dt_snapshots=50.0
  sed -i -e 's/\(dt_snapshots = \).*\(_num\)/\1'${dt_snapshots}'\2/' src/control.f90
  #sed -i -e 's/dump_mask(1:8) = .*/dump_mask(2:5) = .TRUE./' src/control.f90
  cd ..
done
