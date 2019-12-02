#!/usr/bin/env bash

set -e

module load $MPI_LIB_MODULE
module load $COMPILER_MODULE

for folder in "$@"; do
  cd $folder
  if [[ $folder =~ "-isotropic" ]]; then
    ~/lare3d/tools/build.sh -m archie
  elif [[ $folder =~ "-switching" ]]; then
    ~/lare3d/tools/build.sh -s -m archie
  else
    echo "No option for input folder"
  fi
  cd ..
done
