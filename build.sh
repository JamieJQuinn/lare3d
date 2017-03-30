#!/usr/bin/env bash
machine='office'
n_grid_points='false'
braginskii='false'
switching='false'
output='true'
defines=''
mode=''
compiler=''
mpif90='mpif90'

while getopts "bsdom:p:" flag; do
  case "${flag}" in
    # Enable braginskii visc
    b) braginskii='true' ;;
    # Enable Switching visc
    s) switching='true' ;;
    # Disable output
    o) output='false' ;;
    # Enable debug
    d) mode='debug' ;;
    # Set machine
    m) machine="${OPTARG}" ;;
    # Set num grid points
    p) n_grid_points="${OPTARG}" ;;
  esac
done

if [ "$braginskii" == 'true' ]; then
  defines+=' -DBRAGINSKII_VISCOSITY'
fi

if [ "$switching" == 'true' ]; then
  defines+=' -DSWITCHING_VISCOSITY'
fi

if [ "$output" == 'false' ]; then
  defines+=' -DNO_IO'
fi

if [ "$machine" == "euclid" ]; then
  export PATH=/usr/lib64/openmpi/bin:$PATH
  compiler='gfortran'
elif [ "$machine" == "archie" ]; then
  compiler='intel'
  mpif90='mpiifort'
else
  compiler='gfortran'
  echo "Error: No machine specified, building with gfortran"
fi

if [ "$n_grid_points" != 'false' ]; then
  sed -i -e 's/\(n[xyz]_global = \)[0-9]*/\1'${n_grid_points}'/' src/control.f90
fi

make MPIF90=$mpif90 COMPILER=$compiler DEFINE="$defines" MODE="$mode"
