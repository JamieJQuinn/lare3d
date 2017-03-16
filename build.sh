#!/usr/bin/env bash
machine='office'
n_grid_points='false'
braginskii='false'
output='true'
defines=''
mode=''

while getopts "bdom:p:" flag; do
  case "${flag}" in
    # Enable braginskii visc
    b) braginskii='true' ;;
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

if [ "$output" == 'false' ]; then
  defines+=' -DNO_IO'
fi

if [ "$machine" == "euclid" ]; then
  n_proc=10
  export PATH=/usr/lib64/openmpi/bin:$PATH
  compiler='gfortran'
elif [ "$machine" == "office" ]; then
  n_proc=4
  compiler='gfortran'
else 
  echo "Error: No machine specified"
  exit
fi

if [ "$n_grid_points" != 'false' ]; then
  sed -i -e 's/\(n[xyz]_global = \)[0-9]*/\1'${n_grid_points}'/' src/control.f90
fi

make -j $n_proc COMPILER=$compiler DEFINE="$defines" MODE="$mode"
