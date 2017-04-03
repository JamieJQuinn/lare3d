#!/usr/bin/env bash
build_state_file='bin/build_state'
machine='office'
n_grid_points='false'
braginskii='false'
switching='false'
output='true'
mode=''

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

defines=''
if [ "$braginskii" == 'true' ]; then
  defines+=' -DBRAGINSKII_VISCOSITY'
elif [ "$switching" == 'true' ]; then
  defines+=' -DSWITCHING_VISCOSITY'
fi

if [ "$output" == 'false' ]; then
  defines+=' -DNO_IO'
fi

compiler=''
mpif90='mpif90'
if [ "$machine" == "euclid" ]; then
  export PATH=/usr/lib64/openmpi/bin:$PATH
  compiler='gfortran'
elif [ "$machine" == "archie" ]; then
  module load mpi/intel/impi/5.0.3
  module load compilers/intel/17.0.1
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

echo "" > $build_state_file
echo "$defines" >> $build_state_file
echo "$mode" >> $build_state_file
echo "built for $machine" >> $build_state_file
echo "GRID POINTS:" >> $build_state_file
echo "$(grep 'n[xyz]_global' src/control.f90)" >> $build_state_file
echo "$(grep -o 'dt_snapshots = .*' src/control.f90)" >> $build_state_file
echo "$(grep -o 'nsteps = .*' src/control.f90)" >> $build_state_file
echo "$(grep -o 't_end = .*' src/control.f90)" >> $build_state_file
echo "$(grep -o 'data_dir = .*' src/control.f90)" >> $build_state_file
if grep -q 'dump_mask(20)' src/control.f90 ; then
  echo "ISO: YES" >> $build_state_file
fi
if grep -q 'dump_mask(21)' src/control.f90 ; then
  echo "ANISO: YES" >> $build_state_file
fi
echo "$(grep 'dump_mask([0-9]*\:[0-9]*) =' src/control.f90)" >> $build_state_file
