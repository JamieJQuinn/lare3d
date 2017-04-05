#!/usr/bin/env bash
build_state_file='bin/build_state'
machine=''
n_grid_points='false'
braginskii='false'
switching='false'
output='true'
mode=''
defines=''

while getopts "bsdom:p:" flag; do
  case "${flag}" in
    # Enable braginskii visc
    b) defines+=' -DBRAGINSKII_VISCOSITY' ;;
    # Enable Switching visc
    s) defines+=' -DSWITCHING_VISCOSITY' ;;
    # Output continuous viscous heating
    o) defines+=' -DOUTPUT_CONTINUOUS_VISC_HEATING' ;;
    # Enable debug
    d) mode='debug' ;;
    # Set machine
    m) machine="${OPTARG}" ;;
    # Set num grid points
    p) n_grid_points="${OPTARG}" ;;
  esac
done

# Sort out machine specific settings
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
elif [ "$machine" != "" ]; then
  compiler='gfortran'
  echo "Machine specified but no special conditions, using gfortran"
else
  compiler='gfortran'
  echo "Error: No machine specified, building with gfortran"
fi

# Replace grid points in control with inputs
if [ "$n_grid_points" != 'false' ]; then
  sed -i -e 's/\(n[xyz]_global = \)[0-9]*/\1'${n_grid_points}'/' src/control.f90
fi

# Build
make MPIF90=$mpif90 COMPILER=$compiler DEFINE="$defines" MODE="$mode"

# Create state file and print some relevant build facts
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
