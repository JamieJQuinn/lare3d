#!/usr/bin/env bash

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null && pwd  )"

time_end=500.0
dt_snapshots=50.0

nx=500
ny=$nx
nz=2*$nx

#nx=512
#ny=512
#nz=512

visc_mantissa='1'
resist_mantissa='1'

for n in 4; do
  for m in 4; do
    #for switching_value in 0.5; do
    #for nx in 250 300 350 400 450 500; do
    for visc in '-isotropic' '-switching'; do
      viscosity=$visc_mantissa'.0_num*10.0_num**(-'$n'_num)'
      resistivity=$resist_mantissa'.0_num*10.0_num**(-'$m'_num)'
      #folder=v0r5e-$m$visc
      #folder=v1e-${n}r0$visc
      folder=v${visc_mantissa}e-${n}r${resist_mantissa}e-${m}$visc
      #folder=v${visc_mantissa}e-${n}r${resist_mantissa}e-$m$visc
      cp -rf $SCRIPT_DIR/lare3d $folder
      echo $folder":"
      echo "viscosity = "$viscosity
      echo "resistivity = "$resistivity
      cd $folder
      sed -i -e 's/visc3 =.*/visc3 = '$viscosity'/' src/control.f90
      sed -i -e 's/eta_background =.*/eta_background = '$resistivity'/' src/control.f90
      sed -i -e 's/\(t_end = \).*\(_num\)/\1'${time_end}'\2/' src/control.f90
      sed -i -e 's/\(dt_snapshots = \).*\(_num\)/\1'${dt_snapshots}'\2/' src/control.f90
      sed -i -e 's/\(nx_global = \).*/\1'${nx}'/' src/control.f90
      sed -i -e 's/\(ny_global = \).*/\1'${ny}'/' src/control.f90
      sed -i -e 's/\(nz_global = \).*/\1'${nz}'/' src/control.f90
      # Remove centred magnetic field
      #sed -i -e 's/\(dump_mask(22:25) = \).*/\1.false./' src/control.f90
      # Remap kinetic energy
      #sed -i -e 's/\(rke = \).*/\1.false./' src/control.f90
      # Fix switching
      #sed -i -e 's/\(fix_switching = \).*/\1.true./' src/control.f90
      #sed -i -e 's/\(fixed_switching_function = \).*/\1'${switching_value}'_num/' src/control.f90
      # Anomolous resistivity
      sed -i -e 's/j_max =.*/j_max = 5.0_num/' src/control.f90
      sed -i -e 's/eta0 =.*/eta0 = 1.0e-3_num/' src/control.f90
      # Change size of box
      #sed -i -e 's/\(y_max = \).*/\1'2.05_num'/' src/control.f90
      #sed -i -e 's/\(x_max = \).*/\1'2.1_num'/' src/control.f90
      # Shock viscosity?
      #sed -i -e 's/\(visc1 = \).*/\1'0.1_num'/' src/control.f90
      #sed -i -e 's/\(visc2 = \).*/\1'0.5_num'/' src/control.f90
      # Restart?
      #sed -i -e 's/\(initial = \).*/\1IC_RESTART/' src/control.f90
      #sed -i -e 's/\(restart_snapshot = \).*/\12/' src/control.f90
      cd ..
    done
  #done
  done
done
