#!/usr/bin/env bash
machine='office'
braginskii='false'
output='true'
defines=''
mode=''

while getopts "bdom:" flag; do
  case "${flag}" in
		# Enable braginskii visc
    b) braginskii='true' ;;
		# Disable output
		o) output='false' ;;
    # Enable debug
    d) mode='debug' ;;
    # Set machine
    m) machine="${OPTARG}" ;;
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
fi

make -j $n_proc COMPILER=$compiler DEFINE="$defines" MODE="$mode"
