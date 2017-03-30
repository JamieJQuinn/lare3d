#!/usr/bin/env bash
machine=''
mpi_opts=''
n_proc_arg=''

while getopts "n:m:" flag; do
  case "${flag}" in
    # Set number of processes to launch
    n) n_proc_arg="${OPTARG}" ;;
    # Set machine
    m) machine="${OPTARG}" ;;
  esac
done

if [ "$machine" == "euclid" ]; then
  n_proc=15
  export PATH=/usr/lib64/openmpi/bin:$PATH
  mpi_opts+='--mca pml ob1'
elif [ "$machine" == "office" ]; then
  n_proc=4
else 
  if [ "$n_proc_arg" != '' ]; then
    n_proc=$n_proc_arg
  else
    echo "Error: No machine or number of nodes specified"
    exit
  fi
fi


time mpirun -np $n_proc $mpi_opts bin/lare3d
