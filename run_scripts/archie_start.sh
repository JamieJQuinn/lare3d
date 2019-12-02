#!/usr/bin/env bash

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" >/dev/null && pwd  )"

for folder in "$@"
do
  cd $folder
  cp $SCRIPT_DIR/start-lare3d.sh .
  sbatch -J $folder start-lare3d.sh
  cd ..
done
