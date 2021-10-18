#!/bin/bash

# default call for SLURM submit with snakemake
target=$1
log=$2

if [[ -z ${log} ]] ; then
  log=$(basename ${target}).out
fi

# submit using snakemake -> check for command
command -v snakemake >/dev/null 2>&1 || { echo "snakemake not available, exiting." ; exit 1 ; }

# actual call
nohup nice snakemake --profile ./profiles/slurm \
    ${target} &> ${log} &

echo "Running snakemake. Log-file: ${log}"
