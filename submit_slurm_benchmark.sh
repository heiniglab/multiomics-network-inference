#!/bin/bash

# default call for SLURM submit with snakemake
target=gather_benchmark_results
log=gather_benchmark_results.log

# submit using snakemake -> check for command
command -v snakemake >/dev/null 2>&1 || { echo "snakemake not available, exiting." ; exit 1 ; }

# actual call
nohup nice snakemake --profile ./profiles/slurm_singlenode \
    ${target} &> ${log} &

echo "Running snakemake. Log-file: ${log}"
