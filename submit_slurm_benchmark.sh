#!/bin/bash

# default call for SLURM submit with snakemake
target=gather_benchmark_results
log=gather_benchmark_results.log

# submit using snakemake -> check for command
command -v snakemake >/dev/null 2>&1 || { echo "snakemake not available, exiting." ; exit 1 ; }

# actual call
nohup nice snakemake -u configs/slurm.json --jobs=200 -k --local-cores=1 \
    --latency-wait 40 --cluster "sbatch --nice=10000 -w icb-neu-001 -t {cluster.time} -c {cluster.cpu} --mem-per-cpu {cluster.mem} \
    -p {cluster.partition}" \
    ${target} &> ${log} &

echo "Running snakemake. Log-file: ${log}"
