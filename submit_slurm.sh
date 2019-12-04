#!/bin/bash

# default call for SLURM submit with snakemake
target=$1
log=$2

if [[ -z ${log} ]] ; then
  log=$(basename ${target}).out
fi

# actual call
nohup nice snakemake -u configs/slurm.json --jobs=150 --local-cores=1 \
    --cluster "sbatch -t {cluster.time} -c {cluster.cpu} --mem-per-cpu {cluster.mem} \
    -p {cluster.partition}" \
    ${target} &> ${log} &

echo "Running snakemake. Log-file: ${log}"
