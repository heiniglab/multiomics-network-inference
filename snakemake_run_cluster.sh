#!/bin/bash

target_rule=$1
log_file=$2
max_jobs=$3

# set default for optional 'max_jobs' param
if [ -z "$max_jobs" ]; then
  max_jobs=200
fi

nohup nice snakemake -u configs/slurm.json --jobs=${max_jobs} -k --local-cores=1 \
  --cluster "sbatch -t {cluster.time} -c {cluster.cpu} --mem-per-cpu {cluster.mem} -p {cluster.partition}" \
  ${target_rule} &> ${log_file} &
