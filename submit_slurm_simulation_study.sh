#!/bin/bash

# default call for SLURM submit with snakemake
# targetting the cohort study

echo "Checking for snakemake ------------------------------"

command -v snakemake >/dev/null 2>&1 || { echo "snakemake not available, exiting." ; exit 1 ; }


echo "Executing snakemake ---------------------------------"

target=all_simulation
logFile=./${target}.log

nohup nice snakemake --profile=./profiles/slurm \
  ${target} &> ${logFile} &

echo "Snakemake running -----------------------------------"
echo "Log-file: ${logFile}"
