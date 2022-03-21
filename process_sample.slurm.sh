#!/bin/bash
#SBATCH -A 100humans
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -o cluster_logs/slurm-%x-%j-%N.out

SAMPLE=$1

# set umask to avoid locking each other out of directories
umask 002

# add lockfile to directory to prevent multiple simultaneous jobs
LOCKFILE=samples/${SAMPLE}/process_sample.lock
lockfile -r 0 ${LOCKFILE} || exit 1
trap "rm -f ${LOCKFILE}; exit" SIGINT SIGTERM ERR EXIT

# set singularity temporary directory
source workflow/variables.env
export SINGULARITY_TMPDIR="$TEMP"
export SINGULARITY_BIND="$TEMP"

# execute snakemake
snakemake \
    --config sample=${SAMPLE} \
    --nolock \
    --local-cores 4 \
    --profile workflow/profiles/slurm \
    --default-resources "tmpdir='${TEMP}'" \
    --snakefile workflow/process_sample.smk
