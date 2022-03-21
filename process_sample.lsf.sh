#!/bin/bash
#BSUB -cwd
#BSUB -L /bin/bash
#BSUB -q default
#BSUB -n 4
#BSUB -o ./cluster_logs/lsf-$LSB_JOBNAME-$LSB_JOBID-$HOSTNAME.out
#BSUB -e ./cluster_logs/lsf-$LSB_JOBNAME-$LSB_JOBID-$HOSTNAME.err

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
    --profile workflow/profiles/lsf \
    --default-resources "tmpdir='${TEMP}'" \
    --snakefile workflow/process_sample.smk
