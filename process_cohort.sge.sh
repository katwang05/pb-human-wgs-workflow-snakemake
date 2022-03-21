#!/bin/bash
#$ -A 100humans
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default
#$ -pe smp 4
#$ -o ./cluster_logs/sge-$JOB_NAME-$JOB_ID-$HOSTNAME.out
#$ -e ./cluster_logs/sge-$JOB_NAME-$JOB_ID-$HOSTNAME.err

COHORT=$1

# set umask to avoid locking each other out of directories
umask 002

# add lockfile to directory to prevent multiple simultaneous jobs
mkdir -p cohorts/${COHORT}/
LOCKFILE=cohorts/${COHORT}/process_cohort.lock
lockfile -r 0 ${LOCKFILE} || exit 1
trap "rm -f ${LOCKFILE}; exit" SIGINT SIGTERM ERR EXIT

# set singularity temporary directory
source workflow/variables.env
export SINGULARITY_TMPDIR="$TEMP"
export SINGULARITY_BIND="$TEMP"

# execute snakemake
snakemake \
    --config cohort=${COHORT} \
    --nolock \
    --local-cores 4 \
    --profile workflow/profiles/sge \
    --default-resources "tmpdir='${TEMP}'" \
    --snakefile workflow/process_cohort.smk
