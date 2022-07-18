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

# USAGE: qsub workflow/process_cohort.sge.sh <cohort_id>

COHORT=$1

# set umask to avoid locking each other out of directories
umask 002

# add lockfile to directory to prevent multiple simultaneous jobs
mkdir -p cohorts/${COHORT}/
LOCKFILE=cohorts/${COHORT}/process_cohort.lock
lockfile -r 0 ${LOCKFILE} || exit 1
trap "rm -f ${LOCKFILE}; exit" SIGINT SIGTERM ERR EXIT

# get variables from workflow/variables.env
source workflow/variables.env

# execute snakemake
snakemake \
    --config "cohort='$COHORT'" \
    --nolock \
    --profile workflow/profiles/sge \
    --snakefile workflow/process_cohort.smk
