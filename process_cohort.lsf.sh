#!/bin/bash
#BSUB -P 100humans
#BSUB -cwd
#BSUB -L /bin/bash
#BSUB -q default
#BSUB -n 4
#BSUB -o ./cluster_logs/lsf-$LSB_JOBNAME-$LSB_JOBID-$HOSTNAME.out
#BSUB -e ./cluster_logs/lsf-$LSB_JOBNAME-$LSB_JOBID-$HOSTNAME.err

# USAGE: bsub workflow/process_cohort.lsf.sh <cohort_id>

COHORT=$1

# set umask to avoid locking each other out of directories
umask 002

# add lockfile to directory to prevent multiple simultaneous jobs
COHORTDIR="cohorts/${COHORT}"
LOCKFILE="$COHORTDIR/process_cohort.lock"

if [ -f "$LOCKFILE" ]; then
    echo "lockfile $LOCKFILE already exists. Remove lockfile and try again." && exit 1
else
    mkdir -p "$COHORTDIR" || exit 1
    touch "$LOCKFILE" || exit 1
fi
trap 'rm -f ${LOCKFILE}; exit' SIGINT SIGTERM ERR EXIT

# get variables from workflow/variables.env
source workflow/variables.env

# execute snakemake
snakemake \
    --config "cohort='$COHORT'" \
    --nolock \
    --profile workflow/profiles/lsf \
    --snakefile workflow/process_cohort.smk
