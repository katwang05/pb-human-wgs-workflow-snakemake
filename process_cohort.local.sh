#!/bin/bash

# USAGE: bash workflow/process_cohort.local.sh <cohort_id>

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

# make logs directory if it doesn't exist
mkdir -p logs

# execute snakemake
snakemake \
    --config "cohort='$COHORT'" \
    --nolock \
    --profile workflow/profiles/local \
    --snakefile workflow/process_cohort.smk \
    2>&1 | tee "logs/process_cohort.${COHORT}.$(date -d 'today' +'%Y%m%d%H%M').log"
