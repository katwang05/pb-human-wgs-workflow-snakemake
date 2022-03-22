#!/bin/bash

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

# make logs directory if it doesn't exist
mkdir -p logs

# execute snakemake
snakemake \
    --config cohort=${COHORT} \
    --nolock \
    --snakefile workflow/process_cohort.smk \
    2>&1 | tee "logs/process_cohort.${COHORT}.$(date -d 'today' +'%Y%m%d%H%M').log"
