#!/bin/bash

# USAGE: bash workflow/process_sample.local.sh <sample_id>

SAMPLE=$1

# set umask to avoid locking each other out of directories
umask 002

# add lockfile to directory to prevent multiple simultaneous jobs
SAMPLEDIR="cohorts/${SAMPLE}"
LOCKFILE="$SAMPLEDIR/process_cohort.lock"

if [ -f "$LOCKFILE" ]; then
    echo "lockfile $LOCKFILE already exists. Remove lockfile and try again." && exit 1
else
    mkdir -p "$SAMPLEDIR" || exit 1
    touch "$LOCKFILE" || exit 1
fi
trap 'rm -f ${LOCKFILE}; exit' SIGINT SIGTERM ERR EXIT

# get variables from workflow/variables.env
source workflow/variables.env

# make logs directory if it doesn't exist
mkdir -p logs

# execute snakemake
snakemake \
    --config "sample='$SAMPLE'" \
    --nolock \
    --profile workflow/profiles/local \
    --snakefile workflow/process_sample.smk \
    2>&1 | tee "logs/process_sample.${SAMPLE}.$(date -d 'today' +'%Y%m%d%H%M').log"
