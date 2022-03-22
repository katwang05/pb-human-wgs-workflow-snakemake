#!/bin/bash

SAMPLE=$1

# set umask to avoid locking each other out of directories
umask 002

# add lockfile to directory to prevent multiple simultaneous jobs
LOCKFILE=samples/${SAMPLE}/process_sample.lock
lockfile -r 0 ${LOCKFILE} || exit 1
trap "rm -f ${LOCKFILE}; exit" SIGINT SIGTERM ERR EXIT

# get variables from workflow/variables.env
source workflow/variables.env

# make logs directory if it doesn't exist
mkdir -p logs

# execute snakemake
snakemake \
    --config sample=${SAMPLE} cpu_only=True\
    --nolock \
    --snakefile workflow/process_sample.smk \
    2>&1 | tee "logs/process_sample.${SAMPLE}.$(date -d 'today' +'%Y%m%d%H%M').log"
