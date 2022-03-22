#!/bin/bash

# set umask to avoid locking each other out of directories
umask 002

# get variables from workflow/variables.env
source workflow/variables.env

# make logs directory if it doesn't exist
mkdir -p logs

# execute snakemake
snakemake \
    --profile workflow/profiles/local \
    --snakefile workflow/process_smrtcells.smk \
    2>&1 | tee "logs/process_smrtcells.$(date -d 'today' +'%Y%m%d%H%M').log"
    