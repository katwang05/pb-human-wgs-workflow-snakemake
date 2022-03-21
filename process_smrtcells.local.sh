#!/bin/bash

# set umask to avoid locking each other out of directories
umask 002

# source variables including temporary directory (TEMP)
source workflow/variables.env

# make logs directory if it doesn't exist
mkdir -p logs

# execute snakemake
snakemake --reason \
    --rerun-incomplete \
    --printshellcmds \
    --keep-going \
    --cores \
    --use-conda --conda-frontend mamba \
    --default-resources "tmpdir='${TEMP}'" \
    --snakefile workflow/process_smrtcells.smk \
    2>&1 | tee "logs/process_smrtcells.$(date -d 'today' +'%Y%m%d%H%M').log"
    