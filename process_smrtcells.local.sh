#!/bin/bash

# set umask to avoid locking each other out of directories
umask 002

mkdir -p logs

# execute snakemake
snakemake --reason \
    --rerun-incomplete \
    --printshellcmds \
    --keep-going \
    --cores \
    --use-conda --conda-frontend mamba \
    --default-resources "tmpdir=system_tmpdir" \
    --snakefile workflow/process_smrtcells.smk \
    2>&1 | tee "logs/process_smrtcells.$(date -d 'today' +'%Y%m%d%H%M%S').log"
