#!/bin/bash
#$ -A 100humans
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default
#$ -o ./cluster_logs/sge-$JOB_NAME-$JOB_ID-$HOSTNAME.out
#$ -e ./cluster_logs/sge-$JOB_NAME-$JOB_ID-$HOSTNAME.err

# set umask to avoid locking each other out of directories
umask 002

# source variables including temporary directory (TEMP)
source workflow/variables.env

# execute snakemake
snakemake \
    --local-cores 1 \
    --profile workflow/profiles/sge \
    --default-resources "tmpdir='${TEMP}'" \
    --snakefile workflow/process_smrtcells.smk
