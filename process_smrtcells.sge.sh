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

# set umask to avoid locking each other out of directories
umask 002

# get variables from workflow/variables.env
source workflow/variables.env

# execute snakemake
snakemake \
    --profile workflow/profiles/sge \
    --snakefile workflow/process_smrtcells.smk
