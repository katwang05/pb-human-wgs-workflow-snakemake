#!/bin/bash
#BSUB -P 100humans
#BSUB -cwd
#BSUB -L /bin/bash
#BSUB -q default
#BSUB -n 4
#BSUB -o ./cluster_logs/lsf-$LSB_JOBNAME-$LSB_JOBID-$HOSTNAME.out
#BSUB -e ./cluster_logs/lsf-$LSB_JOBNAME-$LSB_JOBID-$HOSTNAME.err

# USAGE: bsub workflow/process_smrtcells.lsf.sh

# set umask to avoid locking each other out of directories
umask 002

# get variables from workflow/variables.env
source workflow/variables.env

# execute snakemake
snakemake \
    --profile workflow/profiles/lsf \
    --snakefile workflow/process_smrtcells.smk
