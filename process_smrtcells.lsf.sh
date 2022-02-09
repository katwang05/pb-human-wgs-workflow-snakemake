#!/bin/bash
#BSUB -cwd
#BSUB -L /bin/bash
#BSUB -q default
#BSUB -n 1
#BSUB -o ./cluster_logs/lsf-$LSB_JOBNAME-$LSB_JOBID-$HOSTNAME.out
#BSUB -e ./cluster_logs/lsf-$LSB_JOBNAME-$LSB_JOBID-$HOSTNAME.err

# set umask to avoid locking each other out of directories
umask 002

# execute snakemake
snakemake \
    --local-cores 1 \
    --profile workflow/profiles/lsf \
    --snakefile workflow/process_smrtcells.smk
