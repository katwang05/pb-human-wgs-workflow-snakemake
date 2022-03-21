#!/bin/bash
#SBATCH -A 100humans
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -o cluster_logs/slurm-%x-%j-%N.out

# set umask to avoid locking each other out of directories
umask 002

# source variables including temporary directory (TEMP)
source workflow/variables.env

# execute snakemake
snakemake \
    --local-cores 1 \
    --profile workflow/profiles/slurm \
    --default-resources "tmpdir='${TEMP}'" \
    --snakefile workflow/process_smrtcells.smk
