#!/bin/bash
#SBATCH -A 100humans
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -o cluster_logs/slurm-%x-%j-%N.out

# USAGE: sbatch workflow/process_smrtcells.slurm.sh

# set umask to avoid locking each other out of directories
umask 002

# get variables from workflow/variables.env
source workflow/variables.env

# execute snakemake
snakemake \
    --profile workflow/profiles/slurm \
    --snakefile workflow/process_smrtcells.smk
