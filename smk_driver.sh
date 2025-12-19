#!/bin/bash
#SBATCH --job-name=so_cut_run_driver
#SBATCH --account cutnrunta
#SBATCH --partition=ipop-up
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=masterRun.out
#SBATCH --error=masterRun.err

module load snakemake

snakemake --snakefile snakefile \
    --cores $SLURM_CPUS_PER_TASK \
    --latency-wait 60 \
    --rerun-incomplete \
