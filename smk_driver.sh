#!/bin/bash
#SBATCH --job-name=so_cut_run_driver
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=steve.odette@etu.u-paris.fr
#SBATCH --account cutnrunta
#SBATCH --partition=ipop-up
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=master_%j.out
#SBATCH --error=master_%j.err

module load snakemake

snakemake --snakefile snakefile \
    --cores $SLURM_CPUS_PER_TASK \
    --latency-wait 60 \
    --rerun-incomplete \
