#!/bin/bash
#SBATCH --job-name=so_cut_run
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=snakemake.out
#SBATCH --error=snakemake.err

module load snakemake

snakemake --snakefile snakefile --cores 4 --latency-wait 60 --restart-times 3
