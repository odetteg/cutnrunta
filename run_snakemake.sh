#!/bin/bash
#SBATCH --job-name=so_cut_run
#SBATCH --partition=ipop-up
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=snakemake.out
#SBATCH --error=snakemake.err

module load snakemake

snakemake --snakefile snakefile \
          --cores 4 \
          --latency-wait 60 \
          --restart-times 1 \
          --rerun-incomplete \
          --printshellcmds \
          --reason
