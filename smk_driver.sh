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


module load samtools/1.21


SNAKEMAKE_ENV="${PWD}/.snakemake_env"
if [ ! -d "$SNAKEMAKE_ENV" ]; then
    echo "Creating snakemake environment with updated conda..."

    module load snakemake
    conda create -y -p $SNAKEMAKE_ENV -c conda-forge -c bioconda \
        "conda>=24.7.1" \
        snakemake
    module unload snakemake
fi


export PATH="${SNAKEMAKE_ENV}/bin:$PATH"

export CONDA_PREFIX=$SNAKEMAKE_ENV

# 4. Verify versions
echo "Using conda version:"
conda --version
echo "Using snakemake version:"
snakemake --version
echo "Conda location:"
which conda

export CONDA_PKGS_DIRS=$PWD/.conda_pkgs
mkdir -p $CONDA_PKGS_DIRS
conda config --set channel_priority strict

# Download and install SEACR

if [ ! -d ~/tools ]; then
    mkdir -p ~/tools
fi
if [ ! -d ~/tools/SEACR ]; then
    git clone https://github.com/FredHutch/SEACR.git ~/tools/SEACR
fi
chmod +x ~/tools/SEACR/SEACR_1.3.*


snakemake --snakefile snakefile \
    --cores $SLURM_CPUS_PER_TASK \
    --use-conda \
    --conda-prefix /shared/projects/cutnrunta/so_/.snakemake/conda-envs \
    --conda-create-envs-only

snakemake --snakefile snakefile \
    --cores $SLURM_CPUS_PER_TASK \
    --use-conda \
    --conda-prefix /shared/projects/cutnrunta/so_/.snakemake/conda-envs \
    --latency-wait 60 \
    --rerun-incomplete \
    --printshellcmds