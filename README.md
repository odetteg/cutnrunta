My conda environment: conda activate /shared/projects/cutnrunta/so_/conda/envs/cut_snake

Steps I used to create the env:

mkdir -p /shared/projects/cutnrunta/so_/conda/envs
mkdir -p /shared/projects/cutnrunta/so_/conda/pkgs
export CONDA_ENVS_PATH=/shared/projects/cutnrunta/so_/conda/envs
export CONDA_PKGS_DIRS=/shared/projects/cutnrunta/so_/conda/pkgs
conda create -p /shared/projects/cutnrunta/so_/conda/envs/cut_snake python=3.10
conda activate /shared/projects/cutnrunta/so_/conda/envs/cut_snake
#Overview of the data from a first glance:
Most data were of good quality. However, most failed on per sequence content, adapter sequence,
overrepresented sequence, and sequence duplication.
At the beggining, I ran the multiqc before trimming. I will, however, push the multiqc towards the end
once I am confident with the analysis.
