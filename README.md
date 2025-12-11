My conda environment: conda activate /shared/projects/cutnrunta/so_/conda/envs/cut_snake
Steps I used to create the env:
mkdir -p /shared/projects/cutnrunta/so_/conda/envs
mkdir -p /shared/projects/cutnrunta/so_/conda/pkgs
export CONDA_ENVS_PATH=/shared/projects/cutnrunta/so_/conda/envs
export CONDA_PKGS_DIRS=/shared/projects/cutnrunta/so_/conda/pkgs
conda create -p /shared/projects/cutnrunta/so_/conda/envs/cut_snake python=3.10
