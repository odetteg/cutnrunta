
---

# CUT&RUN Analysis Pipeline

This repository contains a **Snakemake-based pipeline** for processing CUT&RUN sequencing data. The workflow is designed for execution on an HPC cluster. Adjustments can me made based on whatever fine you have.

---

## 🛠 Environment and Dependencies

### Conda Environment Activation

```bash
conda activate /shared/projects/cutnrunta/so_/conda/envs/cut_snake

```

### Environment Setup

The Conda environment and package cache are stored within the project directory. 

```bash

mkdir -p /shared/projects/cutnrunta/so_/conda/envs
mkdir -p /shared/projects/cutnrunta/so_/conda/pkgs

export CONDA_ENVS_PATH=/shared/projects/cutnrunta/so_/conda/envs
export CONDA_PKGS_DIRS=/shared/projects/cutnrunta/so_/conda/pkgs

conda create -p /shared/projects/cutnrunta/so_/conda/envs/cut_snake python=3.10
conda activate /shared/projects/cutnrunta/so_/conda/envs/cut_snake

```

---

## 📂 Repository Structure


```Folder/files structure
.
├── Snakefile
├── smk_driver.sh
├── config/
│   └── config.yaml
├── workflows/
│   └── rules/
│       ├── qc.smk
│       └── idx.smk
├── Ref_genome/
│   └── (reference genome files)
├── data_cp/
│   └── data/
│       └── sequenced_2025_11_13/
│           └── (raw FASTQ files)
└── results/
    ├── trim/
    │   └── fastp/
    ├── qc/
    │   ├── fastp/
    │   └── multiqc/  
    └── fastqc/
        └── trim/
            └── fastp/

```
### Key Directory Definitions

| Variable | Path |
| --- | --- |
| **RAW_DATA_DIR** | `/shared/projects/cutnrunta/so_/data_cp/data/sequenced_2025_11_13` |
| **REF_DIR** | `/shared/projects/cutnrunta/so_/Ref_genome` |
| **RESULTS_DIR** | `/shared/projects/cutnrunta/so_/results` |
| **FASTP_TRIMMED_DIR** | `.../results/trim/fastp` |
| **FASTQC_QC_FASTP_DIR** | `.../results/fastqc/trim/fastp` |

---

## 🚀 Workflow Execution

The pipeline is launched from the project root directory using the Slurm driver script:

```bash
sbatch smk_driver.sh

```

The driver script initializes the environment, defines cluster resources, and invokes Snakemake with the appropriate execution parameters.

### Implemented Workflow Rules

The main `Snakefile` includes modular rule files:

* `qc.smk`: Rules for quality control (read-level assessments and post-processing).
* `idx.smk`: Rules for reference genome indexing and preparation.

```snakefile
include: "workflows/rules/qc.smk"
include: "workflows/rules/idx.smk"

```

---

## 📊 Data & Quality Control Strategy

### Data Overview

The input data consists of **Illumina CUT&RUN sequencing reads**.

* **Initial Assessment:** Overall base quality is high.
* **Known Issues:** Multiple FastQC modules (per-sequence content, adapter content, duplication) may show "Fail" status.

### QC Strategy

1. **Baseline:** An initial MultiQC report is generated on untrimmed data to establish a baseline.
2. **Processing:** Trimming and filtering parameters are applied via `fastp`.
3. **Final Summary:** Once finalized, MultiQC is rerun to provide a consolidated summary of the processed data quality.

---

## 📝 Notes and Future Extensions

The workflow is under active development. Future updates will include:

* Peak calling rules (e.g., SEACR or MACS2).
* Signal visualization (BigWig generation).
* Differential binding analysis.

---
