---

# CUT&RUN Analysis Pipeline

This repository contains a **Snakemake-based pipeline** for processing CUT&RUN sequencing data. The workflow is designed for execution on an HPC cluster and can be adjusted based on your specific requirements.

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
```
.
├── Snakefile
├── smk_driver.sh
├── config/
│   └── config.yaml
├── workflows/
│   ├── rules/
│   └── scripts/
├── envs/
├── constants/
├── Ref_genome/
├── data_cp/
│   └── data/
│       └── sequenced_2025_11_13/
└── results/
    ├── trim/
    ├── qc/
    ├── map/
    ├── bam/
    ├── bigwig/
    ├── deeptools/
    ├── peaks/
    ├── stats/
    └── plots/
```

### Key Directory Definitions

| Variable | Path |
|----------|------|
| **RAW_DATA_DIR** | `/shared/projects/cutnrunta/so_/data_cp/data/sequenced_2025_11_13` |
| **REF_DIR** | `/shared/projects/cutnrunta/so_/Ref_genome` |
| **RESULTS_DIR** | `/shared/projects/cutnrunta/so_/results` |

---

## 🚀 Workflow Execution

The pipeline is launched from the project root directory using the Slurm driver script:
```bash
sbatch smk_driver.sh
```

The driver script initializes the environment, defines cluster resources, and invokes Snakemake with the appropriate execution parameters.

---

## 📊 Pipeline Overview

### Workflow Steps

1. **Quality Control**
   - FastQC on raw reads
   - Adapter trimming and quality filtering with fastp
   - FastQC on trimmed reads

2. **Reference Genome Preparation**
   - Bowtie2 index building
   - Genome indexing (samtools faidx)

3. **Read Mapping**
   - Bowtie2 alignment (paired-end, local mode)
   - Filtering (MAPQ ≥ 10, properly paired reads)
   - Sorting and indexing

4. **Duplicate Handling**
   - Mark/remove duplicates with Picard
   - Generate deduplication metrics

5. **Read Shifting**
   - ATAC-seq style Tn5 offset correction (+4/-5 bp)
   - Produces both unshifted and shifted BAM files

6. **Coverage Analysis**
   - Effective genome size calculation
   - Genome coverage statistics with mosdepth
   - BigWig generation (RPGC-normalized)
   - Average BigWig tracks across replicates

7. **Multi-Sample QC**
   - Fragment length distributions
   - Signal enrichment fingerprints
   - PCA for sample clustering
   - Correlation analysis

8. **Peak Calling**
   - MACS3 peak calling (narrow or broad mode)
   - Performed on both shifted and unshifted alignments
   - Treatment vs control comparisons

9. **Report Generation**
   - MultiQC summary report

---

## 🔬 Peak Calling Configuration

### Switching Between Narrow and Broad Modes

Edit `config/config.yaml`:
```yaml
peaks_call:
  macs3:
    use_macs3: true
    broad: False  # Set to True for broad peaks
    qvalue: 0.05
    broad_cutoff: 0.1
```

**Mode Selection:**
- `broad: False` → Narrow peaks (transcription factors, sharp binding sites)
- `broad: True` → Broad peaks (histone modifications, diffuse signals)

The pipeline automatically generates peaks from both unshifted and shifted BAM files for comparison.

---

## 📈 Output Files

### Key Results

- **QC Reports:** `results/qc/multiqc/multiqc.html`
- **Trimmed Reads:** `results/trim/fastp/`
- **Aligned BAMs:** `results/bam/mark_remove_dups/` and `results/bam/shifted/`
- **BigWig Tracks:** `results/bigwig/unshifted/` and `results/bigwig/shifted/`
- **Peaks:** `results/peaks/macs3/unshifted/` and `results/peaks/macs3/shifted/`
- **Statistics:** `results/stats/`
- **Plots:** `results/plots/`

---

## 📝 Important Notes

### Deduplication
The deduped BAM files are processed by Picard MarkDuplicates. Check `config["picard"]["rmv_duplicates"]` in the config file to confirm whether duplicates were removed or just marked.

### Effective Genome Size
Effective genome size is calculated per sample based on filtered read lengths and mappable regions, providing more accurate normalization than a single genome-wide value.

### Mitochondrial Sequences
Mitochondrial sequence handling is configurable via the `deeptools` parameters in the config file.


### Recommended
I recommend running macs3 with use_control false for cut n run. Cut n run has low noise and often have more read depth in treatment than control. MACS3 will try and overcorrect for this leading to no peaks.


---

## 📚 Important References

### Read Shifting in ATAC-seq
Yan et al. (2020). "From reads to insight: a hitchhiker's guide to ATAC-seq data analysis"  
*Genome Biology* 21:22. https://doi.org/10.1186/s13059-020-1929-3

### CUT&RUN Methodology
Skene & Henikoff (2017). "An efficient targeted nuclease strategy for high-resolution mapping of DNA binding sites"  
*eLife* 6:e21856. https://doi.org/10.7554/eLife.21856

### Comparison of the precicion of various peak calling tools
Nooranikhojasteh, A., Tavallaee, G., & Orouji, E. (2025). Benchmarking Peak Calling Methods for CUT&RUN. Bioinformatics, btaf375. https://doi.org/10.1093/bioinformatics/btaf375

---

## 🔧 Configuration

All pipeline parameters are controlled via `config/config.yaml`, including:
- Input/output directories
- Tool-specific parameters
- Resource allocation (CPU, memory, runtime)
- Peak calling settings

Modify the config file to customize the pipeline without changing the workflow code.

---

## 🚧 Future Development

Planned features include:
- I will include a functionality to remove the mitochondrial sequences. They can be significantly enriched.
- Differential binding analysis
- Motif enrichment analysis
- Additional peak callers (SEACR)
- Automated peak annotation
- A binning size of 10 or 15 gives a good visualization in the bigwigs. It does not, however, affect the downstream analysis. In the IGV, try autoscall and not log scale.

---