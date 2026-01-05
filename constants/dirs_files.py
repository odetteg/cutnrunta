import os
import sys
from pathlib import Path
import yaml

sys.path.append(str(Path(__file__).resolve().parent.parent))

# ============================================================================
# BASE DIRECTORIES
# ============================================================================
BASE_DIR = Path(__file__).resolve().parent.parent
CONFIG_FILE_PATH = BASE_DIR / "config" / "config.yaml"
RAW_DATA_DIR = BASE_DIR / "data_cp/data/sequenced_2025_11_13"
RESULTS_DIR = BASE_DIR / "results"
LOGS_DIR = BASE_DIR / "logs"


CONFIG_FILE_PATH = BASE_DIR / "config" / "config.yaml"

with open(CONFIG_FILE_PATH, "r") as confile:
    config = yaml.safe_load(confile)

# ============================================================================
# REFERENCE GENOME
# ============================================================================
REF_DIR = BASE_DIR / "Ref_genome"
BT9_TA_REF_FA = REF_DIR / "BT9_TA_prefixed.fna"
BT9_TA_GTF = REF_DIR / "BT9_TA_fused.gtf"
GTF = REF_DIR / "BT9_TA_fused.gtf"

# ============================================================================
# QC & TRIMMING (FASTP, FASTQC)
# ============================================================================
# Directories
RAW_FASTQC_DIR = RESULTS_DIR / "fastqc/raw"
FASTQC_QC_FASTP_DIR = RESULTS_DIR / "fastqc/trim/fastp"
FASTP_TRIMMED_DIR = RESULTS_DIR / "trim/fastp"
FASTP_QC_REPORTS_DIR = RESULTS_DIR / "qc/fastp"
QC_DIR = RESULTS_DIR / "qc"

# ============================================================================
# ALIGNMENT & FILTERING (BOWTIE2, SAMTOOLS)
# ============================================================================
# Directories
bt2_raw_map_dir = RESULTS_DIR / "map/raw"
temp_sort_dir = RESULTS_DIR / "map/temp"
samtools_sorted_filtered_dir = RESULTS_DIR / "map/filtered_sorted"

# ============================================================================
# BAM PROCESSING (PICARD, SAMTOOLS)
# ============================================================================
# Directories
mark_remove_dups = RESULTS_DIR / "bam/mark_remove_dups"
merged_bam_dir = RESULTS_DIR / "bam/merged"
unsorted_merged_bam = RESULTS_DIR / "bam/merged/unsorted_merged"

# Files
samtools_txt = BASE_DIR / "ad_files/samtools.stats.txt"
agg_stats_csv = BASE_DIR / "ad_files/aggegrates.samtools.stats.csv"

# ============================================================================
# STATS & METRICS
# ============================================================================
# Directories
stats_dir = RESULTS_DIR / "stats"
sorted_filtered_stats_dir = stats_dir / "map/filtered_sorted_stats"

# ============================================================================
# BIGWIG GENERATION (DEEPTOOLS)
# ============================================================================
# Directories
bw_dir = RESULTS_DIR / "bigwig"
deeptools_dir = RESULTS_DIR / "deeptools"
bamCompare = bw_dir / "bamCompare"

# ============================================================================
# PEAK CALLING
# ============================================================================
# Directories
macs3_dir = RESULTS_DIR / "macs3"
seacr_dir = RESULTS_DIR / "seacr"
unfiltered_seacr_peaks_dir = seacr_dir / "raw"
filtered_seacr_peaks_dir = seacr_dir / "filtered"
use_control = config.get("peaks_call", {}).get("seacr", {}).get("use_control", False)
temp_peaks_dir = "peaks_wcontrol" if use_control else "peaks_nocontrol"
dynamic_unfiltered_seacr_dir = unfiltered_seacr_peaks_dir / temp_peaks_dir
dynamic_filtered_seacr_dir = filtered_seacr_peaks_dir / temp_peaks_dir
dynamic_intersected_peaks = seacr_dir / "intersected" / temp_peaks_dir
# dynamic_consensus_peaks = seacr_dir / "consensus" / temp_peaks_dir
dynamic_greedy_consensus_peaks = seacr_dir / "consensus" / "greedy" / temp_peaks_dir
dynamic_conservative_consensus_peaks = (
    seacr_dir / "consensus" / "conservative" / temp_peaks_dir
)

dynamic_reciprocal_consensus_peaks = (
    seacr_dir / "consensus" / "reciprocal" / temp_peaks_dir
)
analysis_dir = RESULTS_DIR / "analysis"
#===========================================================================
#CHIPseeker ANNOTATION
#===========================================================================
ChIPseeker_dir= RESULTS_DIR / "ChIPseeker"
ChIPseeker_seacr_dir= RESULTS_DIR / "ChIPseeker" / "seacr"
dynamic_chipseeker_seacr_raw_dir= ChIPseeker_seacr_dir / "raw" / temp_peaks_dir

# ============================================================================
# AUXILIARY FILES
# ============================================================================
# Directories
ad_files = BASE_DIR / "ad_files"

# Files
sample_egs_csv = ad_files / "sample_egs.csv"
sample_lens_csv = ad_files / "sample_lens.csv"

# ============================================================================
# PLOTS & VISUALIZATIONS
# ============================================================================
# Directories
plots = RESULTS_DIR / "plots"

# Files
frag_len_dst_pdf = plots / "samtools/aggegrates.samtools.fraglens.pdf"
frag_len_dst_facet_pdf = plots / "samtools/aggegrates.samtools.fraglens_facet.pdf"
