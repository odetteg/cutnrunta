import os
from pathlib import Path
import pandas as pd
import yaml
import sys
import re

sys.path.append(str(Path(__file__).resolve().parent.parent))
from constants.common import *


BASE_DIR = Path(__file__).resolve().parent.parent
CONFIG_FILE_PATH = BASE_DIR / "config" / "config.yaml"

RAW_DATA_DIR: BASE_DIR / "data_cp/data/sequenced_2025_11_13"
REF_DIR: BASE_DIR / "Ref_genome"
RESULTS_DIR: BASE_DIR / "results"
FASTP_TRIMMED_DIR: BASE_DIR / "results/trim/fastp"
FASTP_QC_REPORTS_DIR: BASE_DIR / "results/qc/fastp"
FASTQC_QC_FASTP_DIR: BASE_DIR / "results/fastqc/trim/fastp"
samtools_txt = BASE_DIR / "samtools.stats.txt"
agg_stats_csv = BASE_DIR / "aggegrates.samtools.stats.csv"
frag_len_dst_pdf = BASE_DIR / "aggegrates.samtools.fraglens.pdf"
frag_len_dst_facet_pdf = BASE_DIR / "aggegrates.samtools.fraglens_facet.pdf"
