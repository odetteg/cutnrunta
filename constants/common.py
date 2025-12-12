import os
from pathlib import Path
import pandas as pd
import yaml

# --------------------------
# Load configuration
# --------------------------

BASE_DIR = Path(__file__).resolve().parent.parent
CONFIG_FILE_PATH = BASE_DIR / "config" / "config.yaml"

with open(CONFIG_FILE_PATH, 'r') as yamfile:
    config = yaml.safe_load(yamfile)

RAW_DATA_DIR = Path(config.get("RAW_DATA_DIR"))
REF_DIR = Path(config.get("REF_DIR"))
SO_DIR = Path(config.get("SO_DIR"))
RESULTS_DIR = Path(config.get("RESULTS_DIR"))

# FASTQC and logs directories
FASTQC_DIR = RESULTS_DIR / "fastqc"
LOGS_DIR = SO_DIR / "logs"

# --------------------------
# Load samples
# --------------------------
samples_file = BASE_DIR / "samples.txt"
samples_df = pd.read_csv(samples_file, header=None, names=["filename"])

# Extract the sample IDs. The sample path just without the fastq exts
sample_ids = []

for filename in samples_df["filename"]:
    p = Path(filename)
    sample_id = p.name
    for _ in range(len(p.suffixes)):
        sample_id = Path(sample_id).stem
    sample_ids.append(sample_id)

# --------------------------
# Sanity check
# --------------------------

if not sample_ids:
    raise ValueError("No sample IDs were loaded from samples.txt")

