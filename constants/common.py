import os
from pathlib import Path
import pandas as pd
import yaml
import re

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
RAW_FASTQC_DIR = RESULTS_DIR / "fastqc/raw"
LOGS_DIR = SO_DIR / "logs"

# --------------------------
# Load samples
# --------------------------
samples_file = BASE_DIR / "samples.txt"
samples_df = pd.read_csv(samples_file, header=None, names=["filename"])

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


def build_sample_metadata(samples_df):
    """
    Build a metadata table from a DataFrame of FASTQ filenames.
    Keeps full FASTQ filenames as forward/reverse reads,
    and adds a base_sample_id without the _R1_001/_R2_001 suffix.
    """
    df = samples_df.copy()
    df["basename"] = df["filename"].apply(lambda x: Path(x).name)
    def detect_read(filename):
        if "_R1" in filename:
            return "R1"
        elif "_R2" in filename:
            return "R2"
        else:
            return "Unknown"

    df["read"] = df["basename"].apply(detect_read)

 
    df["pair_id"] = df["basename"].apply(lambda x: re.sub(r"_R[12]_\d+.*$", "", x))

    sample_metadata = (
        df.pivot(index="pair_id", columns="read", values="basename")
        .rename_axis(None, axis=1)  
        .reset_index()
        .rename(columns={
            "pair_id": "sample_id",
            "R1": "fwd_read",
            "R2": "rev_read"
        })
    )

    sample_metadata = sample_metadata.set_index("sample_id")
    

    return sample_metadata

sample_metadata = build_sample_metadata(samples_df)


def get_paired_rds(sample_id, sample_metadata=sample_metadata, raw_data_dir=RAW_DATA_DIR):
    """
    Return full paths to the forward and reverse FASTQ files for a sample.
    Returns (None, None) if sample_id is not found.
    """
    if sample_id not in sample_metadata.index:
        print(f"Warning: sample_id '{sample_id}' not found in metadata.")
        return None, None

    fwd = sample_metadata.loc[sample_id, "fwd_read"]
    rev = sample_metadata.loc[sample_id, "rev_read"]

    if pd.isna(fwd) or pd.isna(rev):
        print(f"Warning: forward or reverse read missing for sample '{sample_id}'.")
        return None, None

    return (
        os.path.join(raw_data_dir, fwd),
        os.path.join(raw_data_dir, rev)
    )

base_ids = sample_metadata.index.tolist()

for b in base_ids:
    fwd, rev = get_paired_rds(b)
    print(b, fwd, rev)
