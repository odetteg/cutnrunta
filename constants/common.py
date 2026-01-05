import os
from pathlib import Path
import pandas as pd
import yaml
import re
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))
from constants.dirs_files import *


# --------------------------
# Load configuration
# --------------------------

BASE_DIR = Path(__file__).resolve().parent.parent
CONFIG_FILE_PATH = BASE_DIR / "config" / "config.yaml"

with open(CONFIG_FILE_PATH, "r") as yamfile:
    config = yaml.safe_load(yamfile)


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
        .rename(columns={"pair_id": "sample_id", "R1": "fwd_read", "R2": "rev_read"})
    )

    sample_metadata = sample_metadata.set_index("sample_id")

    return sample_metadata


sample_metadata = build_sample_metadata(samples_df)


def get_paired_rds(
    sample_id, sample_metadata=sample_metadata, raw_data_dir=RAW_DATA_DIR
):
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

    return (os.path.join(raw_data_dir, fwd), os.path.join(raw_data_dir, rev))


base_ids = sample_metadata.index.tolist()


def bw_input(wildcards):
    """
    Returns named input files as dictionary for bigwig rule.
    """
    bw_dict_ = {
        "bam": mark_remove_dups / f"{wildcards.base_id}.deduped.bam",
        "bam_bai": mark_remove_dups / f"{wildcards.base_id}.deduped.bam.bai",
        "sample_egs_csv": ad_files / "sample_egs.csv",
    }

    return bw_dict_


def parse_base_id(base_id):
    """
    Parse a single base_id into cell_line, antibody, replicate, and rest_info.
    """
    cell_line, rest = base_id.split("-", 1)

    match = re.match(r"([A-Za-z0-9]+)(?:-(\d+))?_(.*)", rest)
    if match:
        antibody = match.group(1)
        replicate = match.group(2)
        rest_info = match.group(3)
    else:
        antibody = rest
        replicate = None
        rest_info = ""

    return pd.Series([cell_line, antibody, replicate, rest_info])


baseid_df = pd.DataFrame({"base_id": base_ids})

baseid_df[["cell_line", "antibody", "replicate", "rest_info"]] = baseid_df[
    "base_id"
].apply(parse_base_id)


baseid_df["bw_path"] = baseid_df["base_id"].apply(lambda x: str(bw_dir / f"{x}.bw"))

grouped = (
    baseid_df.groupby(["cell_line", "antibody"])["bw_path"].apply(list).reset_index()
)

cell_lines = baseid_df["cell_line"].unique().tolist()
antibodies = baseid_df["antibody"].unique().tolist()


def parse_study_design(baseid_df, control_patterns):
    """
    Generate study design DataFrame with control vs treatment labels.

    Args:
        baseid_df (pd.DataFrame): Must contain 'base_id' and 'antibody' columns.
        control_patterns (list of str): List of substrings identifying controls.

    Returns:
        pd.DataFrame: Copy of baseid_df with added 'condition' column.
    """
    control_abs = {
        s for s in baseid_df["base_id"] if any(p in s for p in control_patterns)
    }
    df = baseid_df.copy()
    df["condition"] = df["base_id"].apply(
        lambda x: "control" if x in control_abs else "treatment"
    )
    return df


control_patterns = config["controls"]["antibodies"]
design_df = parse_study_design(baseid_df=baseid_df, control_patterns=control_patterns)


def get_ctrl_bam(sample_id):
    cell_line = design_df.loc[design_df.base_id == sample_id, "cell_line"].values[0]
    control_base = design_df.query('condition=="control" and cell_line==@cell_line')[
        "base_id"
    ].values[0]
    return mark_remove_dups / f"{control_base}.deduped.bam"


def get_ctrl_bai(sample_id):
    cell_line = design_df.loc[design_df.base_id == sample_id, "cell_line"].values[0]
    control_base = design_df.query('condition=="control" and cell_line==@cell_line')[
        "base_id"
    ].values[0]
    return mark_remove_dups / f"{control_base}.deduped.bam.bai"


def get_ctrl_bg(sample_id):
    """Get the control bedgraph file for a given treatment sample."""
    try:
        cell_line = design_df.loc[design_df.base_id == sample_id, "cell_line"].values[0]
        control_base = design_df.query(
            'condition == "control" and cell_line == @cell_line'
        )["base_id"].values[0]
        return str(RESULTS_DIR / f"bedtools/{control_base}.fragments.bedgraph")
    except IndexError:
        raise ValueError(f"Could not find control for sample '{sample_id}'")


def get_peak_ids(df):
    return df.loc[df["condition"] == "treatment", "base_id"].tolist()


peak_ids = get_peak_ids(design_df)
temp_peak_idf = design_df.copy()
temp_peak_idf = temp_peak_idf[temp_peak_idf["condition"] == "treatment"]

# May change this later because of the possible combination of replicates
rep1_id = str(config.get("reps", {}).get("rep1"))
rep2_id = str(config.get("reps", {}).get("rep2"))


temp_peak_idf["gen_abd"] = temp_peak_idf["cell_line"] + "-" + temp_peak_idf["antibody"]
design_df["gen_abd"] = design_df["cell_line"] + "-" + design_df["antibody"]

# rep_df = (
#     temp_peak_idf.pivot(index="gen_abd", columns="replicate", values="base_id")
#     .rename(columns={1: "rep1", 2: "rep2"})
#     .reset_index()
# )
rep_df = (
    temp_peak_idf.pivot(index="gen_abd", columns="replicate", values="base_id")
    .reset_index()
)

gen_abd_ids: list = rep_df["gen_abd"].tolist()
ctr_abd_ids = design_df.query("condition == 'control'")["gen_abd"].tolist()
tt_base_ids = design_df.query("condition == 'treatment'")["base_id"].tolist()

def get_reps(gen_abd_id: str, rep_df=rep_df) -> tuple:
    row = rep_df.loc[rep_df["gen_abd"] == gen_abd_id]
    if row.empty:
        raise ValueError(f"Warning!: {gen_abd_id} not found in rep_df")
    return tuple(row[["1", "2"]].iloc[0])  



def get_control_id(gen_abd_id: str) -> str:
    """Get control sample ID for given treatment sample."""
    cell_line = design_df.loc[design_df.gen_abd == gen_abd_id, "cell_line"].values[0]
    control_id = design_df.query('condition=="control" and cell_line==@cell_line')[
        "gen_abd"
    ].values[0]
    return control_id


def get_ctrl_merged_bam(gen_abd_id: str) -> Path:
    """Get merged BAM path for control sample."""
    return merged_bam_dir / f"control/{get_control_id(gen_abd_id)}.merged.bam"


def get_ctrl_merged_bai(gen_abd_id: str) -> Path:
    """Get merged BAI path for control sample."""
    return merged_bam_dir / f"control/{get_control_id(gen_abd_id)}.merged.bam.bai"


# This is a patchwork for now to get things going. Might refactor if a think of a better solution
def get_merge_control_id(ctr_abd_id: str) -> str:
    """Get control base id for given treatment antibody ID"""
    cell_line = design_df.loc[design_df.gen_abd == ctr_abd_id, "cell_line"].values[0]
    merge_control_id = design_df.query(
        'condition=="control" and cell_line==@cell_line'
    )["base_id"].values[0]
    return merge_control_id


def get_ctrl_deduped_bam(ctr_abd_id: str) -> Path:
    """Get deduped BAM path for control sample."""
    return mark_remove_dups / f"{get_merge_control_id(ctr_abd_id)}.deduped.bam"


def get_ctrl_deduped_bai(ctr_abd_id: str) -> Path:
    """Get deduped BAM path for control sample."""
    return mark_remove_dups / f"{get_merge_control_id(ctr_abd_id)}.deduped.bam.bai"




# print(get_ctrl_merged_bam("TBL3-H2K119Ub"))
# print(design_df.query("condition == 'treatment'")["base_id"].tolist())
# print(design_df)
# print(get_reps("TBL3-H2K119Ub"))
# print(design_df.query("condition == 'control'")["gen_abd"].tolist())
# print(design_df["base_id"].str.contains("TBL3-H2K119Ub"))
# print(design_df.query("condition == 'treatment'")["base_id"])
# print(gen_abd_ids)
# print(f"ctr_abd_ids: {ctr_abd_ids}")
# print(get_ctrl_merged_bam("TBL3-H2K119Ub"))
# print(design_df.query("condition == 'treatment'")["base_id"])
# print(design_df.query("condition == 'treatment'")["gen_abd"].tolist())
# print(gen_abd_ids)
# print(design_df)
# for id in gen_abd_ids:
#     print(f"{id}: {get_ctrl_merged_bam(id)}")
# print(tt_base_ids)
# print("ID overlap:", set(tt_base_ids) & set(gen_abd_ids))
