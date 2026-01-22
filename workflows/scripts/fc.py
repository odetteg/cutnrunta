import os
import pandas as pd
import numpy as np
import sys
import re
from constants.common import *

# ==============================================================================
# Get inputs from Snakemake
# ==============================================================================

raw_counts = snakemake.input.raw_counts
chrom_sizes_file = snakemake.input.chrom_sizes
chromosome_occupancy = snakemake.output.chromosome_occupancy
log = snakemake.log[0] if snakemake.log else None
output_file = str(snakemake.output.chromosome_occupancy)

# Get config
config = snakemake.config

# Get threshold settings
use_control = config["fold_change"]["use_control"]
calculate_foldchange = config["fold_change"]["calculate_foldchange"]

if calculate_foldchange and use_control:
    threshold = config["fold_change"]["threshold"]
else:
    threshold = config["fold_change"]["threshold_wo_ctrl"]

print("=" * 70)
print("Calculating Chromosome Occupancy")
print("=" * 70)
print(f"Threshold: {threshold}")
print(f"Use control: {use_control}")
print(f"Calculate fold change: {calculate_foldchange}")

# ==============================================================================
# Select samples based on output file
# ==============================================================================

if "per_rep" in output_file:
    print(f"Analysis type: Per-replicate (detected from output: {output_file})")
    tt_abds = tt_base_ids
    get_ctr_id = get_control_base_id
elif "avg" in output_file:
    print(f"Analysis type: Averaged (detected from output: {output_file})")
    tt_abds = tt_abd_ids
    get_ctr_id = get_control_id
else:
    print(f"WARNING: Could not detect analysis type from output: {output_file}")
    print(f"Defaulting to averaged analysis")
    tt_abds = tt_abd_ids
    get_ctr_id = get_control_id

print(f"Treatment samples: {tt_abds}")

# ==============================================================================
# Load data
# ==============================================================================

print("\n[1] Loading data...")
counts = pd.read_csv(raw_counts, sep="\t")
chrom_sizes = pd.read_csv(
    chrom_sizes_file, sep="\t", names=["chromosome", "size"], header=None
)

# Clean column names
counts.columns = counts.columns.str.replace("'", "")

# Rename chromosome column
if "#chr" in counts.columns:
    counts = counts.rename(columns={"#chr": "chromosome"})
elif "chr" in counts.columns:
    counts = counts.rename(columns={"chr": "chromosome"})
else:
    first_col = counts.columns[0]
    counts = counts.rename(columns={first_col: "chromosome"})

# ==============================================================================
# Filter valid chromosomes
# ==============================================================================

print("\n[2] Filtering chromosomes...")
valid_chroms = set(counts["chromosome"].unique()) & set(
    chrom_sizes["chromosome"].unique()
)
print(f"    Valid chromosomes: {len(valid_chroms)}")

counts = counts[counts["chromosome"].isin(valid_chroms)]
chrom_sizes = chrom_sizes[chrom_sizes["chromosome"].isin(valid_chroms)]

# ==============================================================================
# Calculate occupancy & fold-change per bin
# ==============================================================================

print("\n[3] Calculating occupancy and per-bin fold-change...")

results = []
all_bins = []  # store per-bin fold-change info

# Process each treatment sample
for tt_abd in tt_abds:
    print(f"\n    Processing: {tt_abd}")
    ctr_abd_id = get_ctr_id(tt_abd)
    print(f"      Control: {ctr_abd_id}")

    # Check if sample exists in data
    if tt_abd not in counts.columns:
        print(f"      WARNING: Sample '{tt_abd}' not found in data. Skipping...")
        continue

    # --- Adaptive Pseudo-counts ---
    # Treatment
    non_zero_tt = counts.loc[counts[tt_abd] > 0, tt_abd]
    p_tt = non_zero_tt.median() * 0.1 if not non_zero_tt.empty else 0.001

    # Control
    p_ctrl = 0.001
    if ctr_abd_id and ctr_abd_id in counts.columns:
        non_zero_ctrl = counts.loc[counts[ctr_abd_id] > 0, ctr_abd_id]
        if not non_zero_ctrl.empty:
            p_ctrl = non_zero_ctrl.median() * 0.1

    print(f"      Pseudo-counts used -> Treatment: {p_tt:.6f}, Control: {p_ctrl:.6f}")

    # Process each chromosome
    for chrom in sorted(valid_chroms):
        chrom_data = counts[counts["chromosome"] == chrom].copy()
        if len(chrom_data) == 0:
            continue

        # Calculate fold-change per bin if needed
        if calculate_foldchange and use_control and ctr_abd_id in counts.columns:
            chrom_data["fold_change"] = (chrom_data[tt_abd] + p_tt) / (
                chrom_data[ctr_abd_id] + p_ctrl
            )
            chrom_data["covered_bins"] = chrom_data["fold_change"] >= threshold
        else:
            chrom_data["fold_change"] = np.nan
            chrom_data["covered_bins"] = chrom_data[tt_abd] >= threshold

        # Calculate bin size
        chrom_data["binSize"] = chrom_data["end"] - chrom_data["start"]

        # Calculate coverage metrics
        n_covered_bins = chrom_data["covered_bins"].sum()
        covered_bases = chrom_data.loc[chrom_data["covered_bins"], "binSize"].sum()

        # Get chromosome size
        chrom_size_vals = chrom_sizes[chrom_sizes["chromosome"] == chrom]["size"].values
        if len(chrom_size_vals) == 0:
            continue
        chrom_size = chrom_size_vals[0]

        # Percentage occupancy
        pct_occupancy = (covered_bases / chrom_size) * 100

        # Store chromosome-level results
        results.append(
            {
                "chromosome": chrom,
                "condition": tt_abd,
                "pct_occupancy": pct_occupancy,
                "covered_bases": covered_bases,
                "total_bases": chrom_size,
                "n_bins_covered": n_covered_bins,
                "n_bins_total": len(chrom_data),
            }
        )

        # Append per-bin data for saving
        chrom_data["condition"] = tt_abd
        chrom_data["chromosome_size"] = chrom_size
        all_bins.append(chrom_data)

# ==============================================================================
# Create DataFrames and extract information
# ==============================================================================

print("\n[4] Creating results dataframe...")
results_df = pd.DataFrame(results)
if len(results_df) == 0:
    print("ERROR: No results generated!")
    sys.exit(1)

# Extract mark/state
try:
    marks = set()
    states = set()
    for sample in tt_abds:
        parts = sample.split("-")
        if len(parts) >= 2:
            states.add(parts[0])
            marks.add(parts[1])
    mark_pattern = "(" + "|".join(map(re.escape, marks)) + ")"
    state_pattern = "(" + "|".join(map(re.escape, states)) + ")"

    results_df["mark"] = results_df["condition"].str.extract(mark_pattern)
    results_df["state"] = results_df["condition"].str.extract(state_pattern)
except Exception as e:
    print(f"    Warning: Could not extract mark/state info: {e}")
    results_df["mark"] = "unknown"
    results_df["state"] = "unknown"

# ==============================================================================
# Summary statistics & Save
# ==============================================================================

summary = results_df.groupby("condition")["pct_occupancy"].agg(
    ["mean", "median", "min", "max", "std"]
)
bin_coverage = results_df.groupby("condition")[["n_bins_covered", "n_bins_total"]].sum()
bin_coverage["pct_bins_covered"] = (
    bin_coverage["n_bins_covered"] / bin_coverage["n_bins_total"]
) * 100


results_df.to_csv(chromosome_occupancy, index=False)


summary.to_csv(snakemake.output.occupancy_summary_stats)
bin_coverage.to_csv(snakemake.output.bin_coverage_stats)
if all_bins:
    per_bin_df = pd.concat(all_bins)
    per_bin_df.to_csv(snakemake.output.per_bin_fc, index=False)

if log:
    with open(log, "w") as f:
        f.write("Chromosome occupancy calculation completed successfully\n")
        f.write(summary.to_string())

print("\nDONE!")
