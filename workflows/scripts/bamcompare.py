from snakemake.shell import shell
import pandas as pd


base_id = snakemake.wildcards.base_id


# Input
tt_bam = snakemake.input.tt_bam
tt_bai = snakemake.input.tt_bai
ctrl_bam = snakemake.input.ctrl_bam
ctrl_bai = snakemake.input.ctrl_bai
egs = snakemake.input.egs

# Output
bam_compare_bw = snakemake.output.bam_compare_bw

# Threads
threads = snakemake.threads

# Log
std_out = snakemake.log.std_out
std_err = snakemake.log.std_err

# Params
normalization = snakemake.params.normalization
scale_factor_method = snakemake.params.scale_factor_method
operation = snakemake.params.operation
binsize = snakemake.params.binsize
extra = snakemake.params.extra

# Read effective genome size for this sample. The egs is a csv file with macthed base len and effective genome size

egs_df = pd.read_csv(egs)

mask = egs_df["sample"].str.contains(base_id)

if mask.any():
    effective_genome_size = int(egs_df.loc[mask, "effective_genome_size"].values[0])
else:
    # Keeping this here for now. Migh resolve later.
    with open(std_err, "a") as f:
        f.write(
            f"Warning: Could not correctly map {base_id} to its effective genome size. Using default: 2432564026\n"
        )
    effective_genome_size = 2432564026

cmd = f"""
bamCompare -b1 {tt_bam} \
    -b2 {ctrl_bam} \
    -o {bam_compare_bw} \
    --numberOfProcessors {threads} \
    --effectiveGenomeSize {effective_genome_size} \
    --normalizeUsing {normalization} \
    --scaleFactorsMethod {scale_factor_method} \
    --operation {operation} \
    --binSize {binsize} \
    {extra} \
    1> {std_out} 2>>{std_err}
"""

shell(cmd)
