from snakemake.shell import shell
import pandas as pd
import json

bam = snakemake.input.bam
sample_lens_csv_path = snakemake.input.sample_egs_csv
bw = snakemake.output.bw_
egs_json = snakemake.input.egs_json
binsize = snakemake.params.binsize
normalization = snakemake.params.normalization
extras = snakemake.params.extras
threads = snakemake.threads
base_id = snakemake.wildcards.base_id
config = snakemake.config

with open(egs_json, "r") as f:
    egs = json.load(f)

sample_lens_df = pd.read_csv(sample_lens_csv_path)

if base_id not in sample_lens_df["sample"].values:
    raise ValueError(f"Sample {base_id} not found in {sample_lens_csv_path}")

mode_len = int(
    sample_lens_df.loc[sample_lens_df["sample"] == base_id, "mode_len"].values[0]
)
available_egs = list(egs.keys())
if str(mode_len) not in available_egs:
    mode_len = int(min(available_egs, key=lambda x: abs(int(x) - mode_len)))

effective_genome_size = int(egs[str(mode_len)])

if config.get("peaks_call", {}).get("lanceotron", {}).get("use_lanceotron", True):
    effectiveGenomeSize_params = f"--effectiveGenomeSize {effective_genome_size} "

    binsize = config.get("peaks_call", {}).get("lanceotron", {}).get("binsize", binsize)
    normalization = (
        config.get("peaks_call", {})
        .get("lanceotron", {})
        .get("normalization", normalization)
    )
else:
    effectiveGenomeSize_params = ""

cmd = (
    f"bamCoverage --numberOfProcessors {threads} "
    f"{effectiveGenomeSize_params}"
    f"--bam {bam} --binSize {binsize} "
    f"--normalizeUsing {normalization} --outFileName {bw} {extras}"
)

shell(cmd)
