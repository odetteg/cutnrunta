from snakemake.shell import shell
from pathlib import Path
import csv
import pandas as pd
import shutil

# Inputs
tt_bam = snakemake.input.tt_bam
ctrl_bam = snakemake.input.ctrl_bam
egs = snakemake.input.egs

# Outputs
peak = snakemake.output.peak

# Params
mode_param = snakemake.params.mode
use_control = snakemake.params.use_control
qvalue = snakemake.params.qvalue
extra = snakemake.params.extra
outdir = snakemake.params.outdir
threads = snakemake.threads
std_out = snakemake.log.std_out
sample_id = snakemake.wildcards.sample_id

# Handle broad mode
if mode_param == "broad":
    broad_cutoff = snakemake.params.broad_cutoff
    mode_flag = f"--broad --broad-cutoff {broad_cutoff}"
else:
    mode_flag = ""

# Read the sample metadata
df = pd.read_csv(egs)

# Get effective genome size
effective_genome_size = int(
    df.loc[df["sample"] == sample_id, "effective_genome_size"].values[0]
)

# screen if use there is a need to use control. MACS3 can yield no signals for cutnrun data with low background signal
if use_control:
    control_flag = f"--control {ctrl_bam}"
else:
    control_flag = ""

cmd = f"""
macs3 callpeak \
    --treatment {tt_bam} \
    {control_flag} \
    --gsize {effective_genome_size} \
    --qvalue {qvalue} \
    --format BAMPE \
    --keep-dup all \
    --outdir {outdir} \
    --name {sample_id} \
    {mode_flag} \
    {extra} \
    > {std_out} 2>&1
"""

shell(cmd)


mode_dir = Path(outdir).parent
aux_dir = mode_dir / "aux_files"
aux_dir.mkdir(parents=True, exist_ok=True)


xls_source = Path(outdir) / f"{sample_id}_peaks.xls"
xls_dest = aux_dir / f"{sample_id}_peaks.xls"
if xls_source.exists():
    shutil.move(str(xls_source), str(xls_dest))


if mode_param == "broad":
    gapped_source = Path(outdir) / f"{sample_id}_peaks.gappedPeak"
    gapped_dest = aux_dir / f"{sample_id}_peaks.gappedPeak"
    if gapped_source.exists():
        shutil.move(str(gapped_source), str(gapped_dest))
else:
    bed_source = Path(outdir) / f"{sample_id}_summits.bed"
    bed_dest = aux_dir / f"{sample_id}_summits.bed"
    if bed_source.exists():
        shutil.move(str(bed_source), str(bed_dest))
