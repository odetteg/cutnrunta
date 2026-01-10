from snakemake.shell import shell
from constants.dirs_files import *

sample_id = snakemake.wildcards.sample_id

tt_bam = snakemake.input.tt_bam
ctrl_bam = snakemake.input.get("ctrl_bam", None)
mode = snakemake.params.mode
minreads = snakemake.params.minreads
pvalue = snakemake.params.pvalue
extra = snakemake.params.extra
slide = snakemake.params.slide
minwidth = snakemake.params.minwidth
step = snakemake.params.step
log = snakemake.log.std_out
mdist = snakemake.params.mdistance
chromsizes = snakemake.input.chromsizes

if mode == "broad":
    mdist = mdist
    mode_param = f"--broad --mdist {mdist}"
    peaks = f"{gopeaks_dir}/peaks/broad/{sample_id}"
else:
    mode_param = f"--step {step} --slide {slide}"
    peaks = f"{gopeaks_dir}/peaks/narrow/{sample_id}"

cmd_ = f"""
gopeaks -b {tt_bam} \
    {'-c ' + ctrl_bam if ctrl_bam else ''} \
    -o {peaks} \
    -s {chromsizes} \
    {mode_param} \
    --minreads {minreads} \
    --pval {pvalue} \
    --minwidth {minwidth} \
    {extra} \
    2> {log}
"""

shell(cmd_)
