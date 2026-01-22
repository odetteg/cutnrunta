import pandas as pd
import subprocess
import os
from pathlib import Path
import sys


sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))
from constants.dirs_files import *


bam = snakemake.input["bam"]
base_id = snakemake.wildcards.base_id


output_path = snakemake.output.sample_len_parts
os.makedirs(os.path.dirname(output_path), exist_ok=True)


cmd_ = (
    "samtools view -F 4 {bam} | head -n 1000000 | "
    "cut -f 10 | awk '{{print length}}' | "
    "sort | uniq -c | sort -nr | head -n 1"
)


def parse_mode(result):
    """Safely extract the modal read length from subprocess output."""
    out = result.stdout.strip()
    if not out:
        return 0
    return int(out.split()[-1])


res = subprocess.run(
    cmd_.format(bam=bam),
    shell=True,
    capture_output=True,
    text=True,
    check=True,
)
mode_len = parse_mode(res)

lens_df = pd.DataFrame([{"sample": base_id, "mode_len": mode_len}])
lens_df.to_csv(output_path, index=False)
