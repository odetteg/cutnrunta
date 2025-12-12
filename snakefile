import os
import sys
from pathlib import Path
from constants.common import *
configfile: 'config/config.yaml'
include: "workflows/rules/qc.smk"

# ----------------------------
# All global variables
# ----------------------------
EXTS=["fastq.gz"]


# --------------------------
# Le grand
# --------------------------
rule all:
    input:
        expand(
            str(FASTQC_DIR) + "/{sample_id}_fastqc.html",
            sample_id=sample_ids,
        ),
        expand(
            str(FASTQC_DIR) + "/{sample_id}_fastqc.zip",
            sample_id=sample_ids,
        )