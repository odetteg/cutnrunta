import os
from pathlib import Path
from constants.common import *

# --------------------------
# Rule for every other common thing
# --------------------------
# rule comm_:
#     shell:
#         """
#         mkdir -p {output.results_dir}
#         """

# --------------------------
# FASTQC rule
# --------------------------
rule fastqc:
    input:
        fq = lambda wc: RAW_DATA_DIR / f"{wc.sample_id}.fastq.gz"
    output:
        html = FASTQC_DIR / "{sample_id}_fastqc.html",
        zip  = FASTQC_DIR / "{sample_id}_fastqc.zip"
    log:
        out = LOGS_DIR / "{sample_id}_fastqc_raw_out.log",
        err = LOGS_DIR / "{sample_id}_fastqc_raw_err.log"
    resources:
        cpus=4
    shell:
        """
        module load fastqc/0.11.9
        mkdir -p {FASTQC_DIR}
        mkdir -p {LOGS_DIR}
        fastqc -o {FASTQC_DIR} -t {resources.cpus} -f fastq {input.fq} \
        1> {log.out} 2> {log.err}
        """
