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
# rule fastp_trim:
#     input:


rule multiqc:
    input:
        expand(str(FASTQC_DIR) + "/{sample_id}_fastqc.zip", sample_id=sample_ids)
    output:
        multqc_rpt = RESULTS_DIR / "qc/multiqc/multiqc.html"
    params:
        mqc_dir = RESULTS_DIR / "qc/multiqc",
    threads: config["resources"]["multiqc"]["cpu"]
    log:
        out = "logs/multiqc/multiqc.out.log",
        err = "logs/multiqc/multiqc.err.log"
    shell:
        """
        module load -f multiqc/1.29
        mkdir -p {params.mqc_dir}
        mkdir -p logs/multiqc
        multiqc --force --outdir {params.mqc_dir} -n multiqc.html {input} \
        1>{log.out} 2>{log.err}
        """
