import os
import sys
from pathlib import Path
from constants.common import *
configfile: 'config/config.yaml'
include: "workflows/rules/qc.smk"

# ----------------------------
# All global variables
# ----------------------------



# --------------------------
# Le grand
# --------------------------
rule all:
    input:
        fastqc_rpts = expand(
            str(FASTQC_DIR) + "/{sample_id}_fastqc.{ext}",
            ext=["html", "zip"],
            sample_id=sample_ids,
        ),
        fastp_trimmed_fqc_rpt = expand(
                    config["FASQC_QC_FASTP_DIR"] + "/{base_id}.fastp.trimmed.{read}.fastqc.{ext}",
                    base_id=base_ids, read=["R1", "R2"], ext=["html", "zip"]
        ),
        fastp_trimms = expand(
            config["FASTP_TRIMMED_DIR"] + "/{base_id}.fastp.trimmed.{read}.fastq.gz",
            read=["R1", "R2"],
            base_id=base_ids
        ),

        fastp_rpts = expand(
            config["FASTP_QC_REPORTS_DIR"] + "/{base_id}.fastp.{ext}",
            ext=["html", "json"],
            base_id=base_ids,
        ),
        multiqc_rpt = RESULTS_DIR / "qc/multiqc/multiqc.html"
