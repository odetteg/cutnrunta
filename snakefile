import os
import sys
from pathlib import Path
from constants.common import *


configfile: "config/config.yaml"


include: "workflows/rules/qc.smk"
include: "workflows/rules/idx.smk"
include: "workflows/rules/map.smk"
include: "workflows/rules/multiqc.smk"


rule all:
    input:
        raw_fastqc_rpts=expand(
            str(RAW_FASTQC_DIR) + "/{sample_id}_fastqc.{ext}",
            ext=["html", "zip"],
            sample_id=sample_ids,
        ),
        fastp_trimmed_fqc_rpt=expand(
            config["FASTQC_QC_FASTP_DIR"]
            + "/{base_id}.fastp.trimmed.{read}_fastqc.{ext}",
            base_id=base_ids,
            read=["R1", "R2"],
            ext=["html", "zip"],
        ),
        fastp_trimms=expand(
            config["FASTP_TRIMMED_DIR"] + "/{base_id}.fastp.trimmed.{read}.fastq.gz",
            read=["R1", "R2"],
            base_id=base_ids,
        ),
        fastp_rpts=expand(
            config["FASTP_QC_REPORTS_DIR"] + "/{base_id}.fastp.{ext}",
            ext=["html", "json"],
            base_id=base_ids,
        ),
        bt2_build_index=multiext(
            config["BT9_TA_REF_FA"],
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
        fai=config["BT9_TA_REF_FA"] + ".fai",
        raw_bt2_mapped=expand(
            config["dir_names"]["bt2_raw_map_dir"]
            + "/{base_id}.raw.unsorted.unfiltered.bam",
            base_id=base_ids,
        ),
        sorted_filtered_bam=expand(
            config["dir_names"]["samtools_sorted_filtered_dir"]
            + "/{base_id}.sorted.filtered.bam",
            base_id=base_ids,
        ),
        multiqc_rpt=RESULTS_DIR / "qc/multiqc/multiqc.html",
