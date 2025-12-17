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

rule fastp_trim:
    input:
        read1 = lambda wc: get_paired_rds(wc.base_id)[0],
        read2 = lambda wc: get_paired_rds(wc.base_id)[1]
    output:
        trim_r1 = config["FASTP_TRIMMED_DIR"] + "/{base_id}.fastp.trimmed.R1.fastq.gz",
        trim_r2 = config["FASTP_TRIMMED_DIR"] + "/{base_id}.fastp.trimmed.R2.fastq.gz",
        html = config["FASTP_QC_REPORTS_DIR"] + "/{base_id}.fastp.html",
        json = config["FASTP_QC_REPORTS_DIR"] + "/{base_id}.fastp.json"
    params:
        fastp_qc_rpts = config["FASTP_QC_REPORTS_DIR"],
        fastp_trimmed_dir = config["FASTP_TRIMMED_DIR"],
        extras = config["fastp_trim"]["extra"]
    resources:
        cpu = config["resources"]["fastp_trim"]["cpu"]
    log:
        std_out = "logs/fastp/{base_id}.trimmed.fastp.out.log",
        std_err = "logs/fastp/{base_id}.trimmed.fastp.err.log"
    shell:
        """
        mkdir -p {params.fastp_qc_rpts}
        mkdir -p {params.fastp_trimmed_dir}
        module load fastp/0.23.1
        fastp \
            -i {input.read1} -I {input.read2} \
            -o {output.trim_r1} -O {output.trim_r2} \
            -h {output.html} -j {output.json} \
            {params.extras} -w {resources.cpu} \
            1> {log.std_out} 2> {log.std_err}
        """
        
rule post_trim_fqc:
    input:
        fastp_trimmed = lambda wc: os.path.join(
            config["FASTP_TRIMMED_DIR"], f"{wc.base_id}.fastp.trimmed.{wc.read}.fastq.gz"
        )
    output:
        html = config["FASQC_QC_FASTP_DIR"] + "/{base_id}.fastp.trimmed.{read}.fastqc.html",
        zip  = config["FASQC_QC_FASTP_DIR"] + "/{base_id}.fastp.trimmed.{read}.fastqc.zip"
    log:
        std_out = str(LOGS_DIR) + "/{base_id}.fastqc.fastp.trimmed.{read}.out.log",
        std_err =str(LOGS_DIR) + "/{base_id}.fastqc.fastp.trimmed.{read}.err.log"
    params:
        out_dir = config["FASQC_QC_FASTP_DIR"]
    resources:
        cpus = 4
    shell:
        """
        module load fastqc/0.11.9
        mkdir -p {params.out_dir}
        fastqc -o {params.out_dir} -t {resources.cpus} -f fastq {input.fastp_trimmed} \
        1> {log.std_out} 2> {log.std_err}
        """

rule multiqc:
    input:
        fastqc_raw_fqs = expand(
            FASTQC_DIR / "{sample_id}_fastqc.{ext}",
            sample_id=sample_ids,
            ext=["html", "zip"]
        ),
        fastp_rpts = expand(
            config["FASTP_QC_REPORTS_DIR"] + "/{base_id}.fastp.{ext}",
            ext=["html", "json"],
            base_id=base_ids),
        fqc_fastp_trimm_rpts = expand(config["FASQC_QC_FASTP_DIR"] + "/{base_id}.fastp.trimmed.{read}.fastqc.{ext}",
        base_id=base_ids,
        read=["R1", "R2"],
        ext=["html", "zip"]
        ),
    output:
        multqc_rpt = RESULTS_DIR / "qc/multiqc/multiqc.html"
    params:
        mqc_dir = RESULTS_DIR / "qc/multiqc"
    threads:
        config["resources"]["multiqc"]["cpu"]
    log:
        std_out = "logs/multiqc/multiqc.out.log",
        std_err = "logs/multiqc/multiqc.err.log"
    shell:
        """
        module load -f multiqc/1.29
        mkdir -p {params.mqc_dir} logs/multiqc
        multiqc --force \
                --outdir {params.mqc_dir} \
                -n multiqc.html \
                {input} \
                1> {log.std_out} 2> {log.std_err}
        """
