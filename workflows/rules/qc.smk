import os
from pathlib import Path
from constants.common import *
from constants.dirs_files import *

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
        fq=lambda wc: RAW_DATA_DIR / f"{wc.sample_id}.fastq.gz",
    output:
        html=RAW_FASTQC_DIR / "{sample_id}_fastqc.html",
        zip=RAW_FASTQC_DIR / "{sample_id}_fastqc.zip",
    log:
        out=LOGS_DIR / "fastqc/{sample_id}_fastqc_raw_out.log",
        err=LOGS_DIR / "fastqc/{sample_id}_fastqc_raw_err.log",
    params:
        out_dir=RAW_FASTQC_DIR,
    resources:
        cpus=4,
    shell:
        """
        module load fastqc/0.11.9
        fastqc -o {params.out_dir} -t {resources.cpus} -f fastq {input.fq} \
        > {log.out} 2> {log.err}
        """


rule fastp_trim:
    input:
        read1=lambda wc: get_paired_rds(wc.base_id)[0],
        read2=lambda wc: get_paired_rds(wc.base_id)[1],
    output:
        trim_r1=FASTP_TRIMMED_DIR / "{base_id}.fastp.trimmed.R1.fastq.gz",
        trim_r2=FASTP_TRIMMED_DIR / "{base_id}.fastp.trimmed.R2.fastq.gz",
        html=FASTP_QC_REPORTS_DIR / "{base_id}.fastp.html",
        json=FASTP_QC_REPORTS_DIR / "{base_id}.fastp.json",
    params:
        extras=config["fastp_trim"]["extra"],
        a=config["fastp_trim"]["a"],
        r2_adapter=config["fastp_trim"]["r2_adapter"],
    resources:
        cpu=config["resources"]["fastp_trim"]["cpu"],
    log:
        std_out="logs/fastp/{base_id}.trimmed.fastp.out.log",
        std_err="logs/fastp/{base_id}.trimmed.fastp.err.log",
    shell:
        """
        module load fastp/0.23.1
        fastp \
            -i {input.read1} -I {input.read2} \
            -o {output.trim_r1} -O {output.trim_r2} \
            -h {output.html} -j {output.json} \
            --thread {resources.cpu} \
            {params.extras} \
            > {log.std_out} 2> {log.std_err}
        """


rule post_trim_fqc:
    input:
        fastp_trimmed=lambda wc: os.path.join(
            FASTP_TRIMMED_DIR,
            f"{wc.base_id}.fastp.trimmed.{wc.read}.fastq.gz",
        ),
    output:
        html=FASTQC_QC_FASTP_DIR / "{base_id}.fastp.trimmed.{read}_fastqc.html",
        zip=FASTQC_QC_FASTP_DIR / "{base_id}.fastp.trimmed.{read}_fastqc.zip",
    log:
        std_out=str(LOGS_DIR) + "/{base_id}.fastqc.fastp.trimmed.{read}.out.log",
        std_err=str(LOGS_DIR) + "/{base_id}.fastqc.fastp.trimmed.{read}.err.log",
    params:
        out_dir=FASTQC_QC_FASTP_DIR,
    resources:
        cpus=4,
    shell:
        """
        module load fastqc/0.11.9
        mkdir -p {params.out_dir}
        fastqc -o {params.out_dir} -t {resources.cpus} -f fastq {input.fastp_trimmed} \
        1> {log.std_out} 2> {log.std_err}
        """
