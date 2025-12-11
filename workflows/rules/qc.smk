import sys
from pathlib import Path
from constants.common import *

rule fastqc:
    input:
        expand(
            str(RAW_DATA_DIR) + "/{sample_id}_{meta_id}",
            sample_id=sample_ids,
            meta_id=meta_ids
        )
    output:
        fqs = expand(
            str(fastqc_dir) + "/{sample_id}_{root_id}_fastqc.{ext}",
            sample_id=sample_ids,
            root_id=root_meta_ids,
            ext=["html", "zip"]
        ),
        fqc_out = directory(fastqc_dir)
    message:
        "Running FASTQC on {input}..."
    log:
        out = "logs/{sample_id}_{root_id}_fastqc_raw_out.log",
        err = "logs/{sample_id}_{root_id}_fastqc_raw_err.log"
    shell:
        """
        module load fastqc/0.11.9
        mkdir -p {output.fqc_out}
        mkdir -p logs

        fastqc -o {output.fqc_out} -f fastq {input} \
        1> {log.out} 2> {log.err}
        """
