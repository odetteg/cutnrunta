import os
from pathlib import Path

from constants.dirs_files import *


rule get_default_effective_genome_size:
    input:
        ref_genome=BT9_TA_REF_FA,
    output:
        egs_json=BASE_DIR / "ad_files/egs.json",
    threads: config["resources"]["unique-kmers"]["cpu"]
    params:
        mem_mb=config["resources"]["unique-kmers"]["mem_mb"],
    log:
        std_err="logs/khmer/unique_kmers.err.log",
    conda:
        "../../envs/deeptools.yaml"
    script:
        "../scripts/effective_genome_size.py"


rule egs:
    input:
        bam=mark_remove_dups / "{base_id}.deduped.bam",
        approx_egs=rules.get_default_effective_genome_size.output.egs_json,
    output:
        sample_len_parts=ad_files / "{base_id}.len.csv",
    params:
        sample_id="{base_id}",
    threads: config["resources"]["samtools"]["cpu"]
    script:
        "../scripts/read_len.py"


rule merge_sample_lens:
    input:
        sample_len_parts=expand(ad_files / "{base_id}.len.csv", base_id=base_ids),
        egs_json=rules.get_default_effective_genome_size.output.egs_json,
    output:
        sample_egs_csv=ad_files / "sample_egs.csv",
    run:
        import pandas as pd
        import os
        import json

        os.makedirs(os.path.dirname(output.sample_egs_csv), exist_ok=True)

        # Load the effective genome sizes JSON
        with open(input.egs_json, "r") as f:
            egs = json.load(f)

            # Merge all individual sample length files
        dfs = [pd.read_csv(p) for p in input.sample_len_parts]
        merged = pd.concat(dfs, ignore_index=True)


        # Function to get effective genome size based on fragment length
        def get_effective_genome_size(mode_len):
            available_egs = list(egs.keys())
            mode_len_str = str(int(mode_len))  # Ensure it's a string of int

            # Find exact match or closest read length
            if mode_len_str not in available_egs:
                closest_len = min(
                    available_egs, key=lambda x: abs(int(x) - int(mode_len))
                )
                return int(egs[closest_len])
            return int(egs[mode_len_str])


        merged["effective_genome_size"] = merged["mode_len"].apply(
            get_effective_genome_size
        )

        merged.to_csv(output.sample_egs_csv, index=False)
