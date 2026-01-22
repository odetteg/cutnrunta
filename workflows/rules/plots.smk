import os
from pathlib import Path
from constants.dirs_files import *


rule create_stats_list:
    input:
        stats=expand(
            sorted_filtered_stats_dir  / "{base_id}.stats.txt",
            base_id=base_ids
        ),
    output:
        samtools_txt=samtools_txt,
    params:
        stats_dir = sorted_filtered_stats_dir 
    shell:
        """
        ls {params.stats_dir}/*.stats.txt > {output.samtools_txt}
        """


rule agg_stats:
    input:
        agg_stats_txt=rules.create_stats_list.output.samtools_txt,
    output:
        agg_stats_csv=agg_stats_csv,
    conda: "../../envs/stats.yaml"
    shell:
        """
        which python
        python workflows/scripts/data_pool.py {input.agg_stats_txt} {output.agg_stats_csv}
        """
rule frag_lens:
    input:
        agg_stats_csv = rules.agg_stats.output.agg_stats_csv
    output:
        frag_len_dst_pdf=frag_len_dst_pdf,
        frag_len_dst_facet_pdf=frag_len_dst_facet_pdf
    conda: "../../envs/r_libraries.yaml"
    script:
        "../scripts/frag_size.R"
