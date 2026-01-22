from constants.dirs_files import *


rule ac_lost:
    input:
        bl3_ac=dynamic_greedy_consensus_peaks
        / "BL3-H3K27ac_consensus.greedy.stringent.bed",
        tbl3_ac=dynamic_greedy_consensus_peaks
        / "TBL3-H3K27ac_consensus.greedy.stringent.bed",
    output:
        ac_lost=analysis_dir / "h3k27ac" / "ac_lost_in_tbl3.bed",
    conda:
        "../../envs/bedtools.yaml"
    log:
        "logs/bedtools/ac_lost.log",
    shell:
        """
        bedtools intersect -v -a {input.bl3_ac} -b {input.tbl3_ac} > {output.ac_lost} 2> {log}
        """


rule me_gained:
    input:
        bl3_me=dynamic_greedy_consensus_peaks
        / "BL3-H3K27me3_consensus.greedy.stringent.bed",
        tbl3_me=dynamic_greedy_consensus_peaks
        / "TBL3-H3K27me3_consensus.greedy.stringent.bed",
    output:
        me_gained=analysis_dir / "h3k27me3" / "me_gained_in_tbl3.bed",
    conda:
        "../../envs/bedtools.yaml"
    log:
        "logs/bedtools/me_gained.log",
    shell:
        """
        bedtools intersect -v -a {input.tbl3_me} -b {input.bl3_me} > {output.me_gained} 2> {log}
        """


# rule ac_to_me:
#     input:
#         ac_lost=analysis_dir / "h3k27ac" / "ac_lost_in_tbl3.bed",
#         me_gained=analysis_dir / "h3k27me3" / "me_gained_in_tbl3.bed",
#     output:
#         ac_to_me=analysis_dir / "h3k27_transition" / "ac_to_me_in_tbl3.bed",
#     conda:
#         "../../envs/bedtools.yaml"
#     log:
#         "logs/bedtools/ac_to_me.log",
#     shell:
#         """
#         bedtools intersect -a {input.ac_lost} -b {input.me_gained} > {output.ac_to_me} 2> {log}
#         """

rule bl3_me_lost_in_tbl3:
    input:
        bl3_me=dynamic_greedy_consensus_peaks
        / "BL3-H3K27me3_consensus.greedy.stringent.bed",
        tbl3_me=dynamic_greedy_consensus_peaks
        / "TBL3-H3K27me3_consensus.greedy.stringent.bed",
    output:
        me_lost=analysis_dir / "h3k27me3" / "bl3_me_lost_in_tbl3.bed",
    conda:
        "../../envs/bedtools.yaml"
    log:
        "logs/bedtools/me_lost.log",
    shell:
        """
        bedtools intersect -v -a {input.bl3_me} -b {input.tbl3_me} > {output.me_lost} 2> {log}
        """