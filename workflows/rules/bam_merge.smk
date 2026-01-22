from constants.dirs_files import *

rule merge_sort_bam:
    input:
        all_rep1_deduped=lambda wc: f"{mark_remove_dups}/{get_reps(wc.gen_abd_id)[0]}.deduped.bam",
        all_rep2_deduped=lambda wc: f"{mark_remove_dups}/{get_reps(wc.gen_abd_id)[1]}.deduped.bam",
    output:
        merged_unsorted=temp(merged_bam_dir / "treatment/unsorted_{gen_abd_id}.merged.bam"),
        merged_sorted=merged_bam_dir / "treatment/{gen_abd_id}.merged.bam",
    log:
        std_out="logs/samtools/treatment/{gen_abd_id}.merged.out.log",
        std_err="logs/samtools/treatment/{gen_abd_id}.merged.err.log",
    threads: config["resources"]["bam"]["cpu"]
    shell:
        """
        samtools merge -@ {threads} {output.merged_unsorted} {input.all_rep1_deduped} {input.all_rep2_deduped} \
            > {log.std_out} 2> {log.std_err} \
        && samtools sort -@ {threads} {output.merged_unsorted} -o {output.merged_sorted} \
            >> {log.std_out} 2>> {log.std_err}
        """



rule index_merged:
    input:
        sorted_merged_bam=merged_bam_dir / "treatment/{gen_abd_id}.merged.bam",
    output:
        sorted_merged_bam_bai=merged_bam_dir / "treatment/{gen_abd_id}.merged.bam.bai",
    log:
        std_out="logs/samtools/{gen_abd_id}.idx_merged.out.log",
        std_err="logs/samtools/{gen_abd_id}.idx_merged.err.log",
    conda:
        "../../envs/samtools.yaml"
    threads: config["resources"]["bam"]["cpu"]
    shell:
        """
        samtools index -@ {threads} {input.sorted_merged_bam} \
        > {log.std_out} 2> {log.std_err}
        """


rule merge_ctr_bam:
    input:
        ctr_bam=lambda wc: get_ctrl_deduped_bam(wc.ctr_abd_id),
        ctr_bai=lambda wc: get_ctrl_deduped_bai(wc.ctr_abd_id),
    output:
        merged_ctrl_bam = merged_bam_dir / "control/{ctr_abd_id}.merged.bam",
        merged_ctrl_bai = merged_bam_dir / "control/{ctr_abd_id}.merged.bam.bai"
    log:
        std_out="logs/samtools/{ctr_abd_id}.merge_ctr.out.log",
        std_err="logs/samtools/{ctr_abd_id}.merge_ctr.err.log",
    shell:
        """
        cp {input.ctr_bam} {output.merged_ctrl_bam} 2> {log.std_err}
        cp {input.ctr_bai} {output.merged_ctrl_bai} 2>> {log.std_err}
        """
