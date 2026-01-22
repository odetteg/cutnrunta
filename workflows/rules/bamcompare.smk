from constants.dirs_files import *


rule bamCompare:
    input:
        tt_bam=mark_remove_dups / "{base_id}.deduped.bam",
        tt_bai=mark_remove_dups / "{base_id}.deduped.bam.bai",
        ctrl_bam=lambda wc: get_ctrl_bam(wc.base_id),
        ctrl_bai=lambda wc: get_ctrl_bai(wc.base_id),
        egs=rules.merge_sample_lens.output.sample_egs_csv,
    output:
        bam_compare_bw=bamCompare / "per_replicate/{base_id}.bw",
    conda:
        "../../envs/deeptools.yaml"
    threads: config["resources"]["deeptools"]["cpu"]
    log:
        std_out="logs/deeptools/{base_id}.bamcompare.out.log",
        std_err="logs/deeptools/{base_id}.bamcompare.err.log",
    params:
        normalization=config["deeptools"]["bamcompare"]["normalization"],
        scale_factor_method=config["deeptools"]["scale_factor_method"],
        operation=config["deeptools"]["operation"],
        binsize=config["deeptools"]["bamcompare"]["binsize"],
        extra=config["deeptools"]["bamcompare"]["extra"],
    script:
        "../scripts/bamcompare.py"


rule bamCompareMerged:
    input:
        tt_bam=merged_bam_dir / "treatment/{gen_abd_id}.merged.bam",
        tt_bai=merged_bam_dir / "treatment/{gen_abd_id}.merged.bam.bai",
        ctrl_bam=lambda wc: get_ctrl_merged_bam(wc.gen_abd_id),
        ctrl_bai=lambda wc: get_ctrl_merged_bai(wc.gen_abd_id),
        egs=rules.merge_sample_lens.output.sample_egs_csv,
    output:
        bam_compare_bw=bamCompare / "merged/{gen_abd_id}.merged.bw",
    conda:
        "../../envs/deeptools.yaml"
    threads: config["resources"]["deeptools"]["cpu"]
    log:
        std_out="logs/deeptools/{gen_abd_id}.bamcompare.merged.out.log",
        std_err="logs/deeptools/{gen_abd_id}.bamcompare.merged.err.log",
    params:
        normalization=config["deeptools"]["bamcompare"]["normalization"],
        scale_factor_method=config["deeptools"]["scale_factor_method"],
        operation=config["deeptools"]["operation"],
        binsize=config["deeptools"]["bamcompare"]["binsize"],
        extra=config["deeptools"]["bamcompare"]["extra"],
    script:
        "../scripts/bamcompare_merged.py"
