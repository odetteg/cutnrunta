from constants.dirs_files import *


rule fc_per_replicate:
    input:
        raw_counts=rules.multibigwigsummary_per_rep.output.raw_counts,
        chrom_sizes=rules.get_cs.output.cs,
    output:
        chromosome_occupancy=analysis_dir
        / "fold_change/chromosome_occupancy_per_rep.csv",
        occupancy_summary_stats=analysis_dir
        / "fold_change/occupancy_summary_stats_per_rep.csv",
        per_bin_fc=analysis_dir / "fold_change/per_bin_foldchange_per_rep.csv",
        bin_coverage_stats=analysis_dir / "fold_change/bin_coverage_stats_per_rep.csv",
    params:
        binSize=config["fold_change"]["binSize"],
        tt_ids=tt_base_ids,
        ctr_ids=ctrl_abd_base_ids,
    log:
        "logs/deeptools/fold_change/fc_per_rep.log",
    conda:
        "../../envs/deeptools.yaml"
    script:
        "../scripts/fc.py"


rule fc_avg_bw:
    input:
        raw_counts=rules.multibigwigsummary_avg_bw.output.raw_counts,
        chrom_sizes=rules.get_cs.output.cs,
    output:
        chromosome_occupancy=analysis_dir / "fold_change/chromosome_occupancy_avg.csv",
        occupancy_summary_stats=analysis_dir
        / "fold_change/occupancy_summary_stats_avg.csv",
        per_bin_fc=analysis_dir / "fold_change/per_bin_foldchange_avg.csv",
        bin_coverage_stats=analysis_dir / "fold_change/bin_coverage_stats_avg.csv",
    params:
        binSize=config["fold_change"]["binSize"],
        tt_ids=tt_abd_ids,
        ctr_ids=ctr_abd_ids,
    log:
        "logs/deeptools/fold_change/fc_avg.log",
    conda:
        "../../envs/deeptools.yaml"
    script:
        "../scripts/fc.py"
