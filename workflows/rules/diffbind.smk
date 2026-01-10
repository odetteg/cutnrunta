from constants.dirs_files import *


rule diffBindSeacr:
    input:
        relaxed_sample_sheet=ad_files / "relaxed_sample_sheet.csv",
        stringent_sample_sheet=ad_files / "stringent_sample_sheet.csv",
        stringent_peaks=expand(
            dynamic_unfiltered_seacr_dir / "relaxed/{tt_ids}.relaxed.bed",
            tt_ids=tt_ids,
        ),
        relaxed_peaks=expand(
            dynamic_unfiltered_seacr_dir / "stringent/{tt_ids}.stringent.bed",
            tt_ids=tt_ids,
        ),
    output:
        relaxed_pca=analysis_dir / "diffbind" / "pca" / "relaxed_pca.pdf",
        stringent_pca=analysis_dir / "diffbind" / "pca" / "stringent_pca.pdf",
        relaxed_counts_matrix_csv=analysis_dir
        / "diffbind"
        / "relaxed_counts.matrix.csv",
        stringent_counts_matrix_csv=analysis_dir
        / "diffbind"
        / "stringent_counts.matrix.csv",
        relaxed_h3k27ac_results=analysis_dir
        / "diffbind"
        / "relaxed_H3K27ac_BL3_vs_TBL3.csv",
        stringent_h3k27ac_results=analysis_dir
        / "diffbind"
        / "stringent_H3K27ac_BL3_vs_TBL3.csv",
        relaxed_h3k27me3_results=analysis_dir
        / "diffbind"
        / "relaxed_H3K27me3_BL3_vs_TBL3.csv",
        stringent_h3k27me3_results=analysis_dir
        / "diffbind"
        / "stringent_H3K27me3_BL3_vs_TBL3.csv",
        relaxed_h2k119ub_results=analysis_dir
        / "diffbind"
        / "relaxed_H2K119Ub_BL3_vs_TBL3.csv",
        stringent_h2k119ub_results=analysis_dir
        / "diffbind"
        / "stringent_H2K119Ub_BL3_vs_TBL3.csv",
        diffbind_summary_report=analysis_dir / "diffbind" / "diffbind_summary.txt",
        relaxed_h3k27ac_pca=analysis_dir
        / "diffbind"
        / "pca"
        / "relaxed_H3K27ac_pca.pdf",
        stringent_h3k27ac_pca=analysis_dir
        / "diffbind"
        / "pca"
        / "stringent_H3K27ac_pca.pdf",
        relaxed_h3k27me3_pca=analysis_dir
        / "diffbind"
        / "pca"
        / "relaxed_H3K27me3_pca.pdf",
        stringent_h3k27me3_pca=analysis_dir
        / "diffbind"
        / "pca"
        / "stringent_H3K27me3_pca.pdf",
        relaxed_h2k119ub_pca=analysis_dir
        / "diffbind"
        / "pca"
        / "relaxed_H2K119Ub_pca.pdf",
        stringent_h2k119ub_pca=analysis_dir
        / "diffbind"
        / "pca"
        / "stringent_H2K119Ub_pca.pdf",
        relaxed_h3k27ac_diffbind_rds=analysis_dir
        / "diffbind"
        / "rds"
        / "relaxed_H3K27ac_diffbind.rds",
        stringent_h3k27ac_diffbind_rds=analysis_dir
        / "diffbind"
        / "rds"
        / "stringent_H3K27ac_diffbind.rds",
        relaxed_h3k27me3_diffbind_rds=analysis_dir
        / "diffbind"
        / "rds"
        / "relaxed_H3K27me3_diffbind.rds",
        stringent_h3k27me3_diffbind_rds=analysis_dir
        / "diffbind"
        / "rds"
        / "stringent_H3K27me3_diffbind.rds",
        relaxed_h2k119ub_diffbind_rds=analysis_dir
        / "diffbind"
        / "rds"
        / "relaxed_H2K119Ub_diffbind.rds",
        stringent_h2k119ub_diffbind_rds=analysis_dir
        / "diffbind"
        / "rds"
        / "stringent_H2K119Ub_diffbind.rds",
        relaxed_ac_dba=analysis_dir / "diffbind" / "dba" / "relaxed_H3K27ac_dba.rds",
        stringent_ac_dba=analysis_dir / "diffbind" / "dba" / "stringent_H3K27ac_dba.rds",
        relaxed_me3_dba=analysis_dir / "diffbind" / "dba" / "relaxed_H3K27me3_dba.rds",
        stringent_me3_dba=analysis_dir
        / "diffbind"
        / "dba"
        / "stringent_H3K27me3_dba.rds",
        relaxed_ub_dba=analysis_dir / "diffbind" / "dba" / "relaxed_H2K119Ub_dba.rds",
        stringent_ub_dba=analysis_dir
        / "diffbind"
        / "dba"
        / "stringent_H2K119Ub_dba.rds",
        relaxed_ac_counts_dba=analysis_dir
        / "diffbind"
        / "counts"
        / "relaxed_H3K27ac_counts_dba.rds",
        stringent_ac_counts_dba=analysis_dir
        / "diffbind"
        / "counts"
        / "stringent_H3K27ac_counts_dba.rds",
        relaxed_me3_counts_dba=analysis_dir
        / "diffbind"
        / "counts"
        / "relaxed_H3K27me3_counts_dba.rds",
        stringent_me3_counts_dba=analysis_dir
        / "diffbind"
        / "counts"
        / "stringent_H3K27me3_counts_dba.rds",
        relaxed_ub_counts_dba=analysis_dir
        / "diffbind"
        / "counts"
        / "relaxed_H2K119Ub_counts_dba.rds",
        stringent_ub_counts_dba=analysis_dir
        / "diffbind"
        / "counts"
        / "stringent_H2K119Ub_counts_dba.rds",
        # MA plots
        # relaxed_ac_ma_plot=analysis_dir
        # / "diffbind"
        # / "ma_plots"
        # / "relaxed_H3K27ac_ma_plot.pdf",
        # stringent_ac_ma_plot=analysis_dir
        # / "diffbind"
        # / "ma_plots"
        # / "stringent_H3K27ac_ma_plot.pdf",
        # relaxed_me3_ma_plot=analysis_dir
        # / "diffbind"
        # / "ma_plots"
        # / "relaxed_H3K27me3_ma_plot.pdf",
        # stringent_me3_ma_plot=analysis_dir
        # / "diffbind"
        # / "ma_plots"
        # / "stringent_H3K27me3_ma_plot.pdf",
        # relaxed_ub_ma_plot=analysis_dir
        # / "diffbind"
        # / "ma_plots"
        # / "relaxed_H2K119Ub_ma_plot.pdf",
        # stringent_ub_ma_plot=analysis_dir
        # / "diffbind"
        # / "ma_plots"
        # / "stringent_H2K119Ub_ma_plot.pdf",
    log:
        "logs/r/diffbind.log",
    conda:
        "../../envs/diffbind.yaml"
    script:
        "../scripts/diffbind.R"
