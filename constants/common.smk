seacr_modes = ["stringent", "relaxed"]


def targets():
    """
    Returns all target files and dirs for rule all
    """
    TARGETS = [
        # Raw FastQC
        expand(
            str(RAW_FASTQC_DIR) + "/{sample_id}_fastqc.{ext}",
            sample_id=sample_ids,
            ext=["html", "zip"],
        ),
        # FastQC after fastp
        expand(
            FASTQC_QC_FASTP_DIR / "{base_id}.fastp.trimmed.{read}_fastqc.{ext}",
            base_id=base_ids,
            read=["R1", "R2"],
            ext=["html", "zip"],
        ),
        # fastp trimmed FASTQs
        expand(
            FASTP_TRIMMED_DIR / "{base_id}.fastp.trimmed.{read}.fastq.gz",
            base_id=base_ids,
            read=["R1", "R2"],
        ),
        # fastp reports
        expand(
            FASTP_QC_REPORTS_DIR / "{base_id}.fastp.{ext}",
            base_id=base_ids,
            ext=["html", "json"],
        ),
        # Bowtie2 index
        multiext(
            str(BT9_TA_REF_FA),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
        str(BT9_TA_REF_FA) + ".fai",
        # Mapping outputs
        expand(
            bt2_raw_map_dir / "{base_id}.raw.unsorted.unfiltered.bam",
            base_id=base_ids,
        ),
        expand(
            samtools_sorted_filtered_dir / "{base_id}.sorted.filtered.bam",
            base_id=base_ids,
        ),
        # Merged BAMs
        expand(
            merged_bam_dir / "treatment/{gen_abd_id}.merged.bam",
            gen_abd_id=gen_abd_ids,
        ),
        expand(
            merged_bam_dir / "treatment/{gen_abd_id}.merged.bam.bai",
            gen_abd_id=gen_abd_ids,
        ),
        expand(
            merged_bam_dir / "control/{ctr_abd_id}.merged.bam",
            ctr_abd_id=ctr_abd_ids,
        ),
        expand(
            merged_bam_dir / "control/{ctr_abd_id}.merged.bam.bai",
            ctr_abd_id=ctr_abd_ids,
        ),
        # bamCompare
        expand(
            bamCompare / "per_replicate/{base_id}.bw",
            base_id=tt_base_ids,
        ),
        expand(
            bamCompare / "merged/{gen_abd_id}.merged.bw",
            gen_abd_id=gen_abd_ids,
        ),
        # Deduplication
        expand(
            mark_remove_dups / "{base_id}.deduped.bam",
            base_id=base_ids,
        ),
        expand(
            mark_remove_dups / "{base_id}.deduped.bam.bai",
            base_id=base_ids,
        ),
        expand(
            stats_dir / "picard/{base_id}.dedup.metrics.txt",
            base_id=base_ids,
        ),
        # bedtools fragments
        expand(
            RESULTS_DIR / "bedtools/{base_id}.fragments.bedgraph",
            base_id=base_ids,
        ),
        # MultiQC
        RESULTS_DIR / "qc/multiqc/multiqc.html",
        expand(
            bw_dir / "average_bw/{cell_abd}.bw",
            cell_abd=all_abd_base_ids,
        ),
    ]

    # --------------------------
    # Conditional MACS3 peaks
    # --------------------------
    if config["peaks_call"]["macs3"]["use_macs3"]:
        if config["peaks_call"]["macs3"]["broad"]:
            TARGETS.extend(
                expand(
                    macs3_dir / "broad/peaks/{sample_id}_peaks.broadPeak",
                    sample_id=tt_ids,
                )
            )
        else:
            TARGETS.extend(
                expand(
                    macs3_dir / "narrow/peaks/{sample_id}_peaks.narrowPeak",
                    sample_id=tt_ids,
                )
            )

    # --------------------------
    # Conditional SEACR peaks
    # (treatment-only by construction)
    # --------------------------
    if config["peaks_call"]["seacr"]["use_seacr"]:
        TARGETS.extend(
            expand(
                dynamic_unfiltered_seacr_dir / "{mode}/{sample_id}.{mode}.bed",
                mode=seacr_modes,
                sample_id=tt_ids,
            )
        )

        TARGETS.extend(
            expand(
                dynamic_unfiltered_seacr_dir / "qc/prefiltered_counts.{mode}.{ext}",
                mode=seacr_modes,
                ext=["txt", "hist.txt"],
            )
        )

        TARGETS.extend(
            expand(
                dynamic_filtered_seacr_dir / "{mode}/{peak_id}.filtered.{mode}.bed",
                peak_id=peak_ids,
                mode=seacr_modes,
            )
        )

        TARGETS.extend(
            expand(
                dynamic_filtered_seacr_dir / "qc/filtered_counts.{mode}.{ext}",
                mode=seacr_modes,
                ext=["txt", "hist.txt"],
            )
        )

        TARGETS.extend(
            expand(
                dynamic_intersected_peaks / "{gen_abd_id}_rep1_v_rep2.{mode}.bed",
                gen_abd_id=gen_abd_ids,
                mode=seacr_modes,
            )
        )

        TARGETS.extend(
            expand(
                dynamic_greedy_consensus_peaks
                / "{gen_abd_id}_consensus.greedy.{mode}.bed",
                gen_abd_id=gen_abd_ids,
                mode=seacr_modes,
            )
        )

        TARGETS.extend(
            expand(
                dynamic_conservative_consensus_peaks
                / "{gen_abd_id}_consensus.conservative.{mode}.bed",
                gen_abd_id=gen_abd_ids,
                mode=seacr_modes,
            )
        )

        TARGETS.extend(
            expand(
                dynamic_reciprocal_consensus_peaks
                / "{gen_abd_id}_consensus.reciprocal.{mode}.bed",
                gen_abd_id=gen_abd_ids,
                mode=seacr_modes,
            )
        )
    if config["peaks_call"]["gopeaks"]["use_gopeaks"]:
        if config["peaks_call"]["gopeaks"]["use_broad"]:
            TARGETS.extend(
                expand(
                    gopeaks_dir / "peaks/broad/{sample_id}_peaks.bed",
                    sample_id=tt_ids,
                )
            )
        else:
            TARGETS.extend(
                expand(
                    gopeaks_dir / "peaks/narrow/{sample_id}_peaks.bed",
                    sample_id=tt_ids,
                )
            )
        # if config["peaks_call"]["lanceotron"]["use_lanceotron"]:
        #     TARGETS.extend(expand(lanceotron_dir / "{base_id}.done", base_id=tt_ids))

        # DiffBind outputs
        TARGETS.extend(
            [
                analysis_dir / "diffbind" / "pca" / "relaxed_pca.pdf",
                analysis_dir / "diffbind" / "pca" / "stringent_pca.pdf",
                analysis_dir / "diffbind" / "relaxed_counts.matrix.csv",
                analysis_dir / "diffbind" / "stringent_counts.matrix.csv",
                analysis_dir / "diffbind" / "diffbind_summary.txt",
                analysis_dir / "diffbind" / "pca" / "relaxed_H3K27ac_pca.pdf",
                analysis_dir / "diffbind" / "pca" / "stringent_H3K27ac_pca.pdf",
                analysis_dir / "diffbind" / "pca" / "relaxed_H3K27me3_pca.pdf",
                analysis_dir / "diffbind" / "pca" / "stringent_H3K27me3_pca.pdf",
                analysis_dir / "diffbind" / "pca" / "relaxed_H2K119Ub_pca.pdf",
                analysis_dir / "diffbind" / "pca" / "stringent_H2K119Ub_pca.pdf",
                analysis_dir / "diffbind" / "dba" / "relaxed_H3K27ac_dba.rds",
                analysis_dir / "diffbind" / "dba" / "stringent_H3K27ac_dba.rds",
                analysis_dir / "diffbind" / "dba" / "relaxed_H3K27me3_dba.rds",
                analysis_dir / "diffbind" / "dba" / "stringent_H3K27me3_dba.rds",
                analysis_dir / "diffbind" / "dba" / "relaxed_H2K119Ub_dba.rds",
                analysis_dir / "diffbind" / "dba" / "stringent_H2K119Ub_dba.rds",
                analysis_dir / "diffbind" / "counts" / "relaxed_H3K27ac_counts_dba.rds",
                analysis_dir
                / "diffbind"
                / "counts"
                / "stringent_H3K27ac_counts_dba.rds",
                analysis_dir
                / "diffbind"
                / "counts"
                / "relaxed_H3K27me3_counts_dba.rds",
                analysis_dir
                / "diffbind"
                / "counts"
                / "stringent_H3K27me3_counts_dba.rds",
                analysis_dir
                / "diffbind"
                / "counts"
                / "relaxed_H2K119Ub_counts_dba.rds",
                analysis_dir
                / "diffbind"
                / "counts"
                / "stringent_H2K119Ub_counts_dba.rds",
                analysis_dir / "h3k27ac" / "ac_lost_in_tbl3.bed",
                analysis_dir / "h3k27me3" / "bl3_me_lost_in_tbl3.bed",
                analysis_dir / "h3k27me3" / "me_gained_in_tbl3.bed",
                # analysis_dir / "h3k27_transition" / "ac_to_me_in_tbl3.bed",
                # analysis_dir / "h3k27_transition" / "me_to_ac_in_tbl3.bed",
                # analysis_dir / "diffbind" / "ma_plots" / "relaxed_H3K27ac_ma_plot.pdf",
                # analysis_dir
                # / "diffbind"
                # / "ma_plots"
                # / "stringent_H3K27ac_ma_plot.pdf",
                # analysis_dir / "diffbind" / "ma_plots" / "relaxed_H3K27me3_ma_plot.pdf",
                # analysis_dir
                # / "diffbind"
                # / "ma_plots"
                # / "stringent_H3K27me3_ma_plot.pdf",
                # analysis_dir / "diffbind" / "ma_plots" / "relaxed_H2K119Ub_ma_plot.pdf",
                # analysis_dir
                # / "diffbind"
                # / "ma_plots"
                # / "stringent_H2K119Ub_ma_plot.pdf",
                analysis_dir / "diffbind" / "rds" / "relaxed_H3K27ac_diffbind.rds",
                analysis_dir / "diffbind" / "rds" / "stringent_H3K27ac_diffbind.rds",
                analysis_dir / "diffbind" / "rds" / "relaxed_H3K27me3_diffbind.rds",
                analysis_dir / "diffbind" / "rds" / "stringent_H3K27me3_diffbind.rds",
                analysis_dir / "diffbind" / "rds" / "relaxed_H2K119Ub_diffbind.rds",
                analysis_dir / "diffbind" / "rds" / "stringent_H2K119Ub_diffbind.rds",
                analysis_dir / "fold_change/chromosome_occupancy_avg.csv",
                analysis_dir / "fold_change/chromosome_occupancy_per_rep.csv",
            ]
        )
        TARGETS.extend(
            expand(
                homer_dir / "peaks/{base_id}_peaks.txt",
                base_id=tt_ids,
            )
            + expand(homer_dir / "control/{ctr_id}", ctr_id=ctrl_abd_base_ids)
        )

        TARGETS.extend(
            [
                deeptools_dir / "matrices/matrix.gz",
                plots / "deeptools/h3k27me3_ac_h2k11ub.png",
                deeptools_dir / "matrices/scores_per_transcript.npz",
                deeptools_dir / "raw_counts/scores_per_transcript.tab",
                deeptools_dir / "matrices/scores_per_bin_per_rep.npz",
                deeptools_dir / "matrices/scores_per_bin_avg_bw.npz",
                deeptools_dir / "raw_counts/scores_per_bin_per_rep.tab",
                deeptools_dir / "raw_counts/scores_per_bin_avg_bw.tab",
                plots / "deeptools/scatterplot_PearsonCorr_bigwigScores.png",
                plots / "deeptools/PearsonCorr_bigwigScores.tab",
            ]
        )

    return TARGETS
