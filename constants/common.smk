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
        [str(BT9_TA_REF_FA) + ".fai"],
        # Mapping outputs
        expand(
            bt2_raw_map_dir / "{base_id}.raw.unsorted.unfiltered.bam",
            base_id=base_ids,
        ),
        expand(
            samtools_sorted_filtered_dir / "{base_id}.sorted.filtered.bam",
            base_id=base_ids,
        ),
        expand(
            merged_bam_dir / "treatment/{gen_abd_id}.merged.bam", gen_abd_id=gen_abd_ids
        ),
        expand(
            merged_bam_dir / "treatment/{gen_abd_id}.merged.bam.bai",
            gen_abd_id=gen_abd_ids,
        ),
        expand(
            bamCompare / "per_replicate/{base_id}.bw",
            base_id=tt_base_ids,
        ),
        expand(
            bamCompare / "merged/{gen_abd_id}.merged.bw",
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
        [analysis_dir / "h3k27ac" / "ac_lost_in_tbl3.bed"],
        [analysis_dir / "h3k27me3" / "me_gained_in_tbl3.bed"],
        [analysis_dir / "h3k27_transition" / "ac_to_me_in_tbl3.bed"],
        # Aggregated stats
        [agg_stats_csv],
        [samtools_txt],
        # Fragment length plots
        [frag_len_dst_pdf],
        [frag_len_dst_facet_pdf],
        # Samtools stats
        expand(
            sorted_filtered_stats_dir / "{base_id}.stats.txt",
            base_id=base_ids,
        ),
        expand(
            sorted_filtered_stats_dir / "{base_id}.flagstats.txt",
            base_id=base_ids,
        ),
        # Deduplication
        expand(
            mark_remove_dups / "{base_id}.deduped.bam",
            base_id=base_ids,
        ),
        expand(
            stats_dir / "picard/{base_id}.dedup.metrics.txt",
            base_id=base_ids,
        ),
        expand(
            mark_remove_dups / "{base_id}.deduped.bam.bai",
            base_id=base_ids,
        ),
        # idxstats
        expand(
            stats_dir / "samtools/{base_id}.BT9_TA.dedupped.idxstats.txt",
            base_id=base_ids,
        ),
        # Mosdepth
        expand(
            stats_dir / "mosdepth/{base_id}.{ext}",
            base_id=base_ids,
            ext=["mosdepth.global.dist.txt", "regions.bed.gz"],
        ),
        # EGS
        [BASE_DIR / "ad_files/egs.json"],
        [sample_egs_csv],
        # BigWigs
        expand(
            bw_dir / "replicate_bw/{base_id}.bw",
            base_id=base_ids,
        ),
        expand(QC_DIR / "MT_seq/{base_id}.txt", base_id=base_ids),
        expand(
            bw_dir / "average_bw/{cell_line}_{antibody}.bw",
            cell_line=cell_lines,
            antibody=antibodies,
        ),
        # deepTools: bins + PCA
        [deeptools_dir / "scores_per_bin.npz"],
        [deeptools_dir / "pca.tab"],
        expand(bw_dir / "bedgraph_bw/{base_id}.bw", base_id=base_ids),
        # Fragment length tables
        [plots / "deeptools/fragment_lengths.pdf"],
        # Fingerprints
        [plots / "deeptools/fingerprint.pdf"],
        [deeptools_dir / "qc/fingerprints.tsv"],
        # deepTools: multiBigwigSummary
        [deeptools_dir / "qc/scores_per_bin.csv"],
        # Fragment length tables
        [deeptools_dir / "qc/fragment_lengths.tsv"],
        [deeptools_dir / "qc/fragment_lengths_raw.tsv"],
        # bedtools
        expand(RESULTS_DIR / "bedtools/{base_id}.fragments.bedgraph", base_id=base_ids),
        # MultiQC
        [RESULTS_DIR / "qc/multiqc/multiqc.html"],
    ]

    # Conditional MACS3 peaks
    if config["peaks_call"]["macs3"]["use_macs3"]:
        if not config["peaks_call"]["macs3"]["broad"]:
            TARGETS.extend(
                expand(
                    macs3_dir / "narrow/peaks/{sample_id}_peaks.narrowPeak",
                    sample_id=design_df.query("condition == 'treatment'")["base_id"],
                )
            )
        else:
            TARGETS.extend(
                expand(
                    macs3_dir / "broad/peaks/{sample_id}_peaks.broadPeak",
                    sample_id=design_df.query("condition == 'treatment'")["base_id"],
                )
            )

    # Conditional SEACR peaks
    if config["peaks_call"]["seacr"]["use_seacr"]:
        TARGETS.extend(
            expand(
                dynamic_unfiltered_seacr_dir / "{mode}/{sample_id}.{mode}.bed",
                mode=seacr_modes,
                sample_id=design_df.query("condition=='treatment'")["base_id"],
            )
            + expand(
                dynamic_unfiltered_seacr_dir / "qc/prefiltered_counts.{mode}.{ext}",
                mode=seacr_modes,
                ext=["txt", "hist.txt"],
            )
            + expand(
                dynamic_filtered_seacr_dir / "{mode}/{peak_id}.filtered.{mode}.bed",
                peak_id=peak_ids,
                mode=seacr_modes,
            )
            + expand(
                dynamic_filtered_seacr_dir / "qc/filtered_counts.{mode}.{ext}",
                mode=seacr_modes,
                ext=["txt", "hist.txt"],
            )
            + expand(
                dynamic_intersected_peaks / "{gen_abd_id}_rep1_v_rep2.{mode}.bed",
                gen_abd_id=gen_abd_ids,
                mode=seacr_modes,
            )
            + expand(
                dynamic_greedy_consensus_peaks
                / "{gen_abd_id}_consensus.greedy.{mode}.bed",
                gen_abd_id=gen_abd_ids,
                mode=seacr_modes,
            )
            + expand(
                dynamic_conservative_consensus_peaks
                / "{gen_abd_id}_consensus.conservative.{mode}.bed",
                gen_abd_id=gen_abd_ids,
                mode=seacr_modes,
            )
            + expand(
                dynamic_reciprocal_consensus_peaks
                / "{gen_abd_id}_consensus.reciprocal.{mode}.bed",
                gen_abd_id=gen_abd_ids,
                mode=seacr_modes,
            )
        )
    return TARGETS
