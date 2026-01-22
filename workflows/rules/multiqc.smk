from constants.dirs_files import *
rule multiqc:
    input:
        # FastQC on raw FASTQs
        fastqc_raw_fqs=expand(
            RAW_FASTQC_DIR / "{sample_id}_fastqc.{ext}",
            sample_id=sample_ids,
            ext=["html", "zip"],
        ),

        # fastp reports
        fastp_rpts=expand(FASTP_QC_REPORTS_DIR / "{base_id}.fastp.{ext}",
            base_id=base_ids,
            ext=["html", "json"],
        ),

        # FastQC after fastp trimming
        fastqc_fastp_trimmed=expand(FASTQC_QC_FASTP_DIR / "{base_id}.fastp.trimmed.{read}_fastqc.{ext}",
            base_id=base_ids,
            read=["R1", "R2"],
            ext=["html", "zip"],
        ),

        # Bowtie2 mapping logs
        bt2_raw_map_stats=expand(
            "logs/map/{base_id}.bt2.raw.unsorted.unfiltered.log",
            base_id=base_ids,
        ),

        # Samtools stats
        map_filter_stats=expand(
            sorted_filtered_stats_dir  / "{base_id}.{stat}.txt",
            base_id=base_ids,
            stat=["flagstats", "stats"],
        ),
        MT_stats=expand(QC_DIR / "MT_seq/{base_id}.txt", base_id=base_ids),

        # idxstats
        idxstats=expand(
            stats_dir / "samtools/{base_id}.BT9_TA.dedupped.idxstats.txt",
            base_id=base_ids,
        ),

        # Mosdepth
        mosdepth=expand(
            stats_dir / "mosdepth/{base_id}.{ext}",
            base_id=base_ids,
            ext=["mosdepth.global.dist.txt", "regions.bed.gz"],
        ),

        # deepTools: multiBigwigSummary
        bins_per_rep_npz=deeptools_dir / "matrices/scores_per_bin_per_rep.npz",
        bins_avg_npz=deeptools_dir / "matrices/scores_per_bin_avg_bw.npz",
        score_bin_per_rep=deeptools_dir / "raw_counts/scores_per_bin_per_rep.tab",
        score_bin_per_avg = deeptools_dir / "raw_counts/scores_per_bin_avg_bw.tab",
        chromosome_occupancy_avg =   analysis_dir / "fold_change/chromosome_occupancy_avg.csv",
        chromosome_occupancy_per_rep =    analysis_dir / "fold_change/chromosome_occupancy_per_rep.csv",

        # deepTools: PCA
        pca_tab=deeptools_dir / "matrices/pca.tab",

        # Fragment length tables
        frag_table=deeptools_dir / "qc/fragment_lengths.tsv",
        frag_raw=deeptools_dir / "qc/fragment_lengths_raw.tsv",

        # Fingerprints
        raw_tsv=deeptools_dir / "qc/fingerprints.tsv",

    output:
        multqc_rpt=RESULTS_DIR / "qc/multiqc/multiqc.html",
    params:
        mqc_dir=RESULTS_DIR / "qc/multiqc",
        extra=config["multiqc"]["extra"],
        scan_dirs = f"{RESULTS_DIR} {stats_dir} {sorted_filtered_stats_dir} {deeptools_dir}/qc logs/ {analysis_dir}/fold_change"
    threads: config["resources"]["multiqc"]["cpu"]
    log:
        std_out="logs/multiqc/multiqc.out.log",
        std_err="logs/multiqc/multiqc.err.log", 
    shell:
        """
        module load -f multiqc/1.29
        mkdir -p {params.mqc_dir} logs/multiqc
        multiqc --force {params.extra} \
                --outdir {params.mqc_dir} \
                -n multiqc.html \
                {params.scan_dirs} \
                1> {log.std_out} 2> {log.std_err}
        """
