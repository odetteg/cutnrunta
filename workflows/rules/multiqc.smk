rule multiqc:
    input:
        fastqc_raw_fqs=expand(
            RAW_FASTQC_DIR / "{sample_id}_fastqc.{ext}",
            sample_id=sample_ids,
            ext=["html", "zip"],
        ),
        fastp_rpts=expand(
            config["FASTP_QC_REPORTS_DIR"] + "/{base_id}.fastp.{ext}",
            ext=["html", "json"],
            base_id=base_ids,
        ),
        fqc_fastp_trimm_rpts=expand(
            config["FASTQC_QC_FASTP_DIR"]
            + "/{base_id}.fastp.trimmed.{read}_fastqc.{ext}",
            base_id=base_ids,
            read=["R1", "R2"],
            ext=["html", "zip"],
        ),
        bt2_raw_map_stats=expand("logs/map/{base_id}.bt2.raw.unsorted.unfiltered.log",
        base_id=base_ids),
        map_filter_stats=expand(config["dir_names"]["sorted_filtered_stats_dir"] + "/{base_id}.{stat}.txt",
        base_id=base_ids, stat=["flagstats", "stats"])
    output:
        multqc_rpt=RESULTS_DIR / "qc/multiqc/multiqc.html",
    params:
        mqc_dir=RESULTS_DIR / "qc/multiqc",
        extra=config["multiqc"]["extra"],
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
                {input} \
                1> {log.std_out} 2> {log.std_err}
        """
