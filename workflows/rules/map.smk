rule map:
    input:
        idx=multiext(
            config["BT9_TA_REF_FA"],
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
        trimmed_r1=rules.fastp_trim.output.trim_r1,
        trimmed_r2=rules.fastp_trim.output.trim_r2,
    output:
        raw_bt2_mapped=config["dir_names"]["bt2_raw_map_dir"]
        + "/{base_id}.raw.unsorted.unfiltered.bam",
        # stats=config["dir_names"]["bt2_map_dir"] + "/{base_id}.stats",
    params:
        ref_idx=lambda wc, input: input[0].replace(".1.bt2", ""),
        min_len=config["map"]["min_len"],
        max_len=config["map"]["max_len"],
        MAPQ_cutoff=config["map"]["MAPQ_cutoff"],
        extras=config["map"]["extra"],
    threads: config["resources"]["map"]["cpu"]
    resources:
        runtime=config["resources"]["map"]["time"],
    log:
        bt2_std_err="logs/map/{base_id}.bt2.raw.unsorted.unfiltered.log",
        samtools_std_err="logs/map/{base_id}.samtools.raw.unsorted.unfiltered.err.log",
    shell:
        """
        set -e
        module load bowtie2/2.5.4
        module load samtools/1.21

        bowtie2 \
        {params.extras} \
        -I {params.min_len} \
        -X {params.max_len} \
        --threads {threads} \
        -1 {input.trimmed_r1} \
        -2 {input.trimmed_r2} \
        -x {params.ref_idx} \
        --rg-id '{wildcards.base_id}' --rg 'SM:{wildcards.base_id}' 2> {log.bt2_std_err} | \
        samtools view -bhS > {output.raw_bt2_mapped} 2> {log.samtools_std_err}
        """


rule bam_sort:
    input:
        raw_bam=rules.map.output.raw_bt2_mapped,
    output:
        sorted_unfiltered_bam=temp(
            config["dir_names"]["temp_sort_dir"] + "/{base_id}.sorted.unfiltered.bam"
        ),
    log:
        std_out="logs/samtools_sort/{base_id}.samtools.sort.out.log",
        std_err="logs/samtools_sort/{base_id}.samtools.sort.err.log",
    threads: config["resources"]["bam_sort"]["cpu"]
    resources:
        ram=config["resources"]["bam_sort"]["mem_mb"],
    shell:
        """

        module load samtools/1.21

        samtools sort -@ {threads} \
        {input.raw_bam} \
        -o {output.sorted_unfiltered_bam} \
        1> {log.std_out} 2> {log.std_err}
        """


rule bam_filter:
    input:
        unfiltered_sorted_bam=rules.bam_sort.output.sorted_unfiltered_bam,
    output:
        filtered_sorted_bam=config["dir_names"]["samtools_sorted_filtered_dir"]
        + "/{base_id}.sorted.filtered.bam",
    threads: config["resources"]["bam_filter"]["cpu"]
    log:
        std_out="logs/samtools_sort/{base_id}.samtools.sort.filtered.out.log",
        std_err="logs/samtools_sort/{base_id}.samtools.sort.filtered.err.log",
    params:
        extras=config["bam_filter"]["extra"],
        samtools_sorted_filtered_dir=config["dir_names"]["samtools_sorted_filtered_dir"],
    resources:
        # ram=config["resources"]["bam_filter"]["mem_mb"],
    shell:
        """
        module load samtools/1.21
        samtools view -@ {threads} -b {input.unfiltered_sorted_bam} > {output.filtered_sorted_bam} \
        2> {log.std_err}

        """

rule map_stats:
    input:
        filtered_sorted_bam=rules.bam_filter.output.filtered_sorted_bam
    output:
        stats=config["dir_names"]["sorted_filtered_stats_dir"] + "/{base_id}.stats.txt",
        flag_stats=config["dir_names"]["sorted_filtered_stats_dir"] + "/{base_id}.flagstats.txt"
    log:
        std_err = "logs/stats/{base_id}.stats.log"
    shell:
        """
        module load samtools/1.21
        samtools stats {input.filtered_sorted_bam} > {output.stats} 2> {log.std_err} && \
        samtools flagstat {input.filtered_sorted_bam} > {output.flag_stats} 2>> {log.std_err}
        """


# rule samtools_idx:
#     input:
#         filtered_sorted_bam = rules.bam_filter.filtered_sorted_bam
#     output:
#         filtered_sorted__index_bam = config["samtools_sorted_filtered_dir"] + "/{base_id}.sorted.filtered.bam.bai"
#     log:
#         std_out="logs/samtools_sort/{base_id}.samtools.index.out.log",
#         std_err="logs/samtools_sort/{abse_id}.samtools.index.err.log"
#     threads: config["resources"]["samtools_idx"]["cpu"]
#     shell:
#         """
#         """
