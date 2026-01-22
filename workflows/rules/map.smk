"""
CUT&RUN Analysis Pipeline - Snakemake Workflow
===============================================

Preprocessing workflow for CUT&RUN sequencing data.

Workflow Steps:
    1. QC & Trimming: FastQC, fastp
    2. Alignment: Bowtie2 (local, paired-end)
    3. Filtering: MAPQ, proper pairs, conditional MT removal
    4. Deduplication: Picard MarkDuplicates
    5. QC Stats: samtools stats, flagstat, idxstats, mosdepth
    6. [TODO] Peak Calling: SEACR

Key Design Decisions:
    - No ATAC-seq style shifting (inappropriate for MNase-based CUT&RUN)
    - Conditional MT removal based on contamination threshold
    - Duplicate removal configurable via config file

Author: Steve. G. ODETTE
Date: 2025-01-01
"""

import os
from pathlib import Path
from constants.dirs_files import *
CONTROL_ANTIBODIES = set(config["controls"]["antibodies"])

def is_control(base_id):
    """Determine if a sample is a control based on its antibody."""
    return any(ab in base_id for ab in CONTROL_ANTIBODIES)

rule map:
    input:
        idx=multiext(str(BT9_TA_REF_FA),
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
        raw_bt2_mapped=bt2_raw_map_dir / "{base_id}.raw.unsorted.unfiltered.bam",
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
        sorted_unfiltered_bam=temp(temp_sort_dir / "{base_id}.sorted.unfiltered.bam"
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
        wMT_seq_filtered_sorted_bam= temp(samtools_sorted_filtered_dir / "{base_id}.wMT_seq.sorted.filtered.bam"),
    threads: config["resources"]["bam_filter"]["cpu"]
    log:
        std_out="logs/samtools_sort/{base_id}.samtools.sort.filtered.out.log",
        std_err="logs/samtools_sort/{base_id}.samtools.sort.filtered.err.log",
    params:
        extras=config["bam_filter"]["extra"],
        samtools_sorted_filtered_dir=samtools_sorted_filtered_dir,
    resources:
        # ram=config["resources"]["bam_filter"]["mem_mb"],
    shell:
        """
        module load samtools/1.21
        samtools view -@ {threads} -b {input.unfiltered_sorted_bam} {params.extras} \
        > {output.wMT_seq_filtered_sorted_bam} 2> {log.std_err}
        """

rule check_MT_seq:
    input:
        wMT_seq_filtered_sorted_bam=rules.bam_filter.output.wMT_seq_filtered_sorted_bam,
    output:
        MT_stats=QC_DIR / "MT_seq/{base_id}.txt",
        temp_idx=temp(samtools_sorted_filtered_dir / "{base_id}.wMT_seq.sorted.filtered.bam.bai"),
    params:
        MT_seq=lambda wc: str(config["rmv_MT_seqs"]["MT_seq"])
    threads: config["resources"]["bam_filter"]["cpu"]
    shell:
        """
        module load samtools/1.21
        
        samtools index -@ {threads} {input.wMT_seq_filtered_sorted_bam}
        
        samtools idxstats {input.wMT_seq_filtered_sorted_bam} | \
        awk 'BEGIN {{total=0; mt=0}}
            {{total+=$3}} 
            $1 == "{params.MT_seq}" {{mt+=$3}}
            END {{if (total > 0) print "{wildcards.base_id}", total, mt, \
            100*mt/total; else print "{wildcards.base_id}", 0, 0, 0}}' \
        > {output.MT_stats}
        """

rule filter_MT:
    input:
        filtered_sorted_bam=rules.bam_filter.output.wMT_seq_filtered_sorted_bam,
        MT_stats=rules.check_MT_seq.output.MT_stats,
    output:
        filtered_sorted_bam=samtools_sorted_filtered_dir / "{base_id}.sorted.filtered.bam",
    params:
        MT_threshold=config["rmv_MT_seqs"]["mt_threshold"],
        MT_seq=lambda wc: str(config["rmv_MT_seqs"]["MT_seq"])
    threads: config["resources"]["bam_filter"]["cpu"]
    run:
        with open(input.MT_stats, 'r') as f:
            line = f.readline().strip().split()
            mt_pct = float(line[3])
        if mt_pct > params.MT_threshold:
            shell(
                """
                module load samtools/1.21
                
                samtools view -h {input.filtered_sorted_bam} | \
                awk '$1 ~ /^@/ || $3 != "{params.MT_seq}"' | \
                samtools view -b > {output.filtered_sorted_bam}
                """
            )
        else:
            shell("cp {input.filtered_sorted_bam} {output.filtered_sorted_bam}")

rule map_stats:
    input:
        filtered_sorted_bam=rules.filter_MT.output.filtered_sorted_bam
    output:
        stats=sorted_filtered_stats_dir / "{base_id}.stats.txt",
        flag_stats=sorted_filtered_stats_dir / "{base_id}.flagstats.txt"
    log:
        std_err = "logs/stats/{base_id}.stats.log"
    shell:
        """
        module load samtools/1.21
        samtools stats {input.filtered_sorted_bam} > {output.stats} 2> {log.std_err} && \
        samtools flagstat {input.filtered_sorted_bam} > {output.flag_stats} 2>> {log.std_err}
        """
    
rule duplicates:
    input:
        filtered_sorted_bam=rules.filter_MT.output.filtered_sorted_bam
    output:
        deduped = mark_remove_dups / "{base_id}.deduped.bam",
        dedup_stats = stats_dir / "picard/{base_id}.dedup.metrics.txt"
    log:
        "logs/markdedup/{base_id}.mark_rmv_duplicates.err.log"
    params:
        rmv_duplicates=lambda wc: str(
            config["picard"]["rmv_duplicates_control"]
            if is_control(wc.base_id)
            else config["picard"]["rmv_duplicates_tt"]
        ).lower(),
        extra=config["picard"].get("extra", "")
    shell:
        """
        module load picard/2.23.5

        picard MarkDuplicates \
            -I {input.filtered_sorted_bam} \
            -O {output.deduped} \
            -M {output.dedup_stats} \
            --REMOVE_DUPLICATES {params.rmv_duplicates} \
            {params.extra} \
            2> {log}
        """

rule samtools_idx:
    input:
        deduped = rules.duplicates.output.deduped
    output:
       idx_bam = mark_remove_dups / "{base_id}.deduped.bam.bai"
    threads: config["resources"]["samtools_idx"]["cpu"]
    shell:
        """
        module load samtools/1.21
        samtools index -@ {threads} -b {input.deduped} {output.idx_bam}
        """

rule bam_idx_stats:
    input:
        deduped=rules.duplicates.output.deduped,
        idx_bam=rules.samtools_idx.output.idx_bam  
    output:
        idxstats = stats_dir / "samtools/{base_id}.BT9_TA.dedupped.idxstats.txt"
    log:
        std_err="logs/samtools_stats/{base_id}.samtools.idx.err.log",
    shell:
        """
        module load samtools/1.21
        samtools idxstats {input.deduped} > {output.idxstats} 2> {log.std_err}
        """

rule mosedepth:
    input:
        sorted_bam = rules.duplicates.output.deduped,
        bam_idx = rules.samtools_idx.output.idx_bam
    output:
        regions = stats_dir / "mosdepth/{base_id}.regions.bed.gz", 
        global_stats = stats_dir / "mosdepth/{base_id}.mosdepth.global.dist.txt"

    log:
        deduped_std_err="logs/mosdepths/{base_id}.mosdepth.deduped.err.log",
    threads: config["resources"]["mosdepth"]["cpu"]
    params:
        binSize = config["mosdepth"]["binSize"],
        prefix = lambda wildcards, output: output.regions.replace(".regions.bed.gz", "")
    shell:
        """
        module load mosdepth/0.2.6
        mosdepth --threads {threads} -b {params.binSize} {params.prefix} {input.sorted_bam} 2> {log.deduped_std_err}
        """