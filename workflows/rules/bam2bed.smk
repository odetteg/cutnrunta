rule bam2bed:
    input:
        deduped=rules.duplicates.output.deduped,
        chrom_sizes=rules.get_cs.output.cs,
    output:
        temp_bed=temp(RESULTS_DIR / "bedtools/{base_id}.temp.bed"),
        clean_bed=RESULTS_DIR / "bedtools/{base_id}.clean.bed",
        fragments_bed=RESULTS_DIR / "bedtools/{base_id}.fragments.bed",
        bedgraph=RESULTS_DIR / "bedtools/{base_id}.fragments.bedgraph",
    conda:
        "../../envs/bedtools.yaml"
    threads: config["resources"]["samtools"]["cpu"]
    log:
        log="logs/bedtools/{base_id}.bam2bed.log",
    shell:
        """
        (
        # I do samtools sort -n first to sort by name. If not, bedtools will skip a lot of reads with unmatched mates.
        # If not sorted by name first, you'll loose a lot of data
        samtools sort -n -@ {threads} {input.deduped} | bedtools bamtobed -bedpe -i stdin > {output.temp_bed}
        awk '$1==$4 && $6-$2 < 1000 {{print $0}}' {output.temp_bed} > {output.clean_bed}
        cut -f 1,2,6 {output.clean_bed} | sort -k1,1 -k2,2n -k3,3n > {output.fragments_bed}
        bedtools genomecov -bg -i {output.fragments_bed} -g {input.chrom_sizes} > {output.bedgraph}
        ) &> {log.log}
        """


rule bed2wig:
    input:
        bedgraph=rules.bam2bed.output.bedgraph,
        chrom_sizes=rules.get_cs.output.cs,
    output:
        sorted_bedgraph=temp(bw_dir / "bedgraph_bw/{base_id}.sorted.bedgraph"),
        graph_bw=bw_dir / "bedgraph_bw/{base_id}.bw",
    conda:
        "../../envs/bedtools.yaml"
    log:
        "logs/bedtools/{base_id}.bed2wig.log",
    shell:
        """
        sort -k1,1 -k2,2n {input.bedgraph} > {output.sorted_bedgraph}
        bedGraphToBigWig {output.sorted_bedgraph} {input.chrom_sizes} {output.graph_bw} &> {log}
        """
