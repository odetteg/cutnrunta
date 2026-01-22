from constants.dirs_files import *


rule bigwig:
    input:
        unpack(bw_input),
        egs_json=rules.get_default_effective_genome_size.output.egs_json,
    output:
        bw_=bw_dir / "replicate_bw/{base_id}.bw",
    params:
        genome=config["genome"],
        binsize=config["deeptools"]["binsize"],
        smoothLength = config["deeptools"]["smoothLength"],
        normalization=config["deeptools"]["normalization"],
        extras=config["deeptools"]["extra"],
        rmv_MT_seqs=config["deeptools"]["rmv_MT_seqs"],
    threads: config["resources"]["deeptools"]["cpu"]
    conda:
        "../../envs/deeptools.yaml"
    log:
        "logs/deeptools/{base_id}.bigwig.log",
    script:
        "../scripts/bigwig.py"


rule avg_bw:
    input:
        bw_files=lambda wc: [
            bw_dir / f"replicate_bw/{base_id}.bw"
            for base_id in base_ids
            if base_id.startswith(wc.cell_abd)  
        ],
    output:
        wg_=temp(bw_dir / "average_bw/{cell_abd}.wig"),
    params:
        extra=config["deeptools"]["extra"],
    threads: config["resources"]["deeptools"]["cpu"]
    log:
        "logs/deeptools/average_bw/{cell_abd}.log",
    conda:
        "../../envs/deeptools.yaml"
    script:
        "../scripts/average_bw.py"

rule get_cs:
    input:
        fai=rules.samtools_faidx.output.fai,
    output:
        cs=REF_DIR / "BT9_TA.chrom.sizes",
    shell:
        """
        cut -f1,2 {input.fai} > {output.cs}
        """


rule wig2bw:
    input:
        wig=rules.avg_bw.output.wg_,
        cs=rules.get_cs.output.cs,
    output:
        avg_bw=bw_dir / "average_bw/{cell_abd}.bw",
    log:
        "logs/deeptools/wig2bw/{cell_abd}.log",
    conda:
        "../../envs/deeptools.yaml"
    shell:
        """
        wigToBigWig {input.wig} {input.cs} {output.avg_bw}
        """
