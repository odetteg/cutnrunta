rule bt2_build_index:
    input:
        ref_genome=config["BT9_TA_REF_FA"]
    output:
       idx=multiext(
            config["BT9_TA_REF_FA"],
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        std_out="logs/bt2/bt2_build.out.log",
        std_err="logs/bt2/bt2_build.err.log",
    threads: config["resources"]["bt2_build_index"]["cpu"]
    resources:
        ram=config["resources"]["bt2_build_index"]["mem_mb"],
    shell:
        """
        module load bowtie2/2.5.4
        bowtie2-build --threads {threads} {input.ref_genome} {input.ref_genome} > {log.std_out} 2> {log.std_err}

        """


rule samtools_faidx:
    input:
        ref_genome=config["BT9_TA_REF_FA"],
    output:
        fai=config["BT9_TA_REF_FA"] + ".fai",
    log:
        std_err="logs/samtools/samtools_faidx.err.log",
    shell:
        """
        module load samtools/1.21
        samtools faidx {input.ref_genome} 2> {log.std_err}

        """
