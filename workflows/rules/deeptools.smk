from constants.dirs_files import *


rule multibigwigsummary:
    input:
        bw=expand(bw_dir / "replicate_bw/{base_id}.bw", base_id=base_ids),
    output:
        npz=deeptools_dir / "scores_per_bin.npz",
        csv=deeptools_dir / "qc/scores_per_bin.csv",
    threads: config["resources"]["deeptools"]["cpu"]
    log:
        "logs/deeptools/multibigwigsummary.log",
    params:
        extra=config["deeptools"]["multibigwigsummary"]["extra"],
        names=" ".join(base_ids),
    conda:
        "../../envs/deeptools.yaml"
    shell:
        """
        multiBigwigSummary bins \
        --bwfiles {input.bw} \
        --labels {params.names} \
        --outFileName {output.npz} \
        --numberOfProcessors {threads} \
        --outRawCounts {output.csv} \
        {params.extra} \
        > {log} 2>&1
        """


rule PCA:
    input:
        npz=rules.multibigwigsummary.output.npz,
    output:
        pca_tab=deeptools_dir / "pca.tab",
    params:
        extra=config["deeptools"]["plotpca"]["extra"],
    conda:
        "../../envs/deeptools.yaml"
    log:
        "logs/deeptools/plotpca.log",
    shell:
        """
        plotPCA \
        --corData {input.npz} \
        --outFileNameData {output.pca_tab} \
        --transpose \
        {params.extra} \
        > {log} 2>&1
        """


rule BAM_fragment_sizes:
    input:
        deduped_bam=expand(mark_remove_dups / "{base_id}.deduped.bam", base_id=base_ids),
        deduped_bam_bai=expand(
            mark_remove_dups / "{base_id}.deduped.bam.bai", base_id=base_ids
        ),
    output:
        hist=plots / "deeptools/fragment_lengths.pdf",
        table=deeptools_dir / "qc/fragment_lengths.tsv",
        raw=deeptools_dir / "qc/fragment_lengths_raw.tsv",
    params:
        labels=" ".join(base_ids),
        max_len=config["map"]["max_len"],
    threads: config["resources"]["deeptools"]["cpu"]
    log:
        "logs/deeptools/BAM_fragment_sizes.log",
    conda:
        "../../envs/deeptools.yaml"
    shell:
        """
        bamPEFragmentSize \
        --numberOfProcessors {threads} \
        --maxFragmentLength {params.max_len} \
        --bamfiles {input.deduped_bam} \
        --histogram {output.hist} \
        --table {output.table} \
        --outRawFragmentLengths {output.raw} \
        --samplesLabel {params.labels} \
        --plotTitle "Fragment Length Distribution: Nucleosomal Ladder" \
        > {log} 2>&1
        """


rule plotFingerPrint:
    input:
        deduped_bam=expand(mark_remove_dups / "{base_id}.deduped.bam", base_id=base_ids),
        deduped_bam_bai=expand(
            mark_remove_dups / "{base_id}.deduped.bam.bai", base_id=base_ids
        ),
    output:
        finger_print=plots / "deeptools/fingerprint.pdf",
        raw_tsv=deeptools_dir / "qc/fingerprints.tsv",
    params:
        extra=config["deeptools"]["plotfingerprint"]["extra"],
        max_len=config["map"]["max_len"],
        binsize=config["deeptools"]["binsize"],
        labels=" ".join(base_ids),
    threads: config["resources"]["deeptools"]["cpu"]
    conda:
        "../../envs/deeptools.yaml"
    log:
        "logs/deeptools/plotfingerprint.log",
    shell:
        """
        plotFingerprint \
        --bamfiles {input.deduped_bam} \
        --plotFile {output.finger_print} \
        --maxFragmentLength {params.max_len} \
        --labels {params.labels} \
        --binSize {params.binsize} \
        --plotTitle "Fingerprints of different samples" \
        --numberOfProcessors {threads} \
        --outRawCounts {output.raw_tsv} \
        {params.extra} \
        > {log} 2>&1
        """
