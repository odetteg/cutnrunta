from constants.dirs_files import *


rule multibigwigsummary_per_rep:
    input:
        bw=expand(bw_dir / "replicate_bw/{base_id}.bw", base_id=base_ids),
    output:
        npz=deeptools_dir / "matrices/scores_per_bin_per_rep.npz",
        raw_counts=deeptools_dir / "raw_counts/scores_per_bin_per_rep.tab",
    threads: config["resources"]["deeptools"]["cpu"]
    log:
        "logs/deeptools/multibigwigsummary_per_rep.log",
    params:
        extra=config["deeptools"]["multiBigwigSummary"]["extra"],
        binSize=config["deeptools"]["multiBigwigSummary"]["binSize"],
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
        --outRawCounts {output.raw_counts} \
        {params.extra} \
        > {log} 2>&1
        """


rule multibigwigsummary_avg_bw:
    input:
        rep_bw=expand(bw_dir / "average_bw/{gen_abd}.bw", gen_abd=gen_abd_ids),
        avg_bw=expand(bw_dir / "average_bw/{ctr_abd}.bw", ctr_abd=ctr_abd_ids),
    output:
        npz=deeptools_dir / "matrices/scores_per_bin_avg_bw.npz",
        raw_counts=deeptools_dir / "raw_counts/scores_per_bin_avg_bw.tab",
    threads: config["resources"]["deeptools"]["cpu"]
    log:
        "logs/deeptools/multibigwigsummary_avg_bw.log",
    params:
        extra=config["deeptools"]["multiBigwigSummary"]["extra"],
        binSize=config["deeptools"]["multiBigwigSummary"]["binSize"],
        names=" ".join(base_ids),
    conda:
        "../../envs/deeptools.yaml"
    shell:
        """
        multiBigwigSummary bins \
        --bwfiles {input.rep_bw} {input.avg_bw}\
        --smartLabels \
        --outFileName {output.npz} \
        --numberOfProcessors {threads} \
        --outRawCounts {output.raw_counts} \
        {params.extra} \
        > {log} 2>&1
        """


rule PCA:
    input:
        npz=rules.multibigwigsummary_per_rep.output.npz,
    output:
        pca_tab=deeptools_dir / "matrices/pca.tab",
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


rule gtf_to_bed:
    input:
        gtf=config["GTF"],
    output:
        bed=REF_DIR / "genes.bed",
    log:
        "logs/preprocessing/gtf_to_bed.log",
    shell:
        """
        awk 'BEGIN{{OFS="\\t"}} $3=="gene" {{print $1, $4-1, $5, $10, ".", $7}}' {input.gtf} | \\
        tr -d '";' > {output.bed} 2> {log}
        """


rule computeMatrix:
    input:
        avg_bw=expand(
            bw_dir / "average_bw/{tt_abd_ids}.bw",
            tt_abd_ids=sorted(set(tt_abd_ids)),
        ),
    output:
        matrix=deeptools_dir / "matrices/matrix.gz",
    params:
        args=get_compute_matrix_args(),
        extra=config["deeptools"]["computeMatrix"]["extra"],
    threads: config["resources"]["deeptools"]["cpu"]
    conda:
        "../../envs/deeptools.yaml"
    log:
        "logs/deeptools/computeMatrix.log",
    shell:
        """
        computeMatrix {params.args} \
        -S {input.avg_bw} \
        --outFileName {output.matrix} \
        --numberOfProcessors {threads} \
        --smartLabels \
        {params.extra} \
        > {log} 2>&1
        """


rule plotHeatmap:
    input:
        mat=rules.computeMatrix.output.matrix,
    output:
        heatMap=plots / "deeptools/h3k27me3_ac_h2k11ub.png",
        # outFileSortedRegions = plots / "deeptools/outFileSortedRegions."
    conda:
        "../../envs/deeptools.yaml"
    log:
        "logs/deeptools/plotHeatmap.log",
    shell:
        """
        plotHeatmap -m {input.mat} -o {output.heatMap} > {log} 2>&1
        """


rule compute_per_transcript_npz:
    input:
        avg_bw=expand(
            bw_dir / "average_bw/{tt_abd_ids}.bw",
            tt_abd_ids=sorted(set(tt_abd_ids)),
        ),
        bed=REF_DIR / "genes.bed",
    output:
        npz=deeptools_dir / "matrices/scores_per_transcript.npz",
        raw_count=deeptools_dir / "raw_counts/scores_per_transcript.tab",
    params:
        extra=config["deeptools"]["multiBigwigSummary"]["per_transcript_npz"]["extra"],
    conda:
        "../../envs/deeptools.yaml"
    log:
        "logs/deeptools/compute_per_transcript_npz.log",
    shell:
        """
        multiBigwigSummary BED-file \
        --bwfiles {input.avg_bw} \
        --BED {input.bed} \
        --smartLabels \
        {params.extra} \
        -out {output.npz} --outRawCounts {output.raw_count} > {log} 2>&1
        """


rule plotCorrelation:
    input:
        npz=deeptools_dir / "matrices/scores_per_transcript.npz",
    output:
        pearsonCorr=plots / "deeptools/scatterplot_PearsonCorr_bigwigScores.png",
        pearsonCorrScores=plots / "deeptools/PearsonCorr_bigwigScores.tab",
    params:
        corMethod=config["deeptools"]["plotCorrelation"]["corMethod"],
        extra=config["deeptools"]["plotCorrelation"]["extra"],
        plotTitle=config["deeptools"]["plotCorrelation"]["plotTitle"],
        whatToPlot=config["deeptools"]["plotCorrelation"]["whatToPlot"],
    conda:
        "../../envs/deeptools.yaml"
    log:
        "logs/deeptools/plotCorrelation.log",
    shell:
        """
        plotCorrelation \
        -in {input.npz} \
        --corMethod {params.corMethod} {params.extra} \
        --plotTitle "{params.plotTitle}" \
        --whatToPlot {params.whatToPlot} \
        -o {output.pearsonCorr} \
        --outFileCorMatrix {output.pearsonCorrScores} > {log} 2>&1
        """
