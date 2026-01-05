from constants.dirs_files import *

if config["peaks_call"]["macs3"]["use_macs3"]:
    if not config["peaks_call"]["macs3"]["broad"]:

        rule macs3_narrow:
            input:
                tt_bam=mark_remove_dups / "{sample_id}.deduped.bam",
                tt_bai=mark_remove_dups / "{sample_id}.deduped.bam.bai",
                ctrl_bam=lambda wc: get_ctrl_bam(wc.sample_id),
                ctrl_bai=lambda wc: get_ctrl_bai(wc.sample_id),
                egs=rules.merge_sample_lens.output.sample_egs_csv,
            output:
                peak=macs3_dir / "narrow/peaks/{sample_id}_peaks.narrowPeak",
            params:
                mode="narrow",
                use_control=config["peaks_call"]["macs3"]["use_control"],
                qvalue=config["peaks_call"]["macs3"]["qvalue"],
                extra=config["peaks_call"]["macs3"]["extra"],
                outdir=lambda wc, output: os.path.dirname(output[0]),
            threads: config["resources"]["macs3"]["cpu"]
            conda:
                "../../envs/peaks.yaml"
            log:
                std_out="logs/macs3/{sample_id}_peaks.narrowPeak.log",
            script:
                "../scripts/macs3.py"

    else:

        rule macs3_broad:
            input:
                tt_bam=mark_remove_dups / "{sample_id}.deduped.bam",
                tt_bai=mark_remove_dups / "{sample_id}.deduped.bam.bai",
                ctrl_bam=lambda wc: get_ctrl_bam(wc.sample_id),
                ctrl_bai=lambda wc: get_ctrl_bai(wc.sample_id),
                egs=rules.merge_sample_lens.output.sample_egs_csv,
            output:
                peak=macs3_dir / "broad/peaks/{sample_id}_peaks.broadPeak",
            params:
                mode="broad",
                use_control=config["peaks_call"]["macs3"]["use_control"],
                qvalue=config["peaks_call"]["macs3"]["qvalue"],
                extra=config["peaks_call"]["macs3"]["extra"],
                broad_cutoff=config["peaks_call"]["macs3"]["broad_cutoff"],
                outdir=lambda wc, output: os.path.dirname(output[0]),
            threads: config["resources"]["macs3"]["cpu"]
            conda:
                "../../envs/peaks.yaml"
            log:
                std_out="logs/macs3/{sample_id}_peaks.broadPeaks.log",
            script:
                "../scripts/macs3.py"


if config["peaks_call"]["seacr"]["use_seacr"]:
    if config["peaks_call"]["seacr"]["use_control"]:

        rule seacr_relaxed:
            input:
                treatment_bg=lambda wc: RESULTS_DIR
                / f"bedtools/{wc.sample_id}.fragments.bedgraph",
                control_bg=lambda wc: get_ctrl_bg(wc.sample_id),
            output:
                _peaks=unfiltered_seacr_peaks_dir / "peaks_wcontrol/relaxed/{sample_id}.relaxed.bed",
            params:
                out_prefix=lambda wc: str(
                    unfiltered_seacr_peaks_dir / f"peaks_wcontrol/relaxed/{wc.sample_id}"
                ),
                threshold=config["peaks_call"]["seacr"]["peak_threshold"],
                use_control=config["peaks_call"]["seacr"]["use_control"],
                mode="relaxed",
            conda:
                "../../envs/seacr.yaml"
            log:
                "logs/seacr/{sample_id}_relaxed.log",
            script:
                "../scripts/seacr.py"

        rule seacr_stringent:
            input:
                treatment_bg=lambda wc: RESULTS_DIR
                / f"bedtools/{wc.sample_id}.fragments.bedgraph",
                control_bg=lambda wc: get_ctrl_bg(wc.sample_id),
            output:
                _peaks=unfiltered_seacr_peaks_dir / "peaks_wcontrol/stringent/{sample_id}.stringent.bed",
            params:
                out_prefix=lambda wc: str(
                    unfiltered_seacr_peaks_dir / f"peaks_wcontrol/stringent/{wc.sample_id}"
                ),
                threshold=config["peaks_call"]["seacr"]["peak_threshold"],
                use_control=config["peaks_call"]["seacr"]["use_control"],
                mode="stringent",
            conda:
                "../../envs/seacr.yaml"
            log:
                "logs/seacr/{sample_id}_stringent.log",
            script:
                "../scripts/seacr.py"

    else:

        rule seacr_relaxed_noctrl:
            input:
                treatment_bg=lambda wc: RESULTS_DIR
                / f"bedtools/{wc.sample_id}.fragments.bedgraph",
                # control_bg=lambda wc: get_ctrl_bg(wc.sample_id),
            output:
                _peaks=unfiltered_seacr_peaks_dir / "peaks_nocontrol/relaxed/{sample_id}.relaxed.bed",
            params:
                out_prefix=lambda wc: str(
                    unfiltered_seacr_peaks_dir / f"peaks_nocontrol/relaxed/{wc.sample_id}"
                ),
                threshold=config["peaks_call"]["seacr"]["peak_threshold"],
                use_control=config["peaks_call"]["seacr"]["use_control"],
                mode="relaxed",
            conda:
                "../../envs/seacr.yaml"
            log:
                "logs/seacr/{sample_id}_relaxed.log",
            script:
                "../scripts/seacr.py"

        rule seacr_stringent_noctrl:
            input:
                treatment_bg=lambda wc: RESULTS_DIR
                / f"bedtools/{wc.sample_id}.fragments.bedgraph",
                # control_bg=lambda wc: get_ctrl_bg(wc.sample_id),
            output:
                _peaks=unfiltered_seacr_peaks_dir / "peaks_nocontrol/stringent/{sample_id}.stringent.bed",
            params:
                out_prefix=lambda wc: str(
                    unfiltered_seacr_peaks_dir / f"peaks_nocontrol/stringent/{wc.sample_id}"
                ),
                threshold=config["peaks_call"]["seacr"]["peak_threshold"],
                use_control=config["peaks_call"]["seacr"]["use_control"],
                mode="stringent",
            conda:
                "../../envs/seacr.yaml"
            log:
                "logs/seacr/{sample_id}_stringent.log",
            script:
                "../scripts/seacr.py"
