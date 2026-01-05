import os
from constants.dirs_files import *


rule qc_count:
    input:
        stringent_peaks=lambda wc: [
            dynamic_unfiltered_seacr_dir / f"stringent/{pid}.stringent.bed"
            for pid in peak_ids
        ],
        relaxed_peaks=lambda wc: [
            dynamic_unfiltered_seacr_dir / f"relaxed/{pid}.relaxed.bed"
            for pid in peak_ids
        ],
    output:
        prefiltered_counts_stringent=dynamic_unfiltered_seacr_dir
        / "qc/prefiltered_counts.stringent.txt",
        prefiltered_counts_relaxed=dynamic_unfiltered_seacr_dir
        / "qc/prefiltered_counts.relaxed.txt",
        counts_histogram_stringent=dynamic_unfiltered_seacr_dir
        / "qc/prefiltered_counts.stringent.hist.txt",
        counts_histogram_relaxed=dynamic_unfiltered_seacr_dir
        / "qc/prefiltered_counts.relaxed.hist.txt",
    shell:
        """
        > {output.prefiltered_counts_stringent}  
        > {output.counts_histogram_stringent}
        for input_file in {input.stringent_peaks}; do
            echo -n "$(basename $input_file) " >> {output.prefiltered_counts_stringent}
            wc -l < $input_file >> {output.prefiltered_counts_stringent}
            hist_name=$(basename {output.counts_histogram_stringent})
            echo "$(basename $input_file) " >> {output.counts_histogram_stringent}
            awk '$5 ~ /^[0-9]+$/ {{print int($5/5)*5}}' $input_file | sort -n | uniq -c >> {output.counts_histogram_stringent}
            echo "" >> {output.counts_histogram_stringent}
        done
        
        > {output.prefiltered_counts_relaxed}  
        > {output.counts_histogram_relaxed}
        for input_file in {input.relaxed_peaks}; do
            echo -n "$(basename $input_file) " >> {output.prefiltered_counts_relaxed}
            wc -l < $input_file >> {output.prefiltered_counts_relaxed}
            hist_name=$(basename {output.counts_histogram_relaxed})
            echo "$(basename $input_file) " >>  {output.counts_histogram_relaxed}
            awk '$5 ~ /^[0-9]+$/ {{print int($5/5)*5}}' $input_file | sort -n | uniq -c >> {output.counts_histogram_relaxed}
            echo "" >> {output.counts_histogram_relaxed}
        done
        """


rule filter_peaks:
    input:
        stringent_peaks=dynamic_unfiltered_seacr_dir
        / "stringent/{peak_id}.stringent.bed",
        relaxed_peaks=dynamic_unfiltered_seacr_dir / "relaxed/{peak_id}.relaxed.bed",
    output:
        stringent_filtered_peaks=dynamic_filtered_seacr_dir
        / "stringent/{peak_id}.filtered.stringent.bed",
        relaxed_filtered_peaks=dynamic_filtered_seacr_dir
        / "relaxed/{peak_id}.filtered.relaxed.bed",
    params:
        threshold=config["peaks_call"]["seacr"]["filter_threshold"],
    shell:
        """
        awk '$5 >= {params.threshold}' {input.stringent_peaks} > {output.stringent_filtered_peaks}
        awk '$5 >= {params.threshold}' {input.relaxed_peaks} > {output.relaxed_filtered_peaks}
        """


rule post_filtere_qc:
    input:
        filtered_stringent_peaks=lambda wc: [
            dynamic_filtered_seacr_dir / f"stringent/{pid}.filtered.stringent.bed"
            for pid in peak_ids
        ],
        filtered_relaxed_peaks=lambda wc: [
            dynamic_filtered_seacr_dir / f"relaxed/{pid}.filtered.relaxed.bed"
            for pid in peak_ids
        ],
    output:
        filtered_counts_stringent=dynamic_filtered_seacr_dir
        / "qc/filtered_counts.stringent.txt",
        filtered_counts_relaxed=dynamic_filtered_seacr_dir
        / "qc/filtered_counts.relaxed.txt",
        filtered_counts_histogram_stringent=dynamic_filtered_seacr_dir
        / "qc/filtered_counts.stringent.hist.txt",
        filtered_counts_histogram_relaxed=dynamic_filtered_seacr_dir
        / "qc/filtered_counts.relaxed.hist.txt",
    shell:
        """
        > {output.filtered_counts_stringent}  
        > {output.filtered_counts_histogram_stringent}
        for input_file in {input.filtered_stringent_peaks}; do
            echo -n "$(basename $input_file) " >> {output.filtered_counts_stringent}
            wc -l < $input_file >> {output.filtered_counts_stringent}
            hist_name=$(basename {output.filtered_counts_histogram_stringent})
            echo "$(basename $input_file) " >> {output.filtered_counts_histogram_stringent}
            awk '$5 ~ /^[0-9]+$/ {{print int($5/5)*5}}' $input_file | sort -n | uniq -c >> {output.filtered_counts_histogram_stringent}
            echo "" >> {output.filtered_counts_histogram_stringent}
        done
        
        > {output.filtered_counts_relaxed}  
        > {output.filtered_counts_histogram_relaxed}
        for input_file in {input.filtered_relaxed_peaks}; do
            echo -n "$(basename $input_file) " >> {output.filtered_counts_relaxed}
            wc -l < $input_file >> {output.filtered_counts_relaxed}
            hist_name=$(basename {output.filtered_counts_histogram_relaxed})
            echo "$(basename $input_file) " >>  {output.filtered_counts_histogram_relaxed}
            awk '$5 ~ /^[0-9]+$/ {{print int($5/5)*5}}' $input_file | sort -n | uniq -c >> {output.filtered_counts_histogram_relaxed}
            echo "" >> {output.filtered_counts_histogram_relaxed}
        done
        """


rule intersect:
    input:
        all_rep1_stringent=lambda wc: f"{dynamic_filtered_seacr_dir}/stringent/{get_reps(wc.gen_abd_id)[0]}.filtered.stringent.bed",
        all_rep2_stringent=lambda wc: f"{dynamic_filtered_seacr_dir}/stringent/{get_reps(wc.gen_abd_id)[1]}.filtered.stringent.bed",
        all_rep1_relaxed=lambda wc: f"{dynamic_filtered_seacr_dir}/relaxed/{get_reps(wc.gen_abd_id)[0]}.filtered.relaxed.bed",
        all_rep2_relaxed=lambda wc: f"{dynamic_filtered_seacr_dir}/relaxed/{get_reps(wc.gen_abd_id)[1]}.filtered.relaxed.bed",
    output:
        intersect_stringent=dynamic_intersected_peaks
        / "{gen_abd_id}_rep1_v_rep2.stringent.bed",
        intersect_relaxed=dynamic_intersected_peaks
        / "{gen_abd_id}_rep1_v_rep2.relaxed.bed",
    conda:
        "../../envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -a {input.all_rep1_stringent} -b {input.all_rep2_stringent} -wa -wb > {output.intersect_stringent}
        bedtools intersect -a {input.all_rep1_relaxed} -b {input.all_rep2_relaxed} -wa -wb > {output.intersect_relaxed}
        """


# Min, max
rule greedy_consensus:
    input:
        intersect_stringent=rules.intersect.output.intersect_stringent,
        intersect_relaxed=rules.intersect.output.intersect_relaxed,
    output:
        stringent_max_cons=dynamic_greedy_consensus_peaks
        / "{gen_abd_id}_consensus.greedy.stringent.bed",
        relaxed_max_cons=dynamic_greedy_consensus_peaks
        / "{gen_abd_id}_consensus.greedy.relaxed.bed",
    shell:
        """
        awk 'BEGIN{{OFS="\t"}} {{print $1, ($2<$8?$2:$8), ($3>$9?$3:$9), $4, $5, $6}}' {input.intersect_stringent} > {output.stringent_max_cons}
        awk 'BEGIN{{OFS="\t"}} {{print $1, ($2<$8?$2:$8), ($3>$9?$3:$9), $4, $5, $6}}' {input.intersect_relaxed} > {output.relaxed_max_cons}
        """


# max, min
rule conservative_consensus:
    input:
        intersect_stringent=rules.intersect.output.intersect_stringent,
        intersect_relaxed=rules.intersect.output.intersect_relaxed,
    output:
        stringent_min_cons=dynamic_conservative_consensus_peaks
        / "{gen_abd_id}_consensus.conservative.stringent.bed",
        relaxed_min_cons=dynamic_conservative_consensus_peaks
        / "{gen_abd_id}_consensus.conservative.relaxed.bed",
    shell:
        """
        awk 'BEGIN{{OFS="\t"}} {{start=($2>$8?$2:$8); end=($3<$9?$3:$9); if(start<end) print $1, start, end, $4, $5, $6}}' {input.intersect_stringent} > {output.stringent_min_cons}
        awk 'BEGIN{{OFS="\t"}} {{start=($2>$8?$2:$8); end=($3<$9?$3:$9); if(start<end) print $1, start, end, $4, $5, $6}}' {input.intersect_relaxed} > {output.relaxed_min_cons}
        """

# Reciprocal overlap
rule reciprocal_overlap:
    input:
        all_rep1_stringent=lambda wc: f"{dynamic_filtered_seacr_dir}/stringent/{get_reps(wc.gen_abd_id)[0]}.filtered.stringent.bed",
        all_rep2_stringent=lambda wc: f"{dynamic_filtered_seacr_dir}/stringent/{get_reps(wc.gen_abd_id)[1]}.filtered.stringent.bed",
        all_rep1_relaxed=lambda wc: f"{dynamic_filtered_seacr_dir}/relaxed/{get_reps(wc.gen_abd_id)[0]}.filtered.relaxed.bed",
        all_rep2_relaxed=lambda wc: f"{dynamic_filtered_seacr_dir}/relaxed/{get_reps(wc.gen_abd_id)[1]}.filtered.relaxed.bed",
    output:
        overlap_stringent=dynamic_reciprocal_consensus_peaks
        / "{gen_abd_id}_consensus.reciprocal.stringent.bed",
        overlap_relaxed=dynamic_reciprocal_consensus_peaks
        / "{gen_abd_id}_consensus.reciprocal.relaxed.bed",
    conda:
        "../../envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -a {input.all_rep1_stringent} -b {input.all_rep2_stringent} -wa -r -f 0.5 > {output.overlap_stringent}
        bedtools intersect -a {input.all_rep1_relaxed} -b {input.all_rep2_relaxed} -wa -r -f 0.5 > {output.overlap_relaxed}
        """
