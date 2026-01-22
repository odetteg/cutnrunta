from constants.dirs_files import *

if config["peaks_call"]["seacr"]["use_seacr"]:
    if config["peaks_call"]["seacr"]["use_control"]:
        if config["chipseeker_annotation"]["annotate_raw"]:

            rule ChIPseeker_seacr_relaxed_raw:
                input:
                    bed=rules.seacr_relaxed.output._peaks,
                    gtf=GTF,
                output:
                    txt=dynamic_chipseeker_seacr_raw_dir / "{sample_id}.relaxed.txt",
                conda:
                    "../../envs/ChIPseeker.yaml"
                threads: config["resources"]["chipseeker"]["cpu"]
                log:
                    std_out="logs/chipseeker/{sample_id}.seacr.relaxed.out.log",
                    std_err="logs/chipseeker/{sample_id}.seacr.relaxed.err.log",
                params:
                    upstream=config["chipseeker_annotation"]["upstream"],
                    downstream=config["chipseeker_annotation"]["downstream"],
                resources:
                    mem_mb=config["resources"]["chipseeker"]["mem_mb"],
                    time=config["resources"]["chipseeker"]["time"],
                script:
                    "../scripts/ChIPseeker_seacr.R"

            rule ChIPseeker_seacr_stringent_raw:
                input:
                    bed=rules.seacr_stringent.output._peaks,
                    gtf=GTF,
                output:
                    txt=dynamic_chipseeker_seacr_raw_dir / "{sample_id}.stringent.txt",
                conda:
                    "../../envs/ChIPseeker.yaml"
                threads: config["resources"]["chipseeker"]["cpu"]
                log:
                    std_out="logs/chipseeker/{sample_id}.seacr.stringent.out.log",
                    std_err="logs/chipseeker/{sample_id}.seacr.stringent.err.log",
                params:
                    upstream=config["chipseeker_annotation"]["upstream"],
                    downstream=config["chipseeker_annotation"]["downstream"],
                resources:
                    mem_mb=config["resources"]["chipseeker"]["mem_mb"],
                    time=config["resources"]["chipseeker"]["time"],
                script:
                    "../scripts/ChIPseeker_seacr.R"

    if not config["peaks_call"]["seacr"]["use_control"]:
        if config["chipseeker_annotation"]["annotate_raw_nocontrol"]:

            rule ChIPseeker_seacr_relaxed_raw_noctrl:
                input:
                    bed=rules.seacr_relaxed_noctrl.output._peaks,
                    gtf=GTF,
                output:
                    txt=dynamic_chipseeker_seacr_raw_dir / "{sample_id}.relaxed.txt",
                conda:
                    "../../envs/ChIPseeker.yaml"
                threads: config["resources"]["chipseeker"]["cpu"]
                log:
                    std_out="logs/chipseeker/{sample_id}.seacr.relaxed.noctrl.out.log",
                    std_err="logs/chipseeker/{sample_id}.seacr.relaxed.noctrl.err.log",
                params:
                    upstream=config["chipseeker_annotation"]["upstream"],
                    downstream=config["chipseeker_annotation"]["downstream"],
                resources:
                    mem_mb=config["resources"]["chipseeker"]["mem_mb"],
                    time=config["resources"]["chipseeker"]["time"],
                script:
                    "../scripts/ChIPseeker_seacr.R"

            rule ChIPseeker_seacr_stringent_raw_noctrl:
                input:
                    bed=rules.seacr_stringent_noctrl.output._peaks,
                    gtf=GTF,
                output:
                    txt=dynamic_chipseeker_seacr_raw_dir / "{sample_id}.stringent.txt",
                conda:
                    "../../envs/ChIPseeker.yaml"
                threads: config["resources"]["chipseeker"]["cpu"]
                log:
                    std_out="logs/chipseeker/{sample_id}.seacr.stringent.noctrl.out.log",
                    std_err="logs/chipseeker/{sample_id}.seacr.stringent.noctrl.err.log",
                params:
                    upstream=config["chipseeker_annotation"]["upstream"],
                    downstream=config["chipseeker_annotation"]["downstream"],
                resources:
                    mem_mb=config["resources"]["chipseeker"]["mem_mb"],
                    time=config["resources"]["chipseeker"]["time"],
                script:
                    "../scripts/ChIPseeker_seacr.R"
