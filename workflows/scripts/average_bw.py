from snakemake.shell import shell

bw_files = snakemake.input.bw_files
bw_out = snakemake.output.wg_
log = snakemake.log

shell("wiggletools write {bw_out} mean {bw_files} > {log} 2>&1")
