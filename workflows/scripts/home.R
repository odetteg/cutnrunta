library()
tt_bam = snakemake@[input]["tt_bam"]
ctrl_bam = snakemake@[input]["ctrl_bam"]
REF = snakemake@[input]["REF"]
peaks = snakemake@[output]["peaks"]
style = snakemake@[params]["style"]
extra = snakemake@[params]["extra"]
ttag_dir = snakemake@[params]["ttag_dir"]
cttag_dir = snakemake@[params]["cttag_dir"]

# Create tag directory for input:

makeTagDirectory ttag_dir tt_bam -genome REF
# Create tag directory for control:
makeTagDirectory cttag_dir ctrl_bam -genome REF
# Call peaks using HOMER:

findpeaks ttag_dir -style style -i cttag_dir  -o peaks