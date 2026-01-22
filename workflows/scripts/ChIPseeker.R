# Load libraries

library(GenomicFeatures)
library(ChIPseeker)
library(tidyverse)
library(viridisLite)

bed.file <- snakemake@input[["bed"]]
gtf.file <- snakemake@input[["gtf"]]
output_prefix <- snakemake@params[["output_prefix"]]

txdb <- makeTxDbFromGFF(gtf.file, format="gtf")

peakAnno <- annotatePeak(bed.file, tssRegion=c(-3000, 3000),
                        TxDb=txdb)