library("DiffBind")
library("tidyverse")



safe_load <- function(sample_sheet_path, name) {
    cat(sprintf("\n=== Loading sample sheet for %s ===\n", name))
    cat(sprintf("Path: %s\n", sample_sheet_path))

    if (!file.exists(sample_sheet_path)) {
        stop(sprintf("Sample sheet for %s does not exist at path: %s", name, sample_sheet_path))
    }
    
    # Check if all peak files exist
    ss <- read.csv(sample_sheet_path)
    cat(sprintf("Found %d samples in sample sheet\n", nrow(ss)))

    missing_peaks <- ss$Peaks[!file.exists(ss$Peaks)]
    if (length(missing_peaks) > 0) {
        stop(sprintf("The following peak files do not exist:\n%s", paste(missing_peaks, collapse="\n")))
    }

    # Check if peak files are empty
    empty_peaks <- ss$Peaks[file.size(ss$Peaks) == 0]
    if (length(empty_peaks) > 0) {
        stop(sprintf("The following peak files are empty:\n%s", paste(empty_peaks, collapse="\n")))
    }

    cat("All peak files exist and are non-empty.\n")
    
    # Load and return DBA object
    dba_obj <- dba(sampleSheet = sample_sheet_path)
    cat(sprintf("Successfully loaded %s DBA object\n", name))
    
    return(dba_obj)
}

# Load sample sheets and create DBA objects
relaxed.dbObj <- safe_load(snakemake@input[["relaxed_sample_sheet"]], "relaxed")
stringent.dbObj <- safe_load(snakemake@input[["stringent_sample_sheet"]], "stringent")

cat("\n=== Counting reads ===\n")
# Count reads
relaxed.counts <- dba.count(relaxed.dbObj, bUseSummarizeOverlaps=TRUE)
stringent.counts <- dba.count(stringent.dbObj, bUseSummarizeOverlaps=TRUE)
cat("Read counting complete\n")

cat("\n=== Generating PCA plots ===\n")
# Plot PCA 
pdf(snakemake@output[["relaxed_pca"]], width=10, height=8)
dba.plotPCA(
  relaxed.counts,
  label = DBA_CONDITION,
  labelSize = 0.8
)
dev.off()

pdf(snakemake@output[["stringent_pca"]], width=10, height=8)
dba.plotPCA(
  stringent.counts,  
  label = DBA_CONDITION,
  labelSize = 0.8
)
dev.off()

cat("PCA plots complete\n")

cat("\n=== Saving count matrices ===\n")
# Save count matrices
relaxed_matrix <- dba.peakset(relaxed.counts, bRetrieve=TRUE)
write.csv(relaxed_matrix, snakemake@output[["relaxed_counts_matrix_csv"]], row.names=TRUE)
stringent_matrix <- dba.peakset(stringent.counts, bRetrieve=TRUE)
write.csv(stringent_matrix, snakemake@output[["stringent_counts_matrix_csv"]], row.names=TRUE)
cat("Count matrices saved\n")

cat("\n=== Subsetting by antibody ===\n")
# Subsetting the H3K27ac samples
relaxed.ac <- dba(relaxed.dbObj, mask=relaxed.dbObj$masks$H3K27ac)
stringent.ac <- dba(stringent.dbObj, mask=stringent.dbObj$masks$H3K27ac)

# Subset H3K27me3 samples
relaxed.me3 <- dba(relaxed.dbObj, mask=relaxed.dbObj$masks$H3K27me3)
stringent.me3 <- dba(stringent.dbObj, mask=stringent.dbObj$masks$H3K27me3)

# Subset H2K119Ub samples
relaxed.ub <- dba(relaxed.dbObj, mask=relaxed.dbObj$masks$H2K119Ub)
stringent.ub <- dba(stringent.dbObj, mask=stringent.dbObj$masks$H2K119Ub)
cat("Subsetting complete\n")

cat("\n=== Counting reads per antibody ===\n")
# Getting counts for all
relaxed.ac.counts <- dba.count(relaxed.ac, bUseSummarizeOverlaps=TRUE)
stringent.ac.counts <- dba.count(stringent.ac, bUseSummarizeOverlaps=TRUE)
relaxed.me3.counts <- dba.count(relaxed.me3, bUseSummarizeOverlaps=TRUE)
stringent.me3.counts <- dba.count(stringent.me3, bUseSummarizeOverlaps=TRUE)
relaxed.ub.counts <- dba.count(relaxed.ub, bUseSummarizeOverlaps=TRUE)
stringent.ub.counts <- dba.count(stringent.ub, bUseSummarizeOverlaps=TRUE)
cat("Per-antibody counting complete\n")

cat("\n=== Setting up contrasts ===\n")
# Contrast the samples per antibody
relaxed.ac.contrast <- dba.contrast(relaxed.ac.counts, categories=DBA_CONDITION, minMembers=2)
stringent.ac.contrast <- dba.contrast(stringent.ac.counts, categories=DBA_CONDITION, minMembers=2)
relaxed.me3.contrast <- dba.contrast(relaxed.me3.counts, categories=DBA_CONDITION, minMembers=2)
stringent.me3.contrast <- dba.contrast(stringent.me3.counts, categories=DBA_CONDITION, minMembers=2)
relaxed.ub.contrast <- dba.contrast(relaxed.ub.counts, categories=DBA_CONDITION, minMembers=2)
stringent.ub.contrast <- dba.contrast(stringent.ub.counts, categories=DBA_CONDITION, minMembers=2)
cat("Contrasts set up\n")

cat("\n=== Performing differential analysis ===\n")
# Perform differential analysis
relaxed.ac.diff <- dba.analyze(relaxed.ac.contrast)
stringent.ac.diff <- dba.analyze(stringent.ac.contrast)
relaxed.me3.diff <- dba.analyze(relaxed.me3.contrast)
stringent.me3.diff <- dba.analyze(stringent.me3.contrast)
relaxed.ub.diff <- dba.analyze(relaxed.ub.contrast)
stringent.ub.diff <- dba.analyze(stringent.ub.contrast)
cat("Differential analysis complete\n")


saveRDS(relaxed.ac.diff, file=snakemake@output[["relaxed_h3k27ac_diffbind_rds"]])
saveRDS(stringent.ac.diff, file=snakemake@output[["stringent_h3k27ac_diffbind_rds"]])
saveRDS(relaxed.me3.diff, file=snakemake@output[["relaxed_h3k27me3_diffbind_rds"]])
saveRDS(stringent.me3.diff, file=snakemake@output[["stringent_h3k27me3_diffbind_rds"]])
saveRDS(relaxed.ub.diff, file=snakemake@output[["relaxed_h2k119ub_diffbind_rds"]])
saveRDS(stringent.ub.diff, file=snakemake@output[["stringent_h2k119ub_diffbind_rds"]])

saveRDS(relaxed.ac, file=snakemake@output[["relaxed_ac_dba"]])
saveRDS(stringent.ac, file=snakemake@output[["stringent_ac_dba"]])
saveRDS(relaxed.me3, file=snakemake@output[["relaxed_me3_dba"]])
saveRDS(stringent.me3, file=snakemake@output[["stringent_me3_dba"]])
saveRDS(relaxed.ub, file=snakemake@output[["relaxed_ub_dba"]])
saveRDS(stringent.ub, file=snakemake@output[["stringent_ub_dba"]])

saveRDS(relaxed.ac.counts, file=snakemake@output[["relaxed_ac_counts_dba"]])
saveRDS(stringent.ac.counts, file=snakemake@output[["stringent_ac_counts_dba"]])
saveRDS(relaxed.me3.counts, file=snakemake@output[["relaxed_me3_counts_dba"]])
saveRDS(stringent.me3.counts, file=snakemake@output[["stringent_me3_counts_dba"]])
saveRDS(relaxed.ub.counts, file=snakemake@output[["relaxed_ub_counts_dba"]])
saveRDS(stringent.ub.counts, file=snakemake@output[["stringent_ub_counts_dba"]])


# Generating the map plots. MAP plots generaly tells us at what point do we se differentia binding

# pdf(snakemake@output[["relaxed_ac_ma"]], width=10, height=8)
# dba.plotMA(relaxed.ac.diff, contrast=1, method=DBA_DESEQ2)
# dev.off()
# pdf(snakemake@output[["stringent_ac_ma"]], width=10, height=8)
# dba.plotMA(stringent.ac.diff, contrast=1, method=DBA_DESEQ2)
# dev.off()   
# pdf(snakemake@output[["relaxed_me3_ma"]], width=10, height=8)
# dba.plotMA(relaxed.me3.diff, contrast=1, method=DBA_DESEQ2)
# dev.off()
# pdf(snakemake@output[["stringent_me3_ma"]], width=10, height=8)
# dba.plotMA(stringent.me3.diff, contrast=1, method=DBA_DESEQ2)
# dev.off()
# pdf(snakemake@output[["relaxed_ub_ma"]], width=10, height=8)
# dba.plotMA(relaxed.ub.diff, contrast=1, method=DBA_DESEQ2)
# dev.off()
# pdf(snakemake@output[["stringent_ub_ma"]], width=10, height=8)
# dba.plotMA(stringent.ub.diff, contrast=1, method=DBA_DESEQ2)
# dev.off()

cat("\n=== Generating reports and saving in individual CSV files ===\n")

relaxed.ac.report <- dba.report(relaxed.ac.diff, th=1)
stringent.ac.report <- dba.report(stringent.ac.diff, th=1)
relaxed.me3.report <- dba.report(relaxed.me3.diff, th=1)
stringent.me3.report <- dba.report(stringent.me3.diff, th=1)
relaxed.ub.report <- dba.report(relaxed.ub.diff, th=1)
stringent.ub.report <- dba.report(stringent.ub.diff, th=1)

# Write to CSV
write.csv(as.data.frame(relaxed.ac.report), file=snakemake@output[["relaxed_h3k27ac_results"]], row.names=FALSE)
write.csv(as.data.frame(stringent.ac.report), file=snakemake@output[["stringent_h3k27ac_results"]], row.names=FALSE)
write.csv(as.data.frame(relaxed.me3.report), file=snakemake@output[["relaxed_h3k27me3_results"]], row.names=FALSE)
write.csv(as.data.frame(stringent.me3.report), file=snakemake@output[["stringent_h3k27me3_results"]], row.names=FALSE)
write.csv(as.data.frame(relaxed.ub.report), file=snakemake@output[["relaxed_h2k119ub_results"]], row.names=FALSE)
write.csv(as.data.frame(stringent.ub.report), file=snakemake@output[["stringent_h2k119ub_results"]], row.names=FALSE)

cat("\n=== DiffBind analysis complete! ===\n")

cat("\n=== Checking differential binding results ===\n")


# PCA for each antibody
pdf(snakemake@output[["relaxed_h3k27ac_pca"]], width=10, height=8)
dba.plotPCA(
  relaxed.ac.counts,  
  label = DBA_CONDITION,
  labelSize = 0.8
)
dev.off()

pdf(snakemake@output[["stringent_h3k27ac_pca"]], width=10, height=8)
dba.plotPCA(
  stringent.ac.counts,  
  label = DBA_CONDITION,
  labelSize = 0.8
)
dev.off()

pdf(snakemake@output[["relaxed_h3k27me3_pca"]], width=10, height=8)
dba.plotPCA(
  relaxed.me3.counts,  
  label = DBA_CONDITION,
  labelSize = 0.8
)
dev.off()

pdf(snakemake@output[["stringent_h3k27me3_pca"]], width=10, height=8)
dba.plotPCA(
  stringent.me3.counts,  
  label = DBA_CONDITION,
  labelSize = 0.8
)
dev.off()

pdf(snakemake@output[["relaxed_h2k119ub_pca"]], width=10, height=8)
dba.plotPCA(
  relaxed.ub.counts,  
  label = DBA_CONDITION,
  labelSize = 0.8
)
dev.off()

pdf(snakemake@output[["stringent_h2k119ub_pca"]], width=10, height=8)
dba.plotPCA(
  stringent.ub.counts,  
  label = DBA_CONDITION,
  labelSize = 0.8
)
dev.off()


# Getting the MA plots









cat("\nRelaxed H2K119Ub:\n")
cat(sprintf("  Total peaks analyzed: %d\n", nrow(dba.report(relaxed.ub.diff, th=1))))
cat(sprintf("  Significant peaks (FDR<0.05): %d\n", nrow(dba.report(relaxed.ub.diff))))
cat(sprintf("  Peaks at FDR<0.1: %d\n", nrow(dba.report(relaxed.ub.diff, th=0.1))))

cat("\nStringent H2K119Ub:\n")
cat(sprintf("  Total peaks analyzed: %d\n", nrow(dba.report(stringent.ub.diff, th=1))))
cat(sprintf("  Significant peaks (FDR<0.05): %d\n", nrow(dba.report(stringent.ub.diff))))
cat(sprintf("  Peaks at FDR<0.1: %d\n", nrow(dba.report(stringent.ub.diff, th=0.1))))

cat("\nRelaxed H3K27ac:\n")
cat(sprintf("  Total peaks analyzed: %d\n", nrow(dba.report(relaxed.ac.diff, th=1))))
cat(sprintf("  Significant peaks (FDR<0.05): %d\n", nrow(dba.report(relaxed.ac.diff))))

cat("\nStringent H3K27ac:\n")
cat(sprintf("  Total peaks analyzed: %d\n", nrow(dba.report(stringent.ac.diff, th=1))))
cat(sprintf("  Significant peaks (FDR<0.05): %d\n", nrow(dba.report(stringent.ac.diff))))

cat("\nRelaxed H3K27me3:\n")
cat(sprintf("  Total peaks analyzed: %d\n", nrow(dba.report(relaxed.me3.diff, th=1))))
cat(sprintf("  Significant peaks (FDR<0.05): %d\n", nrow(dba.report(relaxed.me3.diff))))

cat("\nStringent H3K27me3:\n")
cat(sprintf("  Total peaks analyzed: %d\n", nrow(dba.report(stringent.me3.diff, th=1))))
cat(sprintf("  Significant peaks (FDR<0.05): %d\n", nrow(dba.report(stringent.me3.diff))))

cat("\n=== Creating differential binding summary report ===\n")

# Create summary report
summary_file <- file.path(snakemake@output[["diffbind_summary_report"]])

sink(summary_file)
cat("=======================================================\n")
cat("DiffBind Differential Binding Analysis Summary\n")
cat("=======================================================\n")
cat(sprintf("Analysis date: %s\n\n", Sys.time()))

# H2K119Ub results
cat("-------------------------------------------------------\n")
cat("H2K119Ub Results\n")
cat("-------------------------------------------------------\n")
cat("\nRelaxed:\n")
cat(sprintf("  Total peaks analyzed: %d\n", nrow(dba.report(relaxed.ub.diff, th=1))))
cat(sprintf("  Significant peaks (FDR<0.05): %d\n", nrow(dba.report(relaxed.ub.diff))))
cat(sprintf("  Peaks at FDR<0.1: %d\n", nrow(dba.report(relaxed.ub.diff, th=0.1))))

cat("\nStringent:\n")
cat(sprintf("  Total peaks analyzed: %d\n", nrow(dba.report(stringent.ub.diff, th=1))))
cat(sprintf("  Significant peaks (FDR<0.05): %d\n", nrow(dba.report(stringent.ub.diff))))
cat(sprintf("  Peaks at FDR<0.1: %d\n", nrow(dba.report(stringent.ub.diff, th=0.1))))

# H3K27ac results
cat("\n-------------------------------------------------------\n")
cat("H3K27ac Results\n")
cat("-------------------------------------------------------\n")
cat("\nRelaxed:\n")
cat(sprintf("  Total peaks analyzed: %d\n", nrow(dba.report(relaxed.ac.diff, th=1))))
cat(sprintf("  Significant peaks (FDR<0.05): %d\n", nrow(dba.report(relaxed.ac.diff))))
cat(sprintf("  Peaks at FDR<0.1: %d\n", nrow(dba.report(relaxed.ac.diff, th=0.1))))

cat("\nStringent:\n")
cat(sprintf("  Total peaks analyzed: %d\n", nrow(dba.report(stringent.ac.diff, th=1))))
cat(sprintf("  Significant peaks (FDR<0.05): %d\n", nrow(dba.report(stringent.ac.diff))))
cat(sprintf("  Peaks at FDR<0.1: %d\n", nrow(dba.report(stringent.ac.diff, th=0.1))))

# H3K27me3 results
cat("\n-------------------------------------------------------\n")
cat("H3K27me3 Results\n")
cat("-------------------------------------------------------\n")
cat("\nRelaxed:\n")
cat(sprintf("  Total peaks analyzed: %d\n", nrow(dba.report(relaxed.me3.diff, th=1))))
cat(sprintf("  Significant peaks (FDR<0.05): %d\n", nrow(dba.report(relaxed.me3.diff))))
cat(sprintf("  Peaks at FDR<0.1: %d\n", nrow(dba.report(relaxed.me3.diff, th=0.1))))

cat("\nStringent:\n")
cat(sprintf("  Total peaks analyzed: %d\n", nrow(dba.report(stringent.me3.diff, th=1))))
cat(sprintf("  Significant peaks (FDR<0.05): %d\n", nrow(dba.report(stringent.me3.diff))))
cat(sprintf("  Peaks at FDR<0.1: %d\n", nrow(dba.report(stringent.me3.diff, th=0.1))))

cat("\n=======================================================\n")
sink()

cat(sprintf("Summary report written to: %s\n", summary_file))