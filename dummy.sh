#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="data_cp/data/sequenced_2025_11_13"

mkdir -p "${BASE_DIR}"

files=(
  "BL3-H3K27ac-1_S4_L001_R2_001.fastq.gz"
  "TBL3-H3K27me3-1_S9_L001_R1_001.fastq.gz"
  "BL3-H3K27ac-2_S5_L001_R1_001.fastq.gz"
  "BL3-H3K27me3-1_S2_L001_R1_001.fastq.gz"
  "BL3-H3K27me3-1_S2_L001_R2_001.fastq.gz"
  "TBL3-H3K27me3-1_S9_L001_R2_001.fastq.gz"
  "BL3-H3K27me3-2_S3_L001_R1_001.fastq.gz"
  "BL3-H2K119Ub-1_S6_L001_R1_001.fastq.gz"
  "TBL3-H3K27ac-1_S11_L001_R2_001.fastq.gz"
  "TBL3-H3K27ac-2_S12_L001_R2_001.fastq.gz"
  "TBL3-H2K119Ub-1_S13_L001_R2_001.fastq.gz"
  "TBL3-H2K119Ub-2_S14_L001_R1_001.fastq.gz"
  "TBL3-H3K27me3-2_S10_L001_R2_001.fastq.gz"
  "BL3-H3K27ac-2_S5_L001_R2_001.fastq.gz"
  "TBL3-H3K27ac-2_S12_L001_R1_001.fastq.gz"
  "BL3-H2K119Ub-2_S7_L001_R2_001.fastq.gz"
  "TBL3-H3K27me3-2_S10_L001_R1_001.fastq.gz"
  "BL3-H2K119Ub-1_S6_L001_R2_001.fastq.gz"
  "BL3-H3K27ac-1_S4_L001_R1_001.fastq.gz"
  "TBL3-H3K27ac-1_S11_L001_R1_001.fastq.gz"
  "BL3-H3K27me3-2_S3_L001_R2_001.fastq.gz"
  "TBL3-H2K119Ub-1_S13_L001_R1_001.fastq.gz"
  "TBL3-IgG_S8_L001_R1_001.fastq.gz"
  "BL3-IgG_S1_L001_R2_001.fastq.gz"
  "TBL3-IgG_S8_L001_R2_001.fastq.gz"
  "BL3-H2K119Ub-2_S7_L001_R1_001.fastq.gz"
  "TBL3-H2K119Ub-2_S14_L001_R2_001.fastq.gz"
  "BL3-IgG_S1_L001_R1_001.fastq.gz"
)

for f in "${files[@]}"; do
  : > "${BASE_DIR}/${f}"
done

echo "Dummy FASTQ files created in ${BASE_DIR}"
