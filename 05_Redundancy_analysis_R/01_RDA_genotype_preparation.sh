
#!/bin/bash

# VCF to PLINK format conversion - Sphyrna lewini
# VCFtools v0.1.16 and PLINK v1.9

# Load required modules
module purge
module load vcftools/0.1.16
module load plink/1.9

# Input and output directories
VCF="results/vcf_filtering/populations_08.vcf.recode.vcf"
OUT_DIR="results/plink_conversion"

# Create output directory
mkdir -p "${OUT_DIR}"

# 1. Convert VCF to TPED/TFAM
vcftools \
  --vcf "${VCF}" \
  --out "${OUT_DIR}/hammerhead_tped" \
  --plink-tped

# 2. Convert TPED/TFAM to PLINK PED/MAP format
plink1.9 \
  --tped "${OUT_DIR}/hammerhead_tped.tped" \
  --tfam "${OUT_DIR}/hammerhead_tped.tfam" \
  --recode \
  --out "${OUT_DIR}/RDA_gen" \
  --allow-extra-chr
