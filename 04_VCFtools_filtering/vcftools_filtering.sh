#!/bin/bash

module purge
module load vcftools/0.1.16

# 1. Missing data per individual
vcftools \
  --vcf populations.snps.vcf \
  --out SLEW_indiv \
  --missing-indv

# 2. Missing data per locus
vcftools \
  --vcf populations.snps.vcf \
  --out SLEW_site \
  --missing-site

# 3. Remove selected individuals
vcftools \
  --vcf populations.snps.vcf \
  --remove individuals_to_remove.txt \
  --recode \
  --out populations_filtered

# 4. Missing data after individual removal
vcftools \
  --vcf populations_filtered.recode.vcf \
  --out SLEW_indiv_filtered \
  --missing-indv

vcftools \
  --vcf populations_filtered.recode.vcf \
  --out SLEW_site_filtered \
  --missing-site

# 5. Final filtering (80% completeness)
vcftools \
  --vcf populations_filtered.recode.vcf \
  --max-missing 0.8 \
  --recode \
  --out populations_08


