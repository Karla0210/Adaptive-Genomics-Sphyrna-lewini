#!/bin/bash

# Missing data per individual
vcftools \
--vcf populations.snps.vcf \
--out SLEW_indiv \
--missing-indv

# Missing data per locus
vcftools \
--vcf populations.snps.vcf \
--out SLEW_site \
--missing-site

# Remove individuals
vcftools \
--vcf populations.snps.vcf \
--remove individuals_to_remove.txt \
--recode \
--out populations_filtered

# Missing data after filtering
vcftools \
--vcf populations_filtered.recode.vcf \
--out SLEW_indiv_filtered \
--missing-indv

vcftools \
--vcf populations_filtered.recode.vcf \
--out SLEW_site_filtered \
--missing-site

# Final filtering
vcftools \
--vcf populations_filtered.recode.vcf \
--max-missing 0.8 \
--recode \
--out populations_08
