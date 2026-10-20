
# 04. VCF filtering

Missing data per individual and per locus were evaluated using VCFtools. A minimum locus completeness threshold of 80% (`--max-missing 0.8`) was applied to generate the final SNP dataset (9,558 SNPs) used for downstream genomic analyses.

## Missing data per individual

```bash
vcftools \
  --vcf populations.snps.vcf \
  --missing-indv \
  --out SLEW_indiv
```

## Missing data per locus

```bash
vcftools \
  --vcf populations.snps.vcf \
  --missing-site \
  --out SLEW_site
```

## Missing-data filtering

```bash
vcftools \
  --vcf populations_filtered.recode.vcf \
  --max-missing 0.8 \
  --recode \
  --out populations_08
```
