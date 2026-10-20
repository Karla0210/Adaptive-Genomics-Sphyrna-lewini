
# 07. Redundancy analysis (RDA)

Genotype data were prepared for redundancy analysis (RDA) to identify SNPs potentially associated with environmental variables. The filtered VCF dataset was converted using VCFtools and PLINK v1.9 for subsequent analyses in R.

## Genotype data preparation

```bash
# 1. Convert VCF to TPED/TFAM
vcftools \
  --vcf populations_08.vcf.recode.vcf \
  --out hammerhead_tped \
  --plink-tped

# 2. Convert TPED/TFAM to PLINK genotype format
plink1.9 \
  --tped hammerhead_tped.tped \
  --tfam hammerhead_tped.tfam \
  --recode \
  --out RDA_gen \
  --allow-extra-chr
```

