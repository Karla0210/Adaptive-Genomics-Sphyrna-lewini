
# 04. VCF filtering

Missing data per individual and per locus were assessed using VCFtools. Selected individuals were removed, and missing data were reassessed before applying an 80% minimum locus completeness threshold (`--max-missing 0.8`). The final dataset contained 9,558 SNPs for downstream genomic analyses.

## VCFtools commands

```bash
#!/bin/bash

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
```

