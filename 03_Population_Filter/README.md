
# 03. Population filtering

Loci were filtered using the `populations` module in Stacks v2. Different combinations of `-p` and `-r` were evaluated to optimize SNP retention. The p8-r8 and p7-r8 datasets were selected for subsequent filtering with VCFtools.

## Stacks populations command

```bash
nohup populations \
  -P stacks_final/ \
  -M metadata/PopMap.txt \
  -O populations/populations_p7_r8/ \
  -p 7 \
  -r 0.8 \
  --min-maf 0.05 \
  --write-single-snp \
  --vcf \
  &> SLEW_P7_R8_log &
```


