
module purge
module load vcftools/0.1.16
module load plink/1.9

VCF="path/to/populations_08.vcf.recode.vcf"
OUT_DIR="path/to/output_directory"

# 1. Convert VCF to TPED/TFAM
vcftools \
  --vcf ${VCF} \
  --out ${OUT_DIR}/hammerhead_tped \
  --plink-tped

# 2. Convert TPED/TFAM to PLINK genotype format
plink1.9 \
  --tped ${OUT_DIR}/hammerhead_tped.tped \
  --tfam ${OUT_DIR}/hammerhead_tped.tfam \
  --recode \
  --out ${OUT_DIR}/RDA_gen \
  --allow-extra-chr