
#!/bin/bash

# Population filtering - Sphyrna lewini
# Stacks v2.0: populations

module purge
module load stacks/2.0

# Input files and directories
STACKS_DIR="results/stacks_final"
POPMAP="metadata/PopMap.txt"
OUTPUT_BASE="results/populations"

# Common parameters
MIN_MAF=0.05

# p8 r0.8
mkdir -p "${OUTPUT_BASE}/populations_p8_r8"

populations \
-P "${STACKS_DIR}" \
-M "${POPMAP}" \
-O "${OUTPUT_BASE}/populations_p8_r8" \
-p 8 \
-r 0.8 \
--min-maf "${MIN_MAF}" \
--write-single-snp \
--vcf \
> "${OUTPUT_BASE}/populations_p8_r8.log" 2>&1


# p7 r0.8
mkdir -p "${OUTPUT_BASE}/populations_p7_r8"

populations \
-P "${STACKS_DIR}" \
-M "${POPMAP}" \
-O "${OUTPUT_BASE}/populations_p7_r8" \
-p 7 \
-r 0.8 \
--min-maf "${MIN_MAF}" \
--write-single-snp \
--vcf \
> "${OUTPUT_BASE}/populations_p7_r8.log" 2>&1


# p7 r0.5
mkdir -p "${OUTPUT_BASE}/populations_p7_r5"

populations \
-P "${STACKS_DIR}" \
-M "${POPMAP}" \
-O "${OUTPUT_BASE}/populations_p7_r5" \
-p 7 \
-r 0.5 \
--min-maf "${MIN_MAF}" \
--write-single-snp \
--vcf \
> "${OUTPUT_BASE}/populations_p7_r5.log" 2>&1


# p5 r0.8
mkdir -p "${OUTPUT_BASE}/populations_p5_r8"

populations \
-P "${STACKS_DIR}" \
-M "${POPMAP}" \
-O "${OUTPUT_BASE}/populations_p5_r8" \
-p 5 \
-r 0.8 \
--min-maf "${MIN_MAF}" \
--write-single-snp \
--vcf \
> "${OUTPUT_BASE}/populations_p5_r8.log" 2>&1


# p5 r0.5
mkdir -p "${OUTPUT_BASE}/populations_p5_r5"

populations \
-P "${STACKS_DIR}" \
-M "${POPMAP}" \
-O "${OUTPUT_BASE}/populations_p5_r5" \
-p 5 \
-r 0.5 \
--min-maf "${MIN_MAF}" \
--write-single-snp \
--vcf \
> "${OUTPUT_BASE}/populations_p5_r5.log" 2>&1


# p3 r0.8
mkdir -p "${OUTPUT_BASE}/populations_p3_r8"

populations \
-P "${STACKS_DIR}" \
-M "${POPMAP}" \
-O "${OUTPUT_BASE}/populations_p3_r8" \
-p 3 \
-r 0.8 \
--min-maf "${MIN_MAF}" \
--write-single-snp \
--vcf \
> "${OUTPUT_BASE}/populations_p3_r8.log" 2>&1


# p3 r0.5
mkdir -p "${OUTPUT_BASE}/populations_p3_r5"

populations \
-P "${STACKS_DIR}" \
-M "${POPMAP}" \
-O "${OUTPUT_BASE}/populations_p3_r5" \
-p 3 \
-r 0.5 \
--min-maf "${MIN_MAF}" \
--write-single-snp \
--vcf \
> "${OUTPUT_BASE}/populations_p3_r5.log" 2>&1



