
#!/bin/bash

# Trimmomatic - Sphyrna lewini
# Single-end read trimming and Nextera adapter removal

RAW_READS="/path/to/files/"
OUT_DIR="clean_reads"
ADAPTERS="/path/to.adapters/NexteraPE-PE.fa"

mkdir -p ${OUT_DIR}

for R1 in ${RAW_READS}/*.F.fastq.gz; do
  S=$(basename "$R1" .F.fastq.gz)

# Trimming reads using Trimmomatic
trimmomatic SE -phred33 \
${RAW_READS}/${S}.F.fastq.gz \
${OUT_DIR}/${S}.F.clean.fastq.gz \
ILLUMINACLIP:${ADAPTERS}:2:30:10 \
CROP:120 \
MINLEN:120

done

# Check read lengths of the cleaned files
for CLEAN in "${OUT_DIR}"/*.F.clean.fastq.gz; do
  [ -f "$CLEAN" ] || continue
  printf 'Read lengths for %s:\n' "$(basename "$CLEAN")"
  zcat "$CLEAN" | \
  awk '(NR%4==2){print length($0)}' | \
  sort -n | uniq -c
done
