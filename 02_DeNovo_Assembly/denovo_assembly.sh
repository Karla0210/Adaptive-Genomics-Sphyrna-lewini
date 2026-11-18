#!/usr/bin/env bash

# De novo assembly - Sphyrna lewini
# Stacks: denovo_map.pl

# Activate Conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate stacks

# Input and output directories
INPUT_DIR="data/clean_reads"
METADATA_FILE="metadata/PopMap.txt"
OUTPUT_DIR="results/stacks_final"

# Run de novo assembly
nohup denovo_map.pl \
  --samples "${INPUT_DIR}/" \
  --popmap "${METADATA_FILE}" \
  -o "${OUTPUT_DIR}" \
  -m 3 \
  -M 5 \
  -n 4 \
  -T 2 \
  --min-samples-per-pop 0.80 \
  &> process_log &
