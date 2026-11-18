
#!/usr/bin/env bash

# FastQC and MultiQC - Sphyrna lewini
# Quality assessment of trimmed reads

# Input and output directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT_DIR="${1:-${SCRIPT_DIR}/clean_reads}"
FASTQC_DIR="${SCRIPT_DIR}/fastqc_results"
MULTIQC_DIR="${SCRIPT_DIR}/multiqc_results"

# Create output directories
mkdir -p "$FASTQC_DIR" "$MULTIQC_DIR"

# Check software versions
fastqc --version
multiqc --version

# Check input files
shopt -s nullglob
FASTQ_FILES=("$INPUT_DIR"/*.fastq.gz)

if (( ${#FASTQ_FILES[@]} == 0 )); then
    echo "Error: No FASTQ files found in $INPUT_DIR" >&2
    exit 1
fi

# Run FastQC on trimmed reads
fastqc "${FASTQ_FILES[@]}" \
    --threads 2 \
    --outdir "$FASTQC_DIR"

# Generate the MultiQC summary report
multiqc "$FASTQC_DIR" \
    --outdir "$MULTIQC_DIR" \
    --filename multiqc_report.html
