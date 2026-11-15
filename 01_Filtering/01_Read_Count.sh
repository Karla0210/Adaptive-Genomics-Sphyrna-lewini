
#!/bin/bash

# Raw read counting
# Count the number of reads per sample from compressed FASTQ files
# stored in the original sequencing archive.

# Input and output files
INPUT_ARCHIVE="hammerhead.tar.gz"
OUTPUT_FILE="raw_reads_count.tsv"

# Count raw reads per sample
tar -tzf "${INPUT_ARCHIVE}" |
grep -E "\.(fastq|fq)\.gz$" |
while read -r file; do
    lines=$(tar -xOzf "${INPUT_ARCHIVE}" "$file" | zcat | wc -l)
    reads=$((lines / 4))
    printf "%s\t%s\n" "$file" "$reads"
done > "${OUTPUT_FILE}"

