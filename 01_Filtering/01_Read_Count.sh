tar -tzf hammerhead.tar.gz |
grep -E "\.(fastq|fq)\.gz$" |
while read -r file; do
    lines=$(tar -xOzf hammerhead.tar.gz "$file" | zcat | wc -l)
    reads=$((lines / 4))
    printf "%s\t%s\n" "$file" "$reads"
done > raw_reads_count.tsv
