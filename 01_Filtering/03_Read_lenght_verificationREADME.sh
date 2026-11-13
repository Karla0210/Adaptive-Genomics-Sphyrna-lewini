zcat BC_7.F.clean.fastq.gz | \
awk '(NR%4==2){print length($0)}' | \
sort -n | uniq -c

