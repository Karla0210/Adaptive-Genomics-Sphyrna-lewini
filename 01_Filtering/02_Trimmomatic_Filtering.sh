trimmomatic SE -phred33 \
TE_17.F.fastq.gz \
TE_17.F.clean.fastq.gz \
ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 \
CROP:120 \
MINLEN:120
