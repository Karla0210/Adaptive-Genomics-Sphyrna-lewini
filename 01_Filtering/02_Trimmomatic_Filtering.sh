

module load spack
module load trimmomatic

RAW_READS="/path/to/files/"
OUT_DIR="trimmed"
ADAPTERS="/path/to.adapters/NexteraPE-PE.fa"

mkdir -p ${OUT_DIR}

for R1 in ${RAW_READS}/*.F.fastq.gz; do
  S=$(basename "$R1" .F.fastq.gz)

trimmomatic SE -phred33 \
${RAW_READS}/${S}.F.fastq.gz \
${OUT_DIR}/${S}.F.clean.fastq.gz \
ILLUMINACLIP:${ADAPTERS}:2:30:10 \
CROP:120 \
MINLEN:120

done
