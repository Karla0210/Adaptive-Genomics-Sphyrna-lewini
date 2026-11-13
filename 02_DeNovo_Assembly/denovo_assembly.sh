

module purge
module load stacks/2.62


INPUT_DIR="/home/2ebioinfo/DATOS_KARLA/clean_reads"
METADATA_FILE="/home/2ebioinfo/DATOS_KARLA/info/PopMap.txt"
OUTPUT_DIR="/home/2ebioinfo/DATOS_KARLA/stacks_final"   

nohup denovo_map.pl \
  --samples ${INPUT_DIR}/ \
  --popmap ${METADATA_FILE} \
  -o ${OUTPUT_DIR} \
  -m 3 \
  -M 5 \
  -n 4 \
  -T 2 \
  --min-samples-per-pop 0.80 \
  &> process_log &