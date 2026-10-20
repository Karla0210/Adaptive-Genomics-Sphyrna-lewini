
# 02. De novo assembly

De novo assembly and SNP discovery were performed using Stacks v2 (`denovo_map.pl`). Different parameter combinations were evaluated, and the selected assembly used `m = 3`, `M = 5`, and `n = 4`.

## Stacks command

```bash
nohup denovo_map.pl \
  --samples clean_reads/ \
  --popmap metadata/PopMap.txt \
  -o stacks_final \
  -m 3 \
  -M 5 \
  -n 4 \
  -T 2 \
  --min-samples-per-pop 0.80 \
  &> process_log &
```

