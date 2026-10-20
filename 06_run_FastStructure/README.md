
# 10. Population structure analysis (fastStructure)

Population genetic structure was evaluated using fastStructure with a logistic prior. Analyses were performed independently for neutral SNPs, consensus outliers (29 SNPs), FST outliers (19 SNPs), and RDA candidates (64 SNPs), testing K = 1–8 with 10 replicate runs per K.

## fastStructure script

```bash
#!/bin/bash

# Input dataset prefix
DATA="sln"

# Test K = 1 to 8
for K in {1..8}
do
    echo "Running K=$K"

    # Ten replicate runs per K
    for REP in {1..10}
    do
        echo "Replicate $REP for K=$K"

        fastStructure \
            --input "$DATA" \
            --output "${DATA}.${REP}.${K}" \
            --format str \
            -K "$K" \
            --prior logistic
    done
done

echo "All analyses completed."
```

## Execution

```bash
chmod +x run_faststructure.sh
./run_faststructure.sh
```

The script was executed separately for each SNP dataset by modifying the input prefix (`DATA`). Output files were subsequently processed in R using the `pophelper` package to visualize population ancestry.
