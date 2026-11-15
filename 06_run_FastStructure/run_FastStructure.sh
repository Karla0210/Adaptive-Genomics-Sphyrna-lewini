
#!/bin/bash

# Population structure analysis - Sphyrna lewini
# fastStructure: K = 1-8, 10 replicates per K

# Input dataset prefix (change for each SNP dataset)
DATA="sln"

# Run fastStructure for K = 1 to 8
for K in {1..8}; do
    echo "Running K=$K"

    # Ten replicate runs per K
    for REP in {1..10}; do
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
