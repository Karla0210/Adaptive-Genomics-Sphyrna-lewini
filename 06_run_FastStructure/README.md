
# 06. Population structure analysis (fastStructure)

Population genetic structure was evaluated using fastStructure (Raj et al., 2014) with a logistic prior. Analyses were performed independently for neutral SNPs (9,529 SNPs), consensus outliers (29 SNPs), differentiation-based outliers (19 SNPs), and environment-associated outliers (64 SNPs), testing K = 1–8 with 10 replicate runs per K.

## Script

`run_faststructure.sh` contains the commands used to perform the population structure analyses.

## Execution

The script was executed separately for each SNP dataset by modifying the input prefix (`DATA`).

Output files were subsequently processed in R using the `pophelper` package to visualize population ancestry.

