#!/usr/bin/env Rscript

# Pairwise relatedness heatmap - Sphyrna lewini
# VCFtools --relatedness2 | RELATEDNESS_PHI
# Packages: tidyverse, reshape2, pheatmap

# 01. Input and output paths
PROJECT_DIR <- "."
RELATEDNESS_FILE <- file.path(PROJECT_DIR, "relatedness_cervus.relatedness2")
POPMAP_FILE <- file.path(PROJECT_DIR, "PopMap.txt")
OUTPUT_DIR <- file.path(PROJECT_DIR, "Relatedness_results")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(RELATEDNESS_FILE), file.exists(POPMAP_FILE))

# 02. Load packages
library(tidyverse)
library(reshape2)
library(pheatmap)

# 03. Read pairwise relatedness and prepare matrix
rel <- read.table(RELATEDNESS_FILE, header = TRUE, stringsAsFactors = FALSE)
mat <- acast(rel, INDV1 ~ INDV2, value.var = "RELATEDNESS_PHI")
mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
diag(mat) <- NA

# 04. Match individuals to sampling localities
clean_names <- function(x) {
  x %>%
    str_remove("\\.F\\.clean$") %>%
    str_remove("\\.clean$") %>%
    str_trim()
}

rownames(mat) <- clean_names(rownames(mat))
colnames(mat) <- clean_names(colnames(mat))

popmap <- read.table(POPMAP_FILE, header = FALSE, stringsAsFactors = FALSE)
colnames(popmap) <- c("Individual", "Region")
popmap$Individual <- clean_names(popmap$Individual)
popmap$Region <- trimws(popmap$Region)

region_order <- c("LP", "LA", "RF", "MZ", "CM", "TP", "CT", "CA")
annotation <- popmap[match(rownames(mat), popmap$Individual), "Region", drop = FALSE]
rownames(annotation) <- rownames(mat)
colnames(annotation) <- "Region"
annotation$Region <- factor(annotation$Region, levels = region_order)

if (any(is.na(annotation$Region))) {
  warning("Some individuals in the matrix were not found in the PopMap:")
  print(rownames(annotation)[is.na(annotation$Region)])
}

ann_colors <- list(Region = c(
  LP = "#C7323E", LA = "#F3724B", RF = "#FCC55A", MZ = "#66BC7D",
  CM = "#327146", TP = "#70C2DA", CT = "#3598C4", CA = "#316AAF"
))

# 05. Set heatmap colors and options
max_abs <- max(abs(mat), na.rm = TRUE)
breaks_centered <- seq(-max_abs, max_abs, length.out = 101)
colors_centered <- colorRampPalette(c("blue", "white", "red"))(100)

plot_relatedness <- function() {
  pheatmap(
    mat,
    color = colors_centered,
    breaks = breaks_centered,
    annotation_row = annotation,
    annotation_col = annotation,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    border_color = NA,
    na_col = "white",
    fontsize = 9,
    fontsize_row = 7,
    fontsize_col = 7,
    main = ""
  )
}

# 06. Export heatmap
pdf(file.path(OUTPUT_DIR, "heatmap_relatedness_centered.pdf"), width = 12, height = 10)
plot_relatedness()
dev.off()

png(file.path(OUTPUT_DIR, "heatmap_relatedness_centered.png"),
    width = 3000, height = 2500, res = 300)
plot_relatedness()
dev.off()

# Optional: display the heatmap in an interactive R session
if (interactive()) plot_relatedness()
