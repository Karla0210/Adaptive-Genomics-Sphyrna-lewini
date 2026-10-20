# =========================================================
# Pairwise Relatedness Heatmap (RELATEDNESS_PHI)
# VCFtools --relatedness2
# =========================================================

setwd("E:\\Filtrado populations y VCF p7 r8/Relatedness/")

library(tidyverse)
library(reshape2)
library(pheatmap)

# 1. Read relatedness file
rel <- read.table(
  "relatedness_cervus.relatedness2",
  header = TRUE,
  stringsAsFactors = FALSE
)

# 2. Create individual-by-individual matrix
mat <- acast(rel, INDV1 ~ INDV2, value.var = "RELATEDNESS_PHI")

# Make the matrix symmetric
mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]

# Remove diagonal values
diag(mat) <- NA

# 3. Clean individual names
clean_names <- function(x) {
  x %>%
    str_remove("\\.F\\.clean$") %>%
    str_remove("\\.clean$") %>%
    str_trim()
}

rownames(mat) <- clean_names(rownames(mat))
colnames(mat) <- clean_names(colnames(mat))

# 4. Read PopMap
popmap <- read.table(
  "PopMap.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(popmap) <- c("Individual", "Region")

popmap$Individual <- clean_names(popmap$Individual)
popmap$Region <- trimws(popmap$Region)

# 5. Locality order as presented in the manuscript
region_order <- c("LP", "LA", "RF", "MZ", "CM", "TP", "CT", "CA")

annotation <- popmap[match(rownames(mat), popmap$Individual), "Region", drop = FALSE]
rownames(annotation) <- rownames(mat)
colnames(annotation) <- "Region"

annotation$Region <- factor(annotation$Region, levels = region_order)

# Check for unmatched individuals
if (any(is.na(annotation$Region))) {
  warning("Some individuals in the matrix were not found in the PopMap:")
  print(rownames(annotation)[is.na(annotation$Region)])
}

# 6. Locality colors matching the map
ann_colors <- list(
  Region = c(
    LP = "#C7323E",
    LA = "#F3724B",
    RF = "#FCC55A",
    MZ = "#66BC7D",
    CM = "#327146",
    TP = "#70C2DA",
    CT = "#3598C4",
    CA = "#316AAF"
  )
)

# 7. Center the color scale at zero
max_abs <- max(abs(mat), na.rm = TRUE)

breaks_centered <- seq(
  -max_abs,
  max_abs,
  length.out = 101
)

colors_centered <- colorRampPalette(
  c("blue", "white", "red")
)(100)

# 8. Display heatmap in RStudio
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

# 9. Save as PDF
pdf("heatmap_relatedness_centered.pdf", width = 12, height = 10)

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

dev.off()

# 10. Save as PNG
png("heatmap_relatedness_centered.png", width = 3000, height = 2500, res = 300)

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

dev.off()
