#!/usr/bin/env Rscript

# DAPC - neutral SNPs | Sphyrna lewini
# Packages: dplyr, adegenet, vcfR, paletteer, dartR

# 01. Input and output paths
PROJECT_DIR <- "."  # Folder containing the VCF and PopMap.txt
VCF_FILE <- file.path(PROJECT_DIR, "Resumen_SNPn2.recode.vcf")
POPMAP_FILE <- file.path(PROJECT_DIR, "PopMap.txt")
OUTPUT_DIR <- file.path(PROJECT_DIR, "Figuras_DAPC")
stopifnot(file.exists(VCF_FILE), file.exists(POPMAP_FILE))
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# 02. Packages and data
# Neutral snp dataset
graphics.off()


library(dplyr)
library(adegenet)
library(vcfR)
library(paletteer)
library(dartR)


# Load neutral snp dataset
read.vcfR(VCF_FILE) -> sln_vcf


# Convert VCF to genind
vcfR2genind(sln_vcf) -> sln_gi


# Convert genind to genlight
gi2gl(sln_gi) -> sln_gl


# Assign individuals to sampling localities
read.table(POPMAP_FILE, header = FALSE, stringsAsFactors = FALSE) -> ind_loc

# Check the genotype and locality rows match before assigning populations.
stopifnot(nInd(sln_gl) == nrow(ind_loc))
stopifnot(identical(as.character(indNames(sln_gl)), as.character(ind_loc$V1)))
sln_gl$pop <- factor(ind_loc$V2, levels = c("LP", "LA", "RF", "MZ", "CM", "TP", "CT", "CA"))


# Cross-validation to determine the optimal number of pcs
tab(
  sln_gl,
  NA.method = "mean"
) -> slntab

slntab


xvalDapc(
  slntab,
  ind_loc$V2,
  n.pca.max = 500,
  training.set = 0.90,
  result = "groupMean",
  center = TRUE,
  scale = FALSE,
  n.pca = NULL,
  n.rep = 50,
  xval.plot = TRUE
) -> valslo


valslo$`Number of PCs Achieving Highest Mean Success`


# Selected number of principal components: 12
# Optional exploratory clustering (may prompt interactively in RStudio).
if (interactive()) {
  find.clusters(
    sln_gl,
    max.n.clust = 8,
    stat = "BIC",
    n.iter = 100000,
    n.start = 10000
  )
}


# Perform dapc
# Define sampling locality colors
paletteenew <- c(
  LP = "#C7323E",
  LA = "#F3724B",
  RF = "#FCC55A",
  MZ = "#66BC7D",
  CM = "#327146",
  TP = "#70C2DA",
  CT = "#3598C4",
  CA = "#316AAF"
)


# Run DAPC
dapc(sln_gl, n.pca = 12, n.da = 7) -> slndapc


# Selected parameters:
# Number of principal components = 12
# Number of discriminant functions = 7


# Verify population order
popNames(sln_gl)


# Calculate the percentage of explained variance
# Percentage of discriminant variance explained
var_explained <- slndapc$eig /
  sum(slndapc$eig) * 100

var_explained


# Variance explained by DA1 and DA2
DA1_porcentaje <- round(
  var_explained[1],
  2
)

DA2_porcentaje <- round(
  var_explained[2],
  2
)


DA1_porcentaje

DA2_porcentaje


# Function to generate the dapc figure
graficar_DAPC <- function() {
  
  # Save current graphical parameters
  parametros_originales <- par(
    no.readonly = TRUE
  )
  
  # Configure internal and external margins
  par(
    mar = c(1.5, 1.5, 1, 1),
    oma = c(5, 5, 0.5, 0.5),
    cex.axis = 1.25,
    las = 1
  )
  
  
  ##########################################################
  # Main DAPC plot
  ##########################################################
  
  scatter(
    slndapc,
    1, 2,
    
    # Sampling locality colors
    col = paletteenew[
      popNames(sln_gl)
    ],
    
    # Remove default axis labels
    xlab = "",
    ylab = "",
    
    # Solid circles
    pch = 16,
    
    # Point size
    cex = 4.5,
    
    # Solid colors
    solid = 1,
    
    # Hide legend
    legend = FALSE,
    
    # Hide lines connecting individuals to centroids
    cstar = FALSE,
    
    # Hide automatic labels
    clabel = FALSE,
    
    # Ellipse size
    cellipse = 1.5,
    
    # Disable default scree plot
    scree.da = FALSE
  )
  
  
  ##########################################################
  # Panel title
  ##########################################################
  
  limites <- par("usr")
  
  text(
    x = limites[2] -
      0.05 * diff(limites[1:2]),
    
    y = limites[4] -
      0.07 * diff(limites[3:4]),
    
    labels = "Neutral SNPs",
    adj = c(1, 1),
    cex = 2.5,
    font = 2,
    col = "black"
  )
  
  
  ##########################################################
  # Axis labels
  ##########################################################
  
  par(xpd = NA)
  
  mtext(
    text = paste0(
      "DA1 (",
      DA1_porcentaje,
      "%)"
    ),
    side = 1,
    outer = TRUE,
    line = 2.5,
    cex = 1.45,
    font = 1
  )
  
  mtext(
    text = paste0(
      "DA2 (",
      DA2_porcentaje,
      "%)"
    ),
    side = 2,
    outer = TRUE,
    line = 2.5,
    cex = 1.45,
    font = 1,
    las = 0
  )
  
  par(xpd = FALSE)
  
  
  ##########################################################
  # DA eigenvalues inset
  # Bottom-left position inside the main panel
  ##########################################################
  
  # Highlight DA1 and DA2 in dark gray
  colores_eigen <- rep(
    "grey85",
    length(slndapc$eig)
  )
  
  colores_eigen[
    1:min(
      2,
      length(slndapc$eig)
    )
  ] <- "grey50"
  
  
  # Create inset plot
  par(
    fig = c(
      0.01,  # left boundary
      0.25,  # right boundary
      0.01,  # bottom boundary
      0.28   # top boundary
    ),
    mar = c(
      0.6,
      0.6,
      1.6,
      0.4
    ),
    new = TRUE,
    xpd = FALSE
  )
  
  
  # Draw eigenvalue barplot
  barplot(
    height = slndapc$eig,
    col = colores_eigen,
    border = "grey35",
    axes = FALSE,
    names.arg = rep(
      "",
      length(slndapc$eig)
    ),
    space = 0.25
  )
  
  
  # Draw inset border
  box(
    col = "grey45",
    lwd = 0.8
  )
  
  
  # Add inset title
  title(
    main = "DA eigenvalues",
    cex.main = 0.9,
    font.main = 1,
    line = 0.2
  )
  
  
  # Restore graphical parameters
  par(parametros_originales)
}


# Display dapc plot
graphics.off()

if (interactive()) graficar_DAPC()

# Create output directory
output_dir <- OUTPUT_DIR


if (!dir.exists(output_dir)) {
  
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# Export dapc as png
output_png <- file.path(
  output_dir,
  "DAPC_neutros.png"
)


png(
  filename = output_png,
  width = 16,
  height = 10,
  units = "in",
  res = 300,
  bg = "white"
)


graficar_DAPC()


dev.off()


# Export dapc as pdf
output_pdf <- file.path(
  output_dir,
  "DAPC_neutros.pdf"
)


pdf(
  file = output_pdf,
  width = 16,
  height = 10,
  family = "sans",
  useDingbats = FALSE,
  bg = "white"
)


graficar_DAPC()


dev.off()


# Display output file paths
output_png

output_pdf


# Dapc density plot for da1
# Function to generate the DA1 density plot
graficar_DAPC_densidad <- function() {
  
  par(
    mar = c(5, 5, 2, 2),
    cex.axis = 1.2,
    cex.lab = 1.3
  )
  
  scatter(
    slndapc,
    1, 1,
    col = paletteenew[
      popNames(sln_gl)
    ],
    cex = 2.5,
    legend = TRUE,
    posi.da = "bottomright",
    cstar = FALSE,
    solid = 1,
    pch = 16
  )
}


# Display density plot
if (interactive()) graficar_DAPC_densidad()


# Export density plot as png
output_png_densidad <- file.path(
  output_dir,
  "DAPC_neutros_densidad_DA1.png"
)


png(
  filename = output_png_densidad,
  width = 16,
  height = 10,
  units = "in",
  res = 300,
  bg = "white"
)


graficar_DAPC_densidad()


dev.off()


# Export density plot as pdf
output_pdf_densidad <- file.path(
  output_dir,
  "DAPC_neutros_densidad_DA1.pdf"
)


pdf(
  file = output_pdf_densidad,
  width = 16,
  height = 10,
  family = "sans",
  useDingbats = FALSE,
  bg = "white"
)


graficar_DAPC_densidad()


dev.off()


# Calculate explained variance
var_explained <- slndapc$eig /
  sum(slndapc$eig) * 100

var_explained



# Reproducibility metadata
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "R_session_info.txt"))
