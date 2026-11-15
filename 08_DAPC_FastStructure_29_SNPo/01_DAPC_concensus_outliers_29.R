#!/usr/bin/env Rscript

# DAPC - consensus outliers (29 SNPs) | Sphyrna lewini

# 01. Packages
library(dplyr)
library(adegenet)
library(vcfR)
library(paletteer)
library(dartR)

# 02. Input and output paths
PROJECT_DIR <- "."
VCF_FILE <- file.path(PROJECT_DIR, "Resumen_SNPo2.recode.vcf")
POPMAP_FILE <- file.path(PROJECT_DIR, "PopMap.txt")
output_dir <- file.path(PROJECT_DIR, "Figuras_DAPC")

stopifnot(file.exists(VCF_FILE), file.exists(POPMAP_FILE))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 03. Load genotype data and assign localities
slo_vcf <- read.vcfR(VCF_FILE)
slo_gi <- vcfR2genind(slo_vcf)
slo_gl <- gi2gl(slo_gi)
ind_loc <- read.table(POPMAP_FILE)
slo_gl$pop <- as.factor(ind_loc$V2)

# 04. Cross-validation and cluster evaluation
slotab <- tab(slo_gl, NA.method = "mean")
valslo <- xvalDapc(
  slotab,
  ind_loc$V2,
  n.pca.max = 500,
  training.set = 0.90,
  result = "groupMean",
  center = TRUE,
  scale = FALSE,
  n.pca = NULL,
  n.rep = 50,
  xval.plot = TRUE
)
print(valslo$`Number of PCs Achieving Highest Mean Success`)

find.clusters(
  slo_gl,
  max.n.clust = 8,
  stat = "BIC",
  n.iter = 100000,
  n.start = 10000
)

# 05. DAPC and plotting settings
# Define specific colors for each locality
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
slodapc <- dapc(slo_gl, n.pca = 6, n.da = 6)

# Verify population order
popNames(slo_gl)


# Axis Percentages


# Percentage of discriminant variation
var_explained <- slodapc$eig / sum(slodapc$eig) * 100

var_explained


# Percentage of DA1 and DA2
DA1_porcentaje <- round(var_explained[1], 2)

DA2_porcentaje <- round(var_explained[2], 2)


DA1_porcentaje

DA2_porcentaje


# Coordinates For Ca, Lp, And Ct


# Extract DA1 and DA2 coordinates
coord_dapc <- as.data.frame(
  slodapc$ind.coord[, 1:2, drop = FALSE]
)


# Rename columns
colnames(coord_dapc) <- c("DA1", "DA2")


# Add locality
coord_dapc$Localidad <- as.character(slodapc$grp)


# Check localities
table(coord_dapc$Localidad)


# Calculate the centroid of each locality
centroides_dapc <- coord_dapc %>%
  group_by(Localidad) %>%
  summarise(
    DA1 = mean(DA1, na.rm = TRUE),
    DA2 = mean(DA2, na.rm = TRUE),
    .groups = "drop"
  )


centroides_dapc


# Extract CA coordinates
etiqueta_CA <- centroides_dapc %>%
  filter(Localidad == "CA")


# Extract LP coordinates
etiqueta_LP <- centroides_dapc %>%
  filter(Localidad == "LP")


# Extract CT coordinates
etiqueta_CT <- centroides_dapc %>%
  filter(Localidad == "CT")

# Plotting Function

graficar_DAPC <- function() {
  
  # Save original graphical parameters
  parametros_originales <- par(no.readonly = TRUE)
  
  # Inner and outer margins
  par(
    mar = c(1.5, 1.5, 1, 1),
    oma = c(5, 5, 0.5, 0.5),
    cex.axis = 1.25,
    las = 1
  )
  
  # Plot DAPC
  scatter(
    slodapc,
    1, 2,
    
    # Locality colors
    col = paletteenew[popNames(slo_gl)],
    
    # Remove automatic titles
    xlab = "",
    ylab = "",
    
    # Solid circles
    pch = 16,
    
    # Circle size
    cex = 4.5,
    
    # Solid colors
    solid = 1,
    
    # No legend
    legend = FALSE,
    
    # No lines to centroids
    cstar = FALSE,
    
    # No automatic labels
    clabel = FALSE,
    
    # Ellipse size
    cellipse = 1.5,
    
    # The scree plot will be added manually
    scree.da = FALSE
  )
  
  # Main plot limits
  limites <- par("usr")
  
  # Title in the upper-right corner
  text(
    x = limites[2] - 0.08 * diff(limites[1:2]),
    y = limites[4] - 0.07 * diff(limites[3:4]),
    labels = "Consensus SNPs",
    adj = c(1, 1),
    cex = 2.5,
    font = 2,
    col = "black"
  )
  
  # Allow text placement outside the panel
  par(xpd = NA)
  
  # X-axis title
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
  
  # Y-axis title
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
  
  # Automatic offset for labels
  desplazamiento_x <- diff(
    range(coord_dapc$DA1, na.rm = TRUE)
  ) * 0.035
  
  desplazamiento_y <- diff(
    range(coord_dapc$DA2, na.rm = TRUE)
  ) * 0.055
  
  # CA label
  text(
    x = etiqueta_CA$DA1,
    y = etiqueta_CA$DA2 + desplazamiento_y,
    labels = "CA",
    cex = 1.6,
    font = 2,
    col = "#316AAF"
  )
  
  # LP label
  text(
    x = etiqueta_LP$DA1 - desplazamiento_x,
    y = etiqueta_LP$DA2 + desplazamiento_y,
    labels = "LP",
    cex = 1.6,
    font = 2,
    col = "#C7323E"
  )
  
  # CT label
  text(
    x = etiqueta_CT$DA1,
    y = etiqueta_CT$DA2 + desplazamiento_y,
    labels = "CT",
    cex = 1.6,
    font = 2,
    col = "#3598C4"
  )
  
  # Restore clipping
  par(xpd = FALSE)
  

  # DA EIGENVALUES INSET

  
  # Bar colors:
  # DA1 and DA2 in dark gray and the remaining bars in light gray
  colores_eigen <- rep(
    "grey85",
    length(slodapc$eig)
  )
  
  colores_eigen[
    1:min(2, length(slodapc$eig))
  ] <- "grey50"
  
  # Create a small inset in the lower-left corner
  par(
    fig = c(
      0.02,   # left
      0.25,   # right
      0.005,  # bottom
      0.17    # top
    ),
    mar = c(0.6, 0.6, 1.6, 0.4),
    new = TRUE,
    xpd = FALSE
  )
  
  # Draw eigenvalue bars
  barplot(
    height = slodapc$eig,
    col = colores_eigen,
    border = "grey35",
    axes = FALSE,
    names.arg = rep(
      "",
      length(slodapc$eig)
    ),
    space = 0.25
  )
  
  # Inset border
  box(
    col = "grey45",
    lwd = 0.8
  )
  
  # Inset title
  title(
    main = "DA eigenvalues",
    cex.main = 0.9,
    font.main = 1,
    line = 0.2
  )
  
  # Restore graphical parameters
  par(parametros_originales)
}

# Display The Plot In R


graphics.off()

graficar_DAPC()


# Output Directory


# Output directory created above


# Save As Png


output_png <- file.path(
  output_dir,
  "DAPC_outliers_29.png"
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


# Save As Pdf


output_pdf <- file.path(
  output_dir,
  "DAPC_outliers_29.pdf"
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


# Display File Paths


output_png

output_pdf

# Dapc Axis 1: Density Plot

# Function to generate the density plot
graficar_DAPC_densidad <- function() {
  
  par(
    mar = c(5, 5, 2, 2),
    cex.axis = 1.2,
    cex.lab = 1.3
  )
  
  scatter(
    slodapc,
    1, 1,
    col = paletteenew[popNames(slo_gl)],
    cex = 2.5,
    legend = TRUE,
    posi.da = "bottomright",
    cstar = FALSE,
    solid = 1,
    pch = 16
  )
}

# Display in R
graficar_DAPC_densidad()


# Save Density Plot As Png

output_png_densidad <- file.path(
  output_dir,
  "DAPC_outliers_29_densidad_DA1.png"
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


# Save Density Plot As Pdf

output_pdf_densidad <- file.path(
  output_dir,
  "DAPC_outliers_29_densidad_DA1.pdf"
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


# Percentage Of Explained Variation

var_explained <- slodapc$eig / sum(slodapc$eig) * 100

var_explained


# Display Paths

output_png_densidad

output_pdf_densidad
