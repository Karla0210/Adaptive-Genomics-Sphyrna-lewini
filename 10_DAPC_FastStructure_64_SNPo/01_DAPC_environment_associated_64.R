#!/usr/bin/env Rscript

# DAPC - environment-associated SNPs (64) | Sphyrna lewini
# Packages: dplyr, adegenet, vcfR, paletteer, dartR

# 01. Input and output paths
PROJECT_DIR <- "."
VCF_FILE <- file.path(PROJECT_DIR, "SNPoutRDA.recode.vcf")
POPMAP_FILE <- file.path(PROJECT_DIR, "PopMap.txt")
OUTPUT_DIR <- file.path(PROJECT_DIR, "Figuras_DAPC")

stopifnot(file.exists(VCF_FILE), file.exists(POPMAP_FILE))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 02. Packages and genotype data
library(dplyr)
library(adegenet)
library(vcfR)
library(paletteer)
library(dartR)

slo_vcf <- read.vcfR(VCF_FILE)
slo_gi <- vcfR2genind(slo_vcf)
slo_gl <- gi2gl(slo_gi)
ind_loc <- read.table(POPMAP_FILE)

if (nrow(ind_loc) != nInd(slo_gl)) {
  stop("PopMap and VCF have different numbers of individuals.")
}
if (ncol(ind_loc) >= 2 && !identical(as.character(ind_loc[[1]]), indNames(slo_gl))) {
  stop("Sample order differs between PopMap and VCF. Check sample IDs before proceeding.")
}
slo_gl$pop <- as.factor(ind_loc$V2)

# 03. Cross-validation and clustering
slotab <- tab(slo_gl, NA.method = "mean")
valslo <- xvalDapc(
  slotab, ind_loc$V2,
  n.pca.max = 500, training.set = 0.90,
  result = "groupMean", center = TRUE, scale = FALSE,
  n.pca = NULL, n.rep = 50, xval.plot = TRUE
)
print(valslo$`Number of PCs Achieving Highest Mean Success`)

find.clusters(slo_gl, max.n.clust = 8, stat = "BIC",
              n.iter = 100000, n.start = 10000)

# 04. DAPC and locality colors
# Locality colors
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
# DAPC: 8 PCs and 7 discriminant functions
slodapc <- dapc(slo_gl, n.pca = 8, n.da = 7)
popNames(slo_gl)
# Calculate Explained Variance
var_explained <- slodapc$eig / sum(slodapc$eig) * 100

var_explained
DA1_porcentaje <- round(var_explained[1], 2)

DA2_porcentaje <- round(var_explained[2], 2)

DA1_porcentaje

DA2_porcentaje
coord_dapc <- as.data.frame(
  slodapc$ind.coord[, 1:2, drop = FALSE]
)
colnames(coord_dapc) <- c("DA1", "DA2")
coord_dapc$Localidad <- as.character(slodapc$grp)
table(coord_dapc$Localidad)
centroides_dapc <- coord_dapc %>%
  group_by(Localidad) %>%
  summarise(
    DA1 = mean(DA1, na.rm = TRUE),
    DA2 = mean(DA2, na.rm = TRUE),
    .groups = "drop"
  )

centroides_dapc
etiqueta_LA <- centroides_dapc %>%
  filter(Localidad == "LA")
etiqueta_CT <- centroides_dapc %>%
  filter(Localidad == "CT")
etiqueta_CA <- centroides_dapc %>%
  filter(Localidad == "CA")
# Function to Generate the DAPC Plot

# 05. Main DAPC plot
graficar_DAPC <- function() {
  parametros_originales <- par(no.readonly = TRUE)
  par(
    mar = c(1.5, 1.5, 1, 1),
    oma = c(5, 5, 0.5, 0.5),
    cex.axis = 1.25,
    las = 1
  )
  # MAIN DAPC PLOT
  
  scatter(
    slodapc,
    1, 2,
    col = paletteenew[popNames(slo_gl)],
    xlab = "",
    ylab = "",
    pch = 16,
    cex = 3.5,
    solid = 1,
    legend = FALSE,
    cstar = FALSE,
    clabel = FALSE,
    cellipse = 1.5,
    scree.da = FALSE
  )
  # SNP DATASET TITLE
  
  limites <- par("usr")
  
  text(
    x = limites[2] - 0.08 * diff(limites[1:2]),
    y = limites[4] - 0.07 * diff(limites[3:4]),
    labels = "Environment-associated SNPs",
    adj = c(1, 1),
    cex = 2.5,
    font = 2,
    col = "black"
  )
  # AXIS TITLES
  
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
  # LOCALITY LABELS
  
  desplazamiento_y <- diff(
    range(coord_dapc$DA2, na.rm = TRUE)
  ) * 0.055
  text(
    x = etiqueta_CA$DA1,
    y = etiqueta_CA$DA2 + desplazamiento_y,
    labels = "CA",
    cex = 1.6,
    font = 2,
    col = "#316AAF"
  )
  text(
    x = etiqueta_CT$DA1,
    y = etiqueta_CT$DA2 + desplazamiento_y,
    labels = "CT",
    cex = 1.6,
    font = 2,
    col = "#3598C4"
  )
  text(
    x = etiqueta_LA$DA1,
    y = etiqueta_LA$DA2 + desplazamiento_y,
    labels = "LA",
    cex = 1.6,
    font = 2,
    col = "#F3724B"
  )
  par(xpd = FALSE)
  # DA EIGENVALUES INSET
  colores_eigen <- rep(
    "grey85",
    length(slodapc$eig)
  )
  
  colores_eigen[
    1:min(2, length(slodapc$eig))
  ] <- "grey50"
  par(
    fig = c(
      0.76,  # izquierda
      0.99,  # derecha
      0.05,  # abajo
      0.26   # arriba
    ),
    mar = c(0.6, 0.6, 1.6, 0.4),
    new = TRUE,
    xpd = FALSE
  )
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
  box(
    col = "grey45",
    lwd = 0.8
  )
  title(
    main = "DA eigenvalues",
    cex.main = 0.9,
    font.main = 1,
    line = 0.2
  )
  par(parametros_originales)
}
# Display DAPC plot
graficar_DAPC()

# 06. Export main figure
output_png <- file.path(
  OUTPUT_DIR,
  "DAPC_outliers_64.png"
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
# Export DAPC as PDF

output_pdf <- file.path(
  OUTPUT_DIR,
  "DAPC_outliers_64.pdf"
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
# Display Output File Paths

output_png

output_pdf
# DAPC Density Plot (DA1)
# 07. DA1 density plot
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
graficar_DAPC_densidad()
# Export Density Plot as PNG

output_png_densidad <- file.path(
  OUTPUT_DIR,
  "DAPC_outliers_64_densidad_DA1.png"
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
# Export Density Plot as PDF

output_pdf_densidad <- file.path(
  OUTPUT_DIR,
  "DAPC_outliers_64_densidad_DA1.pdf"
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
print(c(output_png_densidad, output_pdf_densidad))
