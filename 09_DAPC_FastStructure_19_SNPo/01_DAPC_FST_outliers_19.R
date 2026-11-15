#!/usr/bin/env Rscript

# DAPC analysis - Sphyrna lewini
# FST-based outlier SNPs (19 SNPs)

library(dplyr)
library(adegenet)
library(vcfR)
library(dartR)

# 01. Input and output paths
PROJECT_DIR <- "."
VCF_FILE <- file.path(PROJECT_DIR, "SNPoutFST.recode.vcf")
POPMAP_FILE <- file.path(PROJECT_DIR, "PopMap.txt")
OUTPUT_DIR <- file.path(PROJECT_DIR, "Figuras_DAPC")

stopifnot(file.exists(VCF_FILE), file.exists(POPMAP_FILE))
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 02. Read SNPs and assign sampling localities
slo_vcf <- read.vcfR(VCF_FILE)
slo_gi <- vcfR2genind(slo_vcf)
slo_gl <- gi2gl(slo_gi)
ind_loc <- read.table(POPMAP_FILE, header = FALSE, stringsAsFactors = FALSE)
stopifnot(nInd(slo_gl) == nrow(ind_loc))
# The PopMap rows must match the genotype order.
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

# 04. DAPC model and variance explained
paletteenew <- c(
  LP = "#C7323E", LA = "#F3724B", RF = "#FCC55A", MZ = "#66BC7D",
  CM = "#327146", TP = "#70C2DA", CT = "#3598C4", CA = "#316AAF"
)
# Parameters documented in the original script: 6 PCs and 6 discriminant functions.
slodapc <- dapc(slo_gl, n.pca = 6, n.da = 6)
print(popNames(slo_gl))
var_explained <- slodapc$eig / sum(slodapc$eig) * 100
DA1_porcentaje <- round(var_explained[1], 2)
DA2_porcentaje <- round(var_explained[2], 2)

# 05. Locality centroids and LP annotation
coord_dapc <- as.data.frame(slodapc$ind.coord[, 1:2, drop = FALSE])
colnames(coord_dapc) <- c("DA1", "DA2")
coord_dapc$Localidad <- as.character(slodapc$grp)
centroides_dapc <- coord_dapc %>%
  group_by(Localidad) %>%
  summarise(DA1 = mean(DA1, na.rm = TRUE),
            DA2 = mean(DA2, na.rm = TRUE), .groups = "drop")
etiqueta_LP <- centroides_dapc %>% filter(Localidad == "LP")

# 06. DAPC figures
graficar_DAPC <- function() {
  
  # Guardar parámetros gráficos originales
  parametros_originales <- par(no.readonly = TRUE)
  
  # Márgenes internos y exteriores
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
    
    # Colores de las localidades
    col = paletteenew[popNames(slo_gl)],
    
    # Quitar títulos automáticos
    xlab = "",
    ylab = "",
    
    # Círculos sólidos
    pch = 16,
    
    # Tamaño de los círculos
    cex = 4.5,
    
    # Colores sólidos
    solid = 1,
    
    # Sin leyenda
    legend = FALSE,
    
    # Sin líneas hacia los centroides
    cstar = FALSE,
    
    # Sin etiquetas automáticas
    clabel = FALSE,
    
    # Tamaño de las elipses
    cellipse = 1.5,
    
    # Desactivar el screeplot automático
    scree.da = FALSE
  )
  
  # SNP DATASET TITLE
  
  limites <- par("usr")
  
  text(
    x = limites[2] - 0.05 * diff(limites[1:2]),
    y = limites[4] - 0.07 * diff(limites[3:4]),
    labels = "FST-based SNPs",
    adj = c(1, 1),
    cex = 2.5,
    font = 2,
    col = "black"
  )
  
  # AXIS TITLES
  
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
  
  # LP LABEL
  
  desplazamiento_x <- diff(
    range(coord_dapc$DA1, na.rm = TRUE)
  ) * 0.035
  
  desplazamiento_y <- diff(
    range(coord_dapc$DA2, na.rm = TRUE)
  ) * 0.055
  
  text(
    x = etiqueta_LP$DA1 - desplazamiento_x,
    y = etiqueta_LP$DA2 + desplazamiento_y,
    labels = "LP",
    cex = 1.6,
    font = 2,
    col = "#C7323E"
  )
  
  # DA EIGENVALUES INSET
  
  # Definir los colores de las barras:
  # DA1 y DA2 en gris oscuro; los demás en gris claro
  colores_eigen <- rep(
    "grey85",
    length(slodapc$eig)
  )
  
  colores_eigen[1:min(2, length(slodapc$eig))] <- "grey50"
  
  # Crear una nueva región pequeña sobre la figura principal
  par(
    fig = c(
      0.05,  # límite izquierdo
      0.27,  # límite derecho
      0.015, # límite inferior
      0.19   # límite superior
    ),
    mar = c(0.8, 0.8, 1.8, 0.5),
    new = TRUE,
    xpd = FALSE
  )
  
  # Dibujar eigenvalues
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
  
  # Borde del recuadro
  box(
    col = "grey45",
    lwd = 0.8
  )
  
  # Título del recuadro
  title(
    main = "DA eigenvalues",
    cex.main = 0.9,
    font.main = 1,
    line = 0.3
  )
  
  # Restaurar parámetros gráficos
  par(parametros_originales)
}

# Main DAPC figure (16 x 10 inches, 300 dpi)
png(file.path(OUTPUT_DIR, "DAPC_outliers_19.png"),
    width = 16, height = 10, units = "in", res = 300, bg = "white")
graficar_DAPC()
dev.off()

pdf(file.path(OUTPUT_DIR, "DAPC_outliers_19.pdf"),
    width = 16, height = 10, family = "sans", useDingbats = FALSE, bg = "white")
graficar_DAPC()
dev.off()

# 07. DA1 density figure
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

png(file.path(OUTPUT_DIR, "DAPC_outliers_19_densidad_DA1.png"),
    width = 16, height = 10, units = "in", res = 300, bg = "white")
graficar_DAPC_densidad()
dev.off()

pdf(file.path(OUTPUT_DIR, "DAPC_outliers_19_densidad_DA1.pdf"),
    width = 16, height = 10, family = "sans", useDingbats = FALSE, bg = "white")
graficar_DAPC_densidad()
dev.off()
