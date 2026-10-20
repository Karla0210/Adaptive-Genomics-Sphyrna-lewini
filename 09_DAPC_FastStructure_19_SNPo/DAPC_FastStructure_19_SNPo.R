############################################################
# FST Outlier SNP Dataset (19 SNPs)
############################################################

rm(list=ls())     
graphics.off()    
gc()


# Librerías necesarias 
library(dplyr)
library(adegenet)
library(vcfR)
library(paletteer)
library(dartR)


# Cargar carpeta 
setwd("E:\\Filtrado populations y VCF p7 r8/Filtrados_solo_fst/")


# Cargar datos
# Para outliers
read.vcfR("SNPoutFST.recode.vcf") -> slo_vcf 


# Convertir archivo VCF a genind
vcfR2genind(slo_vcf) -> slo_gi


# Convertir a genlight
gi2gl(slo_gi) -> slo_gl


# Asignar individuos a su localidad
read.table("PopMap.txt") -> ind_loc

ind_loc

slo_gl$pop <- as.factor(ind_loc$V2)


# Validar el DAPC para outliers
tab(slo_gl, NA.method = "mean") -> slotab

slotab


xvalDapc(
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
) -> valslo


valslo$`Number of PCs Achieving Highest Mean Success`


# Componentes que se usarán: 6
find.clusters(
  slo_gl,
  max.n.clust = 8,
  stat = "BIC",
  n.iter = 100000,
  n.start = 10000
)


############################################################
# Perform DAPC
############################################################


# Definir colores específicos para cada localidad
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


# Correr DAPC
dapc(slo_gl) -> slodapc

# Elegir:
# Número de PCs = 6
# Funciones discriminantes = 6


# Verificar orden de poblaciones
popNames(slo_gl)


############################################################
# Calculate Explained Variance
############################################################


# Porcentaje de variación discriminante
var_explained <- slodapc$eig / sum(slodapc$eig) * 100

var_explained


# Porcentaje de DA1 y DA2
DA1_porcentaje <- round(var_explained[1], 2)

DA2_porcentaje <- round(var_explained[2], 2)


DA1_porcentaje

DA2_porcentaje


############################################################
# Calculate Locality Centroids
############################################################


# Extraer coordenadas de DA1 y DA2
coord_dapc <- as.data.frame(
  slodapc$ind.coord[, 1:2, drop = FALSE]
)


# Renombrar columnas
colnames(coord_dapc) <- c("DA1", "DA2")


# Agregar localidad
coord_dapc$Localidad <- as.character(slodapc$grp)


# Revisar localidades
table(coord_dapc$Localidad)


# Calcular el centro de cada localidad
centroides_dapc <- coord_dapc %>%
  group_by(Localidad) %>%
  summarise(
    DA1 = mean(DA1, na.rm = TRUE),
    DA2 = mean(DA2, na.rm = TRUE),
    .groups = "drop"
  )


centroides_dapc


# Extraer coordenadas de LP
etiqueta_LP <- centroides_dapc %>%
  filter(Localidad == "LP")


############################################################
# Function to Generate the DAPC Plot
############################################################

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
  
  ######################################################
  # MAIN DAPC PLOT
  ######################################################
  
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
  
  ######################################################
  # SNP DATASET TITLE
  ######################################################
  
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
  
  ######################################################
  # AXIS TITLES
  ######################################################
  
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
  
  ######################################################
  # LP LABEL
  ######################################################
  
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
  
  ######################################################
  # DA EIGENVALUES INSET
  ######################################################
  
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


############################################################
# Display DAPC Plot
############################################################

graphics.off()

graficar_DAPC()


############################################################
# Display DAPC Plot
############################################################

graphics.off()

graficar_DAPC()


############################################################
# Create Output Directory
############################################################

output_dir <- file.path(
  getwd(),
  "Figuras_DAPC"
)


if(!dir.exists(output_dir)){
  
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


############################################################
# Export DAPC as PNG
############################################################

output_png <- file.path(
  output_dir,
  "DAPC_outliers_19.png"
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


############################################################
# Export DAPC as PDF
############################################################

output_pdf <- file.path(
  output_dir,
  "DAPC_outliers_19.pdf"
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


############################################################
# Display Output File Paths
############################################################

output_png

output_pdf


############################################################
# DAPC Density Plot (DA1)
############################################################

############################################################
# DAPC Density Plot (DA1)
############################################################

# Función para generar el gráfico de densidad
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

# Mostrar en R
graficar_DAPC_densidad()


############################################################
# Export Density Plot as PNG
############################################################

output_png_densidad <- file.path(
  output_dir,
  "DAPC_outliers_19_densidad_DA1.png"
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


############################################################
# Export Density Plot as PDF
############################################################

output_pdf_densidad <- file.path(
  output_dir,
  "DAPC_outliers_19_densidad_DA1.pdf"
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


############################################################
# Calculate Explained Variance
############################################################

var_explained <- slodapc$eig / sum(slodapc$eig) * 100

var_explained


############################################################
# Display Output File Paths
############################################################

output_png_densidad

output_pdf_densidad


############################################################
# fastStructure Analysis - FST Outlier Dataset
############################################################

############################################################
# Load and Align Q Files
############################################################

#load library for use
library(pophelper)
library(gridExtra)
library(pophelper)
library(paletteer)
library(adegenet)
library(vcfR)
library(dartR)
library(ggplot2)

setwd("E://Tesis final/ediciones paper/faststructure_19/")


qlist <- choose.files()

qlist <- readQ(qlist)

qlist <- alignK(qlist)


mergeQ2 <- function(qlist) {
  
  is.qlist(qlist)
  
  if(diff(range(as.integer(tabulateQ(qlist)$ind))) != 0) {
    
    stop(
      "mergeQ: Number of individuals differ between runs."
    )
  }
  
  # Computes mean cell-wise across dataframes
  # x A list of numeric dataframes
  
  mergy <- function(x) {
    
    return(
      list(
        Reduce(`+`, x) / length(x)
      )
    )
  }
  
  # if all runs have same K, merge as is
  if(diff(range(as.integer(tabulateQ(qlist)$k))) == 0) {
    
    labels <- summariseQ(
      tabulateQ(qlist)
    )$k
    
    x <- mergy(qlist)
    
    names(x) <- labels
    
  } else {
    
    # if runs have different K, split them and merge within sublists
    qlist <- sortQ2(qlist)
    
    labels <- summariseQ(
      tabulateQ(
        qlist,
        sorttable = FALSE
      )
    )$k
    
    x <- unlist(
      lapply(
        splitQ(qlist),
        mergy
      ),
      recursive = FALSE
    )
    
    names(x) <- labels
  }
  
  return(
    as.qlist(x)
  )
}


sortQ2 <- function(
    qlist,
    by = "k",
    decreasing = FALSE,
    debug = FALSE
) {
  
  is.qlist(qlist)
  
  if(length(by) == 0) {
    
    stop(
      "sortQ: Argument 'by' must not be length zero."
    )
  }
  
  if(!is.character(by)) {
    
    stop(
      "sortQ: Argument 'by' must be a character."
    )
  }
  
  if(!is.logical(decreasing)) {
    
    stop(
      "sortQ: Argument 'decreasing' must be a logical datatype."
    )
  }
  
  fun1 <- function(x) {
    
    as.matrix(
      unlist(
        attributes(x)
      )
    )
  }
  
  a <- lapply(
    qlist,
    fun1
  )
  
  if(debug) {
    
    print(a)
  }
  
  if(
    any(
      !sapply(
        a,
        function(x) {
          
          any(
            grepl(
              paste0(
                by,
                collapse = "|"
              ),
              rownames(x)
            )
          )
        }
      )
    )
  ) {
    
    stop(
      paste0(
        "One or more of the attributes provided in by (",
        by,
        ") is missing in one or more runs. ",
        "If 'ind' or 'k' is missing, use 'as.qlist()' to add them."
      )
    )
  }
  
  # get df of attributes
  b <- as.data.frame(
    t(
      as.data.frame(
        lapply(
          a,
          function(x, y) {
            
            x[y, ]
          },
          by
        ),
        stringAsFactors = FALSE
      )
    ),
    stringsAsFactors = FALSE
  )
  
  fun2 <- function(x) {
    
    if(
      all(
        !is.na(
          as.numeric(
            as.character(x)
          )
        )
      )
    ) {
      
      return(
        as.numeric(
          as.character(x)
        )
      )
      
    } else {
      
      return(x)
    }
  }
  
  b <- as.data.frame(
    sapply(
      b,
      fun2
    ),
    stringAsFactors = FALSE
  )
  
  if(debug) {
    
    print(str(b))
    
    print(b)
  }
  
  # order
  ord <- do.call(
    order,
    b[
      ,
      by,
      drop = FALSE
    ]
  )
  
  if(decreasing) {
    
    ord <- rev(ord)
  }
  
  # sort qlist
  return(
    qlist[ord]
  )
}


mergeQ2(
  qlist = qlist
) -> SLN_PLOT

length(SLN_PLOT)


############################################################
# Generate fastStructure Plot
############################################################

paletteer_d(
  "ggthemes::Classic_10_Light"
) -> SLN_PAL

SLN_PAL


data.frame(
  loc = c(
    rep("1.LP", 2),
    rep("2.LA", 2),
    rep("3.RF", 3),
    rep("4.MZ", 3),
    rep("5.CM", 11),
    rep("6.TP", 10),
    rep("7.CT", 3),
    rep("8.CA", 4)
  )
) -> Loc_sln


P1 <- plotQ(
  SLN_PLOT[5],
  
  # Información de las localidades
  grplab = Loc_sln,
  selgrp = "loc",
  ordergrp = TRUE,
  
  # Salida
  returnplot = TRUE,
  exportplot = FALSE,
  imgoutput = "sep",
  
  # Tamaño general: números del eje Y y otros textos
  basesize = 14,
  
  # Eje Y
  showyaxis = TRUE,
  showticks = TRUE,
  
  # No mostrar etiquetas individuales
  showindlab = FALSE,
  sharedindlab = FALSE,
  
  # Colores
  clustercol = SLN_PAL,
  
  # Leyenda superior
  showlegend = TRUE,
  legendtextsize = 25,
  legendkeysize = 25,
  
  # Etiquetas laterales K=2, K=3
  showsp = TRUE,
  sppos = "right",
  splabsize = 25,
  splabangle = -90,
  splabcol = "grey20",
  splabface = "plain",
  
  # Nombres de las localidades
  showgrplab = TRUE,
  grplabsize = 5,
  grplabangle = 0,
  grplabjust = 0.5,
  grplabcol = "grey20",
  grplabface = "plain",
  grplabheight = 1,
  grplabpos = 0.30,
  grplabspacer = 0.02,
  
  # Línea inferior
  linealpha = 1,
  linecol = "grey20",
  linesize = 0.8,
  linepos = 0.5,
  
  # Marcas que separan las localidades
  pointalpha = 1,
  pointcol = "grey20",
  pointsize = 4,
  
  # Separaciones verticales en las barras
  showdiv = TRUE,
  divcol = "white",
  divsize = 0.8,
  
  # Espacio entre los paneles K
  panelspacer = 0.5
)


grid.arrange(
  P1$plot[[1]],
  ncol = 1
)


# Export Figure
output_dir <- "E:/Tesis final/ediciones paper/Imagenes FINALES"


if (!dir.exists(output_dir)) {
  
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


output_png <- file.path(
  output_dir,
  "FastStructure_outliers_19_K5.png"
)


ggsave(
  filename = output_png,
  plot = P1$plot[[1]],
  width = 16,
  height = 4,
  units = "in",
  dpi = 300,
  bg = "white"
)


############################################################
# Export fastStructure Plot as PDF
############################################################

output_pdf <- file.path(
  output_dir,
  "FastStructure_outliers_19_K5.pdf"
)


ggsave(
  filename = output_pdf,
  plot = P1$plot[[1]],
  width = 16,
  height = 4,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)
