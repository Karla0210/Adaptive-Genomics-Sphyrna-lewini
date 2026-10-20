############################################################
# Neutral SNP Dataset
############################################################

rm(list = ls())
graphics.off()
gc()


############################################################
# Load required libraries
############################################################

library(dplyr)
library(adegenet)
library(vcfR)
library(paletteer)
library(dartR)


############################################################
# Set working directory
############################################################

setwd("E:\\Filtrado populations y VCF p7 r8/Filtrados_new/")


############################################################
# Load neutral SNP dataset
############################################################

read.vcfR("Resumen_SNPn2.recode.vcf") -> sln_vcf


# Convert VCF to genind
vcfR2genind(sln_vcf) -> sln_gi


# Convert genind to genlight
gi2gl(sln_gi) -> sln_gl


# Assign individuals to sampling localities
read.table("PopMap.txt") -> ind_loc

ind_loc

sln_gl$pop <- as.factor(ind_loc$V2)


############################################################
# Cross-validation to determine the optimal number of PCs
############################################################

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
find.clusters(
  sln_gl,
  max.n.clust = 8,
  stat = "BIC",
  n.iter = 100000,
  n.start = 10000
)


############################################################
# Perform DAPC
############################################################

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
dapc(sln_gl) -> slndapc


# Selected parameters:
# Number of principal components = 12
# Number of discriminant functions = 7


# Verify population order
popNames(sln_gl)


############################################################
# Calculate the percentage of explained variance
############################################################

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


############################################################
# Function to generate the DAPC figure
############################################################

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


############################################################
# Display DAPC plot
############################################################

graphics.off()

graficar_DAPC()


############################################################
# Create output directory
############################################################

output_dir <- file.path(
  getwd(),
  "Figuras_DAPC"
)


if (!dir.exists(output_dir)) {
  
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


############################################################
# Export DAPC as PDF
############################################################

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


############################################################
# Display output file paths
############################################################

output_png

output_pdf


############################################################
# DAPC density plot for DA1
############################################################

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
graficar_DAPC_densidad()


############################################################
# Export density plot as PNG
############################################################

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


############################################################
# Export density plot as PDF
############################################################

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


############################################################
# Calculate explained variance
############################################################

var_explained <- slndapc$eig /
  sum(slndapc$eig) * 100

var_explained


############################################################
# fastStructure Analysis - Neutral SNP Dataset
############################################################

# Load required libraries
library(pophelper)
library(gridExtra)
library(ggplot2)


# Set working directory
setwd(
  "E://Tesis final/ediciones paper/faststructure_neutros/"
)


############################################################
# Load and align fastStructure Q files
############################################################

qlist <- choose.files()

qlist <- readQ(qlist)

qlist <- alignK(qlist)


############################################################
# Function to merge replicate Q matrices
############################################################

mergeQ2 <- function(qlist) {
  
  is.qlist(qlist)
  
  if (
    diff(
      range(
        as.integer(
          tabulateQ(qlist)$ind
        )
      )
    ) != 0
  ) {
    
    stop(
      paste(
        "mergeQ: Number of individuals",
        "differs between runs."
      )
    )
  }
  
  
  # Compute cell-wise means across Q matrices
  mergy <- function(x) {
    
    return(
      list(
        Reduce(`+`, x) /
          length(x)
      )
    )
  }
  
  
  # Merge directly if all runs have the same K
  if (
    diff(
      range(
        as.integer(
          tabulateQ(qlist)$k
        )
      )
    ) == 0
  ) {
    
    labels <- summariseQ(
      tabulateQ(qlist)
    )$k
    
    x <- mergy(qlist)
    
    names(x) <- labels
    
  } else {
    
    # Split by K and merge replicates within each K
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


############################################################
# Function to sort Q matrices
############################################################

sortQ2 <- function(
  qlist,
  by = "k",
  decreasing = FALSE,
  debug = FALSE
) {
  
  is.qlist(qlist)
  
  if (length(by) == 0) {
    
    stop(
      "sortQ: Argument 'by' must not have length zero."
    )
  }
  
  if (!is.character(by)) {
    
    stop(
      "sortQ: Argument 'by' must be a character."
    )
  }
  
  if (!is.logical(decreasing)) {
    
    stop(
      paste(
        "sortQ: Argument 'decreasing'",
        "must be logical."
      )
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
  
  
  if (debug) {
    
    print(a)
  }
  
  
  if (
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
        "One or more attributes provided in 'by' (",
        by,
        ") are missing in one or more runs."
      )
    )
  }
  
  
  # Create data frame of attributes
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
        stringsAsFactors = FALSE
      )
    ),
    stringsAsFactors = FALSE
  )
  
  
  fun2 <- function(x) {
    
    if (
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
    stringsAsFactors = FALSE
  )
  
  
  if (debug) {
    
    print(str(b))
    
    print(b)
  }
  
  
  # Determine sorting order
  ord <- do.call(
    order,
    b[
      ,
      by,
      drop = FALSE
    ]
  )
  
  
  if (decreasing) {
    
    ord <- rev(ord)
  }
  
  
  # Return sorted Q list
  return(
    qlist[ord]
  )
}


############################################################
# Merge replicate runs
############################################################

mergeQ2(
  qlist = qlist
) -> SLN_PLOT

length(SLN_PLOT)


############################################################
# Generate fastStructure ancestry barplots
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
  SLN_PLOT[c(2, 3)],
  
  # Sampling locality information
  grplab = Loc_sln,
  selgrp = "loc",
  ordergrp = TRUE,
  
  # Output options
  returnplot = TRUE,
  exportplot = FALSE,
  imgoutput = "join",
  
  # General text size
  basesize = 14,
  
  # Y-axis
  showyaxis = TRUE,
  showticks = TRUE,
  
  # Hide individual labels
  showindlab = FALSE,
  sharedindlab = FALSE,
  
  # Cluster colors
  clustercol = SLN_PAL,
  
  # Top legend
  showlegend = TRUE,
  legendtextsize = 25,
  legendkeysize = 25,
  
  # K labels
  showsp = TRUE,
  sppos = "right",
  splabsize = 25,
  splabangle = -90,
  splabcol = "grey20",
  splabface = "plain",
  
  # Sampling locality labels
  showgrplab = TRUE,
  grplabsize = 5,
  grplabangle = 0,
  grplabjust = 0.5,
  grplabcol = "grey20",
  grplabface = "plain",
  grplabheight = 1,
  grplabpos = 0.30,
  grplabspacer = 0.02,
  
  # Bottom grouping line
  linealpha = 1,
  linecol = "grey20",
  linesize = 0.8,
  linepos = 0.5,
  
  # Locality separation marks
  pointalpha = 1,
  pointcol = "grey20",
  pointsize = 4,
  
  # Vertical dividers between individuals
  showdiv = TRUE,
  divcol = "white",
  divsize = 0.8,
  
  # Space between K panels
  panelspacer = 0.5
)


grid.arrange(
  P1$plot[[1]],
  ncol = 1
)


############################################################
# Export fastStructure plot
############################################################

output_dir <- paste0(
  "E://Tesis final/ediciones paper/",
  "faststructure_neutros/"
)


output_file <- file.path(
  output_dir,
  "FastStructure_plot_neutros.png"
)


ggsave(
  filename = output_file,
  plot = P1$plot[[1]],
  width = 16,
  height = 10,
  dpi = 300
)
