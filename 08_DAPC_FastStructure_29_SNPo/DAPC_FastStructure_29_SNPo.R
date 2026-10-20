################ 29 OUTLIER SNPs ################

rm(list=ls())     
graphics.off()    
gc()


# Required libraries 
library(dplyr)
library(adegenet)
library(vcfR)
library(paletteer)
library(dartR)


# Set working directory 
setwd("E:\\Filtrado populations y VCF p7 r8/Filtrados_new/")


# Load data
# For outliers
read.vcfR("Resumen_SNPo2.recode.vcf") -> slo_vcf 


# Convert VCF file to genind
vcfR2genind(slo_vcf) -> slo_gi


# Convert to genlight
gi2gl(slo_gi) -> slo_gl


# Assign individuals to their locality
read.table("PopMap.txt") -> ind_loc

ind_loc

slo_gl$pop <- as.factor(ind_loc$V2)


# Validate the DAPC for outliers
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


# Components to be used: 6
find.clusters(
  slo_gl,
  max.n.clust = 8,
  stat = "BIC",
  n.iter = 100000,
  n.start = 10000
)


################ PERFORM DAPC FOR OUTLIERS ################


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
dapc(slo_gl) -> slodapc

# Select:
# Number of PCs = 6
# Discriminant functions = 6


# Verify population order
popNames(slo_gl)


################ AXIS PERCENTAGES ################


# Percentage of discriminant variation
var_explained <- slodapc$eig / sum(slodapc$eig) * 100

var_explained


# Percentage of DA1 and DA2
DA1_porcentaje <- round(var_explained[1], 2)

DA2_porcentaje <- round(var_explained[2], 2)


DA1_porcentaje

DA2_porcentaje


################ COORDINATES FOR CA, LP, AND CT ################


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




################ PLOTTING FUNCTION ################

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
  
  ######################################################
  # DA EIGENVALUES INSET
  ######################################################
  
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

################ DISPLAY THE PLOT IN R ################


graphics.off()

graficar_DAPC()


################ OUTPUT DIRECTORY ################


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


################ SAVE AS PNG ################


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


################ SAVE AS PDF ################


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


################ DISPLAY FILE PATHS ################


output_png

output_pdf




################ DAPC AXIS 1: DENSITY PLOT ################

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


################ SAVE DENSITY PLOT AS PNG ################

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


################ SAVE DENSITY PLOT AS PDF ################

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


################ PERCENTAGE OF EXPLAINED VARIATION ################

var_explained <- slodapc$eig / sum(slodapc$eig) * 100

var_explained


################ DISPLAY PATHS ################

output_png_densidad

output_pdf_densidad

##################### FASTSTRUCTURE 29 ###############################

## RUN STRUCTURE 
# Load required libraries
library(pophelper)
library(gridExtra)
library(pophelper)
library(paletteer)
library(adegenet)
library(vcfR)
library(dartR)
library(ggplot2)
setwd("E://Tesis final/ediciones paper/faststructure_29/")

qlist<-choose.files()
qlist<-readQ(qlist)
qlist<-alignK(qlist)


mergeQ2 <- function(qlist) {
  
  is.qlist(qlist)
  if(diff(range(as.integer(tabulateQ(qlist)$ind)))!=0) stop("mergeQ: Number of individuals differ between runs.")
  
  # Computes mean cell-wise across dataframes
  # x A list of numeric dataframes
  # 
  mergy <- function(x) {
    return(list(Reduce(`+`, x)/length(x)))
  }
  
  # If all runs have the same K, merge as is
  if(diff(range(as.integer(tabulateQ(qlist)$k)))==0) {
    labels <- summariseQ(tabulateQ(qlist))$k
    x <- mergy(qlist)
    names(x) <- labels
  }else{
    # If runs have different K values, split them and merge within sublists
    qlist <- sortQ2(qlist)
    labels <- summariseQ(tabulateQ(qlist,sorttable=FALSE))$k
    x <- unlist(lapply(splitQ(qlist),mergy),recursive=FALSE)
    names(x) <- labels
  }
  
  return(as.qlist(x))
}

sortQ2 <- function(qlist,by="k",decreasing=FALSE,debug=FALSE) {
  
  is.qlist(qlist)
  if(length(by)==0) stop("sortQ: Argument 'by' must not be length zero.")
  if(!is.character(by)) stop("sortQ: Argument 'by' must be a character.")
  if(!is.logical(decreasing)) stop("sortQ: Argument 'decreasing' must be a logical datatype.")
  
  fun1 <- function(x) as.matrix(unlist(attributes(x)))
  a <- lapply(qlist,fun1)
  if(debug) print(a)
  if(any(!sapply(a,function(x) any(grepl(paste0(by,collapse="|"),rownames(x)))))) {
    stop(paste0("One or more of the attributes provided in by (",by,") is missing in one or more runs. If 'ind' or 'k' is missing, use 'as.qlist()' to add them."))
  }
  
  # Get attribute data frame
  b <- as.data.frame(t(as.data.frame(lapply(a,function(x,y) x[y,],by),stringAsFactors=FALSE)),stringsAsFactors=FALSE)
  fun2 <- function(x) if(all(!is.na(as.numeric(as.character(x))))) {return(as.numeric(as.character(x)))}else{return(x)}
  b <- as.data.frame(sapply(b,fun2),stringAsFactors=FALSE)
  
  if(debug) {print(str(b)); print(b)}
  
  # Order
  ord <- do.call(order,b[,by,drop=FALSE])
  if(decreasing) ord <- rev(ord)
  # Sort qlist
  return(qlist[ord])
}

mergeQ2(qlist = qlist)->SLN_PLOT
length(SLN_PLOT)




## Plot
paletteer_d("ggthemes::Classic_10_Light")->SLN_PAL
SLN_PAL
data.frame(loc=c(rep("1.LP",2),rep("2.LA",2),rep("3.RF",3),
                 rep("4.MZ",3),rep("5.CM",11),rep("6.TP",10),rep("7.CT",3),rep("8.CA",4))) -> Loc_sln

P1 <- plotQ(
  SLN_PLOT[5],
  
  # Locality information
  grplab = Loc_sln,
  selgrp = "loc",
  ordergrp = TRUE,
  
  # Output
  returnplot = TRUE,
  exportplot = FALSE,
  imgoutput = "sep",
  
  # General size: Y-axis numbers and other text
  basesize = 14,
  
  # Y-axis
  showyaxis = TRUE,
  showticks = TRUE,
  
  # Do not display individual labels
  showindlab = FALSE,
  sharedindlab = FALSE,
  
  # Colors
  clustercol = SLN_PAL,
  
  # Top legend
  showlegend = TRUE,
  legendtextsize = 25,
  legendkeysize = 25,
  
  # Side labels K=2, K=3
  showsp = TRUE,
  sppos = "right",
  splabsize = 25,
  splabangle = -90,
  splabcol = "grey20",
  splabface = "plain",
  
  # Locality names
  showgrplab = TRUE,
  grplabsize = 5,
  grplabangle = 0,
  grplabjust = 0.5,       # Centers each name within its locality
  grplabcol = "grey20",
  grplabface = "plain",
  grplabheight = 1,
  grplabpos = 0.30,
  grplabspacer = 0.02,
  
  # Bottom line
  linealpha = 1,
  linecol = "grey20",
  linesize = 0.8,
  linepos = 0.5,
  
  # Marks separating the localities
  pointalpha = 1,
  pointcol = "grey20",
  pointsize = 4,
  
  # Vertical divisions between bars
  showdiv = TRUE,
  divcol = "white",
  divsize = 0.8,
  
  # Space between K panels
  panelspacer = 0.5
)

grid.arrange(P1$plot[[1]], ncol = 1)


# Save image 

output_dir <- "E:/Tesis final/ediciones paper/Imagenes FINALES"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_png <- file.path(
  output_dir,
  "FastStructure_outliers_29_K5.png"
)

## PDF
output_pdf <- file.path(
  output_dir,
  "FastStructure_outliers_29_K5.pdf"
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


ggsave(
  filename = output_png,
  plot = P1$plot[[1]],
  width = 16,
  height = 4,
  units = "in",
  dpi = 300,
  bg = "white"
)
```
