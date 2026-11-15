#!/usr/bin/env Rscript

# fastStructure plotting - consensus outliers (29 SNPs) | Sphyrna lewini

# 01. Packages
library(pophelper)
library(gridExtra)
library(paletteer)
library(ggplot2)

# 02. Input and output paths
Q_DIR <- "."
output_dir <- file.path(Q_DIR, "fastStructure_figures")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 03. Read and align fastStructure Q files
q_files <- list.files(Q_DIR, pattern = "\\.meanQ$", full.names = TRUE)
if (length(q_files) == 0) stop("No .meanQ files found in Q_DIR")
qlist <- readQ(q_files)
qlist <- alignK(qlist)

# 04. Merge replicate Q matrices
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

# Check that the selected element is the intended K = 5 result.
print(names(SLN_PLOT))
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
  
  # K panel labels
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

# Output directory created above

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
