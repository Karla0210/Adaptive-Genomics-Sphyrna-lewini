#!/usr/bin/env Rscript

# fastStructure visualization - Sphyrna lewini
# FST-based outlier SNPs (19 SNPs)

library(pophelper)
library(gridExtra)
library(paletteer)
library(ggplot2)

# 01. Input and output paths
PROJECT_DIR <- "."
Q_DIR <- file.path(PROJECT_DIR, "faststructure_19")
OUTPUT_DIR <- file.path(PROJECT_DIR, "faststructure_figures")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 02. Read and align Q matrices
q_files <- list.files(Q_DIR, pattern = "\\.meanQ$", full.names = TRUE)
if (length(q_files) == 0) stop("No .meanQ files found in: ", Q_DIR)
qlist <- readQ(q_files)
qlist <- alignK(qlist)

# 03. Average replicate Q matrices by K
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

SLN_PLOT <- mergeQ2(qlist)
print(names(SLN_PLOT))

# 04. Sampling localities and cluster colors
SLN_PAL <- paletteer_d("ggthemes::Classic_10_Light")
Loc_sln <- data.frame(loc = c(
  rep("1.LP", 2), rep("2.LA", 2), rep("3.RF", 3),
  rep("4.MZ", 3), rep("5.CM", 11), rep("6.TP", 10),
  rep("7.CT", 3), rep("8.CA", 4)
))

# 05. Plot and export selected K
# Original code selected the fifth element. Confirm it corresponds to K = 5.
P1 <- plotQ(
  SLN_PLOT[5],
  
  grplab = Loc_sln,
  selgrp = "loc",
  ordergrp = TRUE,
  
  returnplot = TRUE,
  exportplot = FALSE,
  imgoutput = "sep",
  
  basesize = 14,
  
  showyaxis = TRUE,
  showticks = TRUE,
  
  showindlab = FALSE,
  sharedindlab = FALSE,
  
  clustercol = SLN_PAL,
  
  showlegend = TRUE,
  legendtextsize = 25,
  legendkeysize = 25,
  
  showsp = TRUE,
  sppos = "right",
  splabsize = 25,
  splabangle = -90,
  splabcol = "grey20",
  splabface = "plain",
  
  showgrplab = TRUE,
  grplabsize = 5,
  grplabangle = 0,
  grplabjust = 0.5,
  grplabcol = "grey20",
  grplabface = "plain",
  grplabheight = 1,
  grplabpos = 0.30,
  grplabspacer = 0.02,
  
  linealpha = 1,
  linecol = "grey20",
  linesize = 0.8,
  linepos = 0.5,
  
  pointalpha = 1,
  pointcol = "grey20",
  pointsize = 4,
  
  showdiv = TRUE,
  divcol = "white",
  divsize = 0.8,
  
  panelspacer = 0.5
)

ggsave(file.path(OUTPUT_DIR, "FastStructure_outliers_19_K5.png"),
       plot = P1$plot[[1]], width = 16, height = 4,
       units = "in", dpi = 300, bg = "white")

ggsave(file.path(OUTPUT_DIR, "FastStructure_outliers_19_K5.pdf"),
       plot = P1$plot[[1]], width = 16, height = 4,
       units = "in", device = cairo_pdf, bg = "white")
