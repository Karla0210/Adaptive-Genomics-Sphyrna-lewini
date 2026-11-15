#!/usr/bin/env Rscript

# fastStructure plots - neutral SNPs | Sphyrna lewini
# Packages: pophelper, gridExtra, ggplot2, paletteer

# 01. Input and output paths
INPUT_Q_DIR <- "."  # Folder containing fastStructure .meanQ files
OUTPUT_DIR <- file.path(INPUT_Q_DIR, "plots")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# 02. Packages and data
# Faststructure analysis - neutral snp dataset
library(pophelper)
library(gridExtra)
library(ggplot2)
library(paletteer)


# Load and align faststructure q files
# Replace interactive file selection with reproducible file discovery.
q_files <- list.files(INPUT_Q_DIR, pattern = "\\.meanQ$", full.names = TRUE)
if (length(q_files) == 0L) stop("No .meanQ files found in INPUT_Q_DIR")
qlist <- q_files

qlist <- readQ(qlist)

qlist <- alignK(qlist)


# Function to merge replicate q matrices
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


# Function to sort q matrices
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


# Merge replicate runs
mergeQ2(
  qlist = qlist
) -> SLN_PLOT

length(SLN_PLOT)


# Generate faststructure ancestry barplots
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


# IMPORTANT: Confirm names(SLN_PLOT) and the intended K values.
# Original script used positions 2 and 3; preserved here.
print(names(SLN_PLOT))
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


if (interactive()) {
  grid.arrange(P1$plot[[1]], ncol = 1)
}


# Export faststructure plot
output_dir <- OUTPUT_DIR


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
