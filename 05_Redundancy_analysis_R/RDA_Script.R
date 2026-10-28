library("vegan")
library("psych")
library("ggplot2")
setwd("E:/Filtrado populations y VCF p7 r8/RDA/")
gen <- read.table("RDA_gen.raw", header = TRUE, sep = " ",row.names = 1)
dim(gen)

gen<-gen[,6:9563]
sum(is.na(gen))
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp)) # No NAs

env <- read.table("AMBIENTAL.txt", header = TRUE)
str(env) # Look at the structure of the data frame
env$Localidad <- as.factor(env$Localidad) 
env$ID <- as.character(env$ID)


# Check how many times each location appears before converting it into a factor
table(env$Localidad)
# Assign levels and convert to a factor
env$Localidad <- factor(env$Localidad, levels = c("LP", "LA", "RF", "MZ", "CM", "TP", "CT", "CA"))

# Check how the distribution of locations looks now
table(env$Localidad)


pairs.panels(env[,2:14], scale=T)
View(env)
#To remove variables
pred <- subset(env, select=-c(ID,Localidad))
pred <- subset(env, select=-c(ID,Localidad,windV,ClorA,Longitude))
#New
pairs.panels(pred, scale=T,cex.cor = 5)

###SAVE pairs panel###
library(psych)

png(
  filename = "E:/Tesis final/ediciones paper/pairs_panels.png",
  width = 3000,
  height = 3000,
  res = 300
)

pairs.panels(
  pred,
  scale = TRUE,
  cex.cor = 5
)

dev.off()

gen.imp <- gen.imp[, -(1:5)]  
View(gen.imp)


#Graph the RDA
wolf.rda <- rda(gen.imp ~ ., data = pred, scale = TRUE)
wolf.rda





#Adjust the model based on the number of predictors
RsquareAdj(wolf.rda)


#The eigenvalues of the restricted axes reflect the variance explained by each axis

summary(eigenvals(wolf.rda, model = "constrained"))
screeplot(wolf.rda)

#Check significance of RDA
signif.full <- anova.cca(wolf.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full

#Check significance of each axis
signif.axis <- anova.cca(wolf.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis

#Permutation test for rda under reduced model

vif.cca(wolf.rda)


#Plot RDA
plot(wolf.rda, scaling=3)          # default is axes 1 and 2
text(wolf.rda, scaling=3, display="bp", col="#0868ac", cex=1) 
plot(wolf.rda, choices = c(1, 3), scaling=3)  # axes 1 and 3

#We will color the individual points based on their ecotype, which we can find in the env dataset
levels(env$Localidad) <- c("LP","LA","RF","MZ","CM","TP","CT","CA")
eco <- env$Localidad
bg <- c("darkolivegreen1","mediumpurple1","#33a02c","plum1","orchid3","deeppink3","steelblue","aquamarine2" ) # 6 nice colors for our ecotypes

env$Localidad <- factor(env$Localidad)  
table(env$Localidad)  # To see how many categories it actually has


#Add components to the chart separately
# axes 1 & 2
plot(wolf.rda, type="n", scaling=3)
points(wolf.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(wolf.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # the wolves
text(wolf.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=0.8, pt.bg=bg)

# axes 1 & 3
plot(wolf.rda, type="n", scaling=3, choices=c(1,3))
points(wolf.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
points(wolf.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco], choices=c(1,3))
text(wolf.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("topleft", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

#POSITION IDENTIFICATION
load.rda <- scores(wolf.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes

hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
#SNPS NUMBER
cand1 <- outliers(load.rda[,1],3) # 19
cand1
cand2 <- outliers(load.rda[,2],3) # 28
cand2
cand3 <- outliers(load.rda[,3],3) # 17
cand3
ncand <- length(cand1) + length(cand2) + length(cand3)
ncand  #64

#To see if there are duplicate SNPs
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)


##
View(pred)
foo <- matrix(nrow=(ncand), ncol=9) #  9 columns for 9 predictors
colnames(foo) <- c("Latitude","TA","TM","Salinity","ClorM","PAR","windU","CurrentU","CurrentV")


for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen.imp[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
View(cand)



#Investigate to the candidates
length(cand$snp[duplicated(cand$snp)])  # 0 duplicate detections

#SNPS summary
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,13] <- names(which.max(abs(bar[4:12]))) # gives the variable
  cand[i,14] <- max(abs(bar[4:12]))              # gives the correlation
}

colnames(cand)[13] <- "predictor"
colnames(cand)[14] <- "correlation"

table(cand$predictor) 
##ClorM= 9, CurrentU= 7, CurrentV= 1, Latitude= 9, PAR= 1, Salinity= 20, TA= 4, TM= 3, windU= 10



############################################################
# PREPARE CANDIDATE SNP NAMES FOR PLOTTING
############################################################

# Remove allele suffixes from candidate SNP names
cand$snp_clean <- gsub(
  "_A$|_T$|_C$|_G$",
  "",
  cand$snp
)

# Create a cleaned version of genotype-column names
gen_snp_names_clean <- gsub(
  "_A$|_T$|_C$|_G$",
  "",
  colnames(gen.imp)
)

# Identify candidate SNPs in the genotype matrix
candidate_match <- gen_snp_names_clean %in% cand$snp_clean

# Check the number of candidate SNPs matched
sum(candidate_match)

# Check which candidate SNPs were not found
unmatched_candidates <- setdiff(
  cand$snp_clean,
  gen_snp_names_clean
)

unmatched_candidates


############################################################
# DEFINE SNP COLORS
############################################################

# Background SNPs are light gray;
# candidate SNPs are highlighted in red with a yellow border
snp_border <- ifelse(
  candidate_match,
  "yellow",
  "gray32"
)

snp_fill <- ifelse(
  candidate_match,
  "red",
  "#F1EEF6"
)

snp_size <- ifelse(
  candidate_match,
  1.25,
  0.70
)


############################################################
# DEFINE LOCALITY COLORS
############################################################

env$Localidad <- factor(
  env$Localidad,
  levels = c(
    "LP", "LA", "RF", "MZ",
    "CM", "TP", "CT", "CA"
  )
)

locality_colors <- c(
  LP = "#C7323E",
  LA = "#F3724B",
  RF = "#FCC55A",
  MZ = "#66BC7D",
  CM = "#327146",
  TP = "#70C2DA",
  CT = "#3598C4",
  CA = "#316AAF"
)

individual_colors <- locality_colors[
  as.character(env$Localidad)
]


############################################################
# FUNCTION TO GENERATE THE FINAL RDA PLOT
############################################################

plot_final_RDA <- function() {
  
  par(
    mar = c(5, 5, 3, 8),
    xpd = TRUE,
    cex.axis = 1.2,
    cex.lab = 1.3
  )
  
  # Empty plotting area for RDA axes 1 and 2
  plot(
    wolf.rda,
    type = "n",
    scaling = 3,
    choices = c(1, 2),
    main = "Environment-associated SNPs",
    ylim = c(-3, 6)
  )
  
  # Plot background and candidate SNPs
  points(
    wolf.rda,
    display = "species",
    pch = 21,
    cex = snp_size,
    col = snp_border,
    bg = snp_fill,
    scaling = 3,
    choices = c(1, 2)
  )
  
  # Plot environmental predictors
  text(
    wolf.rda,
    display = "bp",
    col = "#0868AC",
    cex = 1,
    scaling = 3,
    choices = c(1, 2)
  )
  
  # Plot individuals according to sampling locality
  points(
    wolf.rda,
    display = "sites",
    pch = 21,
    cex = 1.3,
    col = "gray32",
    bg = individual_colors,
    scaling = 3,
    choices = c(1, 2)
  )
  
  # Locality legend
  legend(
    "topright",
    inset = c(-0.15, 0),
    legend = names(locality_colors),
    bty = "n",
    col = "gray32",
    pch = 21,
    pt.bg = locality_colors,
    cex = 0.9,
    pt.cex = 1.2
  )
}


############################################################
# DISPLAY THE FINAL RDA PLOT IN R
############################################################

graphics.off()

plot_final_RDA()


############################################################
# OUTPUT DIRECTORY
############################################################

output_dir <- "E:/Tesis final/ediciones paper"

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


############################################################
# SAVE AS PNG
############################################################

output_png <- file.path(
  output_dir,
  "RDA_environment_associated_SNPs.png"
)

png(
  filename = output_png,
  width = 16,
  height = 10,
  units = "in",
  res = 300,
  bg = "white"
)

plot_final_RDA()

dev.off()


############################################################
# SAVE AS PDF
############################################################

output_pdf <- file.path(
  output_dir,
  "RDA_environment_associated_SNPs.pdf"
)

pdf(
  file = output_pdf,
  width = 16,
  height = 10,
  family = "sans",
  useDingbats = FALSE,
  bg = "white"
)

plot_final_RDA()

dev.off()


############################################################
# SAVE CANDIDATE SNP LIST
############################################################

output_candidates <- file.path(
  output_dir,
  "RDA_candidate_SNPs.csv"
)

write.csv(
  cand,
  file = output_candidates,
  row.names = FALSE
)


############################################################
# DISPLAY OUTPUT PATHS
############################################################

output_png
output_pdf
output_candidates
