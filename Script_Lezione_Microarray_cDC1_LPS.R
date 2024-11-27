#------------------------------------------------------------------------------#
# Pre-processing Microarray Data - QC, Normalization, DEGs and Visualization
#------------------------------------------------------------------------------#
# Analysis of datasets from: GSE203450 (Microarray Affymetric Data from DCs)


# Set Working Directory.--------------------------------------------------------
setwd("") # Set the Working Directory on your computer.

# Part 1. Load required libraries.==============================================
print("loading libraries...")

# Install Bioconductor Manager (needed for some packages).----------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install the required libraries.-----------------------------------------------
# BiocManager::install("GEOquery")
library(GEOquery)

# BiocManager::install("limma")
library(limma)

# Part 2. Import and explore data from GEO.=====================================
print("Importing data...")

id <- "GSE203450"
gse <- getGEO(id, GSEMatrix =TRUE, AnnotGPL=TRUE)

# Get some information about the file.----------------------------------------
list(gse)
names(pData(gse[[1]]))      # print the sample info.
length(gse)                 # check how many platforms was used.

# Select the dataset.---------------------------------------------------------
gse <- gse[[1]]

# if more than one dataset is present, you can analyse the other dataset by 
# changing the number inside the [[...]]
# e.g. gse <- gse[[2]]

# We can also use that kind of line.
# if (length(gse) > 1) idx <- grep("GPL6246", attr(gse, "names")) else idx <- 1


# Explore the expression dataframe.
pData(gse)                    # print the sample info.
fData(gse)                    # print the gene annotation.
exprs(gse)                    # print the expression data.


# Part 3. Group membership for all samples.=====================================
print("Grouping the samples...")

# Select one comparison between groups.
cdc1_lps <- "XXXXXX111000XXXXXX"          # (cDC1 LPS vs cDC1 Ctrl)

# Split Samples.----------------------------------------------------------------
sml <- strsplit(cdc1_lps, split ="")[[1]]


# Part 4. Filter out excluded samples and Create the expression matrix.=========
print("Filtering out excluded samples...")
sel <- which(sml != "X")   # excluded samples marked as "X".
sml <- sml[sel]
gse <- gse[ ,sel]

head(gse)

# Create the expression matrix.-----------------------------------------------
ex <- exprs(gse)
ex[which(ex <= 0)] <- NaN

head(ex)
hist(ex)

# Part 5. Log2 transformation and Normalization of the data.====================
print("Transforming and Normalizing data...")

# Box-and-whisker plot (all samples before normalization).--------------------
par(mar=c(7,4,2,1))
title <- paste ("GSE203450", "/", annotation(gse), sep ="")

# pdf("Microarray_data_before_transformation.pdf")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=F, las=2)
# dev.off()

# Log2 transformation.----------------------------------------------------------
exprs(gse) <- log2(ex)       # Log2 transform.

# Box-and-whisker plot (after log2 transformation).-----------------------------
par(mar=c(7,4,2,1))
title <- paste ("GSE203450", "/", annotation(gse), sep ="")

# pdf("Microarray_data_after_transformation.pdf")
boxplot(exprs(gse), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
# dev.off()

# Expression value distribution plot.-------------------------------------------
par(mar=c(4,4,2,1))
title <- paste ("GSE203450", "/", annotation(gse), " value distribution", 
                sep ="")

# pdf("Densities_after_normalization.pdf")
plotDensities(ex, main=title, legend=F)
# dev.off()

# Part 6. Normalization log-ratios with Limma.==================================
print("Normalizing data with Limma...")
exprs(gse) <- normalizeBetweenArrays(exprs(gse))        # normalize data
gse

plotDensities(gse, main=title, legend=F)


# Part 7. Assign samples to groups and set up design matrix.====================

# Define groups to compare.-----------------------------------------------------
groups <- make.names(c("cDC1_LPS","cDC1_ctrl"))

gs <- factor(sml)
levels(gs) <- groups
gse$group <- gs
design <- model.matrix(~group + 0, gse)
colnames(design) <- levels(gs)

gse <- gse[complete.cases(exprs(gse)), ] # skip missing values

# UMAP plot (dimensionality reduction).---------------------------------------
# library("maptools")  # point labels without overlaps
library("umap")
library("car")

ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates

ump <- umap(t(ex), n_neighbors = 3, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=3", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# Stop Point.
################################################################################