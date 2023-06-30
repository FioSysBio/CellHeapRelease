# Sample SRR11537946 BALF severe

library(Seurat)
library(dplyr)
library(patchwork)
library(scales)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(SoupX)
library(conflicted)

# SoupX - Automatic mode using cellranger outputs 
# Obtain better results - Previous basic clustering information
# Using 10X data mapped with cellranger 
# Default clustering produced by cellranger is automatically loaded and used
# Load cellranger outputs, to create SoupChannel object, and to estime soup profile

SRR11537946_severe_BALF_SpX <- load10X("/scratch/inova-covd19/PRJNA608742/Fase2/CRCD-SRR11537946/outs")

# Estime the "rho" parameter that represents the fraction of contamination
# rho = 0 means no contamination,  rho = 1 means 100% of UMIs into a droplet are "soup"

SRR11537946_severe_BALF_SpX <- autoEstCont(SRR11537946_severe_BALF_SpX, doPlot=FALSE)

#Generate the corrected matrix without contamination.

out_SRR11537946_severe_BALF_SpX <- adjustCounts(SRR11537946_severe_BALF_SpX)

# Using corrected matrix from SoupX directly in the Seurat through the CreateSeuratObject function

SRR11537946_severe_BALF <- CreateSeuratObject(counts = out_SRR11537946_severe_BALF_SpX, project = "PRJNA608742_Severe_BALF",min.cells=1)

#Calculate the percentage of mitochondrial RNA reads
SRR11537946_severe_BALF[["percent.mt"]] <- PercentageFeatureSet(SRR11537946_severe_BALF, pattern = "^human----MT-")

# Metrics
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

png("severe_BALF_before_filters_SRR11537946.png")

VlnPlot(SRR11537946_severe_BALF, features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
dev.off()

# Step of quality control filters defined in previous analyzes based on Wauters et al, 2021
# Cell filtering parameters for selecting cells based on number of genes/cell, UMI counts/cell, and percent mitochondrial genes according to Wauters et al, 2021.
SRR11537946_severe_BALF_filtered <- subset(SRR11537946_severe_BALF, subset=nFeature_RNA>150 & nFeature_RNA<3000 & percent.mt<20 & nCount_RNA>301)

png("severe_BALF_after_filters_SRR11537946.png")

VlnPlot(SRR11537946_severe_BALF_filtered, features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
dev.off()

# The next steps involve looking for viral RNA in the samples - it just kepts the commands from the previous phase 3 script
# Export .tsv and .csv files

out_test <- as.matrix(SRR11537946_severe_BALF_filtered@assays$RNA@counts)
write.table(out_test, file="features_severe_BALF_SRR11537946.tsv", quote=FALSE, sep='\\t', col.names = TRUE)
write.csv(SRR11537946_severe_BALF_filtered@meta.data, file="metadata_severe_BALF_SRR11537946.csv")

# Import of seurat filter file
sample <- read.table("features_severe_BALF_SRR11537946.tsv", sep="\\t", header=T, row.names = 1)

# Select rows containing sarscov2
sample_sarscov2 <- sample[grep("virus-v6", row.names(sample)),,drop=FALSE]

# Transpose data frame
sample_sarscov2_transposta <- t(sample_sarscov2)

# Print rows where all columns are zero 
not_infected <- rownames(sample_sarscov2_transposta)[which(rowSums(sample_sarscov2_transposta)==0)]

# Print rows where all columns are different of zero
infected <- rownames(sample_sarscov2_transposta)[which(rowSums(sample_sarscov2_transposta)>0)]

# Dataframe not infected 
write.table(sample[,not_infected,drop=FALSE], file= "features_not_infected_severe_BALF_SRR11537946.tsv", quote=FALSE, sep='\\t', col.names=TRUE)

# Dataframe  infected 
write.table(sample[,infected,drop=FALSE], file= "features_infected_severe_BALF_SRR11537946.tsv", quote=FALSE, sep='\\t', col.names = TRUE)

