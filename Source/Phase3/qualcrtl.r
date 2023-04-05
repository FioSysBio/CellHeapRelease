{\rtf1\ansi\ansicpg1252\cocoartf2636
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fnil\fcharset0 Menlo-Regular;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\paperw11900\paperh16840\margl1440\margr1440\vieww24680\viewh8400\viewkind0
\deftab20
\pard\tx577\tx1155\tx1733\tx2311\tx2889\tx3467\tx4045\tx4623\tx5201\tx5779\tx6357\tx6935\tx7513\tx8091\tx8669\tx9247\tx9825\tx10403\tx10981\tx11559\tx12137\tx12715\tx13293\tx13871\tx14449\tx15027\tx15605\tx16183\tx16761\tx17339\tx17917\tx18495\tx19072\tx19650\tx20228\tx20806\tx21384\tx21962\tx22540\tx23118\tx23696\tx24274\tx24852\tx25430\tx26008\tx26586\tx27164\tx27742\tx28320\tx28898\tx29476\tx30054\tx30632\tx31210\tx31788\tx32366\tx32944\tx33522\tx34100\tx34678\tx35256\tx35834\tx36412\tx36990\tx37567\tx38145\tx38723\tx39301\tx39879\tx40457\tx41035\tx41613\tx42191\tx42769\tx43347\tx43925\tx44503\tx45081\tx45659\tx46237\tx46815\tx47393\tx47971\tx48549\tx49127\tx49705\tx50283\tx50861\tx51439\tx52017\tx52595\tx53173\tx53751\tx54329\tx54907\tx55485\tx56062\tx56640\tx57218\tx57796\pardeftab20\li577\fi-577\partightenfactor0

\f0\fs24 \cf2 \CocoaLigature0 # Sample SRR11537946 BALF severe\
\
library(Seurat)\
library(dplyr)\
library(patchwork)\
library(scales)\
library(cowplot)\
library(ggplot2)\
library(RColorBrewer)\
library(gplots)\
library(SoupX)\
library(conflicted)\
\
# SoupX - Automatic mode using cellranger outputs \
# Obtain better results - Previous basic clustering information\
# Using 10X data mapped with cellranger \
# Default clustering produced by cellranger is automatically loaded and used\
# Load cellranger outputs, to create SoupChannel object, and to estime soup profile\
\
SRR11537946_severe_BALF_SpX <- load10X("/scratch/inova-covd19/PRJNA608742/Fase2/CRCD-SRR11537946/outs")\
\
# Estime the "rho" parameter that represents the fraction of contamination\
# rho = 0 means no contamination,  rho = 1 means 100% of UMIs into a droplet are "soup"\
\
SRR11537946_severe_BALF_SpX <- autoEstCont(SRR11537946_severe_BALF_SpX, doPlot=FALSE)\
\
#Generate the corrected matrix without contamination.\
\
out_SRR11537946_severe_BALF_SpX <- adjustCounts(SRR11537946_severe_BALF_SpX)\
\
# Using corrected matrix from SoupX directly in the Seurat through the CreateSeuratObject function\
\
SRR11537946_severe_BALF <- CreateSeuratObject(counts = out_SRR11537946_severe_BALF_SpX, project = "PRJNA608742_Severe_BALF",min.cells=1)\
\
#Calculate the percentage of mitochondrial RNA reads\
SRR11537946_severe_BALF[["percent.mt"]] <- PercentageFeatureSet(SRR11537946_severe_BALF, pattern = "^human----MT-")\
\
# Metrics\
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")\
\
png("severe_BALF_before_filters_SRR11537946.png")\
\
VlnPlot(SRR11537946_severe_BALF, features = feats, pt.size = 0.1, ncol = 3) + NoLegend()\
dev.off()\
\
# Step of quality control filters defined in previous analyzes based on Wauters et al, 2021\
# Cell filtering parameters for selecting cells based on number of genes/cell, UMI counts/cell, and percent mitochondrial genes according to Wauters et al, 2021.\
SRR11537946_severe_BALF_filtered <- subset(SRR11537946_severe_BALF, subset=nFeature_RNA>150 & nFeature_RNA<3000 & percent.mt<20 & nCount_RNA>301)\
\
png("severe_BALF_after_filters_SRR11537946.png")\
\
VlnPlot(SRR11537946_severe_BALF_filtered, features = feats, pt.size = 0.1, ncol = 3) + NoLegend()\
dev.off()\
\
# The next steps involve looking for viral RNA in the samples - it just kepts the commands from the previous phase 3 script\
# Export .tsv and .csv files\
\
out_test <- as.matrix(SRR11537946_severe_BALF_filtered@assays$RNA@counts)\
write.table(out_test, file="features_severe_BALF_SRR11537946.tsv", quote=FALSE, sep='\\t', col.names = TRUE)\
write.csv(SRR11537946_severe_BALF_filtered@meta.data, file="metadata_severe_BALF_SRR11537946.csv")\
\
# Import of seurat filter file\
sample <- read.table("features_severe_BALF_SRR11537946.tsv", sep="\\t", header=T, row.names = 1)\
\
# Select rows containing sarscov2\
sample_sarscov2 <- sample[grep("virus-v6", row.names(sample)),,drop=FALSE]\
\
# Transpose data frame\
sample_sarscov2_transposta <- t(sample_sarscov2)\
\
# Print rows where all columns are zero \
not_infected <- rownames(sample_sarscov2_transposta)[which(rowSums(sample_sarscov2_transposta)==0)]\
\
# Print rows where all columns are different of zero\
infected <- rownames(sample_sarscov2_transposta)[which(rowSums(sample_sarscov2_transposta)>0)]\
\
# Dataframe not infected \
write.table(sample[,not_infected,drop=FALSE], file= "features_not_infected_severe_BALF_SRR11537946.tsv", quote=FALSE, sep='\\t', col.names=TRUE)\
\
# Dataframe  infected \
write.table(sample[,infected,drop=FALSE], file= "features_infected_severe_BALF_SRR11537946.tsv", quote=FALSE, sep='\\t', col.names = TRUE)\
}