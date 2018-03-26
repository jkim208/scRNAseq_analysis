# prepareSeuratObj.R
# Updated: 26-March-2018
# Setup Robjs for cross-technology analyses of Droplet-based Single-Cell RNA Sequencing protocols
# Seurat v2.2.1
# Filter parameters (2018): min.cells = 3, min.genes = 200

library(dplyr)
library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4')
setwd("~/Slim")

##############################################################################################################
# Template for making a knee plot from cell read counts and determing the number of real cells collected
a=read.table("~/Slim/ddSeq/Seq3-N707_UCH1/Input/Sassi-Seq3-N707_cell_readcounts.txt", 
             header=F, stringsAsFactors=F) ##
x=cumsum(a$V1)
x=x/max(x)
plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]",
     ylab="cumulative fraction of reads", xlim=c(1,500))
abline(v=200) ##
title(main="ddSeq Seq3-N707 UCH1 UCH1 Knee Plot")
##############################################################################################################

# Function prepSeurat: takes raw count data from a droplet scRNAseq experiment and makes a Seurat object
# Input:  - raw count data from HUMAN CELLS/GENES ONLY (mixed populations should be cleaned manually)
#         - name of the technology used
#         - name of the project/sample
# To access saved files from prepSeurat, use the following example:
# ddSeq_UCH1_Seq2_N706 <- readRDS("Robj/ddSeq/ddSeq_UCH1_Seq2_N706_v2.1.Robj")

prepSeurat <- function(DGE.data, technology, name) {
  # Create Seurat Object
  seuratObj <- CreateSeuratObject(raw.data = DGE.data, min.cells = 3, min.genes = 200, project = name)
  seuratObj <- NormalizeData(object = seuratObj,normalization.method = "LogNormalize", scale.factor = 10000)
  seuratObj <- FindVariableGenes(object = seuratObj, mean.function = ExpMean, dispersion.function = LogVMR, 
                                 do.plot = FALSE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  # Keep track of mitochondrial representation 
  mito.genes <- grep(pattern = "[Mm][Tt]-", x = rownames(x = seuratObj@data), value = TRUE)
  percent.mito <- Matrix::colSums(seuratObj@raw.data[mito.genes, ])/Matrix::colSums(seuratObj@raw.data)
  seuratObj <- AddMetaData(object = seuratObj, metadata = percent.mito, col.name = "percent.mito")
  # Save Seurat object to target destination
  destination <- paste('Robj/', technology, '/', name, '.Robj', sep="")
  saveRDS(seuratObj, destination)
  print(paste("Seurat object was saved to:", destination, sep=" "))
}

#######################################################
# ddSeq
#######################################################

## ddSeq_293_Seq1_N701 and ddSeq_293_Seq2_N703 needed manual adjustment (mouse removed)

ddSeq_UCH2_Seq1_N704 <- read.table(file = "ddSeq/Seq1-N704_UCH2/Input/Sassi-Seq1-N704_DGE.txt",
                                       header = TRUE, row.names = 1)
ddSeq_UCH1_Seq2_N706 <- read.table(file = "ddSeq/Seq2-N706_UCH1/Input/Sassi-Seq2-N706_DGE.txt",
                                       header = TRUE, row.names = 1)
ddSeq_CH19_Seq3_N706 <- read.table(file = "ddSeq/Seq3-N706_CH19/Input/Sassi-Seq3-N706_DGE.txt",
                                       header = TRUE, row.names = 1)
ddSeq_UCH1_Seq3_N707 <- read.table(file = "ddSeq/Seq3-N707_UCH1/Input/Sassi-Seq3-N707_DGE.txt",
                                   header = TRUE, row.names = 1)

prepSeurat(ddSeq_UCH2_Seq1_N704, 'ddSeq', 'ddSeq_UCH2_Seq1_N704')
prepSeurat(ddSeq_UCH1_Seq2_N706, 'ddSeq', 'ddSeq_UCH1_Seq2_N706')
prepSeurat(ddSeq_CH19_Seq3_N706, 'ddSeq', 'ddSeq_CH19_Seq3_N706')
prepSeurat(ddSeq_UCH1_Seq3_N707, 'ddSeq', 'ddSeq_UCH1_Seq3_N707')

#######################################################
# DropSeq
#######################################################

## DS9, DS34, DS45 needed manual adjustment (mouse removed)

DropSeq_UCH2_DS49 <- read.table(file = "DropSeq/DS49_UCH2/server/DS49-UCH2-Seq1_DGE.txt",
                                header = TRUE, row.names = 1)
DropSeq_UCH1_DS50 <- read.table(file = "DropSeq/DS50_UCH1/server/DS50-UCH1_DGE.txt",
                                   header = TRUE, row.names = 1)
DropSeq_UCH1_RT_ID_DS52 <- read.table(file = "DropSeq/DS52_UCH1/server/DS52-RT_ID__DGE.txt",
                                header = TRUE, row.names = 1)
DropSeq_UCH1_HS_TB_DS52 <- read.table(file = "DropSeq/DS52_UCH1/server/DS52-HS_TB__DGE.txt",
                                      header = TRUE, row.names = 1)

prepSeurat(DropSeq_UCH2_DS49, 'DropSeq', 'DropSeq_UCH2_DS49')
prepSeurat(DropSeq_UCH1_DS50, 'DropSeq', 'DropSeq_UCH1_DS50')
prepSeurat(DropSeq_UCH1_RT_ID_DS52, 'DropSeq', 'DropSeq_UCH1_RT_ID_DS52')
prepSeurat(DropSeq_UCH1_HS_TB_DS52, 'DropSeq', 'DropSeq_UCH1_HS_TB_DS52')

#######################################################
# 10X (Droplike-protocol)
#######################################################
# 293HEK data was overloaded; cannot use until another experiment is performed
# Variable name cannot start with a number (10).

t10X_UCH1_s1 <- read.table(file = "10X/sample1_UCH1/droplike_pipeline/s1_UCH1_10X_DGE.txt", header = TRUE, row.names = 1)
t10X_UCH2_s2 <- read.table(file = "10X/sample2_UCH2/droplike_pipeline/s2_UCH2_10X_DGE.txt", header = TRUE, row.names = 1)

prepSeurat(t10X_UCH1_s1, '10X', '10X_UCH1_s1_droplike')
prepSeurat(t10X_UCH2_s2, '10X', '10X_UCH2_s2_droplike')

#######################################################
# inDrop
#######################################################
# 293HEK: "~/Slim/inDrop/input/iD_293HEK_lib23_human.Robj"
inDrop_UCH1_lib1 <- read.table(file = "inDrop/input/lib1_DGE.txt", header = TRUE, row.names = 1)
inDrop_UCH2_lib2 <- read.table(file = "inDrop/input/lib2_DGE.txt", header = TRUE, row.names = 1)

prepSeurat(inDrop_UCH1_lib1, 'inDrop', 'inDrop_UCH1_lib1')
prepSeurat(inDrop_UCH2_lib2, 'inDrop', 'inDrop_UCH2_lib2')

#######################################################
# Chordoma Tumor: Nadia-DropSeq
#######################################################
# Prepare Seurat objects of inDrop
SSa40_A_chord_low <- read.table(file = "SSa40_ChordTumor/DGE/SSa40-A_chord_low_DGE.txt", header = TRUE, row.names = 1)
SSa40_B_chord_low <- read.table(file = "SSa40_ChordTumor/DGE/SSa40-B_chord_low_DGE.txt", header = TRUE, row.names = 1)
SSa40_C_chord_high <- read.table(file = "SSa40_ChordTumor/DGE/SSa40-C_chord_high_DGE.txt", header = TRUE, row.names = 1)
SSa40_D_chord_high <- read.table(file = "SSa40_ChordTumor/DGE/SSa40-D_chord_high_DGE.txt", header = TRUE, row.names = 1)

prepSeurat(SSa40_A_chord_low, 'SSa40', 'SSa40_A_chord_low')
prepSeurat(SSa40_B_chord_low, 'SSa40', 'SSa40_B_chord_low')
prepSeurat(SSa40_C_chord_high, 'SSa40', 'SSa40_C_chord_high')
prepSeurat(SSa40_D_chord_high, 'SSa40', 'SSa40_D_chord_high')

##############################################################################################################

#######################################################
# scRNA-seq QC
#######################################################

# Go through each Robj created from prepSeurat and filter outliers. Replace the Robj with the filtered version
seuratObj <- readRDS("Robj/DropSeq/DropSeq_UCH1_HS_TB_DS52.Robj")
VlnPlot(object = seuratObj, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = seuratObj, gene1 = "nUMI", gene2 = "nGene")
GenePlot(object = seuratObj, gene1 = "nUMI", gene2 = "percent.mito")
seuratObj <- FilterCells(object = seuratObj, subset.names = c("nGene", "percent.mito"), 
                         low.thresholds = c(200, -Inf), high.thresholds = c(Inf, 0.1))
seuratObj <- NormalizeData(object = seuratObj,normalization.method = "LogNormalize", scale.factor = 10000)
seuratObj <- ScaleData(object = seuratObj, vars.to.regress = c("nUMI", "percent.mito"))
seuratObj <- FindVariableGenes(object = seuratObj, mean.function = ExpMean, dispersion.function = LogVMR, 
                               do.plot = FALSE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = seuratObj@var.genes)
saveRDS(seuratObj, "Robj/DropSeq/DropSeq_UCH1_HS_TB_DS52.Robj")

