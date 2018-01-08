library(dplyr)
library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4/Seurat2.1')

# Setup Robj for cross-technology analyses
# Seurat v2.1
# Filter parameters (2017): min.cells = 3, min.genes = 200, high.thresholds = c(5000, 0.05),
# x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5
setwd("~/R/Projects/Seurat")

# To access saved files from prepSeurat, use the following example:
# ddSeq_UCH1_Seq2_N706 <- readRDS("Robj/ddSeq/ddSeq_UCH1_Seq2_N706_v2.1.Robj")

prepSeurat <- function(DGE.data, technology, name) {
  # This function takes in data from HUMAN CELLS/GENES ONLY, name of the technology used,
  # and a name for the output file
  rownames(DGE.data) <- sub('^homo_', '', rownames(DGE.data) ) # remove 'homo' label
  seuratObj <- CreateSeuratObject(raw.data = DGE.data, min.cells = 3, min.genes = 200, project = name)
  mito.genes <- grep(pattern = "[Mm][Tt]-", x = rownames(x = seuratObj@data), value = TRUE)
  percent.mito <- Matrix::colSums(seuratObj@raw.data[mito.genes, ])/Matrix::colSums(seuratObj@raw.data)
  seuratObj <- AddMetaData(object = seuratObj, metadata = percent.mito, col.name = "percent.mito")
  destination <- paste('Robj/', technology, '/', name, '.Robj', sep="")
  saveRDS(seuratObj, destination)
  print(paste("Seurat object was saved to:", destination, sep=" "))
}

#############################################################################
# ddSeq
#############################################################################

## ddSeq_293_Seq1_N701 needed manual adjustment (mouse removed)
ddSeq_293_Seq1_N701 <- read.table(file = "ddSeq/Seq1-N701_293/Input/Seq1-N701_Human_DGE.txt",
                                   header = TRUE, row.names = 1)

ddSeq_UCH2_Seq1_N704 <- read.table(file = "ddSeq/Seq1-N704_UCH2/Input/Sassi-Seq1-N704_DGE.txt",
                                       header = TRUE, row.names = 1)
ddSeq_293_Seq2_N703 <- read.table(file = "ddSeq/Seq2-N703_293/Input/Sassi-Seq2-N703_DGE_Human_Selected.txt",
                                       header = TRUE, row.names = 1)
ddSeq_UCH1_Seq2_N706 <- read.table(file = "ddSeq/Seq2-N706_UCH1/Input/Sassi-Seq2-N706_DGE.txt",
                                       header = TRUE, row.names = 1)
ddSeq_CH19_Seq3_N706 <- read.table(file = "ddSeq/Seq3-N706_CH19/Input/Sassi-Seq3-N706_DGE.txt",
                                       header = TRUE, row.names = 1)
ddSeq_UCH1_Seq3_N707 <- read.table(file = "ddSeq/Seq3-N707_UCH1/Input/Sassi-Seq3-N707_DGE.txt",
                                   header = TRUE, row.names = 1)

#prepSeurat(ddSeq_293_Seq1_N701, 'ddSeq', 'ddSeq_293HEK_Seq1_N701')
prepSeurat(ddSeq_UCH2_Seq1_N704, 'ddSeq', 'ddSeq_UCH2_Seq1_N704')
#prepSeurat(ddSeq_293_Seq2_N703, 'ddSeq', 'ddSeq_293HEK_Seq2_N703')
prepSeurat(ddSeq_UCH1_Seq2_N706, 'ddSeq', 'ddSeq_UCH1_Seq2_N706')
prepSeurat(ddSeq_CH19_Seq3_N706, 'ddSeq', 'ddSeq_CH19_Seq3_N706')
prepSeurat(ddSeq_UCH1_Seq3_N707, 'ddSeq', 'ddSeq_UCH1_Seq3_N707')

#############################################################################
# DropSeq
#############################################################################

## DropSeq_293_DS9_OP needed manual adjustment (mouse removed)
DropSeq_DS9_HEK <- read.table(file = "DropSeq/DS_9_HEK_MEF/server/DS9_HEK_MEF_OP_Human_DGE.txt",
                                header = TRUE, row.names = 1)
DropSeq_DS34_HEK <- read.table(file = "DropSeq/DS_34_HEK_MEF/DS34_Human_DGE.txt",
                               header = TRUE, row.names = 1)
DropSeq_DS45_HEK <- read.table(file = "DropSeq/DS45_293ETN/server/DS45-CG_DGE.txt",
                               header = TRUE, row.names = 1)


DropSeq_UCH2_DS49 <- read.table(file = "DropSeq/DS49_UCH2/server/DS49-UCH2-Seq1_DGE.txt",
                                header = TRUE, row.names = 1)
DropSeq_UCH1_DS50 <- read.table(file = "DropSeq/DS50_UCH1/server/DS50-UCH1_DGE.txt",
                                   header = TRUE, row.names = 1)
DropSeq_UCH1_RT_ID_DS52 <- read.table(file = "DropSeq/DS52_UCH1/server/DS52-RT_ID__DGE.txt",
                                header = TRUE, row.names = 1)
DropSeq_UCH1_HS_TB_DS52 <- read.table(file = "DropSeq/DS52_UCH1/server/DS52-HS_TB__DGE.txt",
                                      header = TRUE, row.names = 1)

#prepSeurat(DropSeq_DS9_HEK, 'DropSeq', 'DropSeq_293HEK_DS9')
#prepSeurat(DropSeq_DS34_HEK, 'DropSeq', 'DropSeq_293HEK_DS34')
#prepSeurat(DropSeq_DS45_HEK, 'DropSeq', 'DropSeq_293HEK_DS45')
prepSeurat(DropSeq_UCH2_DS49, 'DropSeq', 'DropSeq_UCH2_DS49')
prepSeurat(DropSeq_UCH1_DS50, 'DropSeq', 'DropSeq_UCH1_DS50')
prepSeurat(DropSeq_UCH1_RT_ID_DS52, 'DropSeq', 'DropSeq_UCH1_RT_ID_DS52')
prepSeurat(DropSeq_UCH1_HS_TB_DS52, 'DropSeq', 'DropSeq_UCH1_HS_TB_DS52')

###########################################################################
# scRNA-seq QC
# Go through each Robj created from prepSeurat and filter outliers. Replace the Robj with the filtered version
seuratObj <- readRDS("Robj/DropSeq/DropSeq_UCH1_HS_TB_DS52.Robj")
VlnPlot(object = seuratObj, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = seuratObj, gene1 = "nUMI", gene2 = "nGene")
GenePlot(object = seuratObj, gene1 = "nUMI", gene2 = "percent.mito")
seuratObj <- FilterCells(object = seuratObj, subset.names = c("nGene", "percent.mito"), 
                         low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.1))
seuratObj <- NormalizeData(object = seuratObj,normalization.method = "LogNormalize", scale.factor = 10000)
seuratObj <- ScaleData(object = seuratObj, vars.to.regress = c("nUMI", "percent.mito"))
seuratObj <- FindVariableGenes(object = seuratObj, mean.function = ExpMean, dispersion.function = LogVMR, 
                               do.plot = FALSE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = seuratObj@var.genes)
saveRDS(seuratObj, "Robj/DropSeq/DropSeq_UCH1_HS_TB_DS52.Robj")

