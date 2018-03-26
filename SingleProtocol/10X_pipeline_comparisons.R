# 10X_pipeline_comparisons.R
# Updated: 26-March-2018
# Comparison between cellranger pipeline and droplike-pipeline in mapping 10X data
# Seurat v2.2.1

library(Seurat)
library(dplyr)
setwd("~/Slim")

##########################################################################################################
##########################################################################################################
# 1 - custom Droplike pipeline
droplike_reads <- read.table(file = "10X/sample1_UCH1/droplike_pipeline/s1_UCH1_10X_cell_readcounts.txt", 
                             header = FALSE, col.names = c("read.count", "barcode"))
x <- cumsum(droplike_reads$read.count)
x <- x/max(x)
plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]",
     ylab="cumulative fraction of reads", xlim=c(1,4000))
abline(v=2334)
title(main="10X UCH1 Droplike-Pipeline Cumulative Distribution Plot")

droplike.seurat <- readRDS(file="Robj/10X/10X_UCH1_s1_droplike.Robj") # object created in prepareSeuratObj.R

VlnPlot(object = droplike.seurat, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = droplike.seurat, gene1 = "nUMI", gene2 = "nGene")
GenePlot(object = droplike.seurat, gene1 = "nUMI", gene2 = "percent.mito")
##############################################
# 2 - cellranger pipeline
cellranger.data <- Read10X(data.dir = "~/R/Slim/10X/sample1_UCH1/filtered_gene_bc_matrices/GRCh38")
cellranger.seurat <- CreateSeuratObject(raw.data = cellranger.data, project = "10X_sample.seurat", min.cells = 3, min.genes = 200)
mito.genes <- grep(pattern = "MT-", x = rownames(x = cellranger.seurat@data), value = TRUE)
percent.mito <- Matrix::colSums(cellranger.seurat@raw.data[mito.genes, ])/Matrix::colSums(cellranger.seurat@raw.data)
cellranger.seurat <- AddMetaData(object = cellranger.seurat, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = cellranger.seurat, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
##########################################################################################################
##########################################################################################################

# Downstream processing from Droplike pipeline to match cellranger
droplike.seurat <- FilterCells(object = droplike.seurat, subset.names = c("nGene","nUMI","percent.mito"), 
                         low.thresholds = c(548, 2364, -Inf), high.thresholds = c(Inf, Inf, Inf))
droplike_reads_filtered <- droplike_reads[droplike_reads$barcode %in% droplike.seurat@cell.names,]

droplike.seurat <- NormalizeData(object = droplike.seurat,normalization.method = "LogNormalize", scale.factor = 10000)
droplike.seurat <- FindVariableGenes(object = droplike.seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                               do.plot = FALSE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(droplike.seurat@cell.names)
saveRDS(droplike.seurat, "Robj/10X/10X_UCH1_s1_droplike.Robj")

# quick summary
# 1. UCH1
#   A. cellRanger
#      - 2725 cells before any filter (54.5% capture)
#      - 2493 cells after upper threshold filters (nGene: 5000, MT: 0.15)
#   B. droplike (start at 5000 cells)
#      - 2723 cells after low threshold filters (nGene: 548, nUMI: 2364)
#      - 2457 cells after upper threshold filters (nGene: 5000, MT: 0.15)
# 2. UCH2
#   A. cellRanger
#      - 2334 cells before any filter (46.7% capture)
#      - 2036 cells after upper threshold filters (nGene: 6000, MT: 0.06)
#   B. droplike (start at 5000 cells)
#      - 2336 cells after low threshold filters (nGene: 1108, nUMI: 7040)
#      - 1990 cells after upper threshold filters (nGene: 6000, MT: 0.06)
#
