# Seq1-N704 ddSeq Downstream Analysis (UCH1 or UCH2???)

library(dplyr)
library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4/Seurat2.1')
setwd("~/R/Projects/Seurat/ddSeq")
## Prepare Seurat Object

# Load the dataset
DGE.data <- read.table(file = "Seq1-N704_UCH2/Input/Sassi-Seq1-N704_DGE.txt", header = TRUE, row.names = 1)

# Create Seurat object with raw data. We could filter based on cell/gene thresholds
DGE <- CreateSeuratObject(raw.data = DGE.data, min.cells = 3, min.genes = 200)

# Keep track of percentage mitochondrial gene content (to regress out later)
mito.genes <- grep(pattern = "_[Mm][Tt]-", x = rownames(x = DGE@data), value = TRUE)
percent.mito <- Matrix::colSums(DGE@raw.data[mito.genes, ])/Matrix::colSums(DGE@raw.data)
DGE <- AddMetaData(object = DGE, metadata = percent.mito, col.name = "percent.mito")

VlnPlot(object = DGE, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# We filter out cells that have unique gene counts over 5,000 or less than 200
DGE <- FilterCells(object = DGE, subset.names = c("nGene", "percent.mito"), 
                   low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.05))
# Normalize data to a total of 10,000 molecules 
DGE <- NormalizeData(object = DGE,normalization.method = "LogNormalize", scale.factor = 10000)
# Focus on variable genes (default parameter settings for a 1e4 molecule normalization)
DGE <- FindVariableGenes(object = DGE, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
# Regress uninteresting signals out of analysis. nUMI: number of detected molecules per cell
DGE <- ScaleData(object = DGE, vars.to.regress = c("nUMI", "percent.mito"))
length(x = DGE@var.genes)

# Run Dimension reduction
DGE <- RunPCA(object = DGE, pc.genes = DGE@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 30)
# visualize genes and cells that define the PCA
VizPCA(object = DGE, pcs.use = 1:2)
PCAPlot(object = DGE, dim.1 = 1, dim.2 = 2)
DGE <- ProjectPCA(object = DGE, do.print = FALSE)
#PCHeatmap(object = DGE, pc.use = 1, cells.use = 200, do.balanced = TRUE, label.columns = FALSE)
#PCHeatmap(object = DGE, pc.use = 1:12, cells.use = 200, do.balanced = TRUE, 
#          label.columns = FALSE, use.full = FALSE)

# How many components will we keep?
DGE <- JackStraw(object = DGE, num.replicate = 100, do.print = FALSE, num.pc = 20)
JackStrawPlot(object = DGE, PCs = 1:20)
PCElbowPlot(object = DGE, num.pc = 20)


# Determine cluster based on components used
DGE <- FindClusters(object = DGE, reduction.type = "pca", dims.use = 1:20, 
                    resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
DGE <- RunTSNE(object = DGE, dims.use = 1:20, do.fast = TRUE, perplexity = 30)
t <- TSNEPlot(object = DGE, do.return = TRUE)
t + ggtitle('Seq1-N704 tSNE v2.1')

DGE.markers <- FindAllMarkers(object = DGE, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)
DGE.markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_logFC)

# Observe expression patterns of known sex-linked genes
FeaturePlot(object = DGE, features.plot = c("DDX3X", "RPS4X", "USP9Y", "ZFY", "DDX3Y", "RPS4Y1"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
