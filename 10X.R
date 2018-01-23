library(Seurat)
library(dplyr)
library(Matrix)

sample1.data <- Read10X(data.dir = "~/R/Projects/10X/sample1/filtered_gene_bc_matrices/GRCh38")
s1_UCH1 <- CreateSeuratObject(raw.data = sample1.data, min.cells = 3, min.genes = 200, project = "10X_s1_UCH1")
mito.genes <- grep(pattern = "^MT-", x = rownames(x = s1_UCH1@data), value = TRUE)
percent.mito <- Matrix::colSums(s1_UCH1@raw.data[mito.genes, ])/Matrix::colSums(s1_UCH1@raw.data)
s1_UCH1 <- AddMetaData(object = s1_UCH1, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = s1_UCH1, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = s1_UCH1, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = s1_UCH1, gene1 = "nUMI", gene2 = "nGene")
s1_UCH1 <- FilterCells(object = s1_UCH1, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(5300, 0.2))
s1_UCH1 <- NormalizeData(object = s1_UCH1, normalization.method = "LogNormalize", scale.factor = 10000)
s1_UCH1 <- FindVariableGenes(object = s1_UCH1, mean.function = ExpMean, dispersion.function = LogVMR, 
                             do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
s1_UCH1 <- ScaleData(object = s1_UCH1, vars.to.regress = c("nUMI", "percent.mito"))
s1_UCH1 <- RunPCA(object = s1_UCH1, pc.genes = s1_UCH1@var.genes, 
                  do.print = TRUE, pcs.print = 1:12, genes.print = 5, pcs.compute = 20)
PCElbowPlot(object = s1_UCH1)
PCAPlot(object = s1_UCH1, dim.1 = 1, dim.2 = 2)
s1_UCH1 <- RunTSNE(object = s1_UCH1, dims.use = 1:15, do.fast = TRUE, perplexity = 30)
t <- TSNEPlot(object = s1_UCH1, do.return = TRUE, do.label = TRUE, pt.size = 2)
t + ggtitle('10X_s1_UCH1 tSNE 15dim 30perplexity Seurat_v2.1') + theme(text = element_text(size=15))

FeaturePlot(object = s1_UCH1, features.plot = "T", cols.use = c("dark grey", "red"), 
            reduction.use = "tsne", dark.theme = TRUE)


