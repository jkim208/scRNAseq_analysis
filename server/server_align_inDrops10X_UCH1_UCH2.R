library(Seurat)
library(Matrix)
library(dplyr)
setwd("~/R/Slim")

s_10X_UCH1 <- readRDS("Robj/Unscaled/10X/10X_UCH1_Droplike.Robj")
s_10X_UCH2 <- readRDS("Robj/Unscaled/10X/10X_UCH2_Droplike.Robj")
s_inDrop_UCH1 <- readRDS("Robj/Unscaled/inDrop/inDrop_lib21_UCH1.Robj")
s_inDrop_UCH2 <- readRDS("Robj/Unscaled/inDrop/inDrop_lib22_UCH2.Robj")

s_10X_UCH1 <- NormalizeData(s_10X_UCH1)
s_10X_UCH1 <- FindVariableGenes(object = s_10X_UCH1, do.plot = TRUE, x.low.cutoff = 0.05, x.high.cutoff = 3, y.cutoff = 0.5)
s_10X_UCH1 <- ScaleData(object = s_10X_UCH1)

s_inDrop_UCH1 <- NormalizeData(s_inDrop_UCH1)
s_inDrop_UCH1 <- FindVariableGenes(object = s_inDrop_UCH1, do.plot = TRUE, x.low.cutoff = 0.05, x.high.cutoff = 3, y.cutoff = 0.5)
s_inDrop_UCH1 <- ScaleData(object = s_inDrop_UCH1)

length(s_inDrop_UCH1@var.genes)
length(s_10X_UCH1@var.genes)
nrow(s_10X_UCH1@hvg.info)
g.10X <- head(rownames(s_10X_UCH1@hvg.info), 1000)
g.10X <- head(rownames(s_10X_UCH1@hvg.info), 1000)
g.inDrop <- head(rownames(s_inDrop_UCH1@hvg.info), 1000)
genes.use <- unique(c(g.10X, g.inDrop))

genes.use <- intersect(genes.use, rownames(s_10X_UCH1@scale.data))
genes.use <- intersect(genes.use, rownames(s_10X_UCH1@scale.data))
genes.use <- intersect(genes.use, rownames(s_inDrop_UCH1@scale.data))
length(genes.use)

s_inDrop10X_UCH1 <- RunCCA(s_10X_UCH1, s_inDrop_UCH1, genes.use = genes.use, num.cc = 30)
p1 <- DimPlot(object = s_inDrop10X_UCH1, reduction.use = "cca", group.by = "orig.ident", pt.size = 0.5, do.return = TRUE) + 
  theme(text = element_text(size = 22))
p2 <- VlnPlot(object = s_inDrop10X_UCH1, features.plot = "CC1", group.by = "orig.ident", do.return = TRUE) + 
  theme(text = element_text(size = 22))
png("~/p3.png", width = 1600, height = 900)
p3 <- MetageneBicorPlot(s_inDrop10X_UCH1, grouping.var = "orig.ident", dims.eval = 1:30, display.progress = FALSE)

s_inDrop10X_UCH1 <- AlignSubspace(s_inDrop10X_UCH1, reduction.type = "cca", grouping.var = "orig.ident", dims.align = 1:20)
p1 <- VlnPlot(object = s_inDrop10X_UCH1, features.plot = "ACC1", group.by = "orig.ident", do.return = TRUE) + theme(text = element_text(size = 22))
p2 <- VlnPlot(object = s_inDrop10X_UCH1, features.plot = "ACC2", group.by = "orig.ident", do.return = TRUE) + theme(text = element_text(size = 22))
plot_grid(p1, p2)

# t-SNE and Clustering
s_inDrop10X_UCH1 <- RunTSNE(s_inDrop10X_UCH1, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)
s_inDrop10X_UCH1 <- FindClusters(s_inDrop10X_UCH1, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20)

# Visualization
p1 <- TSNEPlot(s_inDrop10X_UCH1, do.return = T, pt.size = 1.2, group.by = "orig.ident") + theme(text = element_text(size = 22))
p2 <- TSNEPlot(s_inDrop10X_UCH1, pt.size = 1.2, do.return=T, do.label=TRUE,label.size=12) + theme(text = element_text(size = 22))
plot_grid(p1, p2)

FeaturePlot(s_inDrop10X_UCH1, features.plot=c("T"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90",do.return=TRUE)[[1]]
roc_4 <- FindMarkers(object = s_inDrop10X_UCH1, ident.1 = 4, min.pct = 0, test.use = "roc", only.pos = TRUE)
write.table(roc_4, file="~/align_UCH1_inDrops10X_rocDE_4-all.csv", quote=FALSE, sep="\t")

