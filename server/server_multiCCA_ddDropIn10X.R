# Seurat Multi-align
require(Seurat, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4/")
library(dplyr, lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4/")
setwd("~/R")

# UCH1 only
ddSeq_chor <- readRDS("Robj/Unscaled/ddSeq/ddSeq_UCH1.Robj")
DropSeq_chor <- readRDS("Robj/Unscaled/DropSeq/DropSeq_UCH1_lowCutoff.Robj")
inDrop_chor <- readRDS("Robj/Unscaled/inDrop/inDrop_UCH1.Robj")
t10X_chor <- readRDS("Robj/Unscaled/10X/10X_UCH1_Droplike.Robj")

ddSeq_chor <- ScaleData(ddSeq_chor)
DropSeq_chor <- ScaleData(DropSeq_chor)
inDrop_chor <- ScaleData(inDrop_chor)
t10X_chor <- ScaleData(t10X_chor)

ddSeq_chor@meta.data$protocol <- "ddSeq"
DropSeq_chor@meta.data$protocol <- "DropSeq"
inDrop_chor@meta.data$protocol <- "inDrop"
t10X_chor@meta.data$protocol <- "10X"

# Determine genes to use for CCA, must be highly variable in at least 2 datasets
sc.list <- list(ddSeq_chor, DropSeq_chor, inDrop_chor, t10X_chor)
genes.use <- c()
for (i in 1:length(sc.list)) {
  genes.use <- c(genes.use, head(rownames(sc.list[[i]]@hvg.info), 1000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(sc.list)) {
  genes.use <- genes.use[genes.use %in% rownames(sc.list[[i]]@scale.data)]
}

sc.cca <- RunMultiCCA(object.list = sc.list, genes.use = var_genes, num.ccs = 15)
#sc.cca <- readRDS("Robj/ddDropIn_cca.Robj")

p1 <- DimPlot(object = sc.cca, reduction.use = "cca", group.by = "protocol", pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = sc.cca, features.plot = "CC1", group.by = "protocol", do.return = TRUE)
plot_grid(p1, p2) 

#MetageneBicorPlot(sc.cca, grouping.var = "protocol", dims.eval = 1:15)

DimHeatmap(object = sc.cca, reduction.type = "cca", cells.use = 1000, dim.use = 1:15, 
           do.balanced = TRUE)
d <- 14 # adjust based on Heatmap
sc.cca <- CalcVarExpRatio(object = sc.cca, reduction.type = "pca", grouping.var = "protocol", dims.use = 1:d)
sc.cca_backup <- sc.cca
length(sc.cca@cell.names) # 5756
sc.cca <- SubsetData(object = sc.cca, subset.name = "var.ratio.pca", accept.low = 0.5) # Will discard some cells
length(sc.cca@cell.names) # 5693
sc.cca <- AlignSubspace(object = sc.cca, reduction.type = "cca", grouping.var = "protocol", dims.align = 1:d)
sc.cca <- RunTSNE(object = sc.cca, reduction.use = "cca.aligned", dims.use = 1:d, do.fast = TRUE)
TSNEPlot(object = sc.cca, group.by = "protocol", do.return = TRUE, pt.size = 2)

sc.cca <- FindClusters(sc.cca, reduction.type = "cca.aligned",
                                    dims.use = 1:d, save.SNN = T, resolution = 0.4)
sc.cca <- SetAllIdent(sc.cca, id="orig.ident")


