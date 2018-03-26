# inDrop10X_chord.R
# Updated: 26-March-2018
# Functional Analysis of scRNAseq data on Chordoma cells UCH1/UCH2 generated from inDrops and 10X protocols
# Includes standard merged analysis and alignment analysis via Seurat
# Seurat v2.2.1

library(Seurat)
library(Matrix)
library(dplyr)
setwd("~/Slim")

# Load in seurat objects to merge
s_10X_UCH1 <- readRDS("Robj/10X/10X_UCH1_Droplike.Robj")
s_10X_UCH2 <- readRDS("Robj/10X/10X_UCH2_Droplike.Robj")
s_inDrop_UCH1 <- readRDS("Robj/inDrop/inDrop_lib21_UCH1.Robj")
s_inDrop_UCH2 <- readRDS("Robj/inDrop/inDrop_lib22_UCH2.Robj")

s_inDrop_10X_chord <- MergeSeurat(object1 = s_10X_UCH1, object2 = s_10X_UCH2, 
                        add.cell.id1 = "10X_UCH1", add.cell.id2 = "10X_UCH2")
s_inDrop_10X_chord <- MergeSeurat(object1 = s_inDrop_10X_chord, object2 = s_inDrop_UCH1, 
                        add.cell.id2 = "inDrop_UCH1")
s_inDrop_10X_chord <- MergeSeurat(object1 = s_inDrop_10X_chord, object2 = s_inDrop_UCH2, 
                        add.cell.id2 = "inDrop_UCH2", project="chord_inDrop_10X")

s_inDrop_10X_chord <- NormalizeData(object = s_inDrop_10X_chord, normalization.method = "LogNormalize", scale.factor = 10000)
s_inDrop_10X_chord <- FindVariableGenes(object = s_inDrop_10X_chord, mean.function = ExpMean, dispersion.function = LogVMR, 
                              do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

#saveRDS(s_inDrop_10X_chord, "Robj/chord_inDrop_10X.Robj")
############################################################################################
s_inDrop_10X_chord <- readRDS("Robj/chord_inDrop_10X.Robj")

table(s_inDrop_10X_chord@meta.data$orig.ident) # total is 10054
VlnPlot(object = s_inDrop_10X_chord, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
s_inDrop_10X_chord <- ScaleData(object = s_inDrop_10X_chord)
s_inDrop_10X_chord <- RunPCA(object = s_inDrop_10X_chord, pc.genes = s_inDrop_10X_chord@var.genes, pcs.compute = 20)
PCAPlot(object = s_inDrop_10X_chord, dim.1 = 1, dim.2 = 2, do.return=TRUE) + 
  theme(text=element_text(size=15)) + ggtitle("inDrop-10X UCH1/UCH2 PCA")

# PC analysis: repeat for each pc and direction
pcGenes <- round(head(sort(s_inDrop_10X_chord@dr$pca@gene.loadings[,1], decreasing = TRUE), n=10), digits = 4)
write.table(pcGenes, file="saveData/inDrop10X_chord_PC/pc1_up.txt", quote=FALSE, sep="\t", col.names=FALSE)

PCElbowPlot(s_inDrop_10X_chord)

s_inDrop_10X_chord <- RunTSNE(object = s_inDrop_10X_chord, dims.use = 1:14, do.fast = TRUE, perplexity = 30)

s_inDrop_10X_chord <- FindClusters(s_inDrop_10X_chord, resolution = 0.6, dims.use = 1:14)
t <- TSNEPlot(object = s_inDrop_10X_chord, do.return = TRUE, pt.size = 1.2, do.label=TRUE,label.size=8)
t + ggtitle('inDrop_UCH1 tSNE 14dim 30perplexity v2.2.1') + theme(text = element_text(size=18))

markers <- FindMarkers(object = s_inDrop_10X_chord, ident.2 = c(3,5), ident.1=c(0,1,2,5),min.pct = 0.25)
head(markers, n=10)
write.table(markers, file="~/Desktop/inDrop_UCH1_clus3_4_DE.csv", quote=FALSE, sep="\t")
# when avg_logFC > 0, ident.1 is OVEREXPRESSED as opposed to ident.2
write.table(rownames(markers[markers$avg_logFC>0,]), file="~/Desktop/inDrop_UCH1_clus3_4_Up_list.txt", 
            quote=FALSE, col.names = FALSE,row.names = FALSE,sep="\t")
write.table(rownames(markers[markers$avg_logFC<0,]), file="~/Desktop/inDrop_UCH1_clus3_4_Down_list.txt", 
            quote=FALSE, col.names = FALSE,row.names = FALSE,sep="\t")

FeaturePlot(object = s_inDrop_10X_chord, features.plot = "T", cols.use = c("dark grey", "red"), reduction.use = "tsne", 
            dark.theme = TRUE, min.cutoff = "q10", max.cutoff = "q90", do.return = TRUE)[[1]] + theme(text = element_text(size=22))

png("./plots/f1.png", width=1600, height=900)
FeaturePlot(s_inDrop_10X_chord, features.plot=c("T","nGene","nUMI","MKI67","CKS2","CRYL1","DNA2","HJURP","SUOX"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(s_inDrop_10X_chord, features.plot=c("T","CCNB1","CCNB2","MCM2","MCM3","MCM4","MCM6","MCM7","MCM10"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(s_inDrop_10X_chord, features.plot=c("ENPP1","CSPG4","COL2A1","ACAN","SLC6A12","KCNK2"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(s_inDrop_10X_chord, features.plot=c("T","SUZ12","TCF7L2","IL6R","TGFB1","MMP14","ITGB1","EGR1","SOX9"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(s_inDrop_10X_chord, features.plot=c("T","LUM","SOX5","SOX6","NFAT5"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(s_inDrop_10X_chord, features.plot=c("T","CA3","TAF13","SYCP2L","KRT19","PDE1A","RRN3P3","CSPG4","WWP2"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(s_inDrop_10X_chord, features.plot=c("T","SAMD5","ENPP1","CD109","TRIL","MRGPRX3","C1QTNF3","EGFLAM","FMOD"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(s_inDrop_10X_chord, features.plot=c("T","HPGDS","MMP19","MIA","EMP2","RAB3B","PIP5K1A","ATF5","HOXA7"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(s_inDrop_10X_chord, features.plot=c("T","GDA","DSG2","GPM6A","AP3B2","FN1","ITGB1","TGFBR1","RELA"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")

VlnPlot(object = s_inDrop_10X_chord, features.plot = c("T"))
mean(s_inDrop_10X_chord@data["T", grep("10X", s_inDrop_10X_chord@cell.names)] > 0) # 85% of 10X UCH1 cells had >0 T expr

############################################################################################
# align 10X and inDrop
s_10X_chord <- MergeSeurat(object1 = s_10X_UCH1, object2 = s_10X_UCH2, project = "10X_chord",
                                  add.cell.id1 = "10X_UCH1", add.cell.id2 = "10X_UCH2")
s_inDrop_chord <- MergeSeurat(object1 = s_inDrop_UCH1, object2 = s_inDrop_UCH2, project = "inDrop_chord",
                                  add.cell.id1 = "inDrop_UCH1", add.cell.id2 = "inDrop_UCH2")

s_10X_chord <- NormalizeData(s_10X_chord)
s_10X_chord <- FindVariableGenes(object = s_10X_chord, do.plot = TRUE, x.low.cutoff = 0.1, x.high.cutoff = 3, y.cutoff = 0.5)
s_10X_chord <- ScaleData(object = s_10X_chord)

s_inDrop_chord <- NormalizeData(s_inDrop_chord)
s_inDrop_chord <- FindVariableGenes(object = s_inDrop_chord, do.plot = TRUE, x.low.cutoff = 0.1, x.high.cutoff = 3, y.cutoff = 0.5)
s_inDrop_chord <- ScaleData(object = s_inDrop_chord)

saveRDS(s_10X_chord, file = '~/Slim/Robj/s_10X_chord.Robj')
saveRDS(s_inDrop_chord, file = '~/Slim/Robj/s_inDrop_chord.Robj')
s_10X_chord <- readRDS(file='~/Slim/Robj/s_10X_chord.Robj')
s_inDrop_chord <- readRDS(file='~/Slim/Robj/s_inDrop_chord.Robj')

# Alignment relies on variable genes to determine source of variation 
g.10X <- head(rownames(s_10X_chord@hvg.info), 1000)
g.inDrop <- head(rownames(s_inDrop_chord@hvg.info), 1000)
genes.use <- unique(c(g.10X, g.inDrop))
genes.use <- intersect(genes.use, rownames(s_10X_chord@scale.data))
genes.use <- intersect(genes.use, rownames(s_inDrop_chord@scale.data))

s_inDrop10X_chord <- RunCCA(s_10X_chord, s_inDrop_chord, genes.use = genes.use, num.cc = 30)
p1 <- DimPlot(object = s_inDrop10X_chord, reduction.use = "cca", group.by = "tech", pt.size = 0.5, do.return = TRUE) + 
  theme(text = element_text(size = 22))
p2 <- VlnPlot(object = s_inDrop10X_chord, features.plot = "CC1", group.by = "tech", do.return = TRUE) + 
  theme(text = element_text(size = 22))
plot_grid(p1, p2)

p3 <- MetageneBicorPlot(s_inDrop10X_chord, grouping.var = "tech", dims.eval = 1:30, display.progress = FALSE)
DimHeatmap(object = s_inDrop10X_chord, reduction.type = "cca", cells.use = 500, dim.use = 1:9, do.balanced = TRUE)
s_inDrop10X_chord <- AlignSubspace(s_inDrop10X_chord, reduction.type = "cca", grouping.var = "tech", dims.align = 1:22)
p1 <- VlnPlot(object = s_inDrop10X_chord, features.plot = "ACC1", group.by = "tech", do.return = TRUE) + 
  theme(text = element_text(size = 22))
p2 <- VlnPlot(object = s_inDrop10X_chord, features.plot = "ACC2", group.by = "tech", do.return = TRUE) + 
  theme(text = element_text(size = 22))
plot_grid(p1, p2)

# t-SNE and Clustering
s_inDrop10X_chord <- RunTSNE(s_inDrop10X_chord, reduction.use = "cca.aligned", dims.use = 1:22, do.fast = T)
s_inDrop10X_chord <- FindClusters(s_inDrop10X_chord, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:22)

# Visualization
p1 <- TSNEPlot(s_inDrop10X_chord, do.return = T, pt.size = 0.5, group.by = "tech")
p2 <- TSNEPlot(s_inDrop10X_chord, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)



#FeaturePlot(s_inDrop10X_chord, features.plot=c("CCNB2","MCM2","MCM3","MCM4","MCM6","MCM7","MCM10","AURKA","AURKB"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
roc_highT_v_lowT<- FindMarkers(object = SSa40_merge, ident.1 = c(0:2,4,6:10), ident.2= c(3,5), min.pct = 0, test.use = "roc", only.pos = TRUE)

write.table(roc_highT_v_lowT, file="~/Desktop/align_chord_rocDE_highT_v_lowT", quote=FALSE, sep="/t")


