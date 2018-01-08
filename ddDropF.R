library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4/Seurat2.1')
library(Matrix)
library(dplyr)
############################################################################################
# Start from Fluidigm and Merged ddSeq-DropSeq Seurat objects
Fluidigm <- readRDS(file = '~/R/Projects/Seurat/Fluidigm/Fluidigm_human.Robj')
mega_seurat <- readRDS(file = '~/R/Projects/Seurat/Robj/ddDrop_mega_seurat_v2.1.Robj')
mega_seurat <- MergeSeurat(object1 = mega_seurat, object2 = Fluidigm, 
                       add.cell.id2 = "flu")

# Normalize data to a total of 10,000 molecules 
mega_seurat <- NormalizeData(object = mega_seurat,normalization.method = "LogNormalize", scale.factor = 10000)
# Focus on variable genes (default parameter settings for a 1e4 molecule normalization)
mega_seurat <- FindVariableGenes(object = mega_seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                                 do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = mega_seurat@var.genes)
# Regress uninteresting signals out of analysis. nUMI: number of detected molecules per cell
mega_seurat <- ScaleData(object = mega_seurat, vars.to.regress = c("nUMI", "percent.mito"))

# Run Dimension reduction
mega_seurat <- RunPCA(object = mega_seurat, pc.genes = mega_seurat@var.genes, 
                      do.print = TRUE, pcs.print = 1:12, genes.print = 5, pcs.compute = 20)
PCElbowPlot(object = mega_seurat)
#mega_seurat <- JackStraw(object = mega_seurat, num.replicate = 50, do.print = FALSE)
#JackStrawPlot(object = mega_seurat, PCs = 1:20)
saveRDS(mega_seurat, file = '~/R/Projects/Seurat/Robj/ddFDrop_backup.Robj')
mega_seurat <- readRDS(file='~/R/Projects/Seurat/Robj/ddFDrop_backup.Robj')



mega_seurat <- FindClusters(object = mega_seurat, reduction.type = "pca", dims.use = 1:8, 
                            resolution = 0.4, save.SNN = TRUE, force.recalc = TRUE)
mega_seurat <- RunTSNE(object = mega_seurat, dims.use = 1:8, 
                       do.fast = TRUE, perplexity = 40)
t <- TSNEPlot(object = mega_seurat, do.return = TRUE, do.label=TRUE)
t + ggtitle('ddDropF tSNE 19dim 40perplexity 7.5k varGenes v2.1')

mega_seurat <- StashIdent(object = mega_seurat, save.name = "findClusters")
# Next, switch the identity class of all cells to reflect replicate ID
mega_seurat <- SetAllIdent(object = mega_seurat, id = "orig.ident")
mega_seurat <- StashIdent(object = mega_seurat, save.name = "adj.ident")
mega_seurat@meta.data$adj.ident[c(grep("COL0[1-9]", mega_seurat@meta.data$adj.ident), # first 10 columns.
                                  grep("COL10", mega_seurat@meta.data$adj.ident) ) ] <-rep("Fluidigm-HEK", 112)
mega_seurat@meta.data$adj.ident[c(grep("COL1[1-9]", mega_seurat@meta.data$adj.ident), # last 10 columns.
                                  grep("COL20", mega_seurat@meta.data$adj.ident) ) ] <-rep("Fluidigm-uch1", 118)
mega_seurat <- SetAllIdent(object = mega_seurat, id = "adj.ident")

FeaturePlot(object = mega_seurat, features.plot = "T", cols.use = c("dark grey", "red"), 
            reduction.use = "tsne", dark.theme = TRUE)

# manually add fluidigm IDs
mega_seurat@meta.data$celltype[2076:2203] <- "HEK"
mega_seurat@meta.data$celltype[2204:2475] <- "uch1"
levels(mega_seurat@meta.data$tech) <- c("ddSeq", "DropSeq", "Fluidigm")
mega_seurat@meta.data$tech[2076:2475] <- "Fluidigm"

mega_seurat@meta.data$class <- as.factor(paste(mega_seurat@meta.data$tech, 
                                     mega_seurat@meta.data$celltype, sep='-')) 

############################################################################################

#saveRDS(mega_seurat, file = '~/R/Projects/Seurat/Robj/dd_f_drop.Robj')
#mega_seurat <- readRDS(file='~/R/Projects/Seurat/Robj/dd_f_drop.Robj')

######################
# Explain DS52 far subcluster
DS52_mini <- TSNEPlot(object = mega_seurat, do.identify = TRUE)
DS52_big <- TSNEPlot(object = mega_seurat, do.identify = TRUE) #careful of the nearby fluidigm

mega_seurat <- SetIdent(object = mega_seurat, cells.use = DS52_mini, ident.use = "DS52_mini")
mega_seurat <- SetIdent(object = mega_seurat, cells.use = DS52_big, ident.use = "DS52_big")
DS52.mark <- FindMarkers(object = mega_seurat, ident.1 = 'DS52_mini', ident.2 = 'DS52_big', min.pct = 0.25)
DS52.mark$genes <- rownames(DS52.mark)
FeaturePlot(object = mega_seurat, features.plot = c("MALAT1", "CSPG4", "OLFML2A", "SLC4A11", "PALLD", "APOE"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
FeaturePlot(object = mega_seurat, features.plot = c("MT-ATP8", "MLEC", "MT-CO1", "DNAJC8", "OMD", "GLIPR1"), cols.use = c("dark grey", "red"), 
            reduction.use = "tsne", dark.theme = TRUE)
DS52.geneSet <- arrange(DS52.mark[DS52.mark$p_val_adj<0.05,], p_val_adj)$genes
write.table(DS52.geneSet, col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/DS52_top_markers.txt')
saveRDS(DS52.mark, file = '~/R/Projects/Seurat/Meeting4/DE_clusters/DS52_markers_DF.Robj')
# Generally, the difference between these distant subclusters is that the mini cluster has more instances of
# expressing genes at low levels but over the entire cluster. While the big cluster has DE genes that are
# more heterogeneous and spread out at higher expression levels. 

### readRDS from previous save state. Don't want to mess up identity meta data ###
######################

###################################################################################
# ddSeq vs. DropSeq
mega_seurat <- SetAllIdent(object = mega_seurat, id = "tech")
markers <- FindMarkers(object = mega_seurat, ident.1 = 'ddSeq', 
                       ident.2 = 'DropSeq', min.pct = 0.25)
FeaturePlot(object = mega_seurat, features.plot = c("MALAT1","MT-RNR2","RPS29","RPL39"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
FeaturePlot(object = mega_seurat, features.plot = c("ENO1","CAPNS1","IRAK1","SQSTM1"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
write.table(rownames(markers[markers$p_val_adj<0.05,]), 
            col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/DE_ddSeq_DropSeq.txt')
###################################################################################
# ddSeq vs. Fluidigm
markers <- FindMarkers(object = mega_seurat, ident.1 = 'ddSeq', 
                       ident.2 = 'Fluidigm', min.pct = 0.25)
FeaturePlot(object = mega_seurat, features.plot = c("GAPDH","FTH1","CA3","RPL30","TUBA1B", "UBB"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
write.table(rownames(markers[markers$p_val_adj<0.05 & abs(markers$avg_logFC)>0.5,]), 
            col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/ddSeq-Fluidigm.txt')
###################################################################################
# DropSeq vs. Fluidigm
markers <- FindMarkers(object = mega_seurat, ident.1 = 'DropSeq', 
                       ident.2 = 'Fluidigm', min.pct = 0.25)
FeaturePlot(object = mega_seurat, features.plot = c("TMSB10","OST4","PPIB","TXN","HSPB1", "RPLP1"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
FeaturePlot(object = mega_seurat, features.plot = c("MT-CYB","MT-ATP6","MT-ND4","MT-CO2"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
write.table(rownames(markers[markers$p_val_adj<0.05 & abs(markers$avg_logFC)>0.5,]), 
            col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/DropSeq-Fluidigm.txt')

###################################################################################
# HEK vs UCH1
mega_seurat <- SetAllIdent(object = mega_seurat, id = "celltype")
markers <- FindMarkers(object = mega_seurat, ident.1 = 'HEK', 
                       ident.2 = 'uch1', min.pct = 0.25)
FeaturePlot(object = mega_seurat, features.plot = c("HSPA1A","LDHB","RPL5","ZNF711"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
FeaturePlot(object = mega_seurat, features.plot = c("TMSB4X", "KRT19", "FN1", "HLA-B"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
write.table(rownames(markers[markers$p_val_adj<0.05 & abs(markers$avg_logFC)>0.5,]), 
            col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/HEK_vs_UCH1.txt')
###################################################################################
# HEK vs UCH2
markers <- FindMarkers(object = mega_seurat, ident.1 = 'HEK', 
                       ident.2 = 'uch2', min.pct = 0.25)
FeaturePlot(object = mega_seurat, features.plot = c("NEFL", "SAT1", "PTMA", "HSPD1"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
FeaturePlot(object = mega_seurat, features.plot = c("TMSB4X", "KRT19", "FN1", "HLA-B"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
write.table(rownames(markers[markers$p_val_adj<0.05 & abs(markers$avg_logFC)>0.5,]), 
            col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/HEK_vs_UCH2.txt')
###################################################################################
# UCH1 vs UCH2
markers <- FindMarkers(object = mega_seurat, ident.1 = 'uch1', 
                       ident.2 = 'uch2', min.pct = 0.25)
FeaturePlot(object = mega_seurat, features.plot = c("CTSC", "ADIRF", "CA2", "ERC1"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
FeaturePlot(object = mega_seurat, features.plot = c("KRT19", "ANXA2", "CCDC80", "CAV1"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
write.table(rownames(markers[markers$p_val_adj<0.05 & abs(markers$avg_logFC)>0.5,]), 
            col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/UCH1_vs_UCH2.txt')


###################################################################################
###################################################################################
# 293HEK in ddSeq and DropSeq
mega_seurat <- SetAllIdent(object = mega_seurat, id = "class")
markers <- FindMarkers(object = mega_seurat, ident.1 = 'ddSeq-HEK', 
                       ident.2 = 'DropSeq-HEK', min.pct = 0.25)
FeaturePlot(object = mega_seurat, features.plot = "TSPYL1", cols.use = c("dark grey", "red"), 
            reduction.use = "tsne", dark.theme = TRUE)
# There appears to be a consistent low-level (1-2) expression of this Y gene in ddSeq 293HEK
FeaturePlot(object = mega_seurat, features.plot = "RPS4Y1", cols.use = c("dark grey", "red"), 
            reduction.use = "tsne", dark.theme = TRUE)
# However this Y gene disappears in all 293HEK and UCH2 samples (as expected)
write.table(rownames(markers[markers$p_val_adj<0.05,]), 
            col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/ddSeqvsDropSeq_HEK.txt')
###################################################################################






###################################################################################



###################################################################################

###################################################################################
