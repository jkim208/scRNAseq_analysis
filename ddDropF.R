library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4/Seurat2.1')
library(Matrix)
library(dplyr)
setwd("~/R/Projects/Seurat")
#####################################################################################
# Read in files
#####################################################################################
# Load filtered/normalized Seurat Objects from DropSeq/Robj and ddSeq/Robj
# DS50_UCH1, DS52_UCH1_HS_TB, DS52_UCH1_RT_ID
# Seq2_N706, Seq3_N707
ddSeq_293HEK_Seq1_N701 <- readRDS("Robj/ddSeq/ddSeq_293HEK_Seq1_N701.Robj")
ddSeq_UCH2_Seq1_N704 <- readRDS("Robj/ddSeq/ddSeq_UCH2_Seq1_N704.Robj")
ddSeq_UCH1_Seq2_N706 <- readRDS("Robj/ddSeq/ddSeq_UCH1_Seq2_N706.Robj")
ddSeq_UCH1_Seq3_N707 <- readRDS("Robj/ddSeq/ddSeq_UCH1_Seq3_N707.Robj")

DropSeq_293HEK_OP_DS9 <- readRDS("Robj/DropSeq/DropSeq_293HEK_OP_DS9.Robj")
DropSeq_293HEK_NO_DS9 <- readRDS("Robj/DropSeq/DropSeq_293HEK_NO_DS9.Robj")
DropSeq_293HEK_DS34 <- readRDS("Robj/DropSeq/DropSeq_293HEK_OP_DS34.Robj")
DropSeq_293HEK_DS45 <- readRDS("Robj/DropSeq/DropSeq_293HEK_OP_DS45.Robj")
DropSeq_UCH2_DS49 <- readRDS("Robj/DropSeq/DropSeq_UCH2_DS49.Robj")
DropSeq_UCH1_DS50 <- readRDS("Robj/DropSeq/DropSeq_UCH1_DS50.Robj")
DropSeq_UCH1_DS52_RT_ID <- readRDS("Robj/DropSeq/DropSeq_UCH1_RT_ID_DS52.Robj")
DropSeq_UCH1_DS52_HS_TB <- readRDS("Robj/DropSeq/DropSeq_UCH1_HS_TB_DS52.Robj")

Fluidigm <- readRDS(file = '~/R/Projects/Seurat/Fluidigm/Fluidigm_human.Robj')
#####################################################
# Merge and compare all samples
ddSeq <- MergeSeurat(object1 = ddSeq_293HEK_Seq1_N701, object2 = ddSeq_UCH2_Seq1_N704, 
                     add.cell.id1 = "S1_N701", add.cell.id2 = "S1_N704")
ddSeq <- MergeSeurat(object1 = ddSeq, object2 = ddSeq_UCH1_Seq2_N706, 
                     add.cell.id2 = "S2_N706")
ddSeq <- MergeSeurat(object1 = ddSeq, object2 = ddSeq_UCH1_Seq3_N707, 
                     add.cell.id2 = "S3_N707", project = "ddSeq")

DropSeq <- MergeSeurat(object1 = DropSeq_293HEK_OP_DS9, object2 = DropSeq_293HEK_NO_DS9, 
                       add.cell.id1 = "DS9_OP", add.cell.id2 = "DS9_NO")
DropSeq <- MergeSeurat(object1 = DropSeq, object2 = DropSeq_293HEK_DS34, 
                       add.cell.id2 = "DS34")
DropSeq <- MergeSeurat(object1 = DropSeq, object2 = DropSeq_293HEK_DS45, 
                       add.cell.id2 = "DS45")
DropSeq <- MergeSeurat(object1 = DropSeq, object2 = DropSeq_UCH2_DS49, 
                       add.cell.id2 = "DS49")
DropSeq <- MergeSeurat(object1 = DropSeq, object2 = DropSeq_UCH1_DS50, 
                       add.cell.id2 = "DS50")
DropSeq <- MergeSeurat(object1 = DropSeq, object2 = DropSeq_UCH1_DS52_RT_ID, 
                       add.cell.id2 = "DS52_RTID")
DropSeq <- MergeSeurat(object1 = DropSeq, object2 = DropSeq_UCH1_DS52_HS_TB, 
                       add.cell.id2 = "DS52_HSTB", project = "DropSeq")

table(ddSeq@meta.data$orig.ident)
table(DropSeq@meta.data$orig.ident)

# You can focus on one technology or analyze all of them
ddDrop <- MergeSeurat(object1 = ddSeq, object2 = DropSeq)
table(ddDrop@meta.data$orig.ident)
saveRDS(ddDrop, file = '~/R/Projects/Seurat/Robj/ddDrop.Robj')
#ddDrop <- readRDS(file='~/R/Projects/Seurat/Robj/ddDrop.Robj')

ddDropF <- MergeSeurat(object1 = ddDrop, object2 = Fluidigm, 
                       add.cell.id2 = "flu")
table(ddDropF@meta.data$orig.ident)


# Normalize data to a total of 10,000 molecules 
ddDropF <- NormalizeData(object = ddDropF, normalization.method = "LogNormalize", scale.factor = 10000)
# Focus on variable genes (default parameter settings for a 1e4 molecule normalization)
ddDropF <- FindVariableGenes(object = ddDropF, mean.function = ExpMean, dispersion.function = LogVMR, 
                                 do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = ddDropF@var.genes)
# Regress uninteresting signals out of analysis. nUMI: number of detected molecules per cell
ddDropF <- ScaleData(object = ddDropF, vars.to.regress = c("nUMI", "percent.mito"))

# Run Dimension reduction
ddDropF <- RunPCA(object = ddDropF, pc.genes = ddDropF@var.genes, 
                      do.print = TRUE, pcs.print = 1:12, genes.print = 5, pcs.compute = 20)
PCElbowPlot(object = ddDropF)
PCAPlot(object = ddDropF, dim.1 = 1, dim.2 = 2)
ddDropF <- ProjectPCA(object = ddDropF, do.print = FALSE)
#PCHeatmap(object = ddDropF, pc.use = 1, cells.use = 200, do.balanced = TRUE, label.columns = FALSE)
#PCHeatmap(object = ddDropF, pc.use = 1:12, cells.use = 200, do.balanced = TRUE, 
#          label.columns = FALSE, use.full = FALSE)
#ddDropF <- JackStraw(object = ddDropF, num.replicate = 50, do.print = FALSE)
#JackStrawPlot(object = ddDropF, PCs = 1:20)

ddDropF <- FindClusters(object = ddDropF, reduction.type = "pca", dims.use = 1:11, 
                            resolution = 0.4, save.SNN = TRUE, force.recalc = TRUE)
ddDropF <- RunTSNE(object = ddDropF, dims.use = 1:11, do.fast = TRUE, perplexity = 30)
t <- TSNEPlot(object = ddDropF, do.return = TRUE, do.label = TRUE, pt.size = 2)
t + ggtitle('ddDropF tSNE 11dim 30perplexity 7.8k varGenes v2.1') + theme(text = element_text(size=15))

FeaturePlot(object = ddDropF, features.plot = "T", cols.use = c("dark grey", "red"), 
            reduction.use = "tsne", dark.theme = TRUE)

############################################################################################
# Explore ddDropF
saveRDS(ddDropF, file = '~/R/Projects/Seurat/Robj/ddFDrop.Robj')
#ddDropF <- readRDS(file='~/R/Projects/Seurat/Robj/ddFDrop.Robj')

ddDropF <- SetAllIdent(object = ddDropF, id = "orig.ident")
VlnPlot(object = ddDropF, features.plot = c("nGene", "nUMI", "percent.mito"), 
        ident.include = c("293HEK", "UCH1", "UCH2")) 
VlnPlot(object = ddDropF, features.plot = c("nGene", "nUMI", "percent.mito"), 
        ident.include = c("ddSeq_293HEK_Seq1_N701", "DropSeq_293HEK_NO_DS9", "DropSeq_293HEK_OP_DS34", "DropSeq_293HEK_OP_DS45", "DropSeq_293HEK_OP_DS9", "Fluidigm_293HEK"), x.lab.rot = 1) 
VlnPlot(object = ddDropF, features.plot = c("nGene", "nUMI", "percent.mito"), 
        ident.include = c("ddSeq_UCH2_Seq1_N704", "DropSeq_UCH2_DS49"), x.lab.rot = 1) 


############################################################################################
############################################################################################
############################################################################################
# Explain DS52 far subcluster
DS52_cells <- rownames(ddDropF@meta.data[grep("_DS52", ddDropF@meta.data$orig.ident),])
DS52_mini <- TSNEPlot(object = ddDropF, do.identify = TRUE)
DS52_big_rows <- !DS52_cells %in% DS52_mini
DS52_big <- DS52_cells[DS52_big_rows]
ddDropF <- SetIdent(object = ddDropF, cells.use = DS52_mini, ident.use = "DS52_mini")
ddDropF <- SetIdent(object = ddDropF, cells.use = DS52_big, ident.use = "DS52_big")
VlnPlot(object = ddDropF, features.plot = c("nGene", "nUMI"), 
        ident.include = c("DS52_big", "DS52_mini")) 
DS52.mark <- FindMarkers(object = ddDropF, ident.1 = 'DS52_big', ident.2 = 'DS52_mini', min.pct = 0.25)
FeaturePlot(object = ddDropF, features.plot = "T", cols.use = c("dark grey", "red"), 
            reduction.use = "tsne", dark.theme = TRUE)
############################################################################################
# Explain Seq1 far subcluster
Seq1_N701_mini <- TSNEPlot(object = ddDropF, do.identify = TRUE)
Seq1_N701_big <- TSNEPlot(object = ddDropF, do.identify = TRUE) 
ddDropF <- SetIdent(object = ddDropF, cells.use = Seq1_N701_mini, ident.use = "Seq1_N701_mini")
ddDropF <- SetIdent(object = ddDropF, cells.use = Seq1_N701_big, ident.use = "Seq1_N701_big")
VlnPlot(object = ddDropF, features.plot = c("nGene", "nUMI", "percent.mito"), 
        ident.include = c("Seq1_N701_big", "Seq1_N701_mini"), x.lab.rot = 1) 
Seq1.mark <- FindMarkers(object = ddDropF, ident.1 = 'Seq1_N701_big', ident.2 = 'Seq1_N701_mini', min.pct = 0.25)
FeaturePlot(object = ddDropF, features.plot = "T", cols.use = c("dark grey", "red"), 
            reduction.use = "tsne", dark.theme = TRUE)
Seq1.mark$genes <- rownames(Seq1.mark)
Seq1.geneSet <- arrange(Seq1.mark[Seq1.mark$p_val_adj<0.05,], p_val_adj)$genes
write.table(Seq1.geneSet, col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/Seq1_top_markers.txt')
saveRDS(Seq1.mark, file = '~/R/Projects/Seurat/Meeting4/DE_clusters/Seq1_markers_DF.Robj')
# Generally, the difference between these distant subclusters is that the mini cluster has more instances of
# expressing genes at low levels but over the entire cluster. While the big cluster has DE genes that are
# more heterogeneous and spread out at higher expression levels. 

### readRDS from previous save state. Don't want to mess up identity meta data ###
######################

###################################################################################
# ddSeq vs. DropSeq
ddDropF <- SetAllIdent(object = ddDropF, id = "tech")
markers <- FindMarkers(object = ddDropF, ident.1 = 'ddSeq', 
                       ident.2 = 'DropSeq', min.pct = 0.25)
FeaturePlot(object = ddDropF, features.plot = c("MALAT1","MT-RNR2","RPS29","RPL39"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
FeaturePlot(object = ddDropF, features.plot = c("ENO1","CAPNS1","IRAK1","SQSTM1"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
write.table(rownames(markers[markers$p_val_adj<0.05,]), 
            col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/DE_ddSeq_DropSeq.txt')
###################################################################################
# ddSeq vs. Fluidigm
markers <- FindMarkers(object = ddDropF, ident.1 = 'ddSeq', 
                       ident.2 = 'Fluidigm', min.pct = 0.25)
FeaturePlot(object = ddDropF, features.plot = c("GAPDH","FTH1","CA3","RPL30","TUBA1B", "UBB"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
write.table(rownames(markers[markers$p_val_adj<0.05 & abs(markers$avg_logFC)>0.5,]), 
            col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/ddSeq-Fluidigm.txt')
###################################################################################
# DropSeq vs. Fluidigm
markers <- FindMarkers(object = ddDropF, ident.1 = 'DropSeq', 
                       ident.2 = 'Fluidigm', min.pct = 0.25)
FeaturePlot(object = ddDropF, features.plot = c("TMSB10","OST4","PPIB","TXN","HSPB1", "RPLP1"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
FeaturePlot(object = ddDropF, features.plot = c("MT-CYB","MT-ATP6","MT-ND4","MT-CO2"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
write.table(rownames(markers[markers$p_val_adj<0.05 & abs(markers$avg_logFC)>0.5,]), 
            col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/DropSeq-Fluidigm.txt')

###################################################################################
# HEK vs UCH1
ddDropF <- SetAllIdent(object = ddDropF, id = "celltype")
markers <- FindMarkers(object = ddDropF, ident.1 = 'HEK', 
                       ident.2 = 'uch1', min.pct = 0.25)
FeaturePlot(object = ddDropF, features.plot = c("HSPA1A","LDHB","RPL5","ZNF711"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
FeaturePlot(object = ddDropF, features.plot = c("TMSB4X", "KRT19", "FN1", "HLA-B"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
write.table(rownames(markers[markers$p_val_adj<0.05 & abs(markers$avg_logFC)>0.5,]), 
            col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/HEK_vs_UCH1.txt')
###################################################################################
# HEK vs UCH2
markers <- FindMarkers(object = ddDropF, ident.1 = 'HEK', 
                       ident.2 = 'uch2', min.pct = 0.25)
FeaturePlot(object = ddDropF, features.plot = c("NEFL", "SAT1", "PTMA", "HSPD1"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
FeaturePlot(object = ddDropF, features.plot = c("TMSB4X", "KRT19", "FN1", "HLA-B"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
write.table(rownames(markers[markers$p_val_adj<0.05 & abs(markers$avg_logFC)>0.5,]), 
            col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/HEK_vs_UCH2.txt')
###################################################################################
# UCH1 vs UCH2
markers <- FindMarkers(object = ddDropF, ident.1 = 'uch1', 
                       ident.2 = 'uch2', min.pct = 0.25)
FeaturePlot(object = ddDropF, features.plot = c("CTSC", "ADIRF", "CA2", "ERC1"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
FeaturePlot(object = ddDropF, features.plot = c("KRT19", "ANXA2", "CCDC80", "CAV1"), 
            cols.use = c("dark grey", "red"), reduction.use = "tsne", dark.theme = TRUE)
write.table(rownames(markers[markers$p_val_adj<0.05 & abs(markers$avg_logFC)>0.5,]), 
            col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/UCH1_vs_UCH2.txt')


###################################################################################
###################################################################################
# 293HEK in ddSeq and DropSeq
ddDropF <- SetAllIdent(object = ddDropF, id = "class")
markers <- FindMarkers(object = ddDropF, ident.1 = 'ddSeq-HEK', 
                       ident.2 = 'DropSeq-HEK', min.pct = 0.25)
FeaturePlot(object = ddDropF, features.plot = "TSPYL1", cols.use = c("dark grey", "red"), 
            reduction.use = "tsne", dark.theme = TRUE)
# There appears to be a consistent low-level (1-2) expression of this Y gene in ddSeq 293HEK
FeaturePlot(object = ddDropF, features.plot = "RPS4Y1", cols.use = c("dark grey", "red"), 
            reduction.use = "tsne", dark.theme = TRUE)
# However this Y gene disappears in all 293HEK and UCH2 samples (as expected)
write.table(rownames(markers[markers$p_val_adj<0.05,]), 
            col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting4/DE_clusters/ddSeqvsDropSeq_HEK.txt')
###################################################################################






###################################################################################



###################################################################################

###################################################################################
