library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4/Seurat2.1')
library(Matrix)
library(dplyr)
set.seed(42)
setwd("~/R/Projects/Seurat")
# Analysis starts from line 67

#####################################################################################
# Read in files of chordoma data (UCH1 and UCH2)
#####################################################################################
# Load filtered/normalized Seurat Objects from DropSeq/Robj and ddSeq/Robj
# DS50_UCH1, DS52_UCH1_HS_TB, DS52_UCH1_RT_ID
# Seq2_N706, Seq3_N707
ddSeq_UCH2_Seq1_N704 <- readRDS("Robj/ddSeq/ddSeq_UCH2_Seq1_N704.Robj")
ddSeq_UCH1_Seq2_N706 <- readRDS("Robj/ddSeq/ddSeq_UCH1_Seq2_N706.Robj")
ddSeq_UCH1_Seq3_N707 <- readRDS("Robj/ddSeq/ddSeq_UCH1_Seq3_N707.Robj")

DropSeq_UCH2_DS49 <- readRDS("Robj/DropSeq/DropSeq_UCH2_DS49.Robj")
DropSeq_UCH1_DS50 <- readRDS("Robj/DropSeq/DropSeq_UCH1_DS50.Robj")
DropSeq_UCH1_DS52_RT_ID <- readRDS("Robj/DropSeq/DropSeq_UCH1_RT_ID_DS52.Robj")
DropSeq_UCH1_DS52_HS_TB <- readRDS("Robj/DropSeq/DropSeq_UCH1_HS_TB_DS52.Robj")

Fluidigm_UCH1 <- readRDS(file = 'Robj/Fluidigm/Fluidigm_UCH1.Robj')
#####################################################
# Merge and compare all samples
ddSeq_chor <- MergeSeurat(object1 = ddSeq_UCH2_Seq1_N704, object2 = ddSeq_UCH1_Seq2_N706, 
                     add.cell.id1 = "S1_N704", add.cell.id2 = "S2_N706")
ddSeq_chor <- MergeSeurat(object1 = ddSeq_chor, object2 = ddSeq_UCH1_Seq3_N707, 
                     add.cell.id2 = "S3_N707", project = "ddSeq_chor")

DropSeq_chor <- MergeSeurat(object1 = DropSeq_UCH2_DS49, object2 = DropSeq_UCH1_DS50, 
                       add.cell.id1 = "DS49", add.cell.id2 = "DS50")
DropSeq_chor <- MergeSeurat(object1 = DropSeq_chor, object2 = DropSeq_UCH1_DS52_RT_ID, 
                       add.cell.id2 = "DS52_RTID")
DropSeq_chor <- MergeSeurat(object1 = DropSeq_chor, object2 = DropSeq_UCH1_DS52_HS_TB, 
                       add.cell.id2 = "DS52_HSTB", project = "DropSeq_chor")

table(ddSeq_chor@meta.data$orig.ident)
table(DropSeq_chor@meta.data$orig.ident)

saveRDS(ddSeq_chor, file = '~/R/Projects/Seurat/Robj/ddDropF_chordoma/ddSeq_chor.Robj')
saveRDS(DropSeq_chor, file = '~/R/Projects/Seurat/Robj/ddDropF_chordoma/DropSeq_chor.Robj')

# You can focus on one technology or analyze all of them
ddDrop_chor <- MergeSeurat(object1 = ddSeq_chor, object2 = DropSeq_chor)
table(ddDrop_chor@meta.data$orig.ident)
saveRDS(ddDrop_chor, file = '~/R/Projects/Seurat/Robj/ddDropF_chordoma/ddDrop_chor.Robj')
#ddDrop_chor <- readRDS(file='~/R/Projects/Seurat/Robj/ddDrop_chor.Robj')

ddDropF_chor <- MergeSeurat(object1 = ddDrop_chor, object2 = Fluidigm_UCH1, 
                       add.cell.id2 = "Flu")
table(ddDropF_chor@meta.data$orig.ident)

# Normalize data to a total of 10,000 molecules 
ddDropF_chor <- NormalizeData(object = ddDropF_chor, normalization.method = "LogNormalize", scale.factor = 10000)
# Focus on variable genes (default parameter settings for a 1e4 molecule normalization)
ddDropF_chor <- FindVariableGenes(object = ddDropF_chor, mean.function = ExpMean, dispersion.function = LogVMR, 
                                  do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = ddDropF_chor@var.genes)
# Regress uninteresting signals out of analysis. nUMI: number of detected molecules per cell
ddDropF_chor <- ScaleData(object = ddDropF_chor, vars.to.regress = c("nUMI", "percent.mito"))

saveRDS(ddDropF_chor, file = '~/R/Projects/Seurat/Robj/ddDropF_chordoma/ddDropF_chor.Robj')
#####################################################
#####################################################
# Start Here if continuing previous analysis
ddDropF_chor <- readRDS(file='~/R/Projects/Seurat/Robj/ddDropF_chordoma/ddDropF_chor.Robj')
VlnPlot(object = ddDropF_chor, features.plot = c("nGene", "nUMI", "percent.mito"))

# Run Dimension reduction
ddDropF_chor <- RunPCA(object = ddDropF_chor, pc.genes = ddDropF_chor@var.genes, 
                  do.print = TRUE, pcs.print = 1:12, genes.print = 5, pcs.compute = 20)
#PCAPlot(object = ddDropF_chor, dim.1 = 1, dim.2 = 2)
#PCHeatmap(object = ddDropF_chor, pc.use = 9:18, cells.use = 200, do.balanced = TRUE, 
#          label.columns = FALSE, use.full = FALSE)
#PCElbowPlot(object = ddDropF_chor)
#ddDropF_chor <- JackStraw(object = ddDropF_chor, num.replicate = 30, do.print = FALSE)
#JackStrawPlot(object = ddDropF_chor, PCs = 1:20)

ddDropF_chor <- SetAllIdent(ddDropF_chor, id="orig.ident")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ddDropF_chor <- RunTSNE(object = ddDropF_chor, dims.use = 1:16, do.fast = TRUE, perplexity = 40)
t <- TSNEPlot(object = ddDropF_chor, do.return = TRUE, do.label = TRUE, pt.size = 2)
t + ggtitle('ddDropF_chor tSNE 16dim 40perplexity 7.8kVarGenes Seurat v2.1') + 
  theme(text = element_text(size=15)) + scale_colour_manual(values=cbPalette)

FeaturePlot(object = ddDropF_chor, features.plot = "T", cols.use = c("dark grey", "red"), 
            reduction.use = "tsne", dark.theme = TRUE)

#####################################################
#####################################################
# Custom analysis of clusters in ddDropF_chor
Fluidigm_sm <- TSNEPlot(object = ddDropF_chor, do.identify = TRUE) # 31 cells
Fluidigm_lg <- TSNEPlot(object = ddDropF_chor, do.identify = TRUE) # 71 cells
ddDropF_chor <- SetIdent(object = ddDropF_chor, cells.use = Fluidigm_sm, ident.use = "Fluidigm_sm")
ddDropF_chor <- SetIdent(object = ddDropF_chor, cells.use = Fluidigm_lg, ident.use = "Fluidigm_lg")
VlnPlot(object = ddDropF_chor, features.plot = c("nGene", "nUMI", "percent.mito"), 
        ident.include = c("Fluidigm_sm", "Fluidigm_lg")) 
Flu.mark <- FindMarkers(object = ddDropF_chor, ident.1 = 'Fluidigm_sm', ident.2 = 'Fluidigm_lg', min.pct = 0.25)
# None of the markers are significant

########################
DS49_sm <- TSNEPlot(object = ddDropF_chor, do.identify = TRUE) # 73 cells
DS49_lg <- TSNEPlot(object = ddDropF_chor, do.identify = TRUE) # 199 cells
ddDropF_chor <- SetIdent(object = ddDropF_chor, cells.use = DS49_sm, ident.use = "DS49_sm")
ddDropF_chor <- SetIdent(object = ddDropF_chor, cells.use = DS49_lg, ident.use = "DS49_lg")
VlnPlot(object = ddDropF_chor, features.plot = c("nGene", "nUMI", "percent.mito"), 
        ident.include = c("DS49_sm", "DS49_lg")) 
DS49.mark <- FindMarkers(object = ddDropF_chor, ident.1 = 'DS49_sm', ident.2 = 'DS49_lg', min.pct = 0.25)
# None of the markers are significant

#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
# Alignment Section
#rm(list=ls())
#####################################################
ddSeq_chor <- readRDS(file='~/R/Projects/Seurat/Robj/ddDropF_chordoma/ddSeq_chor.Robj')
DropSeq_chor <- readRDS(file='~/R/Projects/Seurat/Robj/ddDropF_chordoma/DropSeq_chor.Robj')

# Alignment relies on variable genes to determine source of variation 
hvg.ddSeq_chor.group <- rownames(x = head(x = ddSeq_chor@hvg.info, n = 4000)) # take 2000 genes with highest dispersion
hvg.DropSeq_chor.group <- rownames(x = head(x = DropSeq_chor@hvg.info, n = 4000))
hvg.union <- union(x = hvg.ddSeq_chor.group, y = hvg.DropSeq_chor.group)

ddSeq_chor@meta.data[, "protocol"] <- "ddSeq"
DropSeq_chor@meta.data[, "protocol"] <- "DropSeq"

#####################################################
align_dd_drop <- RunCCA(object = ddSeq_chor, object2 = DropSeq_chor, genes.use = hvg.union)
p1 <- DimPlot(object = align_dd_drop, reduction.use = "cca", group.by = "protocol", pt.size = 2,
              do.return = TRUE) + ggtitle('ddSeq_chor vs DropSeq_chor CCA All Cells')
p2 <- VlnPlot(object = align_dd_drop, features.plot = "CC1", group.by = "protocol", do.return = TRUE)
plot_grid(p1 + theme(text = element_text(size=20)), p2 + theme(text = element_text(size=20))) 
#####################################################
DimHeatmap(object = align_dd_drop, reduction.type = "cca", cells.use = 200, dim.use = 1:12, 
           do.balanced = TRUE)
align_dd_drop <- CalcVarExpRatio(object = align_dd_drop, reduction.type = "pca", grouping.var = "protocol", 
                                 dims.use = 1:8)
# We will NOT remove cells (as far as we know, these cells are the same)
discard <- SubsetData(object = align_dd_drop, subset.name = "var.ratio.pca", accept.high = 0.5)
#align_dd_drop <- SubsetData(object = align_dd_drop, subset.name = "var.ratio.pca", accept.low = 0.5)
align_dd_drop <- AlignSubspace(object = align_dd_drop, reduction.type = "cca", grouping.var = "protocol", 
                               dims.align = 1:8)
p1 <- VlnPlot(object = align_dd_drop, features.plot = "ACC1", group.by = "protocol", 
              do.return = TRUE)
p2 <- VlnPlot(object = align_dd_drop, features.plot = "ACC2", group.by = "protocol", 
              do.return = TRUE)
plot_grid(p1, p2)
#####################################################
align_dd_drop <- RunTSNE(object = align_dd_drop, reduction.use = "cca.aligned", 
                         dims.use = 1:8, do.fast = TRUE)
p1 <- TSNEPlot(object = align_dd_drop, group.by = "cellType", do.return = TRUE, pt.size = 2)
p1 + ggtitle('Aligned ddDrop (Group by cellType) tSNE v2.1') + theme(text=element_text(size=15))

align_dd_drop@meta.data[ ,"cellType"] <- NULL
align_dd_drop@meta.data[grep("UCH1", align_dd_drop@meta.data$orig.ident),"cellType"] <- "UCH1"
align_dd_drop@meta.data[grep("UCH2", align_dd_drop@meta.data$orig.ident),"cellType"] <- "UCH2"

# Group 1 is UCH1 dominant, but there are rogue UCH2 cells worth investigating
g1 <- TSNEPlot(object = align_dd_drop, group.by = "cellType", do.identify=TRUE, pt.size = 1) # grab UCH1
g1_uch2 <- g1[grep("DS49|S1_N704", g1)]
g1_uch1 <- g1[!g1 %in% g1_uch2]
g2 <- rownames(align_dd_drop@meta.data)[!rownames(align_dd_drop@meta.data) %in% g1]
align_dd_drop <- SetIdent(object = align_dd_drop, cells.use = g1_uch1, ident.use = "g1_uch1")
align_dd_drop <- SetIdent(object = align_dd_drop, cells.use = g1_uch2, ident.use = "g1_uch2")
align_dd_drop <- SetIdent(object = align_dd_drop, cells.use = g2, ident.use = "g2")
VlnPlot(object = align_dd_drop, features.plot = c("nGene", "nUMI", "percent.mito"), 
        ident.include = c("g1_uch1", "g1_uch2", "g2")) 
g1_uch2.mark <- FindMarkers(object = align_dd_drop, ident.1 = 'g1_uch2', ident.2 = 'g2', min.pct = 0.25)
FeaturePlot(object = align_dd_drop, features.plot = "T", cols.use = c("dark grey", "red"), 
            reduction.use = "tsne", dark.theme = TRUE)

g.mark <- FindMarkers(object = align_dd_drop, ident.1 = 'g1_uch1', ident.2 = 'g2', min.pct = 0.25)
# A (-) logFC means that Group 2 has higher avg expression than Group 1

saveRDS(align_dd_drop, file = '~/R/Projects/Seurat/Robj/ddDropF_chordoma/align_ddDrop_chord.Robj')
#align_dd_drop <- readRDS(file = '~/R/Projects/Seurat/Robj/ddDropF_chordoma/align_ddDrop_chord.Robj')

#####################################################
# Align with Fluidigm
align_dd_drop <- FindVariableGenes(object = align_dd_drop, mean.function = ExpMean, dispersion.function = LogVMR, 
                                   do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
hvg.align_dd_drop.group <- rownames(x = head(x = align_dd_drop@hvg.info, n = 2000)) # take 2000 genes with highest dispersion
hvg.Fluidigm.group <- rownames(x = head(x = Fluidigm_UCH1@hvg.info, n = 2000))
hvg.union <- union(x = hvg.align_dd_drop.group, y = hvg.Fluidigm.group)

Fluidigm_UCH1@meta.data[, "protocol"] <- "Fluidigm"
align_dd_drop@meta.data[, "protocol2"] <- "dd_Drop"
Fluidigm_UCH1@meta.data[, "protocol2"] <- "Fluidigm2"
#####################################################
a_ddDrop_F <- RunCCA(object = align_dd_drop, object2 = Fluidigm_UCH1, genes.use = hvg.union)
p1 <- DimPlot(object = a_ddDrop_F, reduction.use = "cca", group.by = "protocol", pt.size = 2,
              do.return = TRUE) + ggtitle('ddDrop_F CCA All Cells')
p2 <- VlnPlot(object = a_ddDrop_F, features.plot = "CC1", group.by = "protocol", do.return = TRUE)
plot_grid(p1 + theme(text = element_text(size=20)), p2 + theme(text = element_text(size=20))) 
#####################################################
DimHeatmap(object = a_ddDrop_F, reduction.type = "cca", cells.use = 200, dim.use = 4:15, 
           do.balanced = TRUE)
a_ddDrop_F <- CalcVarExpRatio(object = a_ddDrop_F, reduction.type = "pca", grouping.var = "protocol2", 
                                 dims.use = 1:9)
# We will NOT remove cells (as far as we know, these cells are the same)
discard <- SubsetData(object = a_ddDrop_F, subset.name = "var.ratio.pca", accept.high = 0.5)
#align_dd_drop <- SubsetData(object = align_dd_drop, subset.name = "var.ratio.pca", accept.low = 0.5)
a_ddDrop_F <- AlignSubspace(object = a_ddDrop_F, reduction.type = "cca", grouping.var = "protocol2", 
                               dims.align = 1:9)
p1 <- VlnPlot(object = a_ddDrop_F, features.plot = "ACC1", group.by = "protocol2", 
              do.return = TRUE)
p2 <- VlnPlot(object = a_ddDrop_F, features.plot = "ACC2", group.by = "protocol2", 
              do.return = TRUE)
plot_grid(p1, p2)
#####################################################
a_ddDrop_F <- RunTSNE(object = a_ddDrop_F, reduction.use = "cca.aligned", 
                         dims.use = 1:9, do.fast = TRUE)
p1 <- TSNEPlot(object = a_ddDrop_F, group.by = "orig.ident", do.return = TRUE, pt.size = 2)
p1 + ggtitle('Aligned ddDropF tSNE v2.1') + theme(text=element_text(size=15))

##########################################################################################################


