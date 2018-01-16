# Seurat Alignment of Different scRNAseq Technologies
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

ddSeq <- FindVariableGenes(object = ddSeq, mean.function = ExpMean, dispersion.function = LogVMR, 
                               do.plot = FALSE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
ddSeq <- ScaleData(object = ddSeq, vars.to.regress = c("nUMI", "percent.mito"))

DropSeq <- FindVariableGenes(object = DropSeq, mean.function = ExpMean, dispersion.function = LogVMR, 
                               do.plot = FALSE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
DropSeq <- ScaleData(object = DropSeq, vars.to.regress = c("nUMI", "percent.mito"))

saveRDS(ddSeq, file = '~/R/Projects/Seurat/Robj/ddSeq.Robj')
saveRDS(DropSeq, file = '~/R/Projects/Seurat/Robj/DropSeq.Robj')
#ddSeq <- readRDS(file='~/R/Projects/Seurat/Robj/ddSeq.Robj')
#DropSeq <- readRDS(file='~/R/Projects/Seurat/Robj/DropSeq.Robj')

# Alignment relies on variable genes to determine source of variation 
hvg.ddSeq.group <- rownames(x = head(x = ddSeq@hvg.info, n = 2000)) # take 2000 genes with highest dispersion
hvg.DropSeq.group <- rownames(x = head(x = DropSeq@hvg.info, n = 2000))
hvg.union <- union(x = hvg.ddSeq.group, y = hvg.DropSeq.group)

ddSeq@meta.data[, "protocol"] <- "ddSeq"
DropSeq@meta.data[, "protocol"] <- "DropSeq"

#####################################################

align_dd_drop <- RunCCA(object = ddSeq, object2 = DropSeq, genes.use = hvg.union)
p1 <- DimPlot(object = align_dd_drop, reduction.use = "cca", group.by = "protocol", pt.size = 0.5, 
              do.return = TRUE) + ggtitle('ddSeq vs DropSeq CCA All Cells')
p2 <- VlnPlot(object = align_dd_drop, features.plot = "CC1", group.by = "protocol", do.return = TRUE)
plot_grid(p1, p2) 
#####################################################

DimHeatmap(object = align_dd_drop, reduction.type = "cca", cells.use = 500, dim.use = 1:10, 
           do.balanced = TRUE)
# We leverage CCA to identify cells that cannot be aligned between the two datasets. Briefly, we
# quantify how well the low-dimensional space defined by CCA explains each cell's expression
# profile, and compare this to PCA, which is performed on each dataset independently.
# Cells where the percent variance explained is reduced by a userdefined cutoff in CCA compared 
# to PCA are therefore defined by sources of variance that are not shared between the datasets.
# We use a cutoff of 50% for all examples here to identify these cells, and discard them
# from the alignment workflow

align_dd_drop <- CalcVarExpRatio(object = align_dd_drop, 
                                      reduction.type = "pca", grouping.var = "protocol", 
                                      dims.use = 1:8)
align_dd_drop.preSubset <- align_dd_drop

# Attempt to remove dataset-specific cells (cells that supposedly only appear in 1 condition)
align_dd_drop <- SubsetData(object = align_dd_drop, 
                                 subset.name = "var.ratio.pca", accept.low = 0.50)
align_dd_drop.discard <- SubsetData(align_dd_drop.preSubset, subset.name = "var.ratio.pca", accept.high = 0.5)
ncol(align_dd_drop.discard@data)
#[1] 67
length(which(align_dd_drop.discard@meta.data$protocol=='DropSeq'))
#[1] 42
ncol(align_dd_drop@data)
#[1] 2782
length(which(align_dd_drop.preSubset@meta.data$protocol=='ddSeq'))
#[1] 754
length(which(align_dd_drop.preSubset@meta.data$protocol=='DropSeq'))
#[1] 2095
100*754/(754+2095) #percent of ddSeq (preSubset) to total
#[1] 26.47%
100*67/2782 # percent thrown out cells
#[1] 2.41%
#write.table(rownames(align_dd_drop.discard@meta.data), file ='~/R/Projects/Seurat/Seurat_dropped_cells.txt')

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
align_dd_drop <- FindClusters(object = align_dd_drop, reduction.type = "cca.aligned", 
                                   dims.use = 1:8, save.SNN = TRUE, force.recalc = TRUE)
p1 <- TSNEPlot(object = align_dd_drop, group.by = "protocol", do.return = TRUE, pt.size = 0.7)
p2 <- TSNEPlot(object = align_dd_drop, group.by = "orig.ident", 
               do.return = TRUE, pt.size = 0.9, do.label=TRUE)
p3 <- TSNEPlot(object = align_dd_drop, do.return = TRUE, pt.size = 0.7)
p1 + ggtitle('ddSeq vs DropSeq (Group by Protocol) tSNE v2.1')
p2 + ggtitle('ddSeq vs DropSeq (Group by Run) tSNE v2.1')
p3 + ggtitle('ddSeq vs DropSeq (Group by Original Clusters) tSNE v2.1')

FeaturePlot(object = align_dd_drop, features.plot = c("T"), cols.use = c("grey", "red"), 
            reduction.use = "tsne", do.return = TRUE)

saveRDS(align_dd_drop, "~/R/Projects/Seurat/Robj/align_dd_drop.Robj")


############################################################################
############################################################################
############################################################################
ddDrop <- MergeSeurat(object1 = ddSeq, object2 = DropSeq, project = "ddDrop")
ddDrop <- FindVariableGenes(object = ddDrop, mean.function = ExpMean, dispersion.function = LogVMR, 
                             do.plot = FALSE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
ddDrop <- ScaleData(object = ddDrop, vars.to.regress = c("nUMI", "percent.mito"))

# With Fluidigm
# Alignment relies on variable genes to determine source of variation 
hvg.ddDrop.group <- rownames(x = head(x = ddDrop@hvg.info, n = 2000)) # take 2000 genes with highest dispersion
hvg.Fluidigm.group <- rownames(x = head(x = Fluidigm@hvg.info, n = 2000))
hvg.union <- union(x = hvg.ddDrop.group, y = hvg.Fluidigm.group)

ddDrop@meta.data[, "protocol"] <- "ddDrop"
Fluidigm@meta.data[, "protocol"] <- "Fluidigm"

#####################################################

ddDrop_F <- RunCCA(object = ddDrop, object2 = Fluidigm, genes.use = hvg.union)
p1 <- DimPlot(object = ddDrop_F, reduction.use = "cca", group.by = "protocol", pt.size = 0.5, 
              do.return = TRUE) + ggtitle('ddSeq vs DropSeq CCA All Cells')
p2 <- VlnPlot(object = ddDrop_F, features.plot = "CC1", group.by = "protocol", do.return = TRUE)
plot_grid(p1, p2) 
#####################################################

DimHeatmap(object = ddDrop_F, reduction.type = "cca", cells.use = 500, dim.use = 1:10, 
           do.balanced = TRUE)
# We leverage CCA to identify cells that cannot be aligned between the two datasets. Briefly, we
# quantify how well the low-dimensional space defined by CCA explains each cell's expression
# profile, and compare this to PCA, which is performed on each dataset independently.
# Cells where the percent variance explained is reduced by a userdefined cutoff in CCA compared 
# to PCA are therefore defined by sources of variance that are not shared between the datasets.
# We use a cutoff of 50% for all examples here to identify these cells, and discard them
# from the alignment workflow

ddDrop_F <- CalcVarExpRatio(object = ddDrop_F, 
                                 reduction.type = "pca", grouping.var = "protocol", 
                                 dims.use = 1:8)
ddDrop_F.preSubset <- ddDrop_F

# Attempt to remove dataset-specific cells (cells that supposedly only appear in 1 condition)
ddDrop_F <- SubsetData(object = ddDrop_F, 
                            subset.name = "var.ratio.pca", accept.low = 0.50)
ddDrop_F.discard <- SubsetData(ddDrop_F.preSubset, subset.name = "var.ratio.pca", accept.high = 0.5)
ncol(ddDrop_F.discard@data)
#[1] 67
length(which(ddDrop_F.discard@meta.data$protocol=='DropSeq'))
#[1] 42
ncol(ddDrop_F@data)
#[1] 2782
length(which(ddDrop_F.preSubset@meta.data$protocol=='ddSeq'))
#[1] 754
length(which(ddDrop_F.preSubset@meta.data$protocol=='DropSeq'))
#[1] 2095
100*754/(754+2095) #percent of ddSeq (preSubset) to total
#[1] 26.47%
100*67/2782 # percent thrown out cells
#[1] 2.41%
#write.table(rownames(ddDrop_F.discard@meta.data), file ='~/R/Projects/Seurat/Seurat_dropped_cells.txt')

ddDrop_F <- AlignSubspace(object = ddDrop_F, reduction.type = "cca", grouping.var = "protocol", 
                               dims.align = 1:8)
p1 <- VlnPlot(object = ddDrop_F, features.plot = "ACC1", group.by = "protocol", 
              do.return = TRUE)
p2 <- VlnPlot(object = ddDrop_F, features.plot = "ACC2", group.by = "protocol", 
              do.return = TRUE)
plot_grid(p1, p2)

#####################################################

ddDrop_F <- RunTSNE(object = ddDrop_F, reduction.use = "cca.aligned", 
                         dims.use = 1:8, do.fast = TRUE)
ddDrop_F <- FindClusters(object = ddDrop_F, reduction.type = "cca.aligned", 
                              dims.use = 1:8, save.SNN = TRUE, force.recalc = TRUE)
p1 <- TSNEPlot(object = ddDrop_F, group.by = "protocol", do.return = TRUE, pt.size = 0.7)
p2 <- TSNEPlot(object = ddDrop_F, group.by = "orig.ident", 
               do.return = TRUE, pt.size = 0.9, do.label=TRUE)
p3 <- TSNEPlot(object = ddDrop_F, do.return = TRUE, pt.size = 0.7)
p1 + ggtitle('ddSeq vs DropSeq (Group by Protocol) tSNE v2.1')
p2 + ggtitle('ddSeq vs DropSeq (Group by Run) tSNE v2.1')
p3 + ggtitle('ddSeq vs DropSeq (Group by Original Clusters) tSNE v2.1')

FeaturePlot(object = ddDrop_F, features.plot = c("T"), cols.use = c("grey", "red"), 
            reduction.use = "tsne", do.return = TRUE)

saveRDS(ddDrop_F, "~/R/Projects/Seurat/Robj/ddDrop_F.Robj")







