# Analysis with HEK and UCH1 only (Fluidigm did not have UCH2)
# Object preparation and sample alignment of Fluidigm with DropSeq
# Focus on technology differences

library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4/Seurat2.1')
library(Matrix);library(dplyr)
set.seed(42)
setwd("~/R/Projects/Seurat")

##########################################################################################################
# Seurat Object Preparation
##########################################################################################################
ddSeq_293HEK_Seq1_N701 <- readRDS("Robj/ddSeq/ddSeq_293HEK_Seq1_N701.Robj")
ddSeq_UCH1_Seq2_N706 <- readRDS("Robj/ddSeq/ddSeq_UCH1_Seq2_N706.Robj")
ddSeq_UCH1_Seq3_N707 <- readRDS("Robj/ddSeq/ddSeq_UCH1_Seq3_N707.Robj")

ddSeq <- MergeSeurat(object1 = ddSeq_293HEK_Seq1_N701, object2 = ddSeq_UCH1_Seq2_N706, 
                     add.cell.id1 = "S1_N701", add.cell.id2 = "S2_N706")
ddSeq <- MergeSeurat(object1 = ddSeq, object2 = ddSeq_UCH1_Seq3_N707, 
                     add.cell.id2 = "S3_N707", project = "ddSeq")

DropSeq_293HEK_OP_DS9 <- readRDS("Robj/DropSeq/DropSeq_293HEK_OP_DS9.Robj")
DropSeq_293HEK_NO_DS9 <- readRDS("Robj/DropSeq/DropSeq_293HEK_NO_DS9.Robj")
DropSeq_293HEK_DS34 <- readRDS("Robj/DropSeq/DropSeq_293HEK_OP_DS34.Robj")
DropSeq_293HEK_DS45 <- readRDS("Robj/DropSeq/DropSeq_293HEK_OP_DS45.Robj")
DropSeq_UCH1_DS50 <- readRDS("Robj/DropSeq/DropSeq_UCH1_DS50.Robj")
DropSeq_UCH1_DS52_RT_ID <- readRDS("Robj/DropSeq/DropSeq_UCH1_RT_ID_DS52.Robj")
DropSeq_UCH1_DS52_HS_TB <- readRDS("Robj/DropSeq/DropSeq_UCH1_HS_TB_DS52.Robj")

DropSeq <- MergeSeurat(object1 = DropSeq_293HEK_OP_DS9, object2 = DropSeq_293HEK_NO_DS9, 
                       add.cell.id1 = "DS9_OP", add.cell.id2 = "DS9_NO")
DropSeq <- MergeSeurat(object1 = DropSeq, object2 = DropSeq_293HEK_DS34, add.cell.id2 = "DS34")
DropSeq <- MergeSeurat(object1 = DropSeq, object2 = DropSeq_293HEK_DS45, add.cell.id2 = "DS45")
DropSeq <- MergeSeurat(object1 = DropSeq, object2 = DropSeq_UCH1_DS50,  add.cell.id2 = "DS50")
DropSeq <- MergeSeurat(object1 = DropSeq, object2 = DropSeq_UCH1_DS52_RT_ID, add.cell.id2 = "DS52_RTID")
DropSeq <- MergeSeurat(object1 = DropSeq, object2 = DropSeq_UCH1_DS52_HS_TB, 
                       add.cell.id2 = "DS52_HSTB", project = "DropSeq")

tech <- NormalizeData(object = ddSeq, normalization.method = "LogNormalize", scale.factor = 10000)
tech <- FindVariableGenes(object = tech, mean.function = ExpMean, dispersion.function = LogVMR, 
                                  do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
tech <- ScaleData(object = tech, vars.to.regress = c("nUMI", "percent.mito"))

#saveRDS(tech, file = '~/R/Projects/Seurat/Robj/a_ddDropF/DropSeq_HEK_UCH1.Robj')
#saveRDS(tech, file = '~/R/Projects/Seurat/Robj/a_ddDropF/ddSeq_HEK_UCH1.Robj')

##########################################################################################################
# Alignment Section
### Fluidigm and DropSeq
##########################################################################################################
Fluidigm <- readRDS(file = '~/R/Projects/Seurat/Fluidigm/Fluidigm_human.Robj')
DropSeq <- readRDS(file='~/R/Projects/Seurat/Robj/a_ddDropF/DropSeq_HEK_UCH1.Robj')
ddSeq <- readRDS(file='~/R/Projects/Seurat/Robj/a_ddDropF/ddSeq_HEK_UCH1.Robj')

# Alignment relies on variable genes to determine source of variation 
hvg.Fluidigm.group <- rownames(x = head(x = Fluidigm@hvg.info, n = 2000)) # take 2000 genes with highest dispersion
hvg.DropSeq.group <- rownames(x = head(x = DropSeq@hvg.info, n = 2000))
hvg.union <- union(x = hvg.Fluidigm.group, y = hvg.DropSeq.group)

Fluidigm@meta.data[, "protocol_a1"] <- "Fluidigm"
DropSeq@meta.data[, "protocol_a1"] <- "DropSeq"
#####################################################
a_F_Drop <- RunCCA(object = DropSeq, object2 = Fluidigm, genes.use = hvg.union)
p1 <- DimPlot(object = a_F_Drop, reduction.use = "cca", group.by = "protocol_a1", pt.size = 2,
              do.return = TRUE) + ggtitle('Fluidigm vs DropSeq CCA All Cells')
p2 <- VlnPlot(object = a_F_Drop, features.plot = "CC1", group.by = "protocol_a1", do.return = TRUE)
plot_grid(p1 + theme(text = element_text(size=20)), p2 + theme(text = element_text(size=20))) 
#####################################################
DimHeatmap(object = a_F_Drop, reduction.type = "cca", cells.use = 200, dim.use = 1:12, do.balanced = TRUE)
d <- 5
a_F_Drop <- CalcVarExpRatio(object = a_F_Drop, reduction.type = "pca", grouping.var = "protocol_a1", 
                                 dims.use = 1:d)
# We will NOT remove cells (as far as we know, these cells are the same)
backup_a <- a_F_Drop
#discard <- SubsetData(object = backup_a, subset.name = "var.ratio.pca", accept.high = 0.5)
#rownames(discard@meta.data) # discarded cells for proper alignment. Whats wrong with them? 
a_F_Drop <- SubsetData(object = a_F_Drop, subset.name = "var.ratio.pca", accept.low = 0.5)
a_F_Drop <- AlignSubspace(object = a_F_Drop, reduction.type = "cca", grouping.var = "protocol_a1", 
                               dims.align = 1:d)
p1 <- VlnPlot(object = a_F_Drop, features.plot = "ACC1", group.by = "protocol_a1", do.return = TRUE)
p2 <- VlnPlot(object = a_F_Drop, features.plot = "ACC2", group.by = "protocol_a1", do.return = TRUE)
plot_grid(p1, p2)
#####################################################
a_F_Drop@meta.data[ ,"cellType"] <- NULL
a_F_Drop@meta.data[grep("HEK", a_F_Drop@meta.data$orig.ident),"cellType"] <- "HEK"
a_F_Drop@meta.data[grep("UCH1", a_F_Drop@meta.data$orig.ident),"cellType"] <- "UCH1"
a_F_Drop@meta.data[ ,"class"] <- paste(a_F_Drop@meta.data$protocol_a1, a_F_Drop@meta.data$cellType, sep="-")

a_F_Drop <- RunTSNE(object = a_F_Drop, reduction.use = "cca.aligned", dims.use = 1:d, do.fast = TRUE)
p1 <- TSNEPlot(object = a_F_Drop, group.by = "class", do.return = TRUE, pt.size = 2)
p1 + ggtitle('Aligned Fluidigm and DropSeq (Group by cellType) tSNE v2.1') + theme(text=element_text(size=15))

a_F_Drop <- SetAllIdent(a_F_Drop, id="class") # reset identities
VlnPlot(object = a_F_Drop, features.plot = c("nGene", "nUMI", "percent.mito"), x.lab.rot = TRUE)
g1 <- TSNEPlot(object = a_F_Drop, do.identify = TRUE) # bigger cluster that has UCH1
g2 <- TSNEPlot(object = a_F_Drop, do.identify = TRUE) # smaller cluster that has 293HEK
HEK_in_g1 <- g1[grep("HEK", a_F_Drop@ident[g1])] #cells from g1 (UCH1) that were actually HEK and thus misplaced. All Fluidigm
Flu_in_g2 <- g2[grep("Fluidigm", a_F_Drop@ident[g2])]
# 38/134. 28.35% of g1 Fluidigm cells were misappropriated as UCH1 when they were actually 293HEK
# 8/82. 9.76% of g2 Fluidigm cells were misappropriated as 293HEK when they were actually UCH1
# 46/216. 21.30% of Fluidigm cells were misidentified overall
a_F_Drop <- SetIdent(object = a_F_Drop, cells.use = HEK_in_g1, ident.use = "HEK_in_g1")
a_F_Drop <- SetIdent(object = a_F_Drop, cells.use = Flu_in_g2, ident.use = "Flu_in_g2")
mark <- FindMarkers(object = a_F_Drop, ident.1 = 'HEK_in_g1', ident.2 = 'Flu_in_g2', min.pct = 0.25)


