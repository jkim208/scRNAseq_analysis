# Goal: Examine DS52 subclusters (apparent from ddDropF Analysis)
# Violin plot of the subclusters suggests differences in expression pattern associated with cell proliferation

library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4/Seurat2.1')
library(Matrix)
library(dplyr)

setwd("~/R/Projects/Seurat")
DropSeq_UCH1_DS52_RT_ID <- readRDS("Robj/DropSeq/DropSeq_UCH1_RT_ID_DS52.Robj")
DropSeq_UCH1_DS52_HS_TB <- readRDS("Robj/DropSeq/DropSeq_UCH1_HS_TB_DS52.Robj")
DS52 <- MergeSeurat(object1 = DropSeq_UCH1_DS52_RT_ID, object2 = DropSeq_UCH1_DS52_HS_TB, 
                    add.cell.id1 = "DS52_RTID", add.cell.id2 = "DS52_HSTB", project = "DropSeq")
DS52 <- NormalizeData(object = DS52, normalization.method = "LogNormalize", scale.factor = 10000)
DS52 <- FindVariableGenes(object = DS52, mean.function = ExpMean, dispersion.function = LogVMR, 
                             do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = DS52@var.genes)
DS52 <- ScaleData(object = DS52, vars.to.regress = c("nUMI", "percent.mito"))

DS52 <- RunPCA(object = DS52, pc.genes = DS52@var.genes, 
                  do.print = TRUE, pcs.print = 1:12, genes.print = 5, pcs.compute = 20)
PCElbowPlot(object = DS52)
PCAPlot(object = DS52, dim.1 = 1, dim.2 = 2)
DS52 <- ProjectPCA(object = DS52, do.print = FALSE)
PCHeatmap(object = DS52, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
#DS52 <- JackStraw(object = DS52, num.replicate = 50, do.print = FALSE)
#JackStrawPlot(object = DS52, PCs = 1:20)

DS52 <- FindClusters(object = DS52, reduction.type = "pca", dims.use = 1:5, 
                        resolution = 0.4, save.SNN = TRUE, force.recalc = TRUE)
DS52 <- SetAllIdent(object = DS52, id = "orig.ident")
DS52 <- RunTSNE(object = DS52, dims.use = 1:5, do.fast = TRUE, perplexity = 30)
t <- TSNEPlot(object = DS52, do.return = TRUE, do.label=TRUE, pt.size = 3, label.size = 10)
t + ggtitle('DS52 tSNE 5dim 30perplexity 7.8k varGenes Seurat2.1') + theme(text = element_text(size=20))

FeaturePlot(object = DS52, features.plot = "T", cols.use = c("dark grey", "red"), reduction.use = "tsne", 
            dark.theme = TRUE, pt.size = 3)

DS52 <- SetIdent(object = DS52, cells.use = DS52_mini, ident.use = "DS52_mini")
DS52 <- SetIdent(object = DS52, cells.use = DS52_big, ident.use = "DS52_big")

VlnPlot(object = DS52, features.plot = c("nGene", "nUMI", "percent.mito"), 
        ident.include = c("DS52_big", "DS52_mini"), x.lab.rot = 1) 
df <- DS52@meta.data
df <- df[order(df$nUMI),]
ggplot(df, aes(nUMI, index)) + geom_line() + geom_hline(yintercept = 940)

DS52.mark <- FindMarkers(object = DS52, ident.1 = 'DS52_big', ident.2 = 'DS52_mini', min.pct = 0.25)
write.table(rownames(DS52.mark), col.names=FALSE, quote=FALSE, row.names=FALSE,
            file ='~/R/Projects/Seurat/Meeting5/DS52/DS52_subclusters_top_markers.txt')
saveRDS(DS52, file = '~/R/Projects/Seurat/Robj/DS52_subclusters.Robj')
DS52 <- readRDS(file = '~/R/Projects/Seurat/Robj/DS52_subclusters.Robj')


#########################################################################
# Alignment of DS52 runs
# Alignment relies on variable genes to determine source of variation 
hvg.DS52_RT_ID.group <- rownames(x = head(x = DropSeq_UCH1_DS52_RT_ID@hvg.info, n = 2000)) # take 2000 genes with highest dispersion
hvg.DS52_HS_TB.group <- rownames(x = head(x = DropSeq_UCH1_DS52_HS_TB@hvg.info, n = 2000))
hvg.union <- union(x = hvg.DS52_RT_ID.group, y = hvg.DS52_HS_TB.group)

DropSeq_UCH1_DS52_RT_ID@meta.data[, "protocol"] <- "DropSeq_UCH1_DS52_RT_ID"
DropSeq_UCH1_DS52_HS_TB@meta.data[, "protocol"] <- "DropSeq_UCH1_DS52_HS_TB"

#####################################################

align_dd_drop <- RunCCA(object = DropSeq_UCH1_DS52_RT_ID, 
                        object2 = DropSeq_UCH1_DS52_HS_TB, genes.use = hvg.union)
p1 <- DimPlot(object = align_dd_drop, reduction.use = "cca", group.by = "protocol", pt.size = 2,
              do.return = TRUE) + ggtitle('DS52_RT_ID vs DS52_HS_TB CCA All Cells')
p2 <- VlnPlot(object = align_dd_drop, features.plot = "CC1", group.by = "protocol", do.return = TRUE)
plot_grid(p1 + theme(text = element_text(size=20)), p2 + theme(text = element_text(size=20))) 
#####################################################

DimHeatmap(object = align_dd_drop, reduction.type = "cca", cells.use = 200, dim.use = 1:12, 
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
                                 dims.use = 1:5)

# We will NOT remove cells (as far as we know, these cells are the same)
align_dd_drop <- SubsetData(object = align_dd_drop, subset.name = "var.ratio.pca", accept.low = 0.5)

align_dd_drop <- AlignSubspace(object = align_dd_drop, reduction.type = "cca", grouping.var = "protocol", 
                               dims.align = 1:5)
p1 <- VlnPlot(object = align_dd_drop, features.plot = "ACC1", group.by = "protocol", 
              do.return = TRUE)
p2 <- VlnPlot(object = align_dd_drop, features.plot = "ACC2", group.by = "protocol", 
              do.return = TRUE)
plot_grid(p1, p2)

#####################################################

align_dd_drop <- RunTSNE(object = align_dd_drop, reduction.use = "cca.aligned", 
                         dims.use = 1:5, do.fast = TRUE)
p1 <- TSNEPlot(object = align_dd_drop, group.by = "protocol", do.return = TRUE, pt.size = 3)
p1 + ggtitle('Aligned DropSeq_UCH1_DS52_RT_ID vs DropSeq_UCH1_DS52_HS_TB (Group by Protocol) tSNE v2.1')

saveRDS(align_dd_drop, file = '~/R/Projects/Seurat/Robj/align_DS52.Robj')
#align_dd_drop <- readRDS(file = '~/R/Projects/Seurat/Robj/align_DS52.Robj')

DS52_mini <- sub("DS52_RTID_", "", DS52_mini); DS52_mini <- sub("DS52_HSTB_", "", DS52_mini)
DS52_big <- sub("DS52_RTID_", "", DS52_big); DS52_big <- sub("DS52_HSTB_", "", DS52_big)
DS52_mini <- DS52_mini[DS52_mini %in% rownames(align_dd_drop@meta.data)]
DS52_big <- DS52_big[DS52_big %in% rownames(align_dd_drop@meta.data)]
align_dd_drop <- SetIdent(object = align_dd_drop, cells.use = DS52_mini, ident.use = "DS52_mini")
align_dd_drop <- SetIdent(object = align_dd_drop, cells.use = DS52_big, ident.use = "DS52_big")

TSNEPlot(object = align_dd_drop, do.return = TRUE, pt.size = 3) +
  ggtitle('Aligned DS52 runs tSNE v2.1') + theme(text = element_text(size=20))

FeaturePlot(object = align_dd_drop, features.plot = c("T"), cols.use = c("grey", "red"), 
            reduction.use = "tsne", pt.size = 3, dark.theme = TRUE)

