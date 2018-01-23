# Seurat Alignment of Different scRNAseq Technologies
# Generic labels are used to facilitate multiple analyses
# Attempt to align >2 datasets together
# 23-Jan-2018

library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4/Seurat2.1')
library(Matrix)
library(dplyr)
set.seed(42)
setwd("~/R/Projects/Seurat")

a_ddDrop <- readRDS('~/R/Projects/Seurat/Robj/a_ddDropF/ddDrop_HEK_UCH1.Robj')
a_ddDrop <- FindVariableGenes(object = a_ddDrop, mean.function = ExpMean, dispersion.function = LogVMR, 
                             do.plot = FALSE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
Fluidigm <- readRDS(file = '~/R/Projects/Seurat/Fluidigm/Fluidigm_human.Robj')

a <- a_ddDrop
b <- Fluidigm
##########################
a.group <- rownames(x = head(x = a@hvg.info, n = 2000)) # take 2000 genes with highest dispersion
b.group <- rownames(x = head(x = b@hvg.info, n = 2000))
union <- union(x = a.group, y = b.group)

# when doing multiple alignments, avoid overwriting the protocol group
a@meta.data[, "protocol2"] <- "a_ddDrop" #a@meta.data[, "protocol1"] <- "a_ddDrop"
b@meta.data[, "protocol2"] <- "Fluidigm" #b@meta.data[, "protocol1"] <- "Fluidigm"

ab_align <- RunCCA(object = a, object2 = b, genes.use = union)
p1 <- DimPlot(object = ab_align, reduction.use = "cca", group.by = "protocol", pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = ab_align, features.plot = "CC1", group.by = "protocol", do.return = TRUE)
plot_grid(p1, p2) 

DimHeatmap(object = ab_align, reduction.type = "cca", cells.use = 500, dim.use = 1:12, 
           do.balanced = TRUE)
d <- 8 # adjust based on Heatmap
ab_align <- CalcVarExpRatio(object = ab_align, reduction.type = "pca", grouping.var = "protocol", dims.use = 1:d)
ab_align_backup <- ab_align
ab_align <- SubsetData(object = ab_align, subset.name = "var.ratio.pca", accept.low = 0.50) # Will discard some cells
ab_align <- AlignSubspace(object = ab_align, reduction.type = "cca", grouping.var = "protocol", dims.align = 1:d)
ab_align <- RunTSNE(object = ab_align, reduction.use = "cca.aligned", dims.use = 1:d, do.fast = TRUE)
TSNEPlot(object = ab_align, group.by = "protocol", do.return = TRUE, pt.size = 2)

ab_align@meta.data[ ,"cellType"] <- NULL
ab_align@meta.data[grep("HEK", ab_align@meta.data$orig.ident),"cellType"] <- "HEK"
ab_align@meta.data[grep("UCH1", ab_align@meta.data$orig.ident),"cellType"] <- "UCH1"
ab_align@meta.data[grep("UCH2", ab_align@meta.data$orig.ident),"cellType"] <- "UCH2"
ab_align@meta.data[ ,"class"] <- paste(ab_align@meta.data$protocol, ab_align@meta.data$cellType, sep="-")
TSNEPlot(object = ab_align, group.by = "class", do.return = TRUE, pt.size = 2)

##########################################################################
# Save created Seurat objects
#saveRDS(ab_align, file = '~/R/Projects/Seurat/Robj/a_ddDropF/ddDrop_HEK_UCH1.Robj') # aligned ddSeq and DropSeq
#saveRDS(ab_align, file = '~/R/Projects/Seurat/Robj/a_ddDropF/a_ddDrop_F_HEK_UCH1.Robj') # aligned ddDrop and Fluidigm


