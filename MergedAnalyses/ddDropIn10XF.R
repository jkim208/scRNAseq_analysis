# ddDropIn10XF.R
# Updated: 26-March-2018
# Examine UCH1 samples across all available protocols via merge
# Seurat v2.2.1
library(Seurat)
library(dplyr)
setwd("~/Slim")

############################################################################################
# 1. Load in individual objects and merge based on technology (Skip to #2 if merged object exists)
############################################################################################
s.ddSeq <- readRDS("Robj/ddSeq/ddSeq_UCH1.Robj")
s.10X <- readRDS("Robj/10X/10X_UCH1_Droplike.Robj")
s.DropSeq <- readRDS("Robj/DropSeq/DropSeq_UCH1.Robj")
s.inDrop <- readRDS("Robj/inDrop/inDrop_lib21_UCH1.Robj")
s.Fluidigm <- readRDS("Robj/Fluidigm/Fluidigm_UCH1.Robj")

#########################

techComp <- MergeSeurat(object1 = s.ddSeq, object2 = s.inDrop, 
                       add.cell.id1 = "ddSeq", add.cell.id2 = "inDrop")
techComp <- MergeSeurat(object1 = techComp, object2 = s.DropSeq, 
                       add.cell.id2 = "DropSeq")
techComp <- MergeSeurat(object1 = techComp, object2 = s.10X, 
                       add.cell.id2 = "10X")
techComp <- MergeSeurat(object1 = techComp, object2 = s.Fluidigm, 
                       add.cell.id2 = "Fluidigm", project="techComp")
saveRDS(techComp, file="Robj/techComp.Robj")

############################################################################################
# 2. Read in merged object and begin meta-analysis
############################################################################################
techComp <- readRDS("Robj/techComp.Robj")

techComp <- NormalizeData(object = techComp, normalization.method = "LogNormalize", scale.factor = 10000)
techComp <- FindVariableGenes(object = techComp, mean.function = ExpMean, dispersion.function = LogVMR, 
                             do.plot = FALSE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
backup <- techComp
techComp <- ScaleData(object = techComp) #vars.to.regress = c("nUMI", "percent.mito")
#techComp <- ScaleData(object = techComp, vars.to.regress = c("nUMI", "percent.mito"))
techComp <- RunPCA(object = techComp, pc.genes = techComp@var.genes, pcs.compute = 20, genes.print = 50)
PCAPlot(techComp)
PCHeatmap(object = techComp, pc.use = 1, cells.use = 50, num.genes = 50, do.balanced = FALSE, label.columns = TRUE)
x <- PCHeatmap(object = techComp, pc.use = 1, cells.use = 50, num.genes = 50, do.balanced = FALSE, label.columns = FALSE, do.return = TRUE)
techComp <- SetIdent(techComp, cells.use = colnames(x), ident.use = "Top_PC_genes")
techComp <- SetIdent(techComp, cells.use = techComp@cell.names[!techComp@cell.names %in% colnames(x)], ident.use = "Not_Selected")

techComp <- SetAllIdent(object = techComp, id = "orig.ident")
techComp <- SetAllIdent(object = techComp, id = "tech")

PCElbowPlot(techComp)
techComp <- RunTSNE(object = techComp, dims.use = 1:12, do.fast = TRUE, perplexity = 30)
t <- TSNEPlot(object = techComp, do.return = TRUE, pt.size = 1.5)
t + ggtitle('ddSeq-DropSeq-InDrop UCH1/UCH2 tSNE 12dim 30perplexity v2.1') + theme(text = element_text(size=15))

#####################
# Add tech and cell type labeling to Seurat object 
f1 <- rep(NA, length(techComp@meta.data$orig.ident))
f1[grep("ddSeq_", techComp@meta.data$orig.ident)] <- "ddSeq"
f1[grep("DropSeq_", techComp@meta.data$orig.ident)] <- "DropSeq"
f1[grep("10X_", techComp@meta.data$orig.ident)] <- "10X"
f1[grep("iD_", techComp@meta.data$orig.ident)] <- "inDrop"
f1[grep("Fluidigm", techComp@meta.data$orig.ident)] <- "Fluidigm"
techComp@meta.data$tech <- f1

techComp@meta.data$cellType <- "UCH1"
