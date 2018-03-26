# TCorGenes_cellLines.R
# Updated: 26-March-2018
# Find genes that have similar expression patterns to T in chordoma cells
# Seurat v2.2.1

library(Seurat)
library(Matrix)
library(dplyr)
setwd("~/R/Slim")

############################################################
# 1. Prepare chordoma merged seurat objects (skip if already completed)
##############################
ddSeq2_N706 <- readRDS("Robj/ddSeq/ddSeq_UCH1_Seq2_N706.Robj")
ddSeq3_N707 <- readRDS("Robj/ddSeq/ddSeq_UCH1_Seq3_N707.Robj")
ddSeq1_N704 <- readRDS("Robj/ddSeq/ddSeq_UCH2_Seq1_N704.Robj")

DropSeq_DS50 <- readRDS("Robj/DropSeq/DropSeq_UCH1_DS50.Robj")
DropSeq_DS52HSTB <- readRDS("Robj/DropSeq/DropSeq_UCH1_DS52_HS_TB.Robj")
DropSeq_DS52RTID <- readRDS("Robj/DropSeq/DropSeq_UCH1_DS52_RT_ID.Robj")
DropSeq_DS49 <- readRDS("Robj/DropSeq/DropSeq_UCH2_DS49.Robj")

inDrop_UCH1 <- readRDS("Robj/inDrop/iD_lib21_UCH1.Robj")
inDrop_UCH2 <- readRDS("Robj/inDrop/iD_lib22_UCH2.Robj")

t10X_UCH1 <- readRDS("Robj/10X/10X_UCH1.Robj")
t10X_UCH2 <- readRDS("Robj/10X/10X_UCH2.Robj")
###############################

ddSeq_UCH1 <- MergeSeurat(object1 = ddSeq2_N706, object2 = ddSeq3_N707,
                           add.cell.id1 = "ddSeq2_N706", add.cell.id2 = "ddSeq3_N707")
ddSeq_UCH1 <- MergeSeurat(object1 = ddSeq_UCH1, object2 = DropSeq_DS52RTID,
                            project = "DropSeq_UCH1_adj", add.cell.id2 = "DropSeq_DS52RTID")
#VlnPlot(object = ddSeq_UCH1, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
#sample <- FilterCells(object = sample, subset.names = "nGene", high.thresholds = 1150)
ddSeq_UCH1 <- NormalizeData(object = ddSeq_UCH1, normalization.method = "LogNormalize", scale.factor = 10000)
ddSeq_UCH1 <- FindVariableGenes(object = ddSeq_UCH1, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot=FALSE)
#saveRDS(ddSeq_UCH1, file="Robj/ddSeq/ddSeq_UCH1.Robj")

###############################################################################################
# 2. Apply correlation function to find correlation between T and each gene across all cells in one protocol
##############################
t10X <- readRDS(file="Robj/10X/10X_chor.Robj")
drop <- readRDS(file="Robj/DropSeq/DropSeq_chor.Robj")
dd <- readRDS(file="Robj/ddSeq/ddSeq_chor.Robj")
id <- readRDS(file="Robj/inDrop/inDrop_chor.Robj")

sample <- readRDS("Robj/10X/10X_UCH1.Robj")
###############################
# Gather statistics data from dataset
corT <- apply(sample@data, 1, function (x) round(cor(sample@data["T",], x), digits = 4))
exprCellCount <- apply(sample@data, 1, function (x) sum(x > 0))
avgExpr <- apply(sample@data, 1, function (x) round( mean(x[x > 0]), digits = 4))
sdExpr <- apply(sample@data, 1, function (x) round( sd(x[x > 0]), digits = 4))

corT_df <- data.frame(corT, exprCellCount, avgExpr, sdExpr)

write.table(corT_df, file="saveData/corT/10X_UCH2_corT.csv", quote=FALSE, col.names = TRUE, sep="\t")

# Find genes with similar pattern to T. Consider correlation to T and degree of expression (percent of cells expressed)
t.cells <- sample@cell.names[sample@data["T",] > 0] # T-expressing cells (> 0 expression)
t.data <- as.data.frame(as.matrix(sample@data[, t.cells]))
z <- sort(apply(t.data, 1, function(x) sum(x > 0)), decreasing = TRUE) # cell count for each gene where the geneExp > 0
z <- z / length(t.cells) # list of genes sorted by decreasing number of expressed cells (of that gene)
sum(z > .8) # list of genes with expression presence in cells(which express T) above a threshold
Tcell_expressed_genes <- names(z)[z > .8] ##### save for later. Intersect with next list of genes.

tl <- length(t.cells)
t.filter <- tl + c(-0.2*tl, 0.2*tl) # high and low thresholds for expressed cell count (keep it near T cell count)
corT_DropSeq_UCH2$gene <- rownames(corT_DropSeq_UCH2)
t.df <-corT_DropSeq_UCH2 %>%
  select(gene, everything()) %>% # move gene column to front of data frame
  filter(exprCellCount >= t.filter[1]) %>%
  filter(exprCellCount <= t.filter[2]) %>%
  arrange(desc(corT)) 
head(t.df, n=15)
intersect(head(t.df$gene, n = 100), Tcell_expressed_genes)
head(t.df[t.df$gene %in% intersect(head(t.df$gene, n = 100), Tcell_expressed_genes),], n=15)

##############################################################################
# 3.  Meta-analysis of corT values
### Look at how each gene compares to gene "T" with respect to mean exprs, total count, correlation, etc
corT_10X_UCH1 <- read.table(file="saveData/corT/10X_UCH1_corT.csv", sep="\t")
corT_10X_UCH2 <- read.table(file="saveData/corT/10X_UCH2_corT.csv", sep="\t")
corT_DropSeq_UCH1 <- read.table(file="saveData/corT/DropSeq_UCH1_corT.csv", sep="\t")
corT_DropSeq_UCH2 <- read.table(file="saveData/corT/DropSeq_UCH2_corT.csv", sep="\t")
corT_ddSeq_UCH1 <- read.table(file="saveData/corT/ddSeq_UCH1_corT.csv", sep="\t")
corT_ddSeq_UCH2 <- read.table(file="saveData/corT/ddSeq_UCH2_corT.csv", sep="\t")
corT_inDrop_UCH1 <- read.table(file="saveData/corT/inDrop_UCH1_corT.csv", sep="\t")
corT_inDrop_UCH2 <- read.table(file="saveData/corT/inDrop_UCH2_corT.csv", sep="\t")

# Function corT_dblFilter repeats analysis from section 2
corT_dblFilter <- function (sample, tech) {
  csv <- read.table(file=paste("saveData/corT/", sample, "_corT.csv", sep=""), sep="\t")
  robj <- readRDS(paste("Robj/", tech, "/", sample, ".Robj", sep=""))
  t.cells <- robj@cell.names[robj@data["T",] > 0] # T-expressing cells (> 0 expression)
  t.data <- as.data.frame(as.matrix(robj@data[, t.cells]))
  z <- sort(apply(t.data, 1, function(x) sum(x > 0)), decreasing = TRUE) # cell count for each gene where the geneExp > 0
  z <- z / length(t.cells) # list of genes sorted by decreasing number of expressed cells (of that gene)
  sum(z > .8) # list of genes with expression presence in cells(which express T) above a threshold
  Tcell_expressed_genes <- names(z)[z > .8] ##### save for later. Intersect with next list of genes.
  
  tl <- length(t.cells)
  t.filter <- tl + c(-0.2*tl, 0.2*tl) # high and low thresholds for expressed cell count (keep it near T cell count)
  csv$gene <- rownames(csv)
  t.df <- csv %>%
    select(gene, everything()) %>% # move gene column to front of data frame
    filter(exprCellCount >= t.filter[1]) %>%
    filter(exprCellCount <= t.filter[2]) %>%
    arrange(desc(corT)) 
  #head(t.df, n=15)
  #intersect(head(t.df$gene, n = 100), Tcell_expressed_genes)
  final.df <- t.df[t.df$gene %in% intersect(t.df$gene, Tcell_expressed_genes),]
}

df <- corT_dblFilter("10X_UCH1", "10X")

#############
# Downstream Visualization
topCor <- head(sort(corT, decreasing=TRUE), n=30)
topCor

sample <- ScaleData(object = sample)
sample <- RunPCA(object = sample, pc.genes = sample@var.genes, pcs.compute = 30)
PCElbowPlot(sample, num.pc = 30)

sample <- SetAllIdent(sample, id="orig.ident")

sample <- RunTSNE(object = sample, dims.use = 1:12, do.fast = TRUE, perplexity = 30)
TSNEPlot(sample, do.label=TRUE, do.return=TRUE) + ggtitle("10X_UCH1 tSNE 12d30p")

VlnPlot(sample, features.plot=names(topCor[1:9]))
FeaturePlot(sample, features.plot=names(topCor[1:9]), 
            cols.use=c("dark grey", "red"), reduction.use="tsne", dark.theme = TRUE,
            min.cutoff ="q10", max.cutoff = "q90", pt.size = .8)

sum(corT < .1 & corT > -.1, na.rm = TRUE) / nrow(sample@data) # q1: percent of genes where -0.1 < cor < 0.1
sum(corT > .3 | corT < -.3, na.rm = TRUE) # q2: number of genes with cor > 0.3 or cor < -0.3. Includes T. Should be always > 0.

###################################################################################


