# 10X_cellranger.R
# Updated: 26-March-2018
# Examine scRNA-seq data generated by the 10X platform and mapped by the cellranger pipeline
# Seurat v2.2.1

library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4')
library(Matrix);library(dplyr)
set.seed(42)
setwd("~/Slim")

##################################################################################
##################################################################################
# Prepare Seurat objects of 10X data
UCH1_path <- "./10X/sample1_UCH1/filtered_gene_bc_matrices/GRCh38"
UCH2_path <- "./10X/sample2_UCH2/filtered_gene_bc_matrices/GRCh38"
HEK_path  <- "./10X/sample3_293-3T3/filtered_gene_bc_matrices/hg19"
ERCC_path <- "./10X/sample4_ercc/raw_gene_bc_matrices/ercc92"

sample.data <- Read10X(data.dir = ERCC_path)
sample.seurat <- CreateSeuratObject(raw.data = sample.data, project = "10X_sample.seurat",
                                    min.cells = 3, min.genes = 200)
mito.genes <- grep(pattern = "MT-", x = rownames(x = sample.seurat@data), value = TRUE)
percent.mito <- Matrix::colSums(sample.seurat@raw.data[mito.genes, ])/Matrix::colSums(sample.seurat@raw.data)
sample.seurat <- AddMetaData(object = sample.seurat, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = sample.seurat, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = sample.seurat, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = sample.seurat, gene1 = "nUMI", gene2 = "nGene")
sample.seurat <- FilterCells(object = sample.seurat, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.06))
sample.seurat <- NormalizeData(object = sample.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
sample.seurat <- FindVariableGenes(object = sample.seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                             do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
sample.seurat <- ScaleData(object = sample.seurat, vars.to.regress = c("nUMI", "percent.mito"))

#saveRDS(sample.seurat, file='~/R/Projects/Seurat/Robj/10X/10X_s1_UCH1.Robj')
#saveRDS(sample.seurat, file='~/R/Projects/Seurat/Robj/10X/10X_s2_UCH2.Robj')
#saveRDS(sample.seurat, file='~/R/Projects/Seurat/Robj/10X/10X_s3_293HEK.Robj')
#saveRDS(sample.seurat, file='~/R/Projects/Seurat/Robj/10X/10X_s4_ercc.Robj')
#t10X_s1_UCH1 <- readRDS(file='~/R/Projects/Seurat/Robj/10X/10X_s1_UCH1.Robj')
#t10X_s2_UCH2 <- readRDS(file='~/R/Projects/Seurat/Robj/10X/10X_s2_UCH2.Robj')
#t10X_s3_293HEK <- readRDS(file='~/R/Projects/Seurat/Robj/10X/10X_s3_293HEK.Robj')
#t10X_s4_ercc <- readRDS(file='~/R/Projects/Seurat/Robj/10X/10X_s4_ercc.Robj')

so_10X <- RunPCA(object = so_10X, pc.genes = so_10X@var.genes, 
                  do.print = TRUE, pcs.print = 1:12, genes.print = 5, pcs.compute = 20)
PCAPlot(object = so_10X, dim.1 = 1, dim.2 = 2)
PCElbowPlot(object = so_10X)
so_10X <- RunTSNE(object = so_10X, dims.use = 1:16, do.fast = TRUE, perplexity = 30)
t <- TSNEPlot(object = so_10X, do.return = TRUE, do.label = TRUE, pt.size = 2)
t + ggtitle('10X tSNE Seurat_v2.1') + theme(text = element_text(size=15))

FeaturePlot(object = so_10X, features.plot = "T", cols.use = c("dark grey", "red"), 
            reduction.use = "tsne", dark.theme = TRUE)
so_10X@meta.data$cellType[grep("s1_UCH1_", so_10X@cell.names),] <- "UCH1"

##################################################################################
##################################################################################
# ERCC data analysis

# Visualize and summarize data
filtered_ercc <- Read10X(data.dir = "~/R/Projects/10X/sample4/raw_gene_bc_matrices/ercc92")
mean_per_ercc <- apply(filtered_ercc, 1, mean)
sd_per_ercc <- apply(filtered_ercc, 1, sd)

bc_transcripts <- colSums(as.matrix(filtered_ercc)) # use raw data for knee plot
bc_transcripts <- sort(bc_transcripts,decreasing = TRUE)
x <- cumsum(bc_transcripts)/max(cumsum(bc_transcripts))
plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]",
     ylab="cumulative fraction of reads", xlim=c(1, 100000)) # knee plot. Consider using raw data for this
abline(v=76000) # cellranger auto-filter at 74,235 cells
qplot(bc_transcripts, geom="histogram", binwidth = 100, main = "10X ERCC Histogram") # histogram
summary(bc_transcripts)
sd(bc_transcripts)
2641/100000 # mean number of transcripts divided by theoretical amount of ercc transcripts

t <- read.csv("~/Documents/ercc_spike.txt", sep="\t")
colnames(t) <- c("Re-sort ID", "ERCC ID", "subgroup", "concentration in Mix 1 (attomoles/ul)",   
                 "Molecules / uL", "Molecules / nL", "Diluted to 0.32% Molec / nL")
t$ercc_theoretical <- t$`Diluted to 0.32% Molec / nL`/2.5 # we will assume our molecule total is 50% of DropSeq's
t <- t[order(t$`ERCC ID`),]
t$ercc_observed <- mean_per_ercc
t$log2theoretical <- log2(t$ercc_theoretical)
t$log2observed <- log2(t$ercc_observed)

#plot(x=t$log2theoretical, y=t$log2observed)
ggplot(t,aes(log2theoretical,log2observed))+ geom_point() + 
  geom_smooth(method='lm',formula=y~x)
fit <- lm(log2observed ~ log2theoretical, data = t)
summary(fit)

##################################################################################
##################################################################################
# Merged analysis
so_10X <- MergeSeurat(object1 = t10X_s1_UCH1, object2 = t10X_s2_UCH2, project = "10X",
                     add.cell.id1 = "10X_s1_UCH1", add.cell.id2 = "10X_s2_UCH2")
so_10X <- MergeSeurat(object1 = so_10X, object2 = t10X_s3_293HEK, 
                     add.cell.id2 = "10X_s3_HEK")
so_10X <- NormalizeData(object = so_10X, normalization.method = "LogNormalize", scale.factor = 10000)
so_10X <- FindVariableGenes(object = so_10X, mean.function = ExpMean, dispersion.function = LogVMR, 
                             do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
so_10X <- ScaleData(object = so_10X)

so_10X <- RunPCA(object = so_10X, pc.genes = so_10X@var.genes, 
                 do.print = TRUE, pcs.print = 1:12, genes.print = 5, pcs.compute = 20)
so_10X <- RunTSNE(object = so_10X, dims.use = 1:16, do.fast = TRUE, perplexity = 30)
t <- TSNEPlot(object = so_10X, do.return = TRUE, do.label = TRUE, pt.size = 2)
t + ggtitle('10X tSNE Seurat_v2.1') + theme(text = element_text(size=15))

so_10X@meta.data$cellType <- as.character(so_10X@meta.data$orig.ident)
so_10X@meta.data[grep("s1_UCH1_", so_10X@cell.names), "cellType"] <- "UCH1"
so_10X@meta.data[grep("s2_UCH2_", so_10X@cell.names), "cellType"] <- "UCH2"
so_10X <- SetAllIdent(so_10X, id="cellType")

saveRDS(so_10X, "~/R/Projects/Seurat/Robj/so_10X_chordoma.Robj")
# Chordoma specific meta-analysis
ddSeq_chor <- readRDS(file='~/R/Projects/Seurat/Robj/ddDropF_chordoma/ddSeq_chor.Robj')
DropSeq_chor <- readRDS(file='~/R/Projects/Seurat/Robj/ddDropF_chordoma/DropSeq_chor.Robj')

ddDrop10X_chor <- MergeSeurat(object1 = so_10X, object2 = ddSeq_chor)
ddDrop10X_chor <- MergeSeurat(object1 = ddDrop10X_chor, object2 = DropSeq_chor, project = "ddDrop10X_chor")
saveRDS(ddDrop10X_chor, file='~/R/Projects/Seurat/Robj/ddDrop10X_chor.Robj')
ddDrop10X_chor <- NormalizeData(object = ddDrop10X_chor, normalization.method = "LogNormalize", scale.factor = 10000)
ddDrop10X_chor <- FindVariableGenes(object = ddDrop10X_chor, mean.function = ExpMean, dispersion.function = LogVMR, 
                                   do.plot = FALSE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
ddDrop10X_chor <- ScaleData(object = ddDrop10X_chor)

ddDrop10X_chor <- RunPCA(object = ddDrop10X_chor, pc.genes = ddDrop10X_chor@var.genes, 
                 do.print = TRUE, pcs.print = 1:12, genes.print = 5, pcs.compute = 20)
PCAPlot(object = ddDrop10X_chor, dim.1 = 1, dim.2 = 2)
PCElbowPlot(object = ddDrop10X_chor)
ddDrop10X_chor <- RunTSNE(object = ddDrop10X_chor, dims.use = 1:16, do.fast = TRUE, perplexity = 30)
t <- TSNEPlot(object = ddDrop10X_chor, do.return = TRUE, do.label = TRUE, pt.size = 2)
t + ggtitle('10X tSNE Seurat_v2.1') + theme(text = element_text(size=15))




