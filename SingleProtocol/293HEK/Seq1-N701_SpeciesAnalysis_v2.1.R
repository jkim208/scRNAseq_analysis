library(dplyr)
library(ggplot2)
library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4/Seurat2.1')
library(Matrix)

setwd("~/R/Projects/Seurat/ddSeq")
set.seed(42)

####
# 2018 work
##################################################################################
###########################################################
# Mixed Species Analysis
human <- read.table("Seq1-N701_293/Input/Sassi-Seq1-N701_DGE_Human_Selected_Summary.txt", header=T) ##
human <-human %>% mutate(sample="human")
mouse <- read.table("Seq1-N701_293/Input/Sassi-Seq1-N701_DGE_Mouse_Selected_Summary.txt", header=T) ##
mouse <- mouse %>% mutate(sample="mouse")
combined <- full_join(human,mouse, by="CELL_BARCODE")

# '.x' means human. '.y' means mouse
combined <- combined %>% mutate(RatioH = NUM_TRANSCRIPTS.x/(NUM_TRANSCRIPTS.x+NUM_TRANSCRIPTS.y))
# arbitrary definition of what constitutes a mix
combined <- mutate(combined, Species = 
                     ifelse(RatioH > 0.8 , "Human", ifelse(RatioH > 0.2 & RatioH < 0.8, "Mixed", "Mouse")))

combined <- combined %>% mutate (x = ifelse(Species=="Human", NUM_TRANSCRIPTS.x, ifelse (Species=="Mixed", NUM_TRANSCRIPTS.x, 0)))
combined <- combined %>% mutate (y = ifelse (Species=="Mixed", NUM_TRANSCRIPTS.y, ifelse (Species=="Mouse", NUM_TRANSCRIPTS.y, 0)))
combined %>% group_by(Species) %>% summarise(number=n())

cp_Palette <- c("#56B4E9", "#009E73", "firebrick2")

p <- ggplot(combined, aes (x,y , colour=(Species) )) + 
  geom_point() + xlab('Human Transcripts') + ylab('Mouse Transcripts') + 
  ggtitle('ddSeq: Seq 1 - N701',subtitle =  '293/3T3 Mixed-Species plot with top 390 cells') + 
  scale_colour_manual(values=cp_Palette) ##
p

humanCells <- as.character(combined[grep('Human', combined$Species),]$CELL_BARCODE)

DGE.data <- read.table(file = "Seq1-N701_293/Input/Sassi-Seq1-N701_DGE.txt", header = TRUE, row.names = 1)
DGE.data <- DGE.data[grepl('^homo_', rownames(DGE.data)),] # keep only human genes
rownames(DGE.data) <- sub('^homo_', '', rownames(DGE.data))  # remove ^homo tags

# Create Seurat object with raw data. We could filter based on cell/gene thresholds
# For some reason, the DGE.data table lacks some of the barcodes in humanCells
DGE <- CreateSeuratObject(raw.data = DGE.data[,intersect(humanCells, colnames(DGE.data))], 
                          min.cells = 3, min.genes = 200, project = "ddSeq_293HEK_Seq1_N701")
# Keep track of percentage mitochondrial gene content (to regress out later)
mito.genes <- grep(pattern = "[Mm][Tt]-", x = rownames(x = DGE@data), value = TRUE)
percent.mito <- Matrix::colSums(DGE@raw.data[mito.genes, ])/Matrix::colSums(DGE@raw.data)
DGE <- AddMetaData(object = DGE, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = DGE, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = DGE, gene1 = "nUMI", gene2 = "nGene")
GenePlot(object = DGE, gene1 = "nUMI", gene2 = "percent.mito")
# We filter out cells that have unique gene counts over 6,000 or less than 200
DGE <- FilterCells(object = DGE, subset.names = c("nGene", "percent.mito"), 
                   low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.01))
# Normalize data to a total of 10,000 molecules 
DGE <- NormalizeData(object = DGE,normalization.method = "LogNormalize", scale.factor = 10000)
# Focus on variable genes (default parameter settings for a 1e4 molecule normalization)
DGE <- FindVariableGenes(object = DGE, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = DGE@var.genes)
ddSeq_293HEK_Seq1_N701 <- DGE
saveRDS(ddSeq_293HEK_Seq1_N701, file='~/R/Projects/Seurat/Robj/NotScaled/ddSeq/ddSeq_293HEK_Seq1_N701.Robj')











####
# 2017 work
##################################################################################
###########################################################
# To do mixed species analysis via dropSeq pipeline, you must:
# 1. BAMTagHistogram to get total read counts of cells
# 2. Use R to find the top barcodes by graphing against a knee plot
# 3. FilterBAM to get mouse/human specific BAM files
# 4. DigitalExpression using the new BAM file and list of top barcodes
# 5. Use R to determine which cells belong to which species

##############################################################
# Find the Top Barcodes Up to a Specific Threshold
a=read.table("Seq1-N701_293/Input/Sassi-Seq1-N701_readcounts_1.6g.txt", 
             header=F, stringsAsFactors=F) ##
x=cumsum(a$V1)
x=x/max(x)
plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]",
     ylab="cumulative fraction of reads", xlim=c(1,450))
abline(v=390) ##
topBarcodes <- a[1:390,2] # Up to knee threshold ##
write.table(x=topBarcodes, file = "Sassi-Seq1-N701_top_barcodes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE) ##
### Send output file to server to make species-specific DGE
###############################################################

# Mixed Species Analysis
human <- read.table("Seq1-N701_293/Input/Sassi-Seq1-N701_DGE_Human_Selected_Summary.txt", header=T) ##
human <-human %>% mutate(sample="human")
mouse <- read.table("Seq1-N701_293/Input/Sassi-Seq1-N701_DGE_Mouse_Selected_Summary.txt", header=T) ##
mouse <- mouse %>% mutate(sample="mouse")
combined <- full_join(human,mouse, by="CELL_BARCODE")

# '.x' means human. '.y' means mouse
combined <- combined %>% mutate(RatioH = NUM_TRANSCRIPTS.x/(NUM_TRANSCRIPTS.x+NUM_TRANSCRIPTS.y))
# arbitrary definition of what constitutes a mix
combined <- mutate(combined, Species = 
                     ifelse(RatioH > 0.8 , "Human", ifelse(RatioH > 0.2 & RatioH < 0.8, "Mixed", "Mouse")))

combined <- combined %>% mutate (x = ifelse(Species=="Human", NUM_TRANSCRIPTS.x, ifelse (Species=="Mixed", NUM_TRANSCRIPTS.x, 0)))
combined <- combined %>% mutate (y = ifelse (Species=="Mixed", NUM_TRANSCRIPTS.y, ifelse (Species=="Mouse", NUM_TRANSCRIPTS.y, 0)))
combined %>% group_by(Species) %>% summarise(number=n())

cp_Palette <- c("#56B4E9", "#009E73", "firebrick2")

p <- ggplot(combined, aes (x,y , colour=(Species) )) + 
  geom_point() + xlab('Human Transcripts') + ylab('Mouse Transcripts') + 
  ggtitle('ddSeq: Seq 1 - N701',subtitle =  '293/3T3 Mixed-Species plot with top 390 cells') + 
  scale_colour_manual(values=cp_Palette) ##
p

humanCells <- as.character(combined[grep('Human', combined$Species),]$CELL_BARCODE)
mouseCells <- combined[grep('Mouse', combined$Species),]$CELL_BARCODE
mixedCells <- combined[grep('Mixed', combined$Species),]$CELL_BARCODE
# end of external data work
##################################################################################

# Load the DGE datasets
DGE.data <- read.table(file = "Seq1-N701_293/Input/Sassi-Seq1-N701_DGE.txt", header = TRUE, row.names = 1)
#mus_sum <- colSums(DGE.data[18501:32961,]) # mouse genes selected
#homo_sum <- colSums(DGE.data[1:18500,]) # human genes selected
#df <- data.frame(homo_sum, mus_sum, homo_sum-mus_sum)

# DGE.data <- DGE.data[, !(names(DGE.data) %in% mixedCells )]
# DGE.data <- DGE.data[, !(names(DGE.data) %in% mouseCells )] 
DGE.data <- DGE.data[grepl('^homo_', rownames(DGE.data)),] # keep only human genes
rownames(DGE.data) <- sub('^homo_', '', rownames(DGE.data))  # remove ^homo tags

#####################################
# This next section may not be correct
# Focuses in on genes with both human and mouse versions
# aggregates expression values so that each gene row will represent both species
#####################################
# Remove the human and mouse tags on the genes, and add to the dataframe as a new column
DGE.rownames <- tolower(rownames(DGE.data))
untaggedGenes <- sub('^homo_|^mus_', '', DGE.rownames)  
duplicatedGenes <- duplicated(untaggedGenes)  # look for any instance of 2 of the same genes
sum(duplicatedGenes)  # how many genes are seen in both mice and humans from the data
duplicatedGenesIndices <- grep(TRUE, duplicatedGenes)  # grab the indices of gene copy pairs

# aggregate gene pairs and sum up corresponding expression values
l <- sum(duplicatedGenes)
## each rowPair will have the indices of a gene that has a mouse and human counterpart
rowPair <- grep(untaggedGenes[duplicatedGenesIndices[1]], untaggedGenes, fixed=TRUE)
intersectingGenes <- DGE.data[rowPair[1],] + DGE.data[rowPair[2],]

for (i in 2:l) {
  rowPair <- grep(untaggedGenes[duplicatedGenesIndices[i]], untaggedGenes, fixed=TRUE)
  intersectingGenes <- rbind(intersectingGenes, DGE.data[rowPair[1],] + DGE.data[rowPair[2],])
} 
save(intersectingGenes, file = "~/R/Projects/Seurat/Robj/Seq1-N701_intersectingGenes_v2.1.Robj")
rownames(intersectingGenes) <- sub('^homo_|^mus_', '', rownames(intersectingGenes))
#####################################
# End section
#####################################
DGE <- CreateSeuratObject(raw.data = DGE.data[1:homo_l, ], 
                               min.cells = 3, min.genes = 200)
## Prepare Seurat Object
# Create Seurat object with raw data. We could filter based on cell/gene thresholds
DGE <- CreateSeuratObject(raw.data = DGE.data, min.cells = 3, min.genes = 200)

# Keep track of percentage mitochondrial gene content (to regress out later)
mito.genes <- grep(pattern = "[Mm][Tt]-", x = rownames(x = DGE@data), value = TRUE)
percent.mito <- Matrix::colSums(DGE@raw.data[mito.genes, ])/Matrix::colSums(DGE@raw.data)
DGE <- AddMetaData(object = DGE, metadata = percent.mito, col.name = "percent.mito")

# Track each cells' species identity
cell.mouse <- rep('Mouse Cell', length(mouseCells))
names(cell.mouse) <- mouseCells 
cell.human <- rep('Human Cell', length(humanCells))
names(cell.human) <- humanCells 
cell.mixed <- rep('Mixed Cell', length(mixedCells))
names(cell.mixed) <- mixedCells 
cell.species <- c(cell.mouse, cell.human, cell.mixed)
DGE <- AddMetaData(object = DGE, metadata = cell.species, col.name = "cell.species")

VlnPlot(object = DGE, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# We filter out cells that have unique gene counts over 6,000 or less than 200
DGE <- FilterCells(object = DGE, subset.names = c("nGene", "percent.mito"), 
                   low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.05))
# Normalize data to a total of 10,000 molecules 
DGE <- NormalizeData(object = DGE,normalization.method = "LogNormalize", scale.factor = 10000)
# Focus on variable genes (default parameter settings for a 1e4 molecule normalization)
DGE <- FindVariableGenes(object = DGE, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
# Regress uninteresting signals out of analysis. nUMI: number of detected molecules per cell
DGE <- ScaleData(object = DGE, vars.to.regress = c("nUMI", "percent.mito"))

ddSeq_293_Seq1_N701 <- DGE
save(ddSeq_293_Seq1_N701, file='~/R/Projects/Seurat/Robj/ddSeq/ddSeq_293_Seq1_N701_v2.1.Robj')

length(x = DGE@var.genes)
# Run Dimension reduction
DGE <- RunPCA(object = DGE, pc.genes = DGE@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 30)
# visualize genes and cells that define the PCA
VizPCA(object = DGE, pcs.use = 1:2)
PCAPlot(object = DGE, dim.1 = 1, dim.2 = 2)
DGE <- ProjectPCA(object = DGE, do.print = FALSE)
PCHeatmap(object = DGE, pc.use = 1, cells.use = 200, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = DGE, pc.use = 1:12, cells.use = 200, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
# How many components will we keep?
DGE <- JackStraw(object = DGE, num.replicate = 100, do.print = FALSE)
JackStrawPlot(object = DGE, PCs = 1:20)
PCElbowPlot(object = DGE, num.pc = 30)
# Determine cluster based on components used
DGE <- FindClusters(object = DGE, reduction.type = "pca", dims.use = 1:14, 
                    resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
#PrintFindClustersParams(object = DGE)
DGE <- RunTSNE(object = DGE, dims.use = 1:14, do.fast = TRUE, perplexity = 30)
TSNEPlot(object = DGE, do.return=TRUE)+ ggtitle('ddSeq Seq1-N701 Human tSNE v2.1')
#save(DGE, file = "Seq1-N701_Species_post-tSNE_v2.1.Robj")

# Find biomarkers
DGE.markers <- FindAllMarkers(object = DGE, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)
DGE.markers <- DGE.markers[order(DGE.markers$p_val_adj),]
topDEgenes <- DGE.markers$gene[which(DGE.markers$p_val_adj < 0.05)]
topDEgenes <- sub('^homo_|^mus_', '', topDEgenes)  
write.table(topDEgenes, 
            file = "~/R/Projects/Seurat/Meeting4/ddSeq1-N701_markers.txt", 
            quote=FALSE, row.names=FALSE, col.names=FALSE)

top5 <- DGE.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
DoHeatmap(object = DGE, genes.use = top5$gene, slim.col.label = TRUE, remove.key = TRUE)

DGE.markers <- DGE.markers[order(DGE.markers$p_val_adj, decreasing=FALSE),]
DGE.markers <- DGE.markers[which(DGE.markers$p_val_adj < 0.05),]
c1 <- DGE.markers[DGE.markers$cluster==1,]
c0 <- DGE.markers[DGE.markers$cluster==0,]
VlnPlot(object = DGE, features.plot = c0$gene[1:9])

FeaturePlot(object = DGE, features.plot = c("GAPDH"), no.legend = FALSE, 
            min.cutoff = "q10", max.cutoff = "q90", cols.use = c("grey", "red"))
DotPlot(object = DGE, genes.plot = 'GAPDH', plot.legend = TRUE)
write.csv(DGE.markers$gene, file='test.csv', row.names=FALSE, col.names = FALSE, quote=FALSE)

#####################################################
# For cross-species analysis
### Compare species differences
DGE <- StashIdent(object = DGE, save.name = "CellType")
# Next, switch the identity class of all cells to reflect replicate ID
DGE <- SetAllIdent(object = DGE, id = "cell.species")
t <- TSNEPlot(object = DGE, do.return=TRUE)
t + ggtitle('ddSeq Seq1-N701 v2.1 Species Analysis')
#####################################################

#####################################################
# For human only analysis

chrY <- read.table("~/R/Projects/Genes_chrY.txt")
DGE.markers$gene <- sub('^homo_', '', DGE.markers$gene)  
DGE.markers[DGE.markers$gene %in% chrY[,1],]

#####################################################
x <- matrix(rnorm(1000*6),1000,6)
design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))
# First set of 20 genes are genuinely differentially expressed
index1 <- 1:20
x[index1,4:6] <- x[index1,4:6]+1
camera(x, index1, design)



clus_id <- DGE@meta.data$res.0.6
names(clus_id) <- rownames(DGE@meta.data)
clus_id <- clus_id[order(names(clus_id))]

y <- DGE@data[,colnames(DGE.data) %in% names(clus_id)]
y <- y[,order(colnames(y))]

design <- cbind(Intercept=1,Group= as.numeric(clus_id))


# First set of 20 genes are genuinely differentially expressed
index1 <- which(rownames(y) %in% DGE.markers$gene)

# Second set of 20 genes are not DE
index2 <- sample(18500, 1230, replace = TRUE)


camera(y, index1, design)
camera(y, index2, design)

#########################################################

clus_id <- DGE@meta.data$res.0.6
names(clus_id) <- rownames(DGE@meta.data)
clus_id <- clus_id[order(names(clus_id))]

y <- DGE.data[,colnames(DGE.data) %in% names(clus_id)]
y <- y[,order(colnames(y))]

zero_i <- as.numeric(which(clus_id=="0"))
one_i <- as.numeric(which(clus_id=="1"))

z <- data.frame(Cluster0 = rowMeans(y[, zero_i]))
z$Cluster1 <- unname(rowMeans(y[, one_i]))

design <- cbind(Intercept=1, Group= c(0,1))


# First set of 20 genes are genuinely differentially expressed
index1 <- which(rownames(y) %in% DGE.markers$gene)

# Second set of 20 genes are not DE
index2 <- sample(18500, 1230, replace = TRUE)

#########################################################
#########################################################
z <- read.table("~/R/Projects/cellCycleGenes.csv")
design <- cbind(Intercept=1,Group= as.numeric(DGE@meta.data$res.0.6))

for (i in 1:7) {
  print(i)
  print(camera(DGE@data, which(rownames(DGE@data) %in% z[z$V2 == i,]$V1), design))
}



