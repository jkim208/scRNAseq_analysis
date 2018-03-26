library(dplyr)
library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4/Seurat2.1')
setwd("~/R/Projects/Seurat/DropSeq")
set.seed(42)

# Determine which cells are Human (293HEK) from the mixed species analysis
human <- read.table("DS_9_HEK_MEF/server/HEK_MEF_NO_DGE_Human_Summary.txt", header=T) ##
human <- human %>% mutate(sample="human")
mouse <- read.table("DS_9_HEK_MEF/server/HEK_MEF_NO_DGE_Mouse_Summary.txt", header=T) ##
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

humanCells <- as.character(combined[grep('Human', combined$Species),]$CELL_BARCODE)
DGE.data <- read.table(file = "DS_9_HEK_MEF/server/HEK_MEF_NO_DGE.txt", header = TRUE, row.names = 1)
DGE.data <- DGE.data[grepl('^homo_', rownames(DGE.data)),] # keep only human genes
rownames(DGE.data) <- sub('^homo_', '', rownames(DGE.data))  # remove ^homo tags

# Create Seurat object with raw data. We could filter based on cell/gene thresholds
# For some reason, the DGE.data table lacks some of the barcodes in humanCells
DGE <- CreateSeuratObject(raw.data = DGE.data[,intersect(humanCells, colnames(DGE.data))], 
                          min.cells = 3, min.genes = 200, project = "DropSeq_293HEK_DS9_NO")
# Keep track of percentage mitochondrial gene content (to regress out later)
mito.genes <- grep(pattern = "[Mm][Tt]-", x = rownames(x = DGE@data), value = TRUE)
percent.mito <- Matrix::colSums(DGE@raw.data[mito.genes, ])/Matrix::colSums(DGE@raw.data)
DGE <- AddMetaData(object = DGE, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = DGE, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = DGE, gene1 = "nUMI", gene2 = "nGene")
GenePlot(object = DGE, gene1 = "nUMI", gene2 = "percent.mito")
# We filter out cells that have unique gene counts over 6,000 or less than 200
DGE <- FilterCells(object = DGE, subset.names = c("nGene", "percent.mito"), 
                   low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.2))
# Normalize data to a total of 10,000 molecules 
DGE <- NormalizeData(object = DGE,normalization.method = "LogNormalize", scale.factor = 10000)
# Focus on variable genes (default parameter settings for a 1e4 molecule normalization)
DGE <- FindVariableGenes(object = DGE, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = FALSE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
DropSeq_293HEK_DS9_NO <- DGE
saveRDS(DropSeq_293HEK_DS9_NO, file='~/R/Projects/Seurat/Robj/NotScaled/DropSeq/DropSeq_293HEK_DS9_NO.Robj')






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
a <- read.table("DS_9_HEK_MEF/server/HEK_MEF_NO_cell_readcounts.txt", 
             header=F, stringsAsFactors=F) ##
x <- cumsum(a$V1)
x <- x/max(x)
plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]",
     ylab="cumulative fraction of reads", xlim=c(1,1000))
abline(v=400) ##
topBarcodes <- a[1:400,2] # Up to knee threshold ##
write.table(x=topBarcodes, file = "DS9_HEK_MEF_NO_top_barcodes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE) ##
### Send output file to server to make species-specific DGE
###############################################################
# Mixed Species Analysis
human <- read.table("DS_9_HEK_MEF/server/HEK_MEF_NO_DGE_Human_Summary.txt", header=T) ##
human <-human %>% mutate(sample="human")
mouse <- read.table("DS_9_HEK_MEF/server/HEK_MEF_NO_DGE_Mouse_Summary.txt", header=T) ##
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
  ggtitle('DS_9_HEK_MEF',subtitle =  'Mixed-Species plot with top 200 cells') + 
  scale_colour_manual(values=cp_Palette) ##
p

humanCells <- combined[grep('Human', combined$Species),]$CELL_BARCODE
mouseCells <- combined[grep('Mouse', combined$Species),]$CELL_BARCODE
mixedCells <- combined[grep('Mixed', combined$Species),]$CELL_BARCODE
# end of external data work
##################################################################################


# Load the dataset. 
DGE.data <- read.table(file = "DS_9_HEK_MEF/server/HEK_MEF_NO_DGE.txt", header = TRUE, row.names = 1)
# DGE.data <- DGE.data[, !(names(DGE.data) %in% mixedCells )] # remove mixed
# DGE.data <- DGE.data[, !(names(DGE.data) %in% mouseCells )] # remove mouse cells for human analysis
# DGE.data <- DGE.data[grepl('^homo_', rownames(DGE.data)),] # human genes only

# Create Seurat object with raw data. We could filter based on cell/gene thresholds
DGE <- CreateSeuratObject(raw.data = DGE.data, min.cells = 3, min.genes = 200)

# Keep track of percentage mitochondrial gene content (to regress out later)
mito.genes <- grep(pattern = "_[Mm][Tt]-", x = rownames(x = DGE@data), value = TRUE)
percent.mito <- Matrix::colSums(DGE@raw.data[mito.genes, ])/Matrix::colSums(DGE@raw.data)
DGE <- AddMetaData(object = DGE, metadata = percent.mito, col.name = "percent.mito")

#### Only if looking at both human and mouse ####
# Track each cells' species identity
cell.mouse <- rep('Mouse Cell', length(mouseCells))
names(cell.mouse) <- mouseCells 
cell.human <- rep('Human Cell', length(humanCells))
names(cell.human) <- humanCells 
cell.mixed <- rep('Mixed Cell', length(mixedCells))
names(cell.mixed) <- mixedCells
cell.species <- c(cell.mouse, cell.human, cell.mixed)
DGE <- AddMetaData(object = DGE, metadata = cell.species, col.name = "cell.species")
#################################################

VlnPlot(object = DGE, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# We filter out cells that have unique gene counts over 5,000 or less than 200
DGE <- FilterCells(object = DGE, subset.names = c("nGene", "percent.mito"), 
                   low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.05))
# Normalize data to a total of 10,000 molecules 
DGE <- NormalizeData(object = DGE,normalization.method = "LogNormalize", scale.factor = 10000)
# Focus on variable genes (default parameter settings for a 1e4 molecule normalization)
DGE <- FindVariableGenes(object = DGE, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
# Regress uninteresting signals out of analysis. nUMI: number of detected molecules per cell
DGE <- ScaleData(object = DGE, vars.to.regress = c("nUMI", "percent.mito"))
length(x = DGE@var.genes)

# Run Dimension reduction
DGE <- RunPCA(object = DGE, pc.genes = DGE@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 20)
# visualize genes and cells that define the PCA
VizPCA(object = DGE, pcs.use = 1:2)
PCAPlot(object = DGE, dim.1 = 1, dim.2 = 2)
DGE <- ProjectPCA(object = DGE, do.print = FALSE)
#PCHeatmap(object = DGE, pc.use = 1, cells.use = 200, do.balanced = TRUE, label.columns = FALSE)
#PCHeatmap(object = DGE, pc.use = 1:12, cells.use = 200, do.balanced = TRUE, 
#          label.columns = FALSE, use.full = FALSE)

# How many components will we keep?
DGE <- JackStraw(object = DGE, num.replicate = 100, do.print = FALSE)
JackStrawPlot(object = DGE, PCs = 1:20)
PCElbowPlot(object = DGE)


# Determine cluster based on components used
DGE <- FindClusters(object = DGE, reduction.type = "pca", dims.use = 1:20, 
                    resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
PrintFindClustersParams(object = DGE)

#### DS9 both species: 16 pca dims, 30 perpl, 200-5000 genes, 0.05 mt
#### DS9 human only: 16 pca dims, 10 perpl, 200-5000 genes, 0.05 mt. 
                      # Note: many human cells are dropped @ 0.05mt
DGE <- RunTSNE(object = DGE, dims.use = 1:16, do.fast = TRUE, perplexity = 30)
t <- TSNEPlot(object = DGE, do.return = TRUE)
t + ggtitle('DS9 HEK-MEF-NO tSNE v2.1')
#save(DGE, file = "NP1.6_DGE_filtered.Robj")


DGE.markers <- FindAllMarkers(object = DGE, only.pos = TRUE, min.pct = 0.25, 
                              thresh.use = 0.25)
DGE.markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_logFC)

top6 <- DGE.markers %>% group_by(cluster) %>% top_n(-6, p_val_adj)
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
DoHeatmap(object = DGE, genes.use = top6$gene, slim.col.label = TRUE, remove.key = TRUE)

# Compare species differences
DGE <- StashIdent(object = DGE, save.name = "CellType")
# Next, switch the identity class of all cells to reflect replicate ID
DGE <- SetAllIdent(object = DGE, id = "cell.species")
t <- TSNEPlot(object = DGE, do.return=TRUE)
t + ggtitle('DS9 HEK-MEF-NO tSNE v2.1 Species Analysis')

#################################################################################
# Save object for future cross-analyses
DS9_HEK_MEF_NO <- DGE
save(DS9_HEK_MEF_NO, file='~/R/Projects/Seurat/DropSeq/Robj/DS9_HEK_MEF_NO.Robj')





