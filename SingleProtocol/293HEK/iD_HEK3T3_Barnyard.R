library(dplyr)
library(ggplot2)
library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4/Seurat2.1')
library(Matrix)

setwd("~/R/Projects/inDrop")
set.seed(42)

# Mixed Species Analysis
##############################################################
# Find the Top Barcodes Up to a Specific Threshold
a=read.table("input/lib3_cell_readcounts.txt", 
             header=F, stringsAsFactors=F) ##
x=cumsum(a$V1)
x=x/max(x)
plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]",
     ylab="cumulative fraction of reads", xlim=c(1,5000))
abline(v=2500) ##
topBarcodes <- a[1:3000,2] # Up to knee threshold ##
write.table(x=topBarcodes, file = "input/inDrop_lib23_top.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE) ##



human <- read.table("input/lib3_Human_DGE_Summary.txt", header=T) ##
human <-human %>% mutate(sample="human")
mouse <- read.table("input/lib3_Mouse_DGE_Summary.txt", header=T) ##
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
  ggtitle('inDrop',subtitle =  '293/3T3 Mixed-Species plot with top 2500 cells') + 
  scale_colour_manual(values=cp_Palette) ##
p


##############################################################
# human specific analysis
humanCells <- as.character(combined[grep('Human', combined$Species),]$CELL_BARCODE)

DGE.data <- read.table(file = "input/lib3_DGE.txt", header = TRUE, row.names = 1)
DGE.human <- DGE.data[grepl('^homo_', rownames(DGE.data)),] # keep only human genes
rownames(DGE.human) <- sub('^homo_', '', rownames(DGE.human))  # remove ^homo tags
DGE <- CreateSeuratObject(raw.data = DGE.human[,intersect(humanCells, colnames(DGE.human))], 
                          min.cells = 3, min.genes = 200, project = "iD_293HEK_lib23")
mito.genes <- grep(pattern = "[Mm][Tt]-", x = rownames(x = DGE@data), value = TRUE)
percent.mito <- Matrix::colSums(DGE@raw.data[mito.genes, ])/Matrix::colSums(DGE@raw.data)
DGE <- AddMetaData(object = DGE, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = DGE, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = DGE, gene1 = "nUMI", gene2 = "nGene")
GenePlot(object = DGE, gene1 = "nUMI", gene2 = "percent.mito")
# We filter out cells that have unique gene counts over 3,500 or less than 200
DGE <- FilterCells(object = DGE, subset.names = c("nGene", "percent.mito"), 
                   low.thresholds = c(200, -Inf), high.thresholds = c(3500, 0.15))
DGE <- NormalizeData(object = DGE,normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(DGE, file="input/iD_293HEK_lib23_human.Robj") #######################
##############################################################
##############################################################
# iD_293-3T3 mixed species analysis
mouseCells <- as.character(combined[grep('Mouse', combined$Species),]$CELL_BARCODE)
DGE <- CreateSeuratObject(raw.data = DGE.data[,intersect(c(humanCells, mouseCells), colnames(DGE.data))], 
                          min.cells = 3, min.genes = 200, project = "iD_293_3T3_lib23")
mito.genes <- grep(pattern = "[Mm][Tt]-", x = rownames(x = DGE@data), value = TRUE)
percent.mito <- Matrix::colSums(DGE@raw.data[mito.genes, ])/Matrix::colSums(DGE@raw.data)
DGE <- AddMetaData(object = DGE, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = DGE, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = DGE, gene1 = "nUMI", gene2 = "nGene")
GenePlot(object = DGE, gene1 = "nUMI", gene2 = "percent.mito")
# We filter out cells that have unique gene counts over 3,500 or less than 200
DGE <- FilterCells(object = DGE, subset.names = c("nGene", "percent.mito"), 
                   low.thresholds = c(200, -Inf), high.thresholds = c(3000, 0.15))
DGE <- NormalizeData(object = DGE,normalization.method = "LogNormalize", scale.factor = 10000)
DGE <- FindVariableGenes(object = DGE, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = FALSE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
DGE <- ScaleData(object = DGE) # vars.to.regress = c("nUMI", "percent.mito"))
DGE <- RunPCA(object = DGE, pc.genes = DGE@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 30)
PCAPlot(object = DGE, dim.1 = 1, dim.2 = 2)
PCElbowPlot(object = DGE, num.pc = 30)
#First 540 cells are human
DGE <- SetIdent(DGE, cells.use = DGE@cell.names[1:540], ident.use = "293")
DGE <- SetIdent(DGE, cells.use = DGE@cell.names[541:length(DGE@cell.names)], ident.use = "3T3")

DGE <- RunTSNE(object = DGE, dims.use = 1:10, do.fast = TRUE, perplexity = 30)
TSNEPlot(object = DGE, do.return=TRUE) + 
  ggtitle('inDrop lib23 293-3T3 Mixed_Species No_Regress tSNE v2.1')















