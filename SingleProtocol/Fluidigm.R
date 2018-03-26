# Fluidigm.R
# Updated: 26-March-2018
# Examine Fluidigm data via Seurat
# Seurat v2.2.1

library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4')
library(Matrix)
library(dplyr)

setwd("~/Slim/Fluidigm")
# file below produced from fluidigm_prep.R run on server
fluidigm_raw <- readRDS(file="~/Slim/Robj/Fluidigm/Fluidigm_total_expression_count.Robj")
cell_features <- fluidigm_raw[1:5,]; rownames(cell_features) <- rownames(fluidigm_raw)[1:5]
backup <- fluidigm_raw
fluidigm_raw <- as.data.frame(fluidigm_raw[-(1:5),]); rownames(fluidigm_raw) <- rownames(backup)[6:nrow(backup)]
rm(backup)
#########################################################
# Determine what cells are from mouse and what cells are from human
mus_sum <- colSums(fluidigm_raw[59638:105175,]) # mouse genes selected
homo_sum <- colSums(fluidigm_raw[1:59637,]) # human genes selected
df <- data.frame(homo_sum, mus_sum, colnames(fluidigm_raw))
df$ratioH <- df$homo_sum / (df$homo_sum + df$mus_sum)
df <- mutate(df, Species = ifelse(ratioH > 0.8 , "Human", 
                                  ifelse(ratioH > 0.2 & ratioH < 0.8, "Mixed", "Mouse")))
df <- df %>% mutate (x = ifelse(Species=="Human", homo_sum, ifelse (Species=="Mixed", homo_sum, 0)))
df <- df %>% mutate (y = ifelse (Species=="Mixed", mus_sum, ifelse (Species=="Mouse", mus_sum, 0)))
df %>% group_by(Species) %>% summarise(number=n())
#  Species number
#    <chr>  <int>
#1   Human    529
#2   Mixed    119
#3   Mouse    152
ggplot(df.filter, aes (x,y , colour=(Species) )) + 
  geom_point() + xlab('Human Transcripts') + ylab('Mouse Transcripts') + 
  ggtitle('Fluidigm Species Analysis',subtitle =  '293HEK-3T3 and UCH1')
# sum(which(df$Species == "Human") < 401)
# Cols 11-20 are all human cells
# Of the remaining 400 cells in the first 10 cols, we can confidently say 129 of them are human

# return to df to check species representation
df[1:400,] %>% group_by(Species) %>% summarise(number=n()) # first 100 cells are 293-3T3
# A tibble: 3 x 2
#Species number
#<chr>  <int>
#1   Human    129
#2   Mixed    119
#3   Mouse    152
HEK_cols <- as.character(df$colnames.fluidigm_raw.[which(df$Species[1:400] == "Human")])
UCH1_cols <- as.character(df$colnames.fluidigm_raw.[401:nrow(df)])

#########################################################
# Repeat analysis with human-only genes and cells
# start from fluidigm DEM
homo_l <- length(grep('homo_', rownames(fluidigm_raw))) # index=59367
rownames(fluidigm_raw) <- gsub('homo_', '', rownames(fluidigm_raw))
# Human cells determined earlier (HEK and UCH1). Grab only human genes(rows).
Fluidigm <- CreateSeuratObject(raw.data = fluidigm_raw[1:homo_l,c(HEK_cols, UCH1_cols)], is.expr = 1,
                               min.cells = 3, min.genes = 200, project = 'fluidigm')
mito.genes <- grep(pattern = "[Mm][Tt]-", x = rownames(x = Fluidigm@data), value = TRUE)
percent.mito <- Matrix::colSums(Fluidigm@raw.data[mito.genes, ])/Matrix::colSums(Fluidigm@raw.data)
Fluidigm <- AddMetaData(object = Fluidigm, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = Fluidigm, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
# Unusually high MT-content in Fluidigm. Especially the 293HEK.
Fluidigm <- FilterCells(object = Fluidigm, subset.names = c("nGene", "percent.mito"), 
                         low.thresholds = c(200, -Inf), high.thresholds = c(Inf, Inf))
Fluidigm <- NormalizeData(object = Fluidigm,normalization.method = "LogNormalize", scale.factor = 10000)
#saveRDS(Fluidigm, file="~/Slim/Robj/Fluidigm/Fluidigm_human.Robj")

