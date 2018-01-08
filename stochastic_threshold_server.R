
library(Seurat, lib.loc='~/R/x86_64-pc-linux-gnu-library/3.4')
library(dplyr, lib.loc='~/R/x86_64-pc-linux-gnu-library/3.4')
library(Matrix, lib.loc='~/R/x86_64-pc-linux-gnu-library/3.4')
setwd("~/R/ddDrop")

mega_seurat <- readRDS(file = '~/R/ddDrop/ddDrop_mega_seurat_v2.1.Robj')
x <- mega_seurat@meta.data$orig.ident
names(x) <- rownames(mega_seurat@meta.data)
hek_index <- c(grep('293', x), grep('HEK', x))
uch1_index <- grep('UCH1', x)
uch2_index <- grep('UCH2', x)
ddSeq_index <- grep('ddSeq', x)
Drop_index <- grep('DropSeq', x)


expressionByCellType <- function(cellType) {
  # Percent of cells where this gene is detected at a level above 0
  cellPerc <- apply(mega_seurat@raw.data[,cellType], 
                    1, function(x) 100*sum(x > 0) /
                      ncol(mega_seurat@raw.data[,cellType]))
  na.mega_seurat.raw <- mega_seurat@raw.data
  na.mega_seurat.raw[na.mega_seurat.raw == 0] <- NA # convert 0 values to NA
  #Average number of transcripts detected, when it is detected (e.g. drop all 0 values from this calculation)
  avg_UMI <- apply(mega_seurat@raw.data[,cellType], 
                   1, function(x) sum(x) / 
                     sum(x > 0))
  avg_UMI[is.nan(avg_UMI)] <- 0 # remove NaN values
  #Standard deviation of the transcript numbers in the values where it is detected
  sd_UMI <- apply(na.mega_seurat.raw[,cellType], 
                  1, function(x) sd(x, na.rm=TRUE))
  sd_UMI[is.na(sd_UMI)] <- 0 # remove NA values
  return(data.frame(cellPerc, avg_UMI, sd_UMI))
}


hek_ddS_EP <- expressionByCellType(intersect(hek_index, ddSeq_index))
hek_drop_EP <- expressionByCellType(intersect(hek_index, Drop_index))
uch1_ddS_EP <- expressionByCellType(intersect(uch1_index, ddSeq_index))
uch1_drop_EP <- expressionByCellType(intersect(uch1_index, Drop_index))
uch2_ddS_EP <- expressionByCellType(intersect(uch2_index, ddSeq_index))
uch2_drop_EP <- expressionByCellType(intersect(uch2_index, Drop_index))

df <- rbind(hek_ddS_EP, hek_drop_EP, uch1_ddS_EP, uch1_drop_EP, uch2_ddS_EP, uch2_drop_EP)
df$cellType <- as.factor(c(rep('HEK_ddSeq', nrow(hek_ddS_EP)), rep('HEK_Drop', nrow(hek_drop_EP)), 
                           rep('uch1_ddSeq_EP', nrow(uch1_ddS_EP)), rep('uch1_drop_EP', nrow(uch1_drop_EP)),
                           rep('uch2_ddSeq_EP', nrow(uch2_ddS_EP)), rep('uch2_drop_EP', nrow(uch2_drop_EP)) ) )
ddDrop.df.ep.tech <- df
saveRDS(ddDrop.df.ep.tech, file = '~/R/ddDrop/ddDrop_expression_profile_ByTech.Robj')
rm(ddDrop.df.ep.tech)





