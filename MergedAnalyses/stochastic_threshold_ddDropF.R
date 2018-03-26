# 17-Dec-2017
# Continued from dd_Drop_analysis_v2.1.R
library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4/Seurat2.1')
library(Matrix)
require(gridExtra)
setwd("~/R/Projects/Seurat")
set.seed(42)

mega_seurat <- readRDS(file='~/R/Projects/Seurat/Robj/dd_f_drop.Robj')

#################################
# NOTE TO SELF: DID I FORGET TO USE LOG FOR THESE HISTOGRAMS?
#################################

############
# Create indices to track cell type and protocol
x <- mega_seurat@meta.data$adj.ident
names(x) <- rownames(mega_seurat@meta.data)
hek_index <- c(grep('293', x), grep('HEK', x))
uch1_index <- c(grep('uch1', x), grep('UCH1', x))
uch2_index <- grep('UCH2', x)
ddSeq_index <- grep('ddSeq', x)
Drop_index <- grep('DropSeq', x)
Fluidigm_index <- grep('Fluidigm', x)
############
exp_hist <- function(rna_df){
  x <- rna_df
  x[x == 0] <- NA # Cells with non-expressed gene are dropped for that gene's analysis
  x <- as.matrix(x)
  y <- rowMeans(x, na.rm=TRUE)
  z <-as.data.frame(y)
  colnames(z) <- c('log_gene_expression')
  return(z) # return transformed data for histogram plotting
}
######################################################
# HEK
z <- exp_hist(mega_seurat@raw.data[,hek_index])
p1 <- ggplot(z, aes(log_gene_expression)) +
  geom_histogram(binwidth = 0.05) + ggtitle('293/HEK Histogram bin=0.05')+ 
  geom_vline(xintercept = z[which(rownames(z)=='T'),1] ) +
  geom_vline(xintercept = mean(z$log_gene_expression, na.rm=TRUE)) +
  annotate("text", x = 1.25, y = 1500, label = "Average")
mean(z$log_gene_expression, na.rm=TRUE)

######################################################
# UCH1
z <- exp_hist(mega_seurat@data[,uch1_index])
p2 <- ggplot(z, aes(log_gene_expression)) + 
  geom_histogram(binwidth = 0.05) + ggtitle('UCH1 Histogram bin=0.05') + 
  geom_vline(xintercept = z[which(rownames(z)=='T'),1], color="green" ) +
  geom_vline(xintercept = mean(z$log_gene_expression, na.rm=TRUE),color="red") +
  annotate("text", x = 1.40, y = 850, label = "Average", size=3) +
  annotate("text", x = 2.25, y = 850, label = "T", size=3)
######################################################
# UCH2
z <- exp_hist(mega_seurat@data[,uch2_index])
p3 <- ggplot(z, aes(log_gene_expression)) +
  geom_histogram(binwidth = 0.05) + ggtitle('UCH2 Histogram bin=0.05') +
  geom_vline(xintercept = z[which(rownames(z)=='T'),1], color="green" ) +
  geom_vline(xintercept = mean(z$log_gene_expression, na.rm=TRUE),color="red") +
  annotate("text", x = 1.38, y = 1500, label = "Average", size=3) +
  annotate("text", x = 1.84, y = 1500, label = "T", size=3)

grid.arrange(p1,p2,p3, ncol=2)
######################################################
######################################################
# Compare across technologies

z <- exp_hist(mega_seurat@data[,intersect(hek_index, Drop_index)])
p1 <- ggplot(z, aes(log_gene_expression)) + 
  geom_histogram(binwidth = 0.05) + ggtitle('HEK293 DropSeq Histogram bin=0.05') + 
  geom_vline(xintercept = z[which(rownames(z)=='T'),1], color="green" ) +
  geom_vline(xintercept = mean(z$log_gene_expression, na.rm=TRUE),color="red") +
  annotate("text", x = 1.75, y = 1000, label = "Average", size=3) 

z <- exp_hist(mega_seurat@data[,intersect(hek_index, ddSeq_index)])
p2 <- ggplot(z, aes(log_gene_expression)) + 
  geom_histogram(binwidth = 0.05) + ggtitle('HEK293 ddSeq Histogram bin=0.05') + 
  geom_vline(xintercept = z[which(rownames(z)=='T'),1], color="green" ) +
  geom_vline(xintercept = mean(z$log_gene_expression, na.rm=TRUE),color="red") +
  annotate("text", x = 1.2, y = 2800, label = "Average", size=3)

z <- exp_hist(mega_seurat@data[,intersect(hek_index, Fluidigm_index)])
p3 <- ggplot(z, aes(log_gene_expression)) + 
  geom_histogram(binwidth = 0.05) + ggtitle('HEK Fluidigm Histogram bin=0.05') + 
  geom_vline(xintercept = mean(z$log_gene_expression, na.rm=TRUE),color="red") +
  annotate("text", x = 1.2, y = 1500, label = "Average", size=3)

z <- exp_hist(mega_seurat@data[,intersect(uch1_index, Drop_index)])
p4 <- ggplot(z, aes(log_gene_expression)) + 
  geom_histogram(binwidth = 0.05) + ggtitle('UCH1 DropSeq Histogram bin=0.05') + 
  geom_vline(xintercept = z[which(rownames(z)=='T'),1], color="green" ) +
  geom_vline(xintercept = mean(z$log_gene_expression, na.rm=TRUE),color="red") +
  annotate("text", x = 1.01, y = 750, label = "Average", size=3) +
  annotate("text", x = 2.15, y = 750, label = "T", size=3)

z <- exp_hist(mega_seurat@data[,intersect(uch1_index, ddSeq_index)])
p5 <- ggplot(z, aes(log_gene_expression)) + 
  geom_histogram(binwidth = 0.05) + ggtitle('UCH1 ddSeq Histogram bin=0.05') + 
  geom_vline(xintercept = z[which(rownames(z)=='T'),1], color="green" ) +
  geom_vline(xintercept = mean(z$log_gene_expression, na.rm=TRUE),color="red") +
  annotate("text", x = 2.48, y = 820, label = "Average", size=3) +
  annotate("text", x = 3.2, y = 820, label = "T", size=3)

z <- exp_hist(mega_seurat@data[,intersect(uch1_index, Fluidigm_index)])
p6 <- ggplot(z, aes(log_gene_expression)) + 
  geom_histogram(binwidth = 0.05) + ggtitle('UCH1 Fluidigm Histogram bin=0.05') + 
  geom_vline(xintercept = z[which(rownames(z)=='T'),1], color="green" ) +
  geom_vline(xintercept = mean(z$log_gene_expression, na.rm=TRUE),color="red") +
  annotate("text", x = 0.30, y = 1000, label = "Average", size=3) +
  annotate("text", x = 1.88, y = 1000, label = "T", size=3)

z <- exp_hist(mega_seurat@data[,intersect(uch2_index, Drop_index)])
p7 <- ggplot(z, aes(log_gene_expression)) + 
  geom_histogram(binwidth = 0.05) + ggtitle('UCH2 DropSeq Histogram bin=0.05') + 
  geom_vline(xintercept = z[which(rownames(z)=='T'),1], color="green" ) +
  geom_vline(xintercept = mean(z$log_gene_expression, na.rm=TRUE),color="red") +
  annotate("text", x = 1.6, y = 450, label = "Average", size=3) +
  annotate("text", x = 2.44, y = 450, label = "T", size=3)

z <- exp_hist(mega_seurat@data[,intersect(uch2_index, ddSeq_index)])
p8 <- ggplot(z, aes(log_gene_expression)) + 
  geom_histogram(binwidth = 0.05) + ggtitle('UCH2 ddSeq Histogram bin=0.05') + 
  geom_vline(xintercept = z[which(rownames(z)=='T'),1], color="green" ) +
  geom_vline(xintercept = mean(z$log_gene_expression, na.rm=TRUE),color="red") +
  annotate("text", x = 0.60, y = 2100, label = "Average", size=3) +
  annotate("text", x = 1.85, y = 2100, label = "T", size=3)

grid.arrange(p1,p2,p3, ncol=3)
grid.arrange(p4,p5,p6, ncol=3)
grid.arrange(p7,p8, ncol=3)
####################################################
####################################################
####################################################
# Expression Profiling (Stefan's instructions)

mega_seurat@raw.data[,UCH1]
100*sum(mega_seurat@raw.data[5,UCH1]>0)/ncol(mega_seurat@raw.data[,UCH1])

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


hek_EP <- expressionByCellType(hek_index)
uch1_EP <- expressionByCellType(uch1_index)
uch2_EP <- expressionByCellType(uch2_index)
df <- rbind(hek_EP, uch1_EP, uch2_EP)
df$cellType <- as.factor(c(rep('HEK', nrow(hek_EP)), rep('UCH1', nrow(hek_EP)), rep('UCH2', nrow(hek_EP))))
ddDrop.df.ep <- df
saveRDS(ddDrop.df.ep, file = '~/R/Projects/Seurat/Robj/ddDrop_expression_profile.Robj')
rm(ddDrop.df.ep)

hek_ddS_EP <- expressionByCellType(intersect(hek_index, ddSeq_index))
hek_drop_EP <- expressionByCellType(intersect(hek_index, Drop_index))
hek_fluidigm_EP <- expressionByCellType(intersect(hek_index, Fluidigm_index))
uch1_ddS_EP <- expressionByCellType(intersect(uch1_index, ddSeq_index))
uch1_drop_EP <- expressionByCellType(intersect(uch1_index, Drop_index))
uch1_fluidigm_EP <- expressionByCellType(intersect(uch1_index, Fluidigm_index))
uch2_ddS_EP <- expressionByCellType(intersect(uch2_index, ddSeq_index))
uch2_drop_EP <- expressionByCellType(intersect(uch2_index, Drop_index))

df <- rbind(hek_ddS_EP, hek_drop_EP, hek_fluidigm_EP, 
            uch1_ddS_EP, uch1_drop_EP, uch1_fluidigm_EP, uch2_ddS_EP, uch2_drop_EP)
df$cellType <- as.factor(c(rep('HEK_ddSeq', nrow(hek_ddS_EP)), rep('HEK_Drop', nrow(hek_drop_EP)),
                           rep('HEK_Fluidigm', nrow(hek_fluidigm_EP)), rep('uch1_ddSeq', nrow(uch1_ddS_EP)), 
                           rep('UCH1_drop', nrow(uch1_drop_EP)), rep('UCH1_fluidigm'),nrow(uch1_fluidigm_EP)),
                           rep('UCH2_ddSeq', nrow(uch2_ddS_EP)), rep('UCH2_drop', nrow(uch2_drop_EP)) ) 
ddDrop.df.ep.tech <- df
saveRDS(ddDrop.df.ep.tech, file = '~/R/Projects/Seurat/Robj/ddDrop_expression_profile_ByTech.Robj')
rm(ddDrop.df.ep.tech)
ddDrop.df.ep.tech <- readRDS("~/R/Projects/Seurat/Robj/ddDrop_expression_profile_ByTech.Robj")

ggplot(df) + 
  geom_point(mapping=aes(cellPerc, avg_UMI, size=sd_UMI, color=sd_UMI, alpha=sd_UMI)) +
  scale_colour_gradient(low = "red", high = "black") +
  scale_size(range = c(3,0.5)) +
  scale_alpha(range = c(0.8, 0.2)) +
  ggtitle("Expression Profile") +
  facet_grid(. ~ cellType) #+
#geom_hline(data = vline.df, aes(yintercept = vl))

### Find mean expression excluding non-expressing cells
#hek_EP[hek_EP==0] <- NA
#mean(hek_EP$avg_UMI, na.rm=TRUE)
#uch1_EP[uch1_EP==0] <- NA
#mean(uch1_EP$avg_UMI, na.rm=TRUE)
#uch2_EP[uch2_EP==0] <- NA
#mean(uch2_EP$avg_UMI, na.rm=TRUE)
#meanUMI <- c(mean(hek_EP$avg_UMI, na.rm=TRUE), mean(uch1_EP$avg_UMI, na.rm=TRUE), 
#             mean(uch2_EP$avg_UMI, na.rm=TRUE))
#vline.df <- data.frame(z=levels(df$cellType), vl= meanUMI)




