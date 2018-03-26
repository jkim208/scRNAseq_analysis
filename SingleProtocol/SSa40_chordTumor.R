# prepareSeuratObj.R
# Updated: 26-March-2018
# Setup Robjs for cross-technology analyses of Droplet-based Single-Cell RNA Sequencing protocols
# Seurat v2.2.1
# Filter parameters (2018): min.cells = 3, min.genes = 200

library(dplyr);library(ggplot2);library(Seurat);library(Matrix)
setwd("~/Slim/SSa40_ChordTumor/")
set.seed(42)
##################################################################################
# Examine readcounts and find the number of "real" cells
# Do this for each sample A-D
a=read.table("readcounts/SSa40-D_chord_high_cell_readcounts.txt", header=F, stringsAsFactors=F)
x=cumsum(a$V1)
x=x/max(x)
plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]",
     ylab="cumulative fraction of reads", xlim=c(1,1000))
title("SSa40-D_chord_high_cell_readcounts")
abline(v = 200)

##################################################################################
# Prepare Seurat objects for merging
SSa40_A_chordLow <- readRDS(file='~/R/Slim/Robj/Unscaled/SSa40/SSa40-A_chord_low_DGE.Robj')
SSa40_B_chordLow <- readRDS(file='~/R/Slim/Robj/Unscaled/SSa40/SSa40-B_chord_low_DGE.Robj')
SSa40_C_chordHigh <- readRDS(file='~/R/Slim/Robj/Unscaled/SSa40/SSa40-C_chord_high_DGE.Robj')
SSa40_D_chordHigh <- readRDS(file='~/R/Slim/Robj/Unscaled/SSa40/SSa40-D_chord_high_DGE.Robj')

VlnPlot(object = SSa40_A_chordLow, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

SSa40_merge <- MergeSeurat(object1 = SSa40_A_chordLow, object2 = SSa40_B_chordLow, add.cell.id1 = "A_chordLow",
                      add.cell.id2 = "B_chordLow")
SSa40_merge <- MergeSeurat(object1 = SSa40_merge, object2 = SSa40_C_chordHigh, add.cell.id2 = "C_chordHigh")
SSa40_merge <- MergeSeurat(object1 = SSa40_merge, object2 = SSa40_D_chordHigh, add.cell.id2 = "D_chordHigh",
                           project="SSa40")
SSa40_merge <- NormalizeData(object = SSa40_merge, normalization.method = "LogNormalize", scale.factor = 10000)
SSa40_merge <- FindVariableGenes(object = SSa40_merge, mean.function = ExpMean, dispersion.function = LogVMR, 
                            do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
#saveRDS(SSa40_merge, file='~/Slim/Robj/SSa40/SSa40_merge.Robj')

#SSa40_merge <- readRDS(file='~/Slim/Robj/SSa40/SSa40_merge.Robj')
VlnPlot(object = SSa40_merge, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
SSa40_merge <- ScaleData(object = SSa40_merge) #vars.to.regress
SSa40_merge <- RunPCA(object = SSa40_merge, pc.genes = SSa40_merge@var.genes, 
                 do.print = TRUE, pcs.print = 1:2, genes.print = 5, pcs.compute = 20)
PCAPlot(SSa40_merge)
PCElbowPlot(SSa40_merge)
SSa40_merge <- RunTSNE(object = SSa40_merge, dims.use = 1:10, do.fast = TRUE, perplexity = 30)
#SSa40_merge <- SetAllIdent(SSa40_merge, id="orig.ident")

SSa40_merge <- FindClusters(SSa40_merge, resolution = 0.6, dims.use = 1:10)
SSa40_merge <- SetAllIdent(SSa40_merge, id="res.0.6")

p1 <- TSNEPlot(object = SSa40_merge, do.return = TRUE, pt.size = 1.2, do.label=TRUE,label.size=12)
p1 + ggtitle('Chordoma Tumor tSNE Seurat_v2.2') + theme(text = element_text(size=15))
p2 <- FeaturePlot(SSa40_merge, features.plot=c("T"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.cutoff="q90", do.return=TRUE)
plot_grid(p1, p2[[1]])

wilcox_0 <- FindMarkers(object = SSa40_merge, ident.1 = 0, min.pct = 0, test.use = "wilcox", only.pos = FALSE)
roc_4_all <- FindMarkers(object = SSa40_merge, ident.1 = 4, ident.2 = c(0,1,2,3,5,9), min.pct = 0, test.use = "roc", only.pos = TRUE) #ignores immune cells
roc_lowUMI_all <- FindMarkers(object = SSa40_merge, ident.1 = c(0,1,2), ident.2 = c(3,4,5,9), min.pct = 0, test.use = "roc", only.pos = TRUE) #ignores immune cells
roc_all_lowUMI <- FindMarkers(object = SSa40_merge, ident.2 = c(0,1,2), ident.1 = c(3,4,5,9), min.pct = 0, test.use = "roc", only.pos = TRUE) #ignores immune cells

roc_4_5<- FindMarkers(object = SSa40_merge, ident.1 = 4, ident.2 = 5, min.pct = 0, test.use = "roc", only.pos = TRUE)
roc_9_all<- FindMarkers(object = SSa40_merge, ident.1 = 9, ident.2 = c(0,1,2,3,4,5), min.pct = 0, test.use = "roc", only.pos = TRUE) 
roc_3_45<- FindMarkers(object = SSa40_merge, ident.1 = 3, ident.2 = c(4,5), min.pct = 0, test.use = "roc", only.pos = TRUE) 
roc_45_3<- FindMarkers(object = SSa40_merge, ident.2 = 3, ident.1 = c(4,5), min.pct = 0, test.use = "roc", only.pos = TRUE) 

write.table(rownames(roc_all_lowUMI), file="./chordTumor_rocDE_0-1-2_v_3-4-5-9.txt", quote=FALSE, row.names=FALSE,col.names=FALSE)

wilcox_lowT <- FindMarkers(object = SSa40_merge, ident.1 = lowT, ident.2= c(4,5), min.pct = 0, test.use = "wilcox", only.pos = TRUE)
write.table(wilcox_lowT, file="./chordTumor_wilcoxDE_lowT_manual_v_4_5.csv", quote=FALSE, sep="\t")

FeaturePlot(SSa40_merge, features.plot=c("T","nGene", "nUMI", "MKI67", "CCNB2", "KRT19", "KRT8", "ACAN", "COL2A1"), 
            cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q00", max.c="q50")
FeaturePlot(SSa40_merge, features.plot=c("T", "SOD2","NNMT","FTH1","SQSTM1","RND3","ZFP36L1","THBS1","HMOX1"), 
            cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q05", max.c="q75") # surface markers

FeaturePlot(SSa40_merge, features.plot=c("T","nGene","nUMI","MKI67","CKS2","CRYL1","DNA2","HJURP","SUOX"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T","CCNB1","CCNB2","MCM2","MCM3","MCM4","MCM6","MCM7","MCM10"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("ENPP1","CSPG4","COL2A1","ACAN","SLC6A12","KCNK2"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T","SUZ12","TCF7L2","IL6R","TGFB1","MMP14","ITGB1","EGR1","SOX9"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T","LUM","SOX5","SOX6","NFAT5"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T","CA3","TAF13","SYCP2L","KRT19","PDE1A","RRN3P3","CSPG4","WWP2"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T","SAMD5","ENPP1","CD109","TRIL","MRGPRX3","C1QTNF3","EGFLAM","FMOD"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T","HPGDS","MMP19","MIA","EMP2","RAB3B","PIP5K1A","ATF5","HOXA7"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T","GDA","DSG2","GPM6A","AP3B2","FN1","ITGB1","TGFBR1","RELA"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T","KRT10","KRT15","KRT18","KRT19","KRT8"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T","PTPRC","CD302","CD320","CD36","CD3D","CD3EAP","CD3G","CD3E"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T", "YAP1", "FGFR1", "FGFR2", "MAP2K1", "MAPK1"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")

FeaturePlot(SSa40_merge, features.plot=c("T","SNAI1", "SNAI2", "SDC1", "MUC1", "CDH1"), 
            cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T","COL4A1","LAMA2", "LAMA3", "LAMC1"), 
            cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
VlnPlot(object = SSa40_merge, features.plot = c("T", "SNAI1", "SNAI2", "SDC1", "MUC1", "CDH1"), nCol = 3, point.size.use = 0.5)
VlnPlot(object = SSa40_merge, features.plot = c("T", "COL4A1","LAMA2", "LAMA3", "LAMC1"), nCol = 3, point.size.use = 0.5)

FeaturePlot(SSa40_merge, features.plot=c("T", "ACTA2","CDH1", "CDH2", "S100A4", "VIM"), 
            cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
VlnPlot(object = SSa40_merge, features.plot = c("T", "ACTA2","CDH1", "CDH2", "S100A4", "VIM"), nCol = 3, point.size.use = 0.5)

prolif_list <- read.table(file = "~/Slim/Docs/Expt156_chord/Ramaker_proliferation_genes.txt")
prolif_list <- prolif_list[prolif_list$V1 %in% rownames(SSa40_merge@data),]
prolif_count <- apply(SSa40_merge@data[prolif_list,], 1, function(x) round(sum, digits = 4))
percent_cells_exprs_above_0 <- apply(SSa40_merge@data[prolif_list,], 1, function(x) round(mean(x > 0)*100, digits = 2))
# prolif_df <- data.frame(prolif_list, prolif_count, percent_cells_exprs_above_0)
prolif_df <- arrange(prolif_df, desc(prolif_count))
write.table(prolif_df, file="./prolif_df.tsv", quote=FALSE, sep="\t")
FeaturePlot(SSa40_merge, features.plot=c("T",as.character(prolif_df$prolif_list[c(1:8,20)])), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")

FeaturePlot(SSa40_merge, features.plot=c("T", "KRT19", "COL2A1", "CA3", "CD24", "MIA", "TNFRSF11B", "RAB3B", "C1QTNF3", "ACAN", "HAPLN1","KRT15"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T", "CD3D", "CD3E", "CD3G", "CD28", "CD4"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T", "CD4", "FOXP3", "IL2RA", "CTLA4"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T", "CD19", "MS4A1", "CD79A", "CD79B","BLNK"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T", "MS4A1", "CD79A", "CD79B","BLNK"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T", "CD79B","BLNK"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T", "CD14", "CD68", "CD163", "CSF1R","FCGR3A"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T", "IL3RA", "NRP1"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T", "FCGR3A", "NCAM1","KLRB1","KLRC1","KLRD1"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T", "VWF", "CDH5","SELE"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("T", "FAP", "THY1","COL1A1","COL3A1"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")
FeaturePlot(SSa40_merge, features.plot=c("DDX3X", "RPS4X", "USP9Y", "ZFY", "DDX3Y", "RPS4Y1"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")

c("T", "THY1", "ENG","POU5F1","NES","SOX2","PROM1") %in% rownames(SSa40_merge@data)
FeaturePlot(SSa40_merge, features.plot=c("T", "THY1", "ENG","NES"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q10", max.c="q90")

VlnPlot(object = SSa40_merge, features.plot = c("T", "THY1", "ENG","NES"), nCol = 2, point.size.use = 0.5)
VlnPlot(object = SSa40_merge, features.plot = c("T", "KRT19", "COL2A1","CA3", "CD24", "THBS1"), nCol = 3, point.size.use = 0.5)
VlnPlot(object = SSa40_merge, features.plot = c( "nUMI"), point.size.use = 0.5,y.max=5000)

SSa40_merge <- SetAllIdent(SSa40_merge, id="names")
p1 <- TSNEPlot(object = SSa40_merge, do.return = TRUE, pt.size = 1.4, do.label=TRUE,label.size=12)
p1 + ggtitle('SSa40 tSNE Seurat_v2.2') + theme(text = element_text(size=18))

##################################################################################
##################################################################################
# Subset data to remove low_nUMI cells believed to be in apoptosis
# SSa40_merge <- SetAllIdent(SSa40_merge, id="res.0.6")
SSa40_merge <- SubsetData(SSa40_merge, ident.remove = c(0,1,2))
SSa40_merge <- NormalizeData(object = SSa40_merge, normalization.method = "LogNormalize", scale.factor = 10000)
SSa40_merge <- FindVariableGenes(object = SSa40_merge, mean.function = ExpMean, dispersion.function = LogVMR, 
                                 do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
VlnPlot(object = SSa40_merge, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
SSa40_merge@meta.data$res.unfiltered <- SSa40_merge@meta.data$res.0.6
SSa40_merge <- ScaleData(object = SSa40_merge) #vars.to.regress
SSa40_merge <- RunPCA(object = SSa40_merge, pc.genes = SSa40_merge@var.genes, 
                      do.print = TRUE, pcs.print = 1:2, genes.print = 5, pcs.compute = 20)
PCAPlot(SSa40_merge)
PCElbowPlot(SSa40_merge)
SSa40_merge <- RunTSNE(object = SSa40_merge, dims.use = 1:10, do.fast = TRUE, perplexity = 30)
#SSa40_merge <- SetAllIdent(SSa40_merge, id="orig.ident")
SSa40_merge <- FindClusters(SSa40_merge, resolution = 0.6, dims.use = 1:10,force.recalc = TRUE)
SSa40_merge <- SetAllIdent(SSa40_merge, id="names")
p1 <- TSNEPlot(object = SSa40_merge, do.return = TRUE, pt.size = 3, do.label=TRUE,label.size=16)
p2 <- FeaturePlot(SSa40_merge, features.plot=c("T"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q05", max.c="q95", pt.size = 3, do.return = TRUE)[[1]]
plot_grid(p1+ggtitle('Chordoma Tumor tSNE Seurat_v2.2') + theme(text = element_text(size=15)), p2)
# B is low T. CD is high T. E is weird.
p2 <- FeaturePlot(SSa40_merge, features.plot=c("nUMI"), cols.use=c("dark grey", "red"), dark.theme=TRUE, min.cutoff="q05", max.c="q95", pt.size = 3, do.return = TRUE)[[1]]
plot_grid(p1+ggtitle('Chordoma Tumor tSNE Seurat_v2.2') + theme(text = element_text(size=15)), p2)



png("align_cellLines_chord_genes_vln.png",width=1920,height=1080)
VlnPlot(object = SSa40_merge, features.plot = c("T", "KRT8", "KRT19", "COL11A1", "COL2A1", "CA3", "CD24", "MIA","C1QTNF3", "ACAN", "HAPLN1","RAB3B"), ident.include = c("B","C","D"), point.size.use = 0.5)
FeaturePlot(SSa40_merge, features.plot=c("T", "KRT8", "KRT19", "COL11A1", "COL2A1", "CA3", "CD24", "MIA","C1QTNF3", "ACAN", "HAPLN1","RAB3B"), cols.use=c("dark grey", "red"), cells.use = BCD_cells, dark.theme=TRUE, min.cutoff="q00", max.c="q75", pt.size = 2.0)

markers <- FindMarkers(object = SSa40_merge, ident.2 = "B", ident.1= c("C","D"), min.pct = 0, test.use = "roc", only.pos = TRUE)
markers$diff <- markers$pct.1 - markers$pct.2
write.table(rownames(markers), file="./chordTumor_rocDE_B_v_CD.tsv", quote=FALSE, sep="\t",row.names = FALSE,col.names = FALSE)
head(markers[markers$diff>0.3,], n=12)
VlnPlot(object = SSa40_merge, features.plot = c("T", "KRT19", "AKR1B10","CD24", "ID4", "RP11-149I23.3","RAPGEF5", "TRIL", "LYST","C1QTNF3","SMOC2","SLC4A11"), ident.include = c("B","C","D"), point.size.use = 0.5)
BCD_cells <- SSa40_merge@cell.names[grep("[BCD]",SSa40_merge@ident)]
x <- TSNEPlot(object = SSa40_merge, do.identify=T)
BCD_cells <- BCD_cells[!BCD_cells %in% x]
FeaturePlot(SSa40_merge, features.plot=c("T", "ZFP36L1", "CALD1","UACA", "NFKBIZ", "CADM1","GPX3", "NNMT", "DAB2"), cols.use=c("dark grey", "red"), cells.use = BCD_cells, dark.theme=TRUE, min.cutoff="q00", max.c="q75", pt.size = 2.0)

VlnPlot(object = SSa40_merge, features.plot = c("T","SNAI2","CDH2"), ident.include = c("B","C","D"), point.size.use = 1)
VlnPlot(object = SSa40_merge, features.plot = c("T","UCHL3", "ALG11","PPP2CB"), ident.include = c("B","C","D"), point.size.use = 1)
FeaturePlot(SSa40_merge, features.plot=c("T","UCHL3", "ALG11","PPP2CB"), cols.use=c("dark grey", "red"), cells.use = BCD_cells, dark.theme=TRUE, min.cutoff="q00", max.c="q75", pt.size = 2.0)


