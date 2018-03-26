# SSa39_pita_edgeR.R
# Updated: 26-March-2018
# Bulk RNA-sequencing experiment where control 293HEK, chordoma UCH2 and CH22 are 
# subject to mock and Pitavastatin treatments. We use edgeR to analyze. 
# Since mock replicates are missing, we will use some from the Zika experiment

library(reshape2)
library(ggplot2)
library("pheatmap")
library("RColorBrewer")
source("https://bioconductor.org/biocLite.R")
library(edgeR)
library(knitr)
library("org.Hs.eg.db")
library(dplyr)
setwd("~/Slim")

# Count Analysis
samples<-read.table("./SSa39_Pita/Data/Counts_File_Mapping.txt", header=T)
head(samples)

y<-readDGE(samples, path="./SSa39_Pita/Data/Counts", group=as.factor(samples$label))
tail(y$counts)

# Remove meta tags:
to_remove<-!(startsWith(rownames(y$counts), "__"))
tail(to_remove)
y <- y[to_remove,,keep.lib.sizes=FALSE]
tail(y$counts)

# Drop the technical replicates into one group:
y<-sumTechReps(y, ID=y$samples$group)
# How many million reads passed through all the mapping and counting stages in each library?
#sum(y$samples[startsWith(as.character(y$samples$group), "Seq1"),"lib.size"])/1000000
#sum(y$samples[startsWith(as.character(y$samples$group), "Seq2"),"lib.size"])/1000000

#Show what we would filter out by cutting at a CPM of 1:
plot(cpm(y)[,1], ylim=c(0,100))
abline(h=1)

dim(y)
keep <- rowSums(cpm(y)>1) >= 1
y<-y[keep, , keep.lib.sizes=FALSE] # remove low-expressing genes (<1 cpm)
dim(y)

# We chose Seq1 mocks since they had all 3 cell types and the counts were consistent with the Pita experiment
#y <- y[,c(1:9)] # No replicates from Zika Seqs
#y <- y[,c(1:10,12,20,24:25,32:33,34,41,42)] # Includes mock replicates from Seq1,2,3
y <- y[,c(1:10,12,20)] # Seq1 mocks
#y <- y[,c(1:9,24:25,32:33)] # Seq2 mocks
#y <- y[,c(1:10,12,20,24:25,32:33)] # Seq1,2 mocks
#y <- y[,c(1:10,12,20,34,41,42)] # Seq1,3 mocks
#y <- y[,c(1:9,34,41,42)] # Seq3 mocks

# Calculate the norm factors 
y<-calcNormFactors(y)

# MDS Plot
tmp <- plotMDS(y)
coord.df <- as.data.frame(tmp$cmdscale.out)
rm(tmp)
colnames(coord.df) <- c("x","y")

coord.df$sample <- rownames(coord.df)
coord.df$sample <- gsub("_", "-", coord.df$sample)
coord.df$sample <- gsub("Control", "Mock", coord.df$sample)
coord.df$sample <- gsub("293ETN", "293", coord.df$sample)
sample.info <- sapply(coord.df$sample, function(x) strsplit(as.character(x), "-", fixed=T))
coord.df$cells <- sapply(sample.info, function(x) unlist(x)[2])
coord.df$treatment <- sapply(sample.info, function(x) unlist(x)[3])
coord.df$source <- sapply(sample.info, function(x) unlist(x)[1])

g <- ggplot(data=coord.df, aes(x=x, y=y)) + geom_text(aes(label=sample, col=cells))
# g <- ggplot(data=coord.df, aes(x=x, y=y)) + geom_text(aes(label=sample, col=treatment))
# g <- ggplot(data=coord.df, aes(x=x, y=y)) + geom_text(aes(label=sample, col=source)) 
g + ggtitle("MDS Plot: CH22 samples from Expt66B and Zika-Seq1-2")

kable(y$samples)
SSa39.df <- data.frame(coord.df[,3:5])
SSa39.df$cell_group <- "Chord"
SSa39.df[ grep("293", SSa39.df$sample),"cell_group"] <- "Non-Chord" 
SSa39.df$cell_group <- as.factor(SSa39.df$cell_group)
SSa39.df$cell_group <- relevel(SSa39.df$cell_group, ref="Non-Chord")
SSa39.df$treatment <- as.factor(SSa39.df$treatment)
SSa39.df$treatment <- relevel(SSa39.df$treatment, ref="Mock")
SSa39.df$cells <- as.factor(SSa39.df$cells)
SSa39.df$cells <- relevel(SSa39.df$cells, ref="293")

# Create the design for the experiment
design <- model.matrix(~cells + cells:treatment, data=SSa39.df)
design <- design[,colSums(design) != 0]

y_glm <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y_glm <- estimateGLMTrendedDisp(y_glm, design)
y_glm <- estimateGLMTagwiseDisp(y_glm, design)
fit <- glmFit(y_glm, design) 

# Define contrasts to study
my_contrasts <- makeContrasts(
  # 1. Difference between UCH2 and 293 when being compared from Mock to Pita-treated states
  UCH2.vs.293.Pita = cellsUCH2.treatmentPita - cells293.treatmentPita,
  # 2. Difference between CH22 and 293 when being compared from Mock to Pita-treated states
  CH22.vs.293.Pita = cellsCH22.treatmentPita - cells293.treatmentPita,
  # 3. Difference between CH22 and UCH2 when being compared from Mock to Pita-treated states
  UCH2.vs.CH22.Pita = cellsUCH2.treatmentPita - cellsCH22.treatmentPita,
  
  # 4. Difference between chordoma and 293 when being compared from Mock to Pita-treated states
  chord.Pita = (cellsCH22.treatmentPita + cellsUCH2.treatmentPita)/2 - (cells293.treatmentPita),
  
  # 5:7. How each cell type has changed from Mock to Pita-treated states
  cellsUCH2.treatmentPita,
  cellsCH22.treatmentPita,
  cells293.treatmentPita,
  
  levels=make.names(colnames(design))
)

plotBCV(y_glm)

# Label contrasts to test
toTest <- c("UCH2.vs.293.Pita", "CH22.vs.293.Pita", "UCH2.vs.CH22.Pita", "chord.Pita", 
            "cellsUCH2.treatmentPita", "cellsCH22.treatmentPita", "cells293.treatmentPita")

# Via edgeR's glm, find the DE genes
lrt <- glmLRT(fit, contrast=my_contrasts[,"chord.Pita"])
UCH2.Pita <- topTags(lrt, n=100, p.value=0.05)

################################################################################
# limma: alternative to base edgeR
################################################################################
# Use if you need to avoid Ensembl to Entrez geneID conversion
limma_des <- design
colnames(limma_des) <- make.names(colnames(limma_des))

v <- voom(y, limma_des, plot=TRUE)
head(sort(v$E[,1], decreasing = TRUE))

# DE analysis using linear modelling
vfit <- lmFit(v, limma_des)
vfit <- contrasts.fit(vfit, contrasts=my_contrasts)
efit <- eBayes(vfit)
plotSA(efit)

dt <- decideTests(efit)
summary(dt)
#topTreat(efit, coef=5, n=10000, p.value=0.05) # coef dictates the contrast column (1:4 in this case)
l <- length(toTest)
path <- "~/Slim/SSa39_Pita/DE_analysis/"
for (i in 1:l) {
  # #print(nrow(tags))
  DE <- topTreat(efit, coef=i, n=1000, p.value=0.05)
  txt_table <- paste(path, toTest[i], "_limma.csv", sep = "") # print full "tags" table of DE genes
  txt_list <- paste(path, toTest[i], "_limma.txt", sep = "") # print list of DE genes
  write.table(rownames(DE), file=txt_list, quote=FALSE, row.names=FALSE, col.names=FALSE) 
  write.table(round(DE, 4), file=txt_table, quote=FALSE)
  png(paste(path, toTest[i], "_limma.png", sep = ""), width=1000, height=650, res = 150)
  plotMD(efit, column=i, status=dt[,i])
  abline(h=c(-1, 1), col="blue")
  dev.off()
}

#########################################
# Gene annotation using camera and mSigDB gene sets
readGMT <- function(inputFile){
  con <- file(inputFile, open = "r")
  dataList <- list()
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    myVector <- unlist(strsplit(oneLine, "\t"))
    dataList <- c(dataList, list(myVector[3:length(myVector)]))
    names(dataList)[length(dataList)] <- myVector[1]
  }
  close(con)
  return(dataList)
}

H_File <- "~/Slim/Downloads/GeneSets/h.all.v6.1.symbols.gmt"
C2_File <- "~/Slim/Downloads/GeneSets/c2.all.v6.1.symbols.gmt"

Hum_H <- readGMT(H_File)
Hum_C2 <- readGMT(C2_File)
idx_H <- ids2indices(Hum_H, id=rownames(v))
idx_C2 <- ids2indices(Hum_C2, id=rownames(v))
l <- length(toTest)
for (i in 1:l) {
  cam.H <- camera(v, idx_H, limma_des, contrast=my_contrasts[,toTest[i]], inter.gene.cor=0.01)
  cam.C2 <- camera(v, idx_C2, limma_des, contrast=my_contrasts[,toTest[i]], inter.gene.cor=0.01)
  H.list <- paste(path, toTest[i], "_mSigDB_H.csv", sep = "")
  C2.list <- paste(path, toTest[i], "_mSigDB_C2.csv", sep = "")
  write.table(head(cam.H, n=30), file=H.list, quote=FALSE) 
  write.table(head(cam.C2, n=30), file=C2.list, quote=FALSE) 
}

#y_glm <- estimateGLMCommonDisp(y, design, verbose=TRUE)
#Disp = 0.00879 , BCV = 0.0937 # without mock replicates
#Disp = 0.03727 , BCV = 0.1931 # with Seq1 mock replicates from Zika
#Disp = 0.04809 , BCV = 0.2193 # with Seq1+2 mock replicates from Zika
#Disp = 0.07928 , BCV = 0.2816 # with Seq1+3 mock replicates 
#Disp = 0.03561 , BCV = 0.1887 # with Seq2 mock replicates
#Disp = 0.06641 , BCV = 0.2577 # with Seq3 mock replicates 

########################################################################
# Alignment Metrics
library(reshape2)
library(ggplot2)
my_dir<-"/home/justin/R/Slim/SSa39_Pita/Data/AlignmentMetrics/"
tmp<-list.files(my_dir, pattern="*txt")
my_metrics<-sapply(tmp,function(x) read.table(paste(my_dir, x, sep="/"), sep="\t", nrows=1, header=T, fill=T))
my_metrics <- t(data.frame(my_metrics, check.names=F))
tmp<-strsplit(rownames(my_metrics), "_", fixed=T)
rownames(my_metrics)<-lapply(tmp, function(x) paste(x[1], collapse="_"))
my_metrics <- as.data.frame(apply(my_metrics, 2, unlist))

plot(sort(my_metrics$PCT_USABLE_BASES), main="Percent Usable Bases")
plot(sort(my_metrics$PCT_CODING_BASES), main="Percent Coding Bases")
plot(sort(my_metrics$PCT_UTR_BASES), main="Percent URT Bases")
plot(sort(my_metrics$PCT_INTERGENIC_BASES), main="Percent Intergenic Bases")

my_subset <- my_metrics[,12:15]
my_subset$SAMPLE <- rownames(my_subset)
my_melted_samp <- melt(my_subset, id.vars="SAMPLE")
ggplot(data = my_melted_samp, aes(x = SAMPLE, y=value, fill = variable)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand=c(0,0)) + ggtitle("Exp66B Read Base Distribution") + ylab("Proportion of Bases")

########################################################################
# Heatmap
logcpm <- cpm(y, prior.count=2, log=TRUE)
colnames(logcpm) <- c("293-Mock-1","293-Pita-1","293-Pita-2","CH22-Mock-1","CH22-Pita-1","CH22-Pita-2",
                      "UCH2-Mock-1","UCH2-Pita-1","UCH2-Pita-2","293-Mock-2","CH22-Mock-2","UCH2-Mock-2")
logcpm <- logcpm[,c(1,10,2,3,4,11,5,6,7,12,8,9)]
logcpm.subset  <-logcpm[rownames(head(UCH2.Pita[[1]], 40)), ]
require(gplots)
col.pan <- colorpanel(100,"blue","white","red")
heatmap.2(logcpm.subset, col=col.pan, Colv=TRUE, scale="none",trace="none", dendrogram="both", cexRow=1, 
          cexCol=1.4, density.info="none",margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
dev.off()
