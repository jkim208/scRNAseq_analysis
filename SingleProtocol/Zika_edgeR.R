# Derived from Stefan's Analysis of the Zika dataset

library(reshape2)
library(ggplot2)
library(edgeR)
library(knitr)
source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library("org.Hs.eg.db")

samples<-read.table("~/R/Slim/ChordomaCor/Data/Counts_File_Mapping.txt", header=T)
head(samples)

y<-readDGE(samples, path="~/R/Slim/ChordomaCor/Data/All_Counts", group=samples$label)
tail(y$counts)


#Remove meta tags:
to_remove<-!(startsWith(rownames(y$counts), "__"))
tail(to_remove)
y <- y[to_remove,,keep.lib.sizes=FALSE]
tail(y$counts)


#Drop the technical replicates into one group:
y<-sumTechReps(y, ID=y$samples$group)
#How many million reads passed through all the mapping and counting stages in each library?
sum(y$samples[startsWith(as.character(y$samples$group), "Seq1"),"lib.size"])/1000000
sum(y$samples[startsWith(as.character(y$samples$group), "Seq2"),"lib.size"])/1000000
sum(y$samples[startsWith(as.character(y$samples$group), "Seq3"),"lib.size"])/1000000

#Show what we would filter out by cutting at a CPM of 1:
#plot(cpm(y)[,1], ylim=c(0,100))
#abline(h=1)

dim(y)
keep <- rowSums(cpm(y)>1) >= 1
y<-y[keep, , keep.lib.sizes=FALSE]
dim(y)

####### Removed repeated code
#"Seq3-293-Mock-1" seems far off from the rest of the 293 samples. Looking at the read base distributions, it does appear something went wrong during the library construction - the majority of the library is intronic/intergenic. The other 293 samples in Seq3 also look off from the read base dist. chart, but seem ok in the MDS plot - let's remove only the first sample
y <- y[,colnames(y) != "Seq3-293-Mock-1"]

#Now let's recalculate the norm factors and MDS plot:
y<-calcNormFactors(y)
tmp <- plotMDS(y)
my_coords <- as.data.frame(tmp$cmdscale.out)
colnames(my_coords) <- c("x","y")

my_coords$sample <- rownames(my_coords)
my_coords$Seq <- rep(c("1", "2", "3"), each=12)[1:35]
tmp<-strsplit(my_coords$sample, "-", fixed=T)
my_coords$Cells<-as.character(lapply(tmp, function(x) paste(x[2], collapse="_")))
my_coords$Cells[1:12] <- rep(c("293", "CH22", "KHOS", "MugChord", "UCH1", "UCH2"), each=2)
my_coords$Cells[13:16] <- "293"

ggplot(data=my_coords, aes(x=x, y=y)) + geom_text(aes(label=sample, col=Seq))
ggplot(data=my_coords, aes(x=x, y=y)) + geom_text(aes(label=sample, col=Cells))

kable(y$samples)


my_samples <- data.frame(my_coords[,3:5])
my_samples$Treatment <- rep(c("Mock", "GFP"), times=18)[1:35]
my_samples$Treatment[13:24] <- rep(c("GFP", "Mock"), each=2, times=3)
my_samples$Treatment[25] <- "Mock"
my_samples$Treatment[26:33] <- rep(c("PolyU", "Mock"), each=2, times=2)
my_samples$Treatment[34:35] <- rep(c("PolyU"), each=2)
my_samples$Cell_Group <- "Chordoma"
my_samples[c(1:2,5:6,13:16,25:27),"Cell_Group"] <- "Non-Chord"
my_samples$Treatment <- as.factor(my_samples$Treatment)
my_samples$Treatment <- relevel(my_samples$Treatment, ref="Mock")
my_samples$Cells <- as.factor(my_samples$Cells)
my_samples$Cells <- relevel(my_samples$Cells, ref="293")
my_samples$Cell_Group <- as.factor(my_samples$Cell_Group)
my_samples$Cell_Group <- relevel(my_samples$Cell_Group, ref="Non-Chord")

# blocking on cell factor. Design is focused on how trt affects each cell type. Not just difference between trts.
design <- model.matrix(~Cells + Cells:Treatment, data=my_samples)
#design <- model.matrix(~Cell_Group + Cell_Group:Cells + Cell_Group:Treatment, data=my_samples)
design <- design[,colSums(design) != 0]
#design <- model.matrix(~Cells + Treatment, data=my_samples)

y_glm <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y_glm <- estimateGLMTrendedDisp(y_glm, design)
y_glm <- estimateGLMTagwiseDisp(y_glm, design)
fit <- glmFit(y_glm, design)

my_contrasts <- makeContrasts(
  # Difference between chordoma cells vs non-chordoma cells that are being compared from Mock to GFP-transfected states
  chord.GFP = (CellsCH22.TreatmentGFP + CellsMugChord.TreatmentGFP + CellsUCH1.TreatmentGFP + CellsUCH2.TreatmentGFP)/4 - (Cells293.TreatmentGFP + CellsKHOS.TreatmentGFP)/2,
  # Difference between chordoma cells vs non-chordoma cells that are being compared from Mock to PolyU-transfected states
  chord.PolyU = (CellsUCH1.TreatmentPolyU + CellsUCH2.TreatmentPolyU)/2 - Cells293.TreatmentPolyU,

  chord.PolyU.vs.GFP = ((CellsUCH1.TreatmentPolyU + CellsUCH2.TreatmentPolyU)/2 - Cells293.TreatmentPolyU) - ((CellsCH22.TreatmentGFP + CellsMugChord.TreatmentGFP + CellsUCH1.TreatmentGFP + CellsUCH2.TreatmentGFP)/4 - (Cells293.TreatmentGFP + CellsKHOS.TreatmentGFP)/2),
  chord = (CellsCH22 + CellsMugChord + CellsUCH1 + CellsUCH2)/4,
  levels=make.names(colnames(design)))

plotBCV(y_glm)


lrt <- glmLRT(fit, contrast=my_contrasts[,"chord.GFP"])
chord.GFP <- topTags(lrt, n=100000, p.value=0.05)

plotSmear(lrt, de.tags=rownames(chord.GFP))
abline(h=c(-1, 1), col="blue")

lrt <- glmLRT(fit, contrast=my_contrasts[,"chord.PolyU"])
chord.PolyU <- topTags(lrt, n=100000, p.value=0.05)

plotSmear(lrt, de.tags=rownames(chord.PolyU))
abline(h=c(-1, 1), col="blue")

lrt <- glmLRT(fit, contrast=my_contrasts[,"chord.PolyU.vs.GFP"])
chord.PolyU.vs.GFP <- topTags(lrt, n=100000, p.value=0.05)

plotSmear(lrt, de.tags=rownames(chord.PolyU.vs.GFP))
abline(h=c(-1, 1), col="blue")

lrt <- glmLRT(fit, contrast=my_contrasts[,"chord"])
chord <- topTags(lrt, n=100000, p.value=0.05)

plotSmear(lrt, de.tags=rownames(chord))
abline(h=c(-1, 1), col="blue")


library(Glimma)
limma_des <- design
colnames(limma_des) <- make.names(colnames(limma_des))

v <- voom(y, limma_des, plot=TRUE)
vfit <- lmFit(v, limma_des)
vfit <- contrasts.fit(vfit, contrasts=my_contrasts)
efit <- eBayes(vfit)
plotSA(efit)

dt <- decideTests(efit)
summary(dt)
#        chord.GFP chord.PolyU chord.PolyU.vs.GFP chord
# Down          39           0                  2  7254
# NotSig     22586       22799              22836  8805
# Up           213          39                  0  6779

topTreat(efit, coef=1, n=50) # coef dictates the contrast column (1:4 in this case)
topTreat(efit, coef=2, n=50)
topTreat(efit, coef=3, n=50)
topTreat(efit, coef=4, n=50)

plotMD(efit, column=1, status=dt[,1])

#Gene set testing with Camera using the Broad MSigDB collection
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

in_File <- "~/R/Slim/Downloads/GeneSets/h.all.v6.1.symbols.gmt"

Hum_H <- readGMT(in_File)
idx_H <- ids2indices(Hum_H, id=rownames(v))

cam.H.ChordGFP <- camera(v, idx_H, limma_des, contrast=my_contrasts[,1], inter.gene.cor=0.01)
head(cam.H.ChordGFP, n=30)
cam.H.ChordPolyU <- camera(v, idx_H, limma_des, contrast=my_contrasts[,2], inter.gene.cor=0.01)
head(cam.H.ChordPolyU, n=30)
cam.H.Chord.Diff.GFP.PolyU <- camera(v, idx_H, limma_des, contrast=my_contrasts[,3], inter.gene.cor=0.01)
head(cam.H.Chord.Diff.GFP.PolyU, n=30)
cam.H.Chord <- camera(v, idx_H, limma_des, contrast=my_contrasts[,4], inter.gene.cor=0.01)
head(cam.H.Chord, n=30)

#write.table(cam.H.ChordGFP, file="/Users/Stefan/Documents/CCIB_Work/Zika/Sequencing/All_Analysis/Chordoma_GFP_H_SigNetworks.txt")
#write.table(cam.H.ChordPolyU, file="/Users/Stefan/Documents/CCIB_Work/Zika/Sequencing/All_Analysis/Chordoma_PolyU_H_SigNetworks.txt")
#write.table(cam.H.Chord.Diff.GFP.PolyU, file="/Users/Stefan/Documents/CCIB_Work/Zika/Sequencing/All_Analysis/Chordoma_Diff-GFP-PolyU_H_SigNetworks.txt")
#write.table(cam.H.Chord, file="/Users/Stefan/Documents/CCIB_Work/Zika/Sequencing/All_Analysis/Chordoma_H_SigNetworks.txt")

in_File <- "~/R/Slim/Downloads/GeneSets/c2.all.v6.1.symbols.gmt"

Hum_C2 <- readGMT(in_File)
idx <- ids2indices(Hum_C2, id=rownames(v))
cam.C2.ChordGFP <- camera(v, idx, limma_des, contrast=my_contrasts[,1], inter.gene.cor=0.01)
head(cam.C2.ChordGFP, n=30)
cam.C2.ChordPolyU <- camera(v, idx, limma_des, contrast=my_contrasts[,2], inter.gene.cor=0.01)
head(cam.C2.ChordPolyU, n=30)
cam.C2.Chord.Diff.GFP.PolyU <- camera(v, idx, limma_des, contrast=my_contrasts[,3], inter.gene.cor=0.01)
head(cam.C2.Chord.Diff.GFP.PolyU, n=100)
cam.C2.Chord <- camera(v, idx, limma_des, contrast=my_contrasts[,4], inter.gene.cor=0.01)
head(cam.C2.Chord, n=30)

#write.table(cam.C2.ChordGFP, file="/Users/Stefan/Documents/CCIB_Work/Zika/Sequencing/All_Analysis/Chordoma_GFP_C2_SigNetworks.txt")
#write.table(cam.C2.ChordPolyU, file="/Users/Stefan/Documents/CCIB_Work/Zika/Sequencing/All_Analysis/Chordoma_PolyU_C2_SigNetworks.txt")
#write.table(cam.C2.Chord.Diff.GFP.PolyU, file="/Users/Stefan/Documents/CCIB_Work/Zika/Sequencing/All_Analysis/Chordoma_Diff-GFP-PolyU_C2_SigNetworks.txt")
#write.table(cam.C2.Chord, file="/Users/Stefan/Documents/CCIB_Work/Zika/Sequencing/All_Analysis/Chordoma_C2_SigNetworks.txt")








