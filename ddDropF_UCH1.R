library(Seurat, lib.loc = '~/R/x86_64-pc-linux-gnu-library/3.4/Seurat2.1')
library(Matrix)
library(dplyr)
setwd("~/R/Projects/Seurat")
#####################################################################################
# Read in files
#####################################################################################
# Load filtered/normalized Seurat Objects from DropSeq/Robj and ddSeq/Robj
# DS50_UCH1, DS52_UCH1_HS_TB, DS52_UCH1_RT_ID
# Seq2_N706, Seq3_N707
ddSeq_293HEK_Seq1_N701 <- readRDS("Robj/ddSeq/ddSeq_293HEK_Seq1_N701.Robj")
ddSeq_UCH2_Seq1_N704 <- readRDS("Robj/ddSeq/ddSeq_UCH2_Seq1_N704.Robj")
ddSeq_UCH1_Seq2_N706 <- readRDS("Robj/ddSeq/ddSeq_UCH1_Seq2_N706.Robj")
ddSeq_UCH1_Seq3_N707 <- readRDS("Robj/ddSeq/ddSeq_UCH1_Seq3_N707.Robj")

DropSeq_293HEK_OP_DS9 <- readRDS("Robj/DropSeq/DropSeq_293HEK_OP_DS9.Robj")
DropSeq_293HEK_NO_DS9 <- readRDS("Robj/DropSeq/DropSeq_293HEK_NO_DS9.Robj")
DropSeq_293HEK_DS34 <- readRDS("Robj/DropSeq/DropSeq_293HEK_OP_DS34.Robj")
DropSeq_293HEK_DS45 <- readRDS("Robj/DropSeq/DropSeq_293HEK_OP_DS45.Robj")
DropSeq_UCH2_DS49 <- readRDS("Robj/DropSeq/DropSeq_UCH2_DS49.Robj")
DropSeq_UCH1_DS50 <- readRDS("Robj/DropSeq/DropSeq_UCH1_DS50.Robj")
DropSeq_UCH1_DS52_RT_ID <- readRDS("Robj/DropSeq/DropSeq_UCH1_RT_ID_DS52.Robj")
DropSeq_UCH1_DS52_HS_TB <- readRDS("Robj/DropSeq/DropSeq_UCH1_HS_TB_DS52.Robj")

Fluidigm <- readRDS(file = '~/R/Projects/Seurat/Fluidigm/Fluidigm_UCH1.Robj')
