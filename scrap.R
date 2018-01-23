so <- RunPCA(object = so, pc.genes = so@var.genes, 
             do.print = TRUE, pcs.print = 1:12, genes.print = 5, pcs.compute = 20)
PCElbowPlot(object = so)
so <- SetAllIdent(so, id="orig.ident")
so <- RunTSNE(object = so, dims.use = 1:6, do.fast = TRUE, perplexity = 30)
t <- TSNEPlot(object = so, do.return = TRUE, do.label = TRUE, pt.size = 2)
t + ggtitle('so tSNE Seurat v2.1') + 
  theme(text = element_text(size=15))


> RunPCA(object = so, pc.genes = so@var.genes, 
         +        do.print = TRUE, pcs.print = 1:12, genes.print = 5, pcs.compute = 20)
[1] "PC1"
[1] "SERPINA3" "SQSTM1"   "S100A4"   "SOD2"     "PMEPA1"  
[1] "RNASE1" "KRT18"  "ID4"    "HOPX"   "GDA"   
[1] "PC2"
[1] "SCAND1" "ECI1"   "ASNS"   "KRT18"  "ACTR1A"
[1] "COL5A2" "CSPG4"  "ATP1B1" "CAPNS1" "COPZ2" 
[1] "PC3"
[1] "MPG"      "NDUFA4L2" "GPRC5A"   "SPINK6"   "S100A11" 
[1] "TLN2"     "ATP6V1E2" "DGKQ"     "SFRP1"    "HERC3"   
[1] "PC4"
[1] "NR3C2"   "ID1"     "FAM161B" "SPOCK1"  "MARS2"  
[1] "RP11-1260E13.1" "STEAP1B"        "T"              "ACAN"           "ALPP"          
[1] "PC5"
[1] "C20orf194" "KRT18P3"   "RGMA"      "SSH3"      "TP53I3"   
[1] "CCDC88B" "TIMP3"   "GJA1"    "KITLG"   "ASL"    
> so <- readRDS(file="Robj/ddSeq.Robj")
> so <- RunPCA(object = so, pc.genes = so@var.genes, 
               +              do.print = TRUE, pcs.print = 1:12, genes.print = 5, pcs.compute = 20)
[1] "PC1"
[1] "CCDC80"   "LGALS1"   "CTSC"     "SERPINE1" "ABI3BP"  
[1] "HSPA1A" "XIST"   "SET"    "ELP5"   "ZIC2"  
[1] "PC2"
[1] "RNASE1" "RGS2"   "CYBA"   "LGALS1" "ATP5J2"
[1] "SQSTM1"   "COL5A2"   "NEAT1"    "SERPINA3" "SLC25A37"
[1] "PC3"
[1] "HLA-B"    "S100A4"   "SERPINA3" "HLA-A"    "IFITM3"  
[1] "NEAT1"    "XIST"     "COL8A1"   "HSP90AA1" "HSPG2"   
[1] "PC4"
[1] "GDA"     "ID4"     "OLFML2A" "HOPX"    "SLC4A11"
[1] "COL8A1"   "TGFBI"    "SOD2"     "SERPINA3" "SERPINE1"
[1] "PC5"
[1] "COL1A2" "SCAMP4" "GCDH"   "IDUA"   "BICD1" 
[1] ""
[1] "CDC20" "T"     "FOSL1" "EZR"   "AZIN1"

