

#install.packages("shinyMethyl_1.0.4.tar.gz", repos=NULL)
library(shinyMethyl)

dataDir <- "/dcs01/oncbio/ejfertig/HNSCC450K"

# RGSetHNSC is an RGChannelSet containing the 369 samples from TCGA HNSC:
setwd(dataDir)
load("RGSetHNSC.Rda")
load("cov.tcga.Rda")


#mappings <- read.csv("HNSC.mappings.csv")


# Slide index:
slide <- substr(colnames(RGSet),1,10) 
slide <- model.matrix(~as.factor(slide)-1)
slide <- matrix(as.numeric(slide), nrow(slide),ncol(slide))


# Normalization with FunNorm:
# preprocessFunnormExtended.R contains an extension of functional normalization:
source("/home/student/jfortin/headMeth/scripts/preprocessFunnormExtended.R")
grSet.norm <- preprocessFunnormExtended(RGSet, nPCs = 2, design = slide)

summary.tcga.norm <- shinySummarize(grSet.norm)
summary.tcga.norm@phenotype <- cov.tcga
save(summary.tcga.norm, file = "summary.tcga.norm.rda")

