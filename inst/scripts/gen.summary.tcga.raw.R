
library(shinyMethyl)
dataDir <- "/dcs01/oncbio/ejfertig/HNSCC450K"


# RGSetHNSC is an RGChannelSet containing the 369 samples from TCGA HNSC
setwd(dataDir)
load("RGSetHNSC.Rda")
load("cov.tcga.Rda")
summary.tcga.raw <- shinySummarize(RGSet)
summary.tcga.raw@phenotype <- cov.tcga

save(summary.tcga.raw, file = "summary.tcga.raw.rda")


