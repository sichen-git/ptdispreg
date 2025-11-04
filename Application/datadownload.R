setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#################################################################################
### Download data from TCGA
### This may be skipped - the data file brcaExp.rda is already available within this folder

#BiocManager::install(c("TCGAbiolinks"))
library(TCGAbiolinks)

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")

GDCdownload(query)
exp.brca<-GDCprepare(query = query, save = TRUE, save.filename = "../data/brcaExp.rda")
#################################################################################
citation('TCGAbiolinks')
