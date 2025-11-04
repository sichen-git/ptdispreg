setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(SummarizedExperiment)
library(DESeq2)
library(sva)
library(parallel)
library(xtable)

load(file = "../data/brcaExp.rda")

# Gender selection for female (Remove missing gender)
dim(data)
data_f=data[,!is.na(data$gender)]
data_f=data_f[,data_f$gender=="female"]
table(data_f$gender)
table(table(colData(data_f)$patient))

# Remove missing Tumour Stage (NA)
table(colData(data_f)$paper_pathologic_stage)
data_f=data_f[,!is.na(data_f$paper_pathologic_stage)]
data_f=data_f[,!(data_f$paper_pathologic_stage=="NA")]
dim(data_f)

table(data_f$gender)
table(table(colData(data_f)$patient))

# Remove duplicated patients (keep data for tissue that was formalin-fixed paraffin-embedded (FFPE))
duplicated<-names(which(table(data_f$patient)>1))
data_f2<-data_f[,-which(pmatch(data_f$patient,duplicated,dup=TRUE)>0&!data_f$is_ffpe)]
dim(data_f2)
table(table(colData(data_f2)$patient))
table(colData(data_f2)$paper_pathologic_stage)

table(colData(data_f2)$paper_BRCA_Subtype_PAM50)
table(colData(data_f2)$paper_BRCA_Subtype_PAM50,colData(data_f2)$paper_pathologic_stage)


#################
### New stage variable Stage IV==1 Stage I,II,II==0
table(data_f2$paper_pathologic_stage)
data_f2=data_f2[,grep("Stage_I",data_f2$paper_pathologic_stage)]  ## remove all but stages I-IV
data_f2$stage=0
data_f2[,grep("Stage_IV",data_f2$paper_pathologic_stage)]$stage=1
table(data_f2$paper_pathologic_stage,data_f2$stage)
#################


# Remove genes where more than 5% of samples have zero count
dim(data_f2)
countzeros<-apply(assay(data_f2),1,function(x){sum(x==0)});summary(countzeros)
cut_percent<-dim(data_f2)[2]*0.05;cut_percent
data_f2<-data_f2[(countzeros<cut_percent),]
dim(data_f2)



# Normalize expression values with respect to library size and then use log transformation
dds <- DESeqDataSet(se = data_f2,design = ~ 1,ignoreRank = T)
dds<-estimateSizeFactors(dds)

saveRDS(dds$sizeFactor, "../data/sizefactor.rds")

rawcounts<-counts(dds,normalized=F)
saveRDS(rawcounts, "../data/raw_counts.rds")


dds<-counts(dds,normalized=TRUE)
exprs.vsd<-dds



saveRDS(exprs.vsd, "../data/normalizedcounts.rds")


# Create two (2) surrogate variables
X<-cbind(data_f2$stage)
X<-as.matrix(X)
Z0<-NULL
mod = model.matrix(~X)
mod0 = model.matrix(~1)

svobj.2= sva(exprs.vsd.log,mod,mod0=NULL,vfilter=2000,n.sv=2)
Z<-cbind(Z0,svobj.2$sv)
head(Z);dim(Z)



saveRDS(X, "../data/cov_X.rds")
saveRDS(Z, "../data/cov_Z.rds")