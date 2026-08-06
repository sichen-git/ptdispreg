setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../pt_functions/pt_base.R")
source("../pt_functions/ptreg_fun.R")


sizefactors=readRDS("../data/sizefactor.rds")
rawcounts=readRDS("../data/raw_counts.rds")
X_cov=readRDS("../data/cov_X.rds")
Z_cov=readRDS("../data/cov_Z.rds")

#LY6E-DT gene

geneid1=which(rownames(rawcounts)=="ENSG00000247317.3") #LY6E-DT

Y1=matrix(rawcounts[geneid1,], ncol=1)

#LINC00674 gene 

geneid2=which(str_starts(rownames(rawcounts), "ENSG00000237854")=="TRUE") #LINC00674

Y2=matrix(rawcounts[geneid2,],ncol=1)


timeparams=read.csv("reg_time.csv", sep=",", header=F)

colnames(timeparams)=c("Gene", 
                       "Run",
                       "Runtime",
                       "nlm code",
                       "nll",
                       "a","b0","b1","b2","b3","g0","g1",
                       "a<1",
                       "pd hess",
                       "se a","se b0","se b1","se b2","se b3","se g0","se g1")

require(kableExtra)

timeparams%>%
  kable("latex", escape = F, longtable=F,
        caption="",
        booktabs = T, linesep = "", align = "c") %>%
  collapse_rows(columns = 1) %>%
  kable_styling(latex_options =c("HOLD_position", "repeat_header", "scale_down" ))  %>%
  collapse_rows(columns = 1)

gene_names=unique(timeparams$Gene)
counts=cbind(Y1, Y2)

gene.table=NULL
for(i in 1:2){
  timeparams.tmp=timeparams[timeparams$Gene==gene_names[i],]
  summary.tmp=cbind(gene_names[i],matrix(summary(timeparams.tmp$Runtime), nrow=1), dim(counts)[1], median(counts[, i]), max(counts[, i]))
  gene.table=rbind(gene.table, summary.tmp)
  }
colnames(gene.table)=c("Gene", "Min.", "1st Qu.",  "Median",    "Mean", "3rd Qu.",    "Max.", "sample size", "med. count", "max. count" )

gene.table%>%
  kable("latex", escape = F, longtable=F,
        caption="",
        booktabs = T, linesep = "", align = "c") %>%
  collapse_rows(columns = 1) %>%
  kable_styling(latex_options =c("HOLD_position", "repeat_header", "scale_down" ))  %>%
  collapse_rows(columns = 1)