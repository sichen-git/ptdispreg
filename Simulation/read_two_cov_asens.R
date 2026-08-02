setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

params_mom=read.csv("two_cov_mom.csv", sep=",", header=F)
params_n2=read.csv("two_cov_inian2.csv", sep=",", header=F)
params_n1=read.csv("two_cov_inian1.csv", sep=",", header=F)
params_0=read.csv("two_cov_inia0.csv", sep=",", header=F)
params_05=read.csv("two_cov_inia05.csv", sep=",", header=F)


beta=matrix(c(log(10), 0.8, -1), ncol=1)
gamma=matrix(c(1, 2, 1.5), ncol=1)


a.uniq=unique(params_mom[,2])



results=NULL
fits=NULL
for(j in 1:length(a.uniq)){
  values_mom=params_mom[params_mom[,2]==a.uniq[j],3:19]
  values_n2=params_n2[params_n2[,2]==a.uniq[j],3:19]
  values_n1=params_n1[params_n1[,2]==a.uniq[j],3:19]
  values_0=params_0[params_0[,2]==a.uniq[j],3:19]
  values_05=params_05[params_05[,2]==a.uniq[j],3:19]

  # •	nlm() code 1 or 2 or any other output;
  code_mom=sum(values_mom[,16]==1|values_mom[,16]==2)
  code_n2=sum(values_n2[,16]|values_n2[,16]==2)
  code_n1=sum(values_n1[,16]|values_n1[,16]==2)
  code_0=sum(values_n1[,16]|values_n1[,16]==2)
  code_05=sum(values_05[,16]|values_05[,16]==2)
  
  
  # •	finite negative log-likelihood;
  finnl_mom=sum(values_mom[,15]>-Inf)
  finnl_n2=sum(values_n2[,15]>-Inf)
  finnl_n1=sum(values_n1[,15]>-Inf)
  finnl_0=sum(values_n1[,15]>-Inf)
  finnl_05=sum(values_05[,15]>-Inf)
  # •	finite estimates.
  finpar_mom=sum(rowSums(values_mom[,1:7]<Inf)==7)
  finpar_n2=sum(rowSums(values_n2[,1:7]<Inf)==7)
  finpar_n1=sum(rowSums(values_n1[,1:7]<Inf)==7)
  finpar_0=sum(rowSums(values_n1[,1:7]<Inf)==7)
  finpar_05=sum(rowSums(values_05[,1:7]<Inf)==7)
  
  # •	numerically positive-definite observed-information matrix;
  hess_mom=sum(values_mom[,17])
  hess_n2=sum(values_n2[,17])
  hess_n1=sum(values_n1[,17])
  hess_0=sum(values_n1[,17])
  hess_05=sum(values_05[,17])
  # •	finite standard errors.
  finse_mom=sum(rowSums(values_mom[,8:14]<Inf)==7)
  finse_n2=sum(rowSums(values_n2[,8:14]<Inf)==7)
  finse_n1=sum(rowSums(values_n1[,8:14]<Inf)==7)
  finse_0=sum(rowSums(values_n1[,8:14]<Inf)==7)
  finse_05=sum(rowSums(values_05[,8:14]<Inf)==7)
  # obtain parameter values log-likelihood (or nll) 
  maxll_mom=-min(values_mom[,15])
  maxll_n2=-min(values_n2[,15])
  maxll_n1=-min(values_n1[,15])
  maxll_0=-min(values_n1[,15])
  maxll_05=-min(values_05[,15])
  # and percentage reaching the same optimum, using an explicit log-likelihood tolerance.
  onepern2=sum(rowSums((values_n2[,1:7]-values_mom[,1:7])/values_mom[,1:7]<0.01)==7)
  onepern1=sum(rowSums((values_n1[,1:7]-values_mom[,1:7])/values_mom[,1:7]<0.01)==7)
  oneper0=sum(rowSums((values_0[,1:7]-values_mom[,1:7])/values_mom[,1:7]<0.01)==7)
  oneper05=sum(rowSums((values_05[,1:7]-values_mom[,1:7])/values_mom[,1:7]<0.01)==7)
  ponepern2=sum(rowSums((values_n2[,1:7]-values_mom[,1:7])/values_mom[,1:7]<0.001)==7)
  ponepern1=sum(rowSums((values_n1[,1:7]-values_mom[,1:7])/values_mom[,1:7]<0.001)==7)
  poneper0=sum(rowSums((values_0[,1:7]-values_mom[,1:7])/values_mom[,1:7]<0.001)==7)
  poneper05=sum(rowSums((values_05[,1:7]-values_mom[,1:7])/values_mom[,1:7]<0.001)==7)
  
  #overall admissibility
  
  overall_mom=min(code_mom, finnl_mom, finpar_mom, hess_mom, finse_mom, 1000, 1000)
  overall_n2=min(code_n2, finnl_n2, finpar_n2, hess_n2, finse_n2, onepern2, ponepern2)
  overall_n1=min(code_n1, finnl_n1, finpar_n1, hess_n1, finse_n1, onepern1, ponepern1)
  overall_0=min(code_0, finnl_0, finpar_0, hess_0, finse_0, oneper0, poneper0)
  overall_05=min(code_05, finnl_05, finpar_05, hess_05, finse_05, oneper05, poneper05)
  
  reslt.tmp=cbind(c(rep(a.uniq[j], 5)),
                  c("Moment", "a=-2", "a=-1", "a=0", "a=0.5"),
                  c(code_mom, code_n2, code_n1, code_0, code_05),
                  c(finnl_mom, finnl_n2, finnl_n1, finnl_0, finnl_05),
                  c(finpar_mom, finpar_n2, finpar_n1, finpar_0, finpar_05),
                  c(hess_mom, hess_n2, hess_n1, hess_0, hess_05),
                  c(finse_mom, finse_n2, finse_n1, finse_0, finse_05),
                  c(maxll_mom, maxll_n2, maxll_n1, maxll_0, maxll_05),
                  c(1000, onepern2, onepern1, oneper0, oneper05),
                  c(1000, ponepern2, ponepern1, poneper0, poneper05),
                  c(overall_mom, overall_n2, overall_n1, overall_0, overall_05))
  
  results=rbind(results, reslt.tmp)

  
  diffll_n2=abs(values_n2[,15]-values_mom[,15])<0.000001
  diffll_n1=abs(values_n1[,15]-values_mom[,15])<0.000001
  diffll_0=abs(values_0[,15]-values_mom[,15])<0.000001
  diffll_05=abs(values_05[,15]-values_mom[,15])<0.000001
  diffll=sum(diffll_n2&diffll_n1&diffll_0&diffll_05)
  fitting.tmp=min(c(code_mom, code_n2, code_n1, code_0, code_05),
                  c(finnl_mom, finnl_n2, finnl_n1, finnl_0, finnl_05),
                  c(hess_mom, hess_n2, hess_n1, hess_0, hess_05),
                  c(finse_mom, finse_n2, finse_n1, finse_0, finse_05),
                  c(finpar_mom, finpar_n2, finpar_n1, finpar_0, finpar_05))
  common.tmp=diffll
  
  fits=rbind(fits, c(a.uniq[j], fitting.tmp, common.tmp, 1000-fitting.tmp))
}
fits.table=data.frame("True a"=fits[,1], "Dataset Refitted"=rep(1000, 4), 
                      "Successful moment-based fits"=paste(fits[,2], "/1000", " ", fits[,2]/10, "%",sep=""), 
                      "Common optimum"=paste(fits[,3], "/1000", " ", fits[,3]/10, "%", sep = ""),
                      "Alternative-start"=fits[,4])
overall=apply(fits, 2, sum)

fits.table=rbind(fits.table,c("Overall", 4000, 
                              paste(overall[2], "/4000", " ", overall[2]/40, "%",sep=""),
                              paste(overall[3], "/4000", " ", overall[3]/40, "%",sep=""),
                              overall[4]))


colnames(results)=c("$a$", "initial $a$", 
                    "convergence", 
                    "ll<∞", 
                    "estimates<∞", 
                    "PD hessian", 
                    "SE<∞", 
                    "max ll",
                    "Rel. diff. estimates <0.01 (compared to moment)",
                    "Rel. diff. estimates <0.001 (compared to moment)",
                    "Overall Admissibility"
)

require(kableExtra)

results%>%
  kable("latex", escape = F, longtable=F,
        caption="",
        booktabs = T, linesep = "", align = "c") %>%
  collapse_rows(columns = 1) %>%
  kable_styling(latex_options =c("HOLD_position", "repeat_header", "scale_down" ))  %>%
  collapse_rows(columns = 1)

fits.table%>%
  kable("latex", escape = F, longtable=F,
        caption="",
        booktabs = F, linesep = "", align = "c") %>%
  collapse_rows(columns = 1) %>%
  kable_styling(latex_options =c("HOLD_position" ))  %>%
  collapse_rows(columns = 1)
