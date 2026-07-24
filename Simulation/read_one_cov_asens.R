setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

params_mom=read.csv("one_cov_mom.csv", sep=",", header=F)
params_n2=read.csv("one_cov_inian2.csv", sep=",", header=F)
params_n1=read.csv("one_cov_inian1.csv", sep=",", header=F)
params_0=read.csv("one_cov_inia0.csv", sep=",", header=F)
params_05=read.csv("one_cov_inia05.csv", sep=",", header=F)


beta=matrix(c(log(10), 0.8), ncol=1)
gamma=matrix(c(1, 2), ncol=1)


a.uniq=unique(params_mom[,2])



results=NULL
for(j in 1:length(a.uniq)){
  values_mom=params_mom[params_mom[,2]==a.uniq[j],3:15]
  values_n2=params_n2[params_n2[,2]==a.uniq[j],3:15]
  values_n1=params_n1[params_n1[,2]==a.uniq[j],3:15]
  values_0=params_0[params_0[,2]==a.uniq[j],3:15]
  values_05=params_05[params_05[,2]==a.uniq[j],3:15]

  # •	nlm() code 1 or 2 or any other output;
  code_mom=sum(values_mom[,12])
  code_n2=sum(values_n2[,12])
  code_n1=sum(values_n1[,12])
  code_0=sum(values_n1[,12])
  code_05=sum(values_05[,12])
  
  # •	finite negative log-likelihood;
  finnl_mom=sum(values_mom[,11]>-Inf)
  finnl_n2=sum(values_n2[,11]>-Inf)
  finnl_n1=sum(values_n1[,11]>-Inf)
  finnl_0=sum(values_n1[,11]>-Inf)
  finnl_05=sum(values_05[,11]>-Inf)
  
  
  # •	numerically positive-definite observed-information matrix;
  hess_mom=sum(values_mom[,13])
  hess_n2=sum(values_n2[,13])
  hess_n1=sum(values_n1[,13])
  hess_0=sum(values_n1[,13])
  hess_05=sum(values_05[,13])
  # •	finite standard errors.
  finse_mom=sum(rowSums(values_mom[,6:10]<Inf)==5)
  finse_n2=sum(rowSums(values_n2[,6:10]<Inf)==5)
  finse_n1=sum(rowSums(values_n1[,6:10]<Inf)==5)
  finse_0=sum(rowSums(values_n1[,6:10]<Inf)==5)
  finse_05=sum(rowSums(values_05[,6:10]<Inf)==5)
  # obtain parameter values log-likelihood (or nll) 
  maxll_mom=-min(values_mom[,1])
  maxll_n2=-min(values_n2[,11])
  maxll_n1=-min(values_n1[,11])
  maxll_0=-min(values_n1[,11])
  maxll_05=-min(values_05[,11])
  # and percentage reaching the same optimum, using an explicit log-likelihood tolerance.
  onepern2=sum(rowSums((values_n2[,1:5]-values_mom[,1:5])/values_mom[,1:5]<0.01)==5)
  onepern1=sum(rowSums((values_n1[,1:5]-values_mom[,1:5])/values_mom[,1:5]<0.01)==5)
  oneper0=sum(rowSums((values_0[,1:5]-values_mom[,1:5])/values_mom[,1:5]<0.01)==5)
  oneper05=sum(rowSums((values_05[,1:5]-values_mom[,1:5])/values_mom[,1:5]<0.01)==5)
  ponepern2=sum(rowSums((values_n2[,1:5]-values_mom[,1:5])/values_mom[,1:5]<0.001)==5)
  ponepern1=sum(rowSums((values_n1[,1:5]-values_mom[,1:5])/values_mom[,1:5]<0.001)==5)
  poneper0=sum(rowSums((values_0[,1:5]-values_mom[,1:5])/values_mom[,1:5]<0.001)==5)
  poneper05=sum(rowSums((values_05[,1:5]-values_mom[,1:5])/values_mom[,1:5]<0.001)==5)
  
  reslt.tmp=cbind(c(rep(a.uniq[j], 5)),
    c("Moment", "a=-2", "a=-1", "a=0", "a=0.5"),
        c(code_mom, code_n2, code_n1, code_0, code_05),
        c(finnl_mom, finnl_n2, finnl_n1, finnl_0, finnl_05),
        c(hess_mom, hess_n2, hess_n1, hess_0, hess_05),
       c(finse_mom, finse_n2, finse_n1, finse_0, finse_05),
        c(maxll_mom, maxll_n2, maxll_n1, maxll_0, maxll_05),
        c(1000, onepern2, onepern1, oneper0, oneper05),
        c(1000, ponepern2, ponepern1, poneper0, poneper05))
  
  results=rbind(results, reslt.tmp)
}


colnames(results)=c("a", "initial a", 
                    "convergence", 
                    "ll<∞", 
                    "PD hessian", 
                    "SE<∞", 
                    "max ll",
                    "Rel. diff. estimates <0.01 (compared to moment)",
                    "Rel. diff. estimates <0.001 (compared to moment)"
                    )

require(kableExtra)

results%>%
  kable("latex", escape = F, longtable=F,
        caption="",
        booktabs = T, linesep = "", align = "c") %>%
  collapse_rows(columns = 1) %>%
  kable_styling(latex_options =c("HOLD_position", "repeat_header", "scale_down" ))  %>%
  collapse_rows(columns = 1)
