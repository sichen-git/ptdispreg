setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

params_stress=read.csv("one_cov_stress.csv", sep=",", header=F)


samp=c(50,100)
a=c(-1, 0, 0.5)

#specify true coefficient values
beta=matrix(c(log(10), 0.8), ncol=1)
gamma=matrix(c(1, 2), ncol=1)


results=NULL
for(j in 1:length(a)){
  for(k in 1:length(samp)){
    values_stress=params_stress[params_stress[,2]==a[j]&params_stress[,3]==samp[k],4:15]
    converge=sum(values_stress[, 12]==0)/dim(values_stress)[1]
    finitese=sum(!is.na(values_stress[, 9]))/dim(values_stress)[1]
    converged=values_stress[values_stress[, 12]==0,]
    mean.est=apply(converged[,1:5],2, mean)
    bias.tmp=mean.est-c(a[j], beta, gamma)
    convse=values_stress[values_stress[, 12]==0&!is.na(values_stress[, 9]),]
    ub=convse[,1:5]+1.96*sqrt(convse[,6:10])
    lb=convse[,1:5]-1.96*sqrt(convse[,6:10])
    test=matrix(rep(c(a[j], beta, gamma), 500), ncol=5, byrow=T)
    coverage.tmp=ub>test&lb<test
    coverageprob=apply(coverage.tmp,2, sum)/500
    
    results.tmp=rbind(c(a[j],samp[k], converge, finitese, "bias", round(bias.tmp, digits = 4)),
                      c(a[j],samp[k], converge, finitese, "95\\% coverage", coverageprob))
    results=rbind(results, results.tmp)
}
}

colnames(results)=c("a", "sample size", 
                    "convergence", 
                    "Finite SE", 
                    "", 
                    "hat a", 
                    "hat beta0",
                    "hat beta1",
                    "hat gamma0",
                    "hat gamma1"
                    )

require(kableExtra)

results%>%
  kable("latex", escape = F, longtable=F,
        caption="",
        booktabs = T, linesep = "", align = "c") %>%
  collapse_rows(columns = c(1:4)) %>%
  kable_styling(latex_options =c("HOLD_position", "repeat_header", "scale_down" ))
