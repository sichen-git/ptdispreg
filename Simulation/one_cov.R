# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#about 12h on sharcnet

require(HelpersMG)

source("../../pt_functions/pt_base.R")
source("../../pt_functions/ptreg_fun.R")

samp=1000
a=seq(-2,0.6, by=0.1)

#specify true coefficient values
beta=matrix(c(log(10), 0.8), ncol=1)
gamma=matrix(c(1, 2), ncol=1)


for(j in a)
  for(i in 1:1000){
    #generate X covariates
    set.seed(i)
    x0=rep(1,samp)
    x1=runif(samp, min=0, max=1)
    X=as.matrix(cbind(x0,x1))
    
    #generate U covariates
    u0=rep(1,samp)
    u1=runif(samp, min=0, max=1)
    U=as.matrix(cbind(u0, u1))
    
    #generate  Y response
    Y=matrix(rpt_reg_disp(j, beta, gamma, X, U, offsetb = NULL, offsetg = NULL), ncol=1)
    
    
    
    a.ini=max(a.moment(skewness(Y),mean(Y),var(Y)/mean(Y)),-5)
    beta_ini=c(log(mean(Y), base = exp(1)), 0)
    gamma_ini=c(log((var(Y)/mean(Y)-1), base = exp(1)), 0)
    initial.v=c(a.ini,beta_ini,gamma_ini)
    
    
    model.disp.tmp=nlm(pt_lik, initial.v, y=Y, X=X, U=U, hessian=T)
    
    var.tmp=SEfromHessian(model.disp.tmp$hessian)^2
    
    # collect results
    param.tmp=matrix(c(i, j, model.disp.tmp$estimate, var.tmp), nrow=1)
    # write to file
    write.table(param.tmp, file="one_cov.csv",
                append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
  }







