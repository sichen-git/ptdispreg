setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# about

require(HelpersMG)
require(matrixcalc)

source("../pt_functions/pt_base.R")
source("../pt_functions/ptreg_fun.R")


samp=1000
a=c(-2, -0.5, 0, 0.5)

#specify true coefficient values
beta=matrix(c(log(10), 0.8, -1), ncol=1)
gamma=matrix(c(1, 2, 1.5), ncol=1)


for(j in a){
  for(i in 1:1000){
    #generate X covariates
    set.seed(i)
    x0=rep(1,samp)
    x1=runif(samp, min=0, max=1)
    x2=runif(samp, min=0, max=1)
    X=as.matrix(cbind(x0,x1, x2))
    
    #generate U covariates
    u0=rep(1,samp)
    u1=runif(samp, min=0, max=1)
    u2=runif(samp, min=0, max=1)
    U=as.matrix(cbind(u0, u1,u2))
    
    #generate  Y response
    Y=matrix(rpt_reg_disp(j, beta, gamma, X, U, offsetb = NULL, offsetg = NULL), ncol=1)
    
    #moment based
    a.ini=max(a.moment(skewness(Y),mean(Y),var(Y)/mean(Y)),-5)
    beta_ini=c(log(mean(Y), base = exp(1)), 0, 0)
    gamma_ini=c(log((var(Y)/mean(Y)-1), base = exp(1)), 0, 0)
    initial.v=c(a.ini,beta_ini,gamma_ini)
    
    
    model.disp.tmp=nlm(pt_lik, initial.v, y=Y, X=X, U=U, hessian=T)
    
    se.tmp=SEfromHessian(model.disp.tmp$hessian)
    hesspd.tmp=is.positive.definite(model.disp.tmp$hessian, tol=1e-8)
    
    # collect results
    param.tmp=matrix(c(i, j, model.disp.tmp$estimate, se.tmp, model.disp.tmp$minimum, model.disp.tmp$code,hesspd.tmp), nrow=1)
    # write to file
    write.table(param.tmp, file="two_cov_mom.csv",
                append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
    
    #a=-2
    a.ini=-2
    beta_ini=c(log(mean(Y), base = exp(1)), 0, 0)
    gamma_ini=c(log((var(Y)/mean(Y)-1), base = exp(1)), 0, 0)
    initial.v=c(a.ini,beta_ini,gamma_ini)
    
    
    model.disp.tmp=nlm(pt_lik, initial.v, y=Y, X=X, U=U, hessian=T)
    
    se.tmp=SEfromHessian(model.disp.tmp$hessian)
    hesspd.tmp=is.positive.definite(model.disp.tmp$hessian, tol=1e-8)
    
    # collect results
    param.tmp=matrix(c(i, j, model.disp.tmp$estimate, se.tmp, model.disp.tmp$minimum, model.disp.tmp$code,hesspd.tmp), nrow=1)
    # write to file
    write.table(param.tmp, file="two_cov_inian2.csv",
                append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
    #a=-1
    a.ini=-1
    beta_ini=c(log(mean(Y), base = exp(1)), 0, 0)
    gamma_ini=c(log((var(Y)/mean(Y)-1), base = exp(1)), 0, 0)
    initial.v=c(a.ini,beta_ini,gamma_ini)
    
    
    model.disp.tmp=nlm(pt_lik, initial.v, y=Y, X=X, U=U, hessian=T)
    
    se.tmp=SEfromHessian(model.disp.tmp$hessian)
    hesspd.tmp=is.positive.definite(model.disp.tmp$hessian, tol=1e-8)
    
    # collect results
    param.tmp=matrix(c(i, j, model.disp.tmp$estimate, se.tmp, model.disp.tmp$minimum, model.disp.tmp$code,hesspd.tmp), nrow=1)
    # write to file
    write.table(param.tmp, file="two_cov_inian1.csv",
                append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
    #a=0
    a.ini=0
    beta_ini=c(log(mean(Y), base = exp(1)), 0, 0)
    gamma_ini=c(log((var(Y)/mean(Y)-1), base = exp(1)), 0, 0)
    initial.v=c(a.ini,beta_ini,gamma_ini)
    
    
    model.disp.tmp=nlm(pt_lik, initial.v, y=Y, X=X, U=U, hessian=T)
    
    se.tmp=SEfromHessian(model.disp.tmp$hessian)
    hesspd.tmp=is.positive.definite(model.disp.tmp$hessian, tol=1e-8)
    
    # collect results
    param.tmp=matrix(c(i, j, model.disp.tmp$estimate, se.tmp, model.disp.tmp$minimum, model.disp.tmp$code,hesspd.tmp), nrow=1)
    # write to file
    write.table(param.tmp, file="two_cov_inia0.csv",
                append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
    
    #a=0.5
    a.ini=0.5
    beta_ini=c(log(mean(Y), base = exp(1)), 0, 0)
    gamma_ini=c(log((var(Y)/mean(Y)-1), base = exp(1)), 0, 0)
    initial.v=c(a.ini,beta_ini,gamma_ini)
    
    
    model.disp.tmp=nlm(pt_lik, initial.v, y=Y, X=X, U=U, hessian=T)
    
    se.tmp=SEfromHessian(model.disp.tmp$hessian)
    hesspd.tmp=is.positive.definite(model.disp.tmp$hessian, tol=1e-8)
    
    # collect results
    param.tmp=matrix(c(i, j, model.disp.tmp$estimate, se.tmp, model.disp.tmp$minimum, model.disp.tmp$code,hesspd.tmp), nrow=1)
    # write to file
    write.table(param.tmp, file="two_cov_inia05.csv",
                append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
  }
}
