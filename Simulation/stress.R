setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


require(HelpersMG)
require(numDeriv)

source("../pt_functions/pt_base.R")
source("../pt_functions/ptreg_fun.R")

samp=c(50,100)
a=c(-1, 0, 0.5)

#specify true coefficient values
beta=matrix(c(log(10), 0.8), ncol=1)
gamma=matrix(c(1, 2), ncol=1)


for(j in a){
  for(k in 1:length(samp)){
  for(i in 1:500){
    #generate X covariates
    set.seed(i)
    x0=rep(1,samp[k])
    x1=rep(0,samp[k])
    index.tmp=sample(1:samp[k], samp[k]*0.02)
    x1[index.tmp]=1
    X=as.matrix(cbind(x0,x1))
    
    #generate U covariates
    u0=rep(1,samp[k])
    u1=x1
    U=as.matrix(cbind(u0, u1))
    
    #generate  Y response
    Y=matrix(rpt_reg_disp(j, beta, gamma, X, U, offsetb = NULL, offsetg = NULL), ncol=1)
    
    
    
    a.ini=max(a.moment(skewness(Y),mean(Y),var(Y)/mean(Y)),-5)
    beta_ini=c(log(mean(Y), base = exp(1)), 0)
    gamma_ini=c(log((var(Y)/mean(Y)-1), base = exp(1)), 0)
    initial.v=c(a.ini,beta_ini,gamma_ini)
    
    model.disp.tmp=optim(initial.v, fn=pt_lik, method="Nelder-Mead", y=Y, X=X, U=U, hessian=F, control = list(maxit = 20000))
    hess.tmp=hessian(pt_lik, model.disp.tmp$par, y=Y, X=X, U=U)
    
    if(any(is.na(hess.tmp))){
    var.tmp=rep("NA", dim(hess.tmp)[1])
    }
    else{
    var.tmp=SEfromHessian(hess.tmp)^2
    }
    # collect results
    param.tmp=matrix(c(i, j, samp[k], model.disp.tmp$par, var.tmp, model.disp.tmp$value, model.disp.tmp$convergence), nrow=1)
    # write to file
    write.table(param.tmp, file="one_cov_stress.csv",
                append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
  }


}

}


