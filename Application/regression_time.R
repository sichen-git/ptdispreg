setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


require(HelpersMG)
library(stringr)
library(kableExtra)
library("Cairo")
require(matrixcalc)

source("../pt_functions/pt_base.R")
source("../pt_functions/ptreg_fun.R")


sizefactors=readRDS("../data/sizefactor.rds")
rawcounts=readRDS("../data/raw_counts.rds")
X_cov=readRDS("../data/cov_X.rds")
Z_cov=readRDS("../data/cov_Z.rds")


offsetx=matrix(sizefactors, ncol=1)


#run the genes each once to remove any overhead
#LY6E-DT gene

geneid1=which(rownames(rawcounts)=="ENSG00000247317.3") #LY6E-DT

Y1=matrix(rawcounts[geneid1,], ncol=1)

Y=Y1
X=cbind(1, X_cov, Z_cov)
U=cbind(1, X_cov)

beta_ini=c(log(mean(Y), base= exp(1)), rep(0, dim(X_cov)[2]),rep(0,  dim(Z_cov)[2]))
gamma_ini=c(log(var(Y)/mean(Y)-1, base= exp(1)), rep(0, dim(X_cov)[2]))

a.ini=max(a.moment(skewness(Y),mean(Y),var(Y)/mean(Y)),-5)
initial.v=c(a.ini,beta_ini,gamma_ini)

model.nlm=nlm(pt_lik, initial.v, y=Y, X=X, U=U, offsetb=offsetx, offsetg=offsetx, hessian=T)


SE.nlm =SEfromHessian(model.nlm$hessian)





#LINC00674 gene 

geneid2=which(str_starts(rownames(rawcounts), "ENSG00000237854")=="TRUE") #LINC00674

Y2=matrix(rawcounts[geneid2,],ncol=1)



Y=Y2
X=cbind(1,X_cov,Z_cov)
U=cbind(1,X_cov)


beta_ini=c(log(mean(Y), base= exp(1)), rep(0, dim(X_cov)[2]),rep(0,  dim(Z_cov)[2]))
gamma_ini=c(log(var(Y)/mean(Y)-1, base= exp(1)), rep(0, dim(X_cov)[2]))

a.ini=max(a.moment(skewness(Y),mean(Y),var(Y)/mean(Y)),-5)
initial.v=c(a.ini,beta_ini,gamma_ini)

model.nlm=nlm(pt_lik, initial.v, y=Y, X=X, U=U, offsetb=offsetx, offsetg=offsetx, hessian=T)


SE.nlm =SEfromHessian(model.nlm$hessian)



Y=Y1
X=cbind(1, X_cov, Z_cov)
U=cbind(1, X_cov)
for(i in 1:20){

  start_time <- Sys.time()
  beta_ini=c(log(mean(Y), base= exp(1)), rep(0, dim(X_cov)[2]),rep(0,  dim(Z_cov)[2]))
  gamma_ini=c(log(var(Y)/mean(Y)-1, base= exp(1)), rep(0, dim(X_cov)[2]))
  
  a.ini=max(a.moment(skewness(Y),mean(Y),var(Y)/mean(Y)),-5)
  initial.v=c(a.ini,beta_ini,gamma_ini)
  
  model.nlm=nlm(pt_lik, initial.v, y=Y, X=X, U=U, offsetb=offsetx, offsetg=offsetx, hessian=T)
  
  
  SE.nlm =SEfromHessian(model.nlm$hessian)
  
  
  end_time <- Sys.time()
  runtime.tmp=end_time - start_time
  hesspd.tmp=is.positive.definite(model.nlm$hessian, tol=1e-8)
  
  
  result.tmp=matrix(c("LY6E-DT", i , runtime.tmp, model.nlm$code, 
                      model.nlm$minimum,   
                      model.nlm$estimate, 
                      model.nlm$estimate[1]<1,
                      hesspd.tmp,
                      SE.nlm), nrow=1)   
  
  write.table(result.tmp, file="reg_time.csv",
              append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
}


Y=Y2
X=cbind(1, X_cov, Z_cov)
U=cbind(1, X_cov)
for(i in 1:20){
  
  start_time <- Sys.time()
  beta_ini=c(log(mean(Y), base= exp(1)), rep(0, dim(X_cov)[2]),rep(0,  dim(Z_cov)[2]))
  gamma_ini=c(log(var(Y)/mean(Y)-1, base= exp(1)), rep(0, dim(X_cov)[2]))
  
  a.ini=max(a.moment(skewness(Y),mean(Y),var(Y)/mean(Y)),-5)
  initial.v=c(a.ini,beta_ini,gamma_ini)
  
  model.nlm=nlm(pt_lik, initial.v, y=Y, X=X, U=U, offsetb=offsetx, offsetg=offsetx, hessian=T)
  
  
  SE.nlm =SEfromHessian(model.nlm$hessian)
  
  
  end_time <- Sys.time()
  runtime.tmp=end_time - start_time
  hesspd.tmp=is.positive.definite(model.nlm$hessian, tol=1e-8)
  
  
  result.tmp=matrix(c("LINC00674", i , runtime.tmp, model.nlm$code, 
                      model.nlm$minimum,   
                      model.nlm$estimate, 
                      model.nlm$estimate[1]<1,
                      hesspd.tmp,
                      SE.nlm), nrow=1)   
  
  write.table(result.tmp, file="reg_time.csv",
              append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
}
