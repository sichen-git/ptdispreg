setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


require(HelpersMG)
library(stringr)
library(kableExtra)

source("../Poisson-Tweedie/pt_base.R")
source("../Poisson-Tweedie/ptreg_fun.R")


sizefactors=readRDS("../data/sizefactor.rds")
rawcounts=readRDS("../data/raw_counts.rds")
X_cov=readRDS("../data/cov_X.rds")
Z_cov=readRDS("../data/cov_Z.rds")


offsetx=matrix(sizefactors, ncol=1)


#LY6E-DT gene

geneid1=which(rownames(rawcounts)=="ENSG00000247317.3") #LY6E-DT

Y1=matrix(rawcounts[geneid1,], ncol=1)



Y=Y1
X=cbind(1, X_cov, Z_cov)
U=cbind(1, X_cov)

# beta_ini=glm.fit(X,Y, family = poisson())$coefficients
# gamma_ini=glm.fit(U,Y, family = poisson())$coefficients
beta_ini=c(log(mean(Y), base= exp(1)), rep(0, dim(X_cov)[2]),rep(0,  dim(Z_cov)[2]))
gamma_ini=c(log(var(Y)/mean(Y)-1, base= exp(1)), rep(0, dim(X_cov)[2]))

a.ini=max(a.moment(skewness(Y),mean(Y),var(Y)/mean(Y)),-5)
initial.v=c(a.ini,beta_ini,gamma_ini)

model.nlm=nlm(pt_lik, initial.v, y=Y, X=X, U=U, offsetb=offsetx, offsetg=offsetx, hessian=T)


SE.nlm =SEfromHessian(model.nlm$hessian)


tval.nlm=model.nlm$estimate/SE.nlm
modelpval.nlm=2*pt(-abs(tval.nlm),df=length(Y)-1)

mytable.nlm=cbind(model.nlm$estimate,SE.nlm,tval.nlm,modelpval.nlm)

rownames(mytable.nlm)=c("a","b0","b1","b2","b3","g0","g1")
colnames(mytable.nlm)=c("MLE","SE","t-val","p-val")

mytable.nlm%>%
  kbl( format = "latex",booktabs = T) %>%
  kable_styling(latex_options="scale_down")


qresidptlye=pPTreg(model.nlm$estimate,Y, X, U, offsetb=offsetx, offsetg=offsetx)
qresidptlyeb=pPTreg(model.nlm$estimate,Y+1, X, U, offsetb=offsetx, offsetg=offsetx)
qqnorm(qnorm(qresidptlye),main=expression(paste(italic("LY6E-DT")," / PT Fit")));abline(0,1,col=2)

png(filename = "qqlyemeanoffsetboth.png", 
    width = 550,           # Width in pixels
    height = 500,          # Height in pixels
    units = "px",          # Units for width/height (can be "in", "cm", "mm")
    pointsize = 12,        # Default pointsize for text
    bg = "white",          # Background color
    antialias = "cleartype" # Antialiasing for text (e.g., "cleartype", "gray")
)

# Create your plot
qqnorm(qnorm(qresidptlye),main=expression(paste(italic("LY6E-DT")," / PT Fit")));abline(0,1,col=2)

# Close the device to save the file
dev.off()


#NB


# beta_ini=glm.fit(X,Y, family = poisson())$coefficients
# gamma_ini=glm.fit(U,Y, family = poisson())$coefficients
beta_ini=c(log(mean(Y), base= exp(1)), rep(0, dim(X_cov)[2]),rep(0,  dim(Z_cov)[2]))
gamma_ini=c(log(var(Y)/mean(Y)-1, base= exp(1)), rep(0, dim(X_cov)[2]))


initial.v=c(beta_ini,gamma_ini)

model.nlm.nb=nlm(nb_lik, initial.v, y=Y, X=X, U=U, offsetb=offsetx, offsetg=offsetx, hessian=T)


SE.nlm.nb =SEfromHessian(model.nlm.nb$hessian)


tval.nlm.nb=model.nlm.nb$estimate/SE.nlm.nb
modelpval.nlm.nb=2*pt(-abs(tval.nlm.nb),df=length(Y)-1)

mytable.nlm.nb=cbind(model.nlm.nb$estimate,SE.nlm.nb,tval.nlm.nb,modelpval.nlm.nb)

rownames(mytable.nlm.nb)=c("b0","b1","b2","b3","g0","g1")
colnames(mytable.nlm.nb)=c("MLE","SE","t-val","p-val")


mytable.nlm.nb%>%
  kbl( format = "latex",booktabs = T) %>%
  kable_styling(latex_options="scale_down")



qresidptlyenb=pnbreg(model.nlm.nb$estimate,Y, X, U, offsetb=offsetx, offsetg=offsetx)
qresidptlyenbb=pnbreg(model.nlm.nb$estimate,Y+1, X, U, offsetb=offsetx, offsetg=offsetx)
qqnorm(qnorm(qresidptlyenb),main=expression(paste(italic("LY6E-DT")," / NB Fit")));abline(0,1,col=2)


png(filename = "qqlyemeanoffsetnbboth.png", 
    width = 550,           # Width in pixels
    height = 500,          # Height in pixels
    units = "px",          # Units for width/height (can be "in", "cm", "mm")
    pointsize = 12,        # Default pointsize for text
    bg = "white",          # Background color
    antialias = "cleartype" # Antialiasing for text (e.g., "cleartype", "gray")
)

# Create your plot
qqnorm(qnorm(qresidptlyenb),main=expression(paste(italic("LY6E-DT")," / NB Fit")));abline(0,1,col=2)

# Close the device to save the file
dev.off()







#LINC00674 gene 

geneid2=which(str_starts(rownames(rawcounts), "ENSG00000237854")=="TRUE") #LINC00674

Y2=matrix(rawcounts[geneid2,],ncol=1)



Y=Y2
X=cbind(1,X_cov,Z_cov)
U=cbind(1,X_cov)

# beta_ini=glm.fit(X,Y, family = poisson())$coefficients
# gamma_ini=glm.fit(U,Y, family = poisson())$coefficients
beta_ini=c(log(mean(Y), base= exp(1)), rep(0, dim(X_cov)[2]),rep(0,  dim(Z_cov)[2]))
gamma_ini=c(log(var(Y)/mean(Y)-1, base= exp(1)), rep(0, dim(X_cov)[2]))

a.ini=max(a.moment(skewness(Y),mean(Y),var(Y)/mean(Y)),-5)
initial.v=c(a.ini,beta_ini,gamma_ini)

model.nlm=nlm(pt_lik, initial.v, y=Y, X=X, U=U, offsetb=offsetx, offsetg=offsetx, hessian=T)


SE.nlm =SEfromHessian(model.nlm$hessian)


tval.nlm=model.nlm$estimate/SE.nlm
modelpval.nlm=2*pt(-abs(tval.nlm),df=length(Y)-1)

mytable.nlm=cbind(model.nlm$estimate,SE.nlm,tval.nlm,modelpval.nlm)

rownames(mytable.nlm)=c("a","b0","b1","b2","b3","g0","g1")
colnames(mytable.nlm)=c("MLE","SE","t-val","p-val")

mytable.nlm%>%
  kbl( format = "latex",booktabs = T) %>%
  kable_styling(latex_options="scale_down")


qresidptlin=pPTreg(model.nlm$estimate,Y, X, U, offsetb=offsetx, offsetg=offsetx)
qresidptlinb=pPTreg(model.nlm$estimate,Y+1, X, U, offsetb=offsetx, offsetg=offsetx)
qqnorm(qnorm(qresidptlin),main=expression(paste(italic("LINC00674")," / PT Fit")));abline(0,1,col=2)


png(filename = "qqlinmeanoffsetboth.png", 
    width = 550,           # Width in pixels
    height = 500,          # Height in pixels
    units = "px",          # Units for width/height (can be "in", "cm", "mm")
    pointsize = 12,        # Default pointsize for text
    bg = "white",          # Background color
    antialias = "cleartype" # Antialiasing for text (e.g., "cleartype", "gray")
)

# Create your plot
qqnorm(qnorm(qresidptlin),main=expression(paste(italic("LINC00674")," / PT Fit")));abline(0,1,col=2)

# Close the device to save the file
dev.off()





#NB


# beta_ini=glm.fit(X,Y, family = poisson())$coefficients
# gamma_ini=glm.fit(U,Y, family = poisson())$coefficients
beta_ini=c(log(mean(Y), base= exp(1)), rep(0, dim(X_cov)[2]),rep(0,  dim(Z_cov)[2]))
gamma_ini=c(log(var(Y)/mean(Y)-1, base= exp(1)), rep(0, dim(X_cov)[2]))


initial.v=c(beta_ini,gamma_ini)

model.nlm.nb=nlm(nb_lik, initial.v, y=Y, X=X, U=U, offsetb=offsetx, offsetg=offsetx, hessian=T)


SE.nlm.nb =SEfromHessian(model.nlm.nb$hessian)


tval.nlm.nb=model.nlm.nb$estimate/SE.nlm.nb
modelpval.nlm.nb=2*pt(-abs(tval.nlm.nb),df=length(Y)-1)

mytable.nlm.nb=cbind(model.nlm.nb$estimate,SE.nlm.nb,tval.nlm.nb,modelpval.nlm.nb)

rownames(mytable.nlm.nb)=c("b0","b1","b2","b3","g0","g1")
colnames(mytable.nlm.nb)=c("MLE","SE","t-val","p-val")


mytable.nlm.nb%>%
  kbl( format = "latex",booktabs = T) %>%
  kable_styling(latex_options="scale_down")


qresidptlinnb=pnbreg(model.nlm.nb$estimate,Y, X, U, offsetb=offsetx, offsetg=offsetx)
qresidptlinnbb=pnbreg(model.nlm.nb$estimate,Y+1, X, U, offsetb=offsetx, offsetg=offsetx)
normresid=qnorm(qresidptlinnb)
normresid[is.infinite(normresid)] <- NA

qqnorm(normresid,main=expression(paste(italic("LINC00674")," / NB Fit")));abline(0,1,col=2)


png(filename = "qqlinmeanoffsetnbboth.png", 
    width = 550,           # Width in pixels
    height = 500,          # Height in pixels
    units = "px",          # Units for width/height (can be "in", "cm", "mm")
    pointsize = 12,        # Default pointsize for text
    bg = "white",          # Background color
    antialias = "cleartype" # Antialiasing for text (e.g., "cleartype", "gray")
)

# Create your plot
qqnorm(normresid,main=expression(paste(italic("LINC00674")," / NB Fit")));abline(0,1,col=2)

# Close the device to save the file
dev.off()



set.seed(123)

randlyeres=runif(length(qresidptlye),min=qresidptlye,max=qresidptlyeb)
randlyeresnb=runif(length(qresidptlyenb),min=qresidptlyenb,max=qresidptlyenbb)
randlinres=runif(length(qresidptlin),min=qresidptlin,max=qresidptlinb)
randlinresnb=runif(length(qresidptlinnb),min=qresidptlinnb,max=qresidptlinnbb)


ks.test(qnorm(randlyeres), "pnorm")
ks.test(qnorm(randlyeresnb), "pnorm")

ks.test(qnorm(randlinres), "pnorm")
ks.test(qnorm(randlinresnb), "pnorm")



