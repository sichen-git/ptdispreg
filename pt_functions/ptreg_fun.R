
library(moments)
# for initial value of a



# generating random regression variables
rpt_reg_disp <- function(a, beta, gamma, X, U,  offsetb = NULL, offsetg = NULL) {
  eta <- X%*%beta
  iota <- U%*%gamma
  if(!is.null(offsetb)) {eta = eta + log(offsetb, base = exp(1))}
  mu <- exp(eta)
  d <- exp(iota)+1
  if(!is.null(offsetg)) {d = d * offsetg}
  y <-  mapply(rPT.muphi, n = 1, a=a, mu=mu, d=d)
  return(y)
}




# pt likelihood
pt_lik <- function(theta, y, X, U,  offsetb = NULL, offsetg = NULL) {
  a=theta[1]
  if(a>1){Inf}
  else{
    lb=dim(X)[2]
    lg=dim(U)[2]
    beta=matrix(theta[2:(lb+1)],ncol=1)
    gamma=matrix(theta[(lb+2):(lb+lg+1)],ncol=1)
    eta <- X%*%beta
    iota <- U%*%gamma
    if(!is.null(offsetb)) {eta = eta + log(offsetb, base = exp(1))}
    mu <- exp(eta)
    d <- exp(iota)+1
    if(!is.null(offsetg)) {d = d * offsetg}
    return( -sum(log(mapply(dPT.muphi,x=y, a=a, mu=mu, d=d)) ))
  }
}


# for pt calculating regression CMF
pPTreg <- function(theta, y, X, U,  offsetb = NULL, offsetg = NULL) {
  a=theta[1]
  lb=dim(X)[2]
  lg=dim(U)[2]
  beta=matrix(theta[2:(lb+1)],ncol=1)
  gamma=matrix(theta[(lb+2):(lb+lg+1)],ncol=1)
  eta <- X%*%beta
  iota <- U%*%gamma
  if(!is.null(offsetb)) {eta = eta + log(offsetb, base = exp(1))}
  mu <- exp(eta)
  d <- exp(iota)+1
  if(!is.null(offsetg)) {d = d * offsetg}
  return( (mapply(pPT.muphi,q=y, a=a, mu=mu, d=d)))
}



# nb likelihood
nb_lik <- function(theta, y, X, U,  offsetb = NULL, offsetg = NULL) {
    lb=dim(X)[2]
    lg=dim(U)[2]
    beta=matrix(theta[1:lb],ncol=1)
    gamma=matrix(theta[(lb+1):(lb+lg)],ncol=1)
    eta <- X%*%beta
    iota <- U%*%gamma
    if(!is.null(offsetb)) {eta = eta + log(offsetb, base = exp(1))}
    if(!is.null(offsetg)) {iota = iota + log(offsetg, base = exp(1))}
    mu <- exp(eta)
    d <- exp(iota)
    p=1/d
    size=mu*p/(1-p)
    return( -sum(log(mapply(dnbinom,x=y, size=size, mu  = mu)) ))
}


# for nb calculating regression CMF
pnbreg <- function(theta, y, X, U,  offsetb = NULL, offsetg = NULL) {
  lb=dim(X)[2]
  lg=dim(U)[2]
  beta=matrix(theta[1:lb],ncol=1)
  gamma=matrix(theta[(lb+1):(lb+lg)],ncol=1)
  eta <- X%*%beta
  iota <- U%*%gamma
  if(!is.null(offsetb)) {eta = eta + log(offsetb, base = exp(1))}
  if(!is.null(offsetg)) {iota = iota + log(offsetg, base = exp(1))}
  mu <- exp(eta)
  d <- exp(iota)
  p=1/d
  size=mu*p/(1-p)
  return( (pnbinom(q=y, size=size, mu  = mu)))
}
