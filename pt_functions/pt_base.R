# # poisson Tweedie
# # core
# use if no Rcpp package (much slower)
# dPTzero=function(x=0,a=0.5,b=1,c=0.5){
#   if(x<0){
#     return("error: x must be non-negative integer")
#   }
#   else{
#     p=NULL
# 
#     #define p_0
#     p_0=(1-c)^b
#     if(is.na(a) || length(a)==0){stop("a = ",paste(a))}
#     if(a!=0){ p_0=exp(b*((1-c)^a-1)/a);}
# 
#     #define succesive probabilities
#     if(x>0){
#       p[1]=b*c*p_0
#     }
#     if(x>1){
#       #define r
#       r=rep(NULL,x)
#       r[1]=(1-a)*c
#       for(j in 2:x){
#         r[j]=((j-2+a)/j)*c*r[j-1]
#       }
#       #define p_2,...,P_x
#       for(i in 2:x){
#         p[i]=(1/i)*(b*c*p[i-1]+sum( c(1:(i-1)) *r[(i-1):1]*p[1:(i-1)]))
#       }
#     }
#   }
#   return(c(p_0,p))
# }

library(Rcpp)

# Define and compile a C++ function for PT
cppFunction('NumericVector dPTzero(const int & x,
                      const double & a,
                      const double & b,
                      const double & c){
  if (x < 0) stop("Error: x must be non-negative integer");
  
  double pzero= pow((1-c) , b);
  if (a != 0) pzero=exp(b*( pow ((1-c) , a)-1)/a);
  
  NumericVector p (x+1);
  p[0]=pzero;
  
  if (x > 0){
    
    double pinit=b*c*pzero;
    p[1] = pinit;
    
    if (x > 1) {
      
      NumericVector r( x );
      r[0]=(1-a)*c;
      for(int k=1; k<x; k++){
        r[k]=(((double)k+1)-2+a)/((double)k+1)*c*r[k-1];
      }
      
      double temp;
      for( int i=1; i<x; i++ ){
        temp = 0;
        for( int j=0;j<i;j++ ){
          temp = temp + ( (double)j + 1 )*r[i-1-j]*p[j+1];
        }
        p[i+1] = ( 1/( (double)i + 1)*(b*c*p[i] + temp) );
      }
    }
  }
  return p;
}')



# #test
# dPT.0(8, 1,1,0.5) # correction 2
# dPTzero(8, 1,1,0.5)
# dpois(0:8,lambda=1*0.5)
# 
# test.var=c(-1.00000000,0.04165747,0.99834790)
# dPTzero(40, test.var[1],test.var[2],test.var[3])




dPT.abc=function(x=0,a=0.5,b=1,c=0.5){
  p=dPTzero(max(x),a,b,c)
  return(p[x+1])
}



dPT.muphi=function(x=0,a=0.5,mu=1, d=1.5){
  c=ifelse(a==1,1,(d-1)/(d-a))
  b=mu*(1-c)^(1-a)/c

  p=dPTzero(max(x),a,b,c)
  return(p[x+1])
}

#test
# dPT.abc(4:1, 1) # correction 2
# dpois(4:1,lambda=1*0.5)

# dPT.muphi(x=1,a=-1,mu=6.17869e-13, d=5.58187e+13)



pPT.abc=function(q,a=0.5,b=1,c=0.5){
  if(min(q)<0){
    return("error: q must be non-negative integer")
  }
  else{
    p=dPTzero(max(q),a,b,c)
    out=cumsum(p)
    return(out[q+1])
  }
}

pPT.muphi=function(q,a=0.5,mu=1, d=1.5){
  c=ifelse(a==1,1,(d-1)/(d-a))
  b=mu*(1-c)^(1-a)/c
  if(min(q)<0){
    return("error: q must be non-negative integer")
  }
  else{
    p=dPTzero(max(q),a,b,c)
    out=cumsum(p)
    return(out[q+1])
  }
}

#test
pPT.abc(0:5, 1) 
ppois(0:5,lambda=1*0.5)

qPT.abc=function(p,a=0.5,b=1,c=0.5){
  if(sum(p<0 | p>1)>0){
    return("error: p must be a probability")
  }
  else{
    
    mu=b*c/(1-c)^(1-a)
    k.upper=ceiling(mu/(1-max(p)))
    
    p.dens=dPTzero(k.upper,a,b,c)
    
    p.cumu=cumsum(p.dens)
    
    out.rand=NULL
    for(i in 1:length(p)){
      out.rand=c(out.rand,sum(p.cumu<=p[i]))
    }
    
    return(out.rand)
  }
}

qPT.muphi=function(p,a=0.5,mu=1, d=1.5){

  if(sum(p<0 | p>1)>0){
    return("error: p must be a probability")
  }
  else{

    c=ifelse(a==1,1,(d-1)/(d-a))
    b=mu*(1-c)^(1-a)/c
    
    mu=b*c/(1-c)^(1-a)
    k.upper=ceiling(mu/(1-max(p)))
    
    p.dens=dPTzero(k.upper,a,b,c)
    
    p.cumu=cumsum(p.dens)
    
    out.rand=NULL
    for(i in 1:length(p)){
      out.rand=c(out.rand,sum(p.cumu<=p[i]))
    }
    
    return(out.rand)
  }
}


#test
# qPT.abc(c(0.001579507,0.012636055,0.075816332,0.303265330),1)
# qpois(c(0.001579507,0.012636055,0.075816332,0.303265330),lambda=1*0.5)



rPT.abc=function(n,a=0.5,b=1,c=0.5){
  # Test:
  # set.seed(5)
  # n=5
  # a=0
  # b=30/4
  # c=0.8
  
  #generate n random number from Unif(0,1)
  rand.unif=runif(n)
  
  r.max=max(rand.unif)
  
  p0=0
  j=0
  while(r.max>p0){
    p.dens=dPTzero(j,a,b,c)
    p0=sum(p.dens)
    j=j+1
    
  }
  
  p.cumu=cumsum(p.dens)
  
  countcompare=function(x){
    return(sum(p.cumu<=x))
  }
  
  out.rand=sapply(rand.unif,countcompare)
  
  return(out.rand)
}



rPT.muphi=function(n,a=0.5,mu=1, d=1.5){
  c=ifelse(a==1,1,(d-1)/(d-a))
  b=mu*(1-c)^(1-a)/c
  #generate n random number from Unif(0,1)
  rand.unif=runif(n)
  
  r.max=max(rand.unif)

  
  p0=0
  j=0
  while(r.max>p0){
    p.dens=dPTzero(j,a,b,c)
    p0=sum(p.dens)
    j=j+1
    
  }
  
  p.cumu=cumsum(p.dens)
  
  
  countcompare=function(x){
    return(sum(p.cumu<=x))
  }
  
  out.rand=sapply(rand.unif,countcompare)
  
  return(out.rand)
}

# set.seed(5)
# rPT.abc(15,1)
# # [1] 0 1 2 0 0
# set.seed(5)
# rpois(15,lambda=1*0.5)
# # [1] 0 1 2 0 0

a.skew=function(a, mu, d){
  if(d==1){Inf}
  else{
    1/(mu^0.5*d^1.5)*((d-1)^2*(1+1/(1-a))+(3*d-2))
  }
}

a.moment=function(lmb, mu, d){
  if(d==1){return(1)}
  ifelse(lmb<=a.skew(-Inf,mu,d),
         -2,
         1-((lmb*mu^{1/2}*d^{3/2}-(3*d-2))/(d-1)^2-1)^{-1})
}
