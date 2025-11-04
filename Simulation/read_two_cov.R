setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

params=read.csv("two_cov.csv", sep=",", header=F)

beta=matrix(c(log(10), 0.8, -1), ncol=1)
gamma=matrix(c(1, 2, 1.5), ncol=1)


a.uniq=unique(params[,2])


mean.values=NULL
lengths=NULL
for(j in 1:length(a.uniq)){
  mean.values=rbind(mean.values, apply(params[params[,2]==a.uniq[j],3:16],2,mean))
  lengths=rbind(lengths,dim((params[params[,2]==a.uniq[j],])) )
}

var.valuesemp=NULL
for(j in 1:length(a.uniq)){
  var.valuesemp=rbind(var.valuesemp, apply(params[params[,2]==a.uniq[j],3:9],2,var))
}




#relateive bias

true.vals=cbind(a.uniq,matrix(rep(c(beta,gamma),dim(mean.values)[1]), ncol=length(c(beta,gamma)), byrow = T))

bias.mean=cbind(a.uniq,(mean.values[,1:7]-true.vals)/abs(true.vals)*100)

bias.mean=as.data.frame(bias.mean)

colnames(bias.mean)=c("a","a.bias","beta.0.bias","beta.1.bias","beta.2.bias","gamma.0.bias","gamma.1.bias","gamma.2.bias")





#variance relative difference difference



mean.var=mean.values[,8:14]

rel.bias.var=cbind(a.uniq,(mean.var-var.valuesemp)/var.valuesemp*100)

rel.bias.var=as.data.frame(rel.bias.var)

colnames(rel.bias.var)=c("a","a.bias","beta.0.bias","beta.1.bias","beta.2.bias","gamma.0.bias","gamma.1.bias","gamma.2.bias")






require(kableExtra)

cbind(a.uniq, mean.values[,1:7])%>%
  kable("latex", escape = F, longtable=F,
        caption="",
        booktabs = T, linesep = "", align = "c") %>%
  kable_styling(latex_options =c("HOLD_position", "repeat_header" ))  %>%
  collapse_rows(columns = 1)

cbind(a.uniq,var.valuesemp, mean.var)%>%
  kable("latex", escape = F, longtable=F,
        caption="",
        booktabs = T, linesep = "", align = "c") %>%
  kable_styling(latex_options =c("HOLD_position", "repeat_header", "scale_down" ))  %>%
  collapse_rows(columns = 1)


library("ggplot2")
library("tidyr")
library("ggthemes")
library("paletteer")
library("scales")

# rel.bias.var.short=rel.bias.var[,-2]

pal1=paletteer_c("ggthemes::Red-Gold", 1)
pal2=paletteer_c("ggthemes::Blue", length(beta))
pal3=paletteer_c("ggthemes::Red", length(gamma))
pal=c(pal1,pal2,pal3)





png(filename = "relbiastwocov.png", 
    width = 750,           # Width in pixels
    height = 700,          # Height in pixels
    units = "px",          # Units for width/height (can be "in", "cm", "mm")
    pointsize = 12,        # Default pointsize for text
    bg = "white",          # Background color
    antialias = "cleartype" # Antialiasing for text (e.g., "cleartype", "gray")
)


bias.mean %>%
  pivot_longer( cols = c("a.bias", starts_with("beta"),starts_with("gamma"))) %>%
  ggplot(aes(x=a, y=value, colour=name), na.value = NA) +
  scale_colour_manual(name = "Details", values = pal
                      ,labels = c("a", bquote("\u03B2" [0]),bquote("\u03B2" [1]),bquote("\u03B2" [2]),bquote("\u03B3" [0]),bquote("\u03B3" [1]),bquote("\u03B3" [2]))
  )+
  scale_x_continuous(name="Poisson-Tweedie parameter a")+
  scale_y_continuous(name="Percent", oob = oob_censor_any)+
  geom_point() + theme_light()


# Close the device to save the file
dev.off()



png(filename = "reldiffvartwocov.png", 
    width = 750,           # Width in pixels
    height = 700,          # Height in pixels
    units = "px",          # Units for width/height (can be "in", "cm", "mm")
    pointsize = 12,        # Default pointsize for text
    bg = "white",          # Background color
    antialias = "cleartype" # Antialiasing for text (e.g., "cleartype", "gray")
)


rel.bias.var %>%
  pivot_longer( cols = c("a.bias", starts_with("beta"),starts_with("gamma"))) %>%
  ggplot(aes(x=a, y=value, colour=name)) +
  scale_colour_manual(name = "Details", values = pal
                      ,labels = c("a",bquote("\u03B2" [0]),bquote("\u03B2" [1]),bquote("\u03B2" [2]),bquote("\u03B3" [0]),bquote("\u03B3" [1]),bquote("\u03B3" [2]))
  )+
  scale_x_continuous(name="Poisson-Tweedie parameter a")+
  scale_y_continuous(name="Percent")+
  geom_point()+ theme_light()

# Close the device to save the file
dev.off()




### Coverage ratio

coverage=NULL
for(k in 1:length(a.uniq)){
  param.tmp=params[params[,2]==a.uniq[k],3:16]
  true.val=c(a.uniq[k],beta,gamma)
  coverage.tmp=matrix(0, ncol=7, nrow=dim(param.tmp)[1])
  for(j in 1:7){
    for(i in 1:dim(param.tmp)[1]){
      ci.upper=param.tmp[i,j]+1.96*sqrt(param.tmp[i,j+7])
      ci.lower=param.tmp[i,j]-1.96*sqrt(param.tmp[i,j+7])
      in.coverage=(true.val[j]<ci.upper)*(true.val[j]>ci.lower)
      coverage.tmp[i,j]=in.coverage
    }
  }
  coverage=rbind(coverage,c(a.uniq[k],apply(coverage.tmp,2,mean)))
}

colnames(coverage)=c("a","a.cov","beta.0.cov","beta.1.cov","beta.2.cov","gamma.0.cov","gamma.1.cov","gamma.2.cov")
coverage=as.data.frame(coverage)


pal1=paletteer_c("ggthemes::Red-Gold", 1)
pal2=paletteer_c("ggthemes::Blue", length(beta))
pal3=paletteer_c("ggthemes::Red", length(gamma))
pal=c(pal1,pal2,pal3)

# coverage %>%
#   pivot_longer( cols = c("a.cov", starts_with("beta"),starts_with("gamma"))) %>%
#   ggplot(aes(x=a, y=value, colour=name)) +
#   scale_colour_manual(name = "Details", values = pal
#                       ,labels = c("a",bquote("\u03B2" [0]),bquote("\u03B2" [1]),
#                                   bquote("\u03B2" [2]),bquote("\u03B3" [0]),
#                                   bquote("\u03B3" [1]),bquote("\u03B3" [2]))
#   )+
#   scale_x_continuous(name="Poisson-Tweedie parameter a")+
#   scale_y_continuous(name="Coverage Probability")+
#   geom_point()


png(filename = "covprobtwocov.png", 
    width = 750,           # Width in pixels
    height = 700,          # Height in pixels
    units = "px",          # Units for width/height (can be "in", "cm", "mm")
    pointsize = 12,        # Default pointsize for text
    bg = "white",          # Background color
    antialias = "cleartype" # Antialiasing for text (e.g., "cleartype", "gray")
)

coverage %>%
  pivot_longer( cols = c("a.cov", starts_with("beta"),starts_with("gamma"))) %>%
  ggplot(aes(x=a, y=value, colour=name)) +
  scale_colour_manual(name = "Details", values = pal
                      ,labels = c("a",bquote("\u03B2" [0]),bquote("\u03B2" [1]),bquote("\u03B2" [2]),bquote("\u03B3" [0]),bquote("\u03B3" [1]),bquote("\u03B3" [2]))
  )+
  scale_x_continuous(name="Poisson-Tweedie parameter a")+
  scale_y_continuous(name="Coverage Probability")+
  geom_point() + theme_light()

# Close the device to save the file
dev.off()