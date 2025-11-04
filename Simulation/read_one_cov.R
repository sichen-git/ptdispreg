setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

params=read.csv("one_cov.csv", sep=",", header=F)


beta=matrix(c(log(10), 0.8), ncol=1)
gamma=matrix(c(1, 2), ncol=1)

a.uniq=unique(params[,2])

mean.values=NULL
lengths=NULL
for(j in 1:length(a.uniq)){
  mean.values=rbind(mean.values, apply(params[params[,2]==a.uniq[j],3:12],2,mean))
  lengths=rbind(lengths,dim((params[params[,2]==a.uniq[j],])) )
}

var.valuesemp=NULL
for(j in 1:length(a.uniq)){
  var.valuesemp=rbind(var.valuesemp, apply(params[params[,2]==a.uniq[j],3:7],2,var))
}




#relateive bias

true.vals=cbind(a.uniq,matrix(rep(c(beta,gamma),dim(mean.values)[1]), ncol=length(c(beta,gamma)), byrow = T))

bias.mean=cbind(a.uniq,(mean.values[,1:5]-true.vals)/abs(true.vals)*100)

bias.mean=as.data.frame(bias.mean)

colnames(bias.mean)=c("a","a.bias","beta.0.bias","beta.1.bias","gamma.0.bias","gamma.1.bias")


#variance relative difference difference

mean.var=mean.values[,6:10]

rel.bias.var=cbind(a.uniq,(mean.var-var.valuesemp)/var.valuesemp*100)

rel.bias.var=as.data.frame(rel.bias.var)

colnames(rel.bias.var)=c("a","a.bias","beta.0.bias","beta.1.bias","gamma.0.bias","gamma.1.bias")


require(kableExtra)

cbind(a.uniq, mean.values[,1:5])%>%
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




png(filename = "relbiasonecov.png", 
    width = 750,           # Width in pixels
    height = 700,          # Height in pixels
    units = "px",          # Units for width/height (can be "in", "cm", "mm")
    pointsize = 12,        # Default pointsize for text
    bg = "white",          # Background color
    antialias = "cleartype" # Antialiasing for text (e.g., "cleartype", "gray")
)


bias.mean %>%
  pivot_longer( cols = c("a.bias", starts_with("beta"),starts_with("gamma"))) %>%
  ggplot(aes(x=a, y=value, colour=name)) +
  scale_colour_manual(name = "Details", values = pal
                      ,labels = c("a",bquote("\u03B2" [0]),bquote("\u03B2" [1]),bquote("\u03B3" [0]),bquote("\u03B3" [1]))
  )+
  scale_x_continuous(name="Poisson-Tweedie parameter a")+
  scale_y_continuous(name="Percent", oob = oob_censor_any)+
  geom_point()+ theme_light()


# Close the device to save the file
dev.off()



png(filename = "reldiffvaronecov.png", 
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
                      ,labels = c("a",bquote("\u03B2" [0]),bquote("\u03B2" [1]),bquote("\u03B3" [0]),bquote("\u03B3" [1]))
  )+
  scale_x_continuous(name="Poisson-Tweedie parameter a")+
  scale_y_continuous(name="Percent")+
  geom_point()+ theme_light()

# Close the device to save the file
dev.off()


### Coverage ratio

coverage=NULL
for(k in 1:length(a.uniq)){
  param.tmp=params[params[,2]==a.uniq[k],3:12]
  true.val=c(a.uniq[k],beta,gamma)
  coverage.tmp=matrix(0, ncol=5, nrow=dim(param.tmp)[1])
  for(j in 1:5){
    for(i in 1:dim(param.tmp)[1]){
      ci.upper=param.tmp[i,j]+1.96*sqrt(param.tmp[i,j+5])
      ci.lower=param.tmp[i,j]-1.96*sqrt(param.tmp[i,j+5])
      in.coverage=(true.val[j]<ci.upper)*(true.val[j]>ci.lower)
      coverage.tmp[i,j]=in.coverage
    }
  }
  coverage=rbind(coverage,c(a.uniq[k],apply(coverage.tmp,2,mean)))
}

colnames(coverage)=c("a","a.cov","beta.0.cov","beta.1.cov","gamma.0.cov","gamma.1.cov")
coverage=as.data.frame(coverage)


pal1=paletteer_c("ggthemes::Red-Gold", 1)
pal2=paletteer_c("ggthemes::Blue", length(beta))
pal3=paletteer_c("ggthemes::Red", length(gamma))
pal=c(pal1,pal2,pal3)

# coverage %>%
#   pivot_longer( cols = c("a.cov", starts_with("beta"),starts_with("gamma"))) %>%
#   ggplot(aes(x=a, y=value, colour=name)) +
#   scale_colour_manual(name = "Details", values = pal
#                       ,labels = c("a",bquote("\u03B2" [0]),bquote("\u03B2" [1]),bquote("\u03B3" [0]),bquote("\u03B3" [1]))
#   )+
#   scale_x_continuous(name="Poisson-Tweedie parameter a")+
#   scale_y_continuous(name="Coverage Probability")+
#   geom_point()


png(filename = "covprobonecov.png", 
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
                      ,labels = c("a",bquote("\u03B2" [0]),bquote("\u03B2" [1]),bquote("\u03B3" [0]),bquote("\u03B3" [1]))
  )+
  scale_x_continuous(name="Poisson-Tweedie parameter a")+
  scale_y_continuous(name="Coverage Probability")+
  geom_point()+ theme_light()

# Close the device to save the file
dev.off()
