setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("dplyr")
library("stringr")
library("ggplot2")
library("Cairo")


normalizedcounts=as.data.frame(readRDS("../data/normalizedcounts.rds"))
X_cov=readRDS("../data/cov_X.rds")


geneidx1=which(rownames(normalizedcounts)=="ENSG00000247317.3") #LY6E-DT
geneidx2=which(str_starts(rownames(normalizedcounts), "ENSG00000237854")=="TRUE") #LINC00674
Y1=matrix(as.numeric(normalizedcounts[geneidx1,]), ncol=1)
Y2=matrix(as.numeric(normalizedcounts[geneidx2,]), ncol=1)

genelye=data.frame(expression=Y1,stageiv=X_cov)
genelin=data.frame(expression=Y2,stageiv=X_cov)

genelye$stage <- factor(genelye$stageiv, 
                        levels = c(0, 1), 
                        labels = c("Stage I-III", "Stage IV"))


genelin$stage <- factor(genelin$stageiv, 
                        levels = c(0, 1), 
                        labels = c("Stage I-III", "Stage IV"))

  
genelye$xformexpression=log(Y1+1, base =  exp(1))
genelin$xformexpression=log(Y2+1, base =  exp(1))


lyemeans=aggregate(genelye[, 1], list(genelye$stage), mean)
linmeans=aggregate(genelin[, 1], list(genelin$stage), mean)


lyemeansiv=as.numeric(genelye$stageiv)*lyemeans[2,2]
lyemeansiii=(-as.numeric(genelye$stageiv)+1)*lyemeans[1,2]
lyemeanssum=lyemeansiv+lyemeansiii

genelye$expressionpointvarmean=(genelye$expression-lyemeanssum)^2/lyemeanssum


linmeansiv=as.numeric(genelin$stageiv)*linmeans[2,2]
linmeansiii=(-as.numeric(genelin$stageiv)+1)*linmeans[1,2]
linmeanssum=linmeansiv+linmeansiii

genelin$expressionpointvarmean=(genelin$expression-linmeanssum)^2/linmeanssum

# png(filename = "lyeboxplot.png", 
#     width = 550,           # Width in pixels
#     height = 500,          # Height in pixels
#     units = "px",          # Units for width/height (can be "in", "cm", "mm")
#     pointsize = 12,        # Default pointsize for text
#     bg = "white",          # Background color
#     antialias = "cleartype" # Antialiasing for text (e.g., "cleartype", "gray")
# )

p=ggplot(genelye, aes(x=stage, y=expression))+  
  stat_boxplot(geom ='errorbar') +
  geom_boxplot()+
  xlab("Cancer Stage")+
  ylab("Expression Level")
ggsave(p,filename="lyeboxplot.pdf", 
       width=5, height=4.5, units="in", device=cairo_pdf)

# # Close the device to save the file
# dev.off()

# png(filename = "linboxplot.png", 
#     width = 550,           # Width in pixels
#     height = 500,          # Height in pixels
#     units = "px",          # Units for width/height (can be "in", "cm", "mm")
#     pointsize = 12,        # Default pointsize for text
#     bg = "white",          # Background color
#     antialias = "cleartype" # Antialiasing for text (e.g., "cleartype", "gray")
# )

p=ggplot(genelin, aes(x=stage, y=expression))+ 
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot()+
  xlab("Cancer Stage")+
  ylab("Expression Level")
ggsave(p,filename="linboxplot.pdf", 
       width=5, height=4.5, units="in", device=cairo_pdf)
# # Close the device to save the file
# dev.off()


# png(filename = "lyeboxplotxform.png", 
#     width = 550,           # Width in pixels
#     height = 500,          # Height in pixels
#     units = "px",          # Units for width/height (can be "in", "cm", "mm")
#     pointsize = 12,        # Default pointsize for text
#     bg = "white",          # Background color
#     antialias = "cleartype" # Antialiasing for text (e.g., "cleartype", "gray")
# )


p=ggplot(genelye, aes(x=stage, y=xformexpression))+  
  stat_boxplot(geom ='errorbar') +
  geom_boxplot()+
  xlab("Cancer Stage")+
  ylab("Expression Level ln(Y+1)")
ggsave(p,filename="lyeboxplotxform.pdf", 
       width=5, height=4.5, units="in", device=cairo_pdf)
# # Close the device to save the file
# dev.off()

# png(filename = "linboxplotxform.png", 
#     width = 550,           # Width in pixels
#     height = 500,          # Height in pixels
#     units = "px",          # Units for width/height (can be "in", "cm", "mm")
#     pointsize = 12,        # Default pointsize for text
#     bg = "white",          # Background color
#     antialias = "cleartype" # Antialiasing for text (e.g., "cleartype", "gray")
# )

p=ggplot(genelin, aes(x=stage, y=xformexpression))+ 
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot()+
  xlab("Cancer Stage")+
  ylab("Expression Level ln(Y+1)")
ggsave(p,filename="linboxplotxform.pdf", 
       width=5, height=4.5, units="in", device=cairo_pdf)

# # Close the device to save the file
# dev.off()



# png(filename = "lyeboxdisp.png", 
#     width = 550,           # Width in pixels
#     height = 500,          # Height in pixels
#     units = "px",          # Units for width/height (can be "in", "cm", "mm")
#     pointsize = 12,        # Default pointsize for text
#     bg = "white",          # Background color
#     antialias = "cleartype" # Antialiasing for text (e.g., "cleartype", "gray")
# )

p=ggplot(genelye, aes(x=stage, y=expressionpointvarmean))+  
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_y_continuous(limits = c(20,500))+
  labs(y=expression(paste("Gene Expression (X- ", bar("X"),")"^2,"/X")),x="Cancer Stage")
ggsave(p,filename="lyeboxdisp.pdf", 
       width=5, height=4.5, units="in", device=cairo_pdf)
# # Close the device to save the file
# dev.off()


# png(filename = "linboxdisp.png", 
#     width = 550,           # Width in pixels
#     height = 500,          # Height in pixels
#     units = "px",          # Units for width/height (can be "in", "cm", "mm")
#     pointsize = 12,        # Default pointsize for text
#     bg = "white",          # Background color
#     antialias = "cleartype" # Antialiasing for text (e.g., "cleartype", "gray")
# )

p=ggplot(genelin, aes(x=stage, y=expressionpointvarmean))+  
  stat_boxplot(geom ='errorbar') +
  geom_boxplot()+
  scale_y_continuous(limits = c(0,500))+
  labs(y=expression(paste("Gene Expression (X- ", bar("X"),")"^2,"/X")),x="Cancer Stage")
ggsave(p,filename="linboxdisp.pdf", 
       width=5, height=4.5, units="in", device=cairo_pdf)
# # Close the device to save the file
# dev.off()

