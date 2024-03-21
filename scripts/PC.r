## PC 
rm(list=ls())
library(ggplot2)
library(readxl)
library(RColorBrewer)
dat <- read.delim("data/diop_233ind_3162snp.GT.txt",sep="\t") 
meta_sub = read.csv("data/233ind_meta.csv")

pc <- prcomp(t(dat))
scores <- data.frame(pc$x) #scores for all individuals through all PCs
scores = data.frame(pop = factor(meta_sub$State),scores)
scores$pop = factor(scores$pop,levels=levels(scores$pop)[c(3,4,8,5,6,1,2,7)])
xbar1 = tapply(scores$PC1,scores$pop,mean)
xbar2 = tapply(scores$PC2,scores$pop,mean)
scores = data.frame(scores[,c("pop","PC1","PC2")],
                    xbar1 = xbar1[match(scores$pop,names(xbar1))],
                    xbar2 = xbar2[match(scores$pop,names(xbar2))])
cols <- brewer.pal(n = length(unique(scores$pop)), name = "Set1")
print(head((((pc$sdev)^2)/sum((pc$sdev)^2))*100))

p2 <- ggplot(scores, aes(x=PC1, y=PC2, colour=pop))+
  geom_point(size=2, alpha = 0.5)+
    scale_color_manual(values = cols)+ 
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
    geom_segment(aes(x=PC1, y=PC2, xend=xbar1, yend=xbar2, colour = pop))+
  ggtitle("All Sampling Locations")+
  theme_classic()+
  labs(y= "PC2 (6.8%)", x= "PC1 (16.7%)", colour = "Sampled Populations")


#find sample IDs for outliers
tmp = scores[scores$PC1>10 & scores$PC1<18 & scores$PC2 > -10 & scores$PC2<0,]

pdf("output/PC.pdf")
print(p2)
print(p2 + annotate("text",x=tmp$PC1,y=tmp$PC2,label=rownames(tmp),adj=0))
dev.off()

#