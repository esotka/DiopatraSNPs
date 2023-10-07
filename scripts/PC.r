## PC 
rm(list=ls())
library(ggplot2)
library(readxl)
library(RColorBrewer)
dat <- read.delim("data/diop_90perc_genome.95perc.GT_noNAs.txt",sep="\t") 
meta = read_xlsx("data/Revised Diop Meta.xlsx")
colnames(dat) = meta$'New Code'[match(colnames(dat),meta$plateID)]
dat = dat[,!is.na(colnames(dat))]
meta_sub = meta[match(colnames(dat),meta$'New Code'),]
meta_sub$State[substr(meta_sub$'New Code',1,4)=="FlAM"] = "Amelia Island"
meta_sub$State[substr(meta_sub$'New Code',1,4)=="FlFP"] = "Fort Pierce"
meta_sub$State[substr(meta_sub$'New Code',1,4)=="FlSB"] = "St Theresa Beach"

pc <- prcomp(t(dat))
scores <- data.frame(pc$x) #scores for all individuals through all PCs
scores = data.frame(pop = meta_sub$State,scores)
xbar1 = tapply(scores$PC1,scores$pop,mean)
xbar2 = tapply(scores$PC2,scores$pop,mean)
scores = data.frame(scores[,c("pop","PC1","PC2")],
                    xbar1 = xbar1[match(scores$pop,names(xbar1))],
                    xbar2 = xbar2[match(scores$pop,names(xbar2))])
cols <- brewer.pal(n = length(unique(scores$pop)), name = "Set1")
#(((pc$sdev)^2)/sum((pc$sdev)^2))*100

p2 <- ggplot(scores, aes(x=PC1, y=PC2, colour=pop))+
  geom_point(size=2, alpha = 0.5)+
    scale_color_manual(values = cols)+ 
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
    geom_segment(aes(x=PC1, y=PC2, xend=xbar1, yend=xbar2, colour = pop))+
  ggtitle("All Sampling Locations")+
  theme_classic()+
  labs(y= "PC2 (6.86%)", x= "PC1 (16.77%)", colour = "Sampled Populations")


#find sample IDs for outliers
tmp = scores[scores$PC1>10 & scores$PC1<18 & scores$PC2 > -10 & scores$PC2<0,]

pdf("output/PC.pdf")
print(p2)
print(p2 + annotate("text",x=tmp$PC1,y=tmp$PC2,label=rownames(tmp),adj=0))
dev.off()

#