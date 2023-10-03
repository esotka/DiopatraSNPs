## PC 
rm(list=ls())
library(ggplot2)
library(readxl)
library(RColorBrewer)
dat = read.delim("data/diop_80perc_genome.95perc.GT_noNAs.txt")
meta <- read_xlsx("data/IndMeta_Final.xlsx")
meta_sub = meta[match(colnames(dat),meta$'New Code'),]
meta_sub$State[substr(meta_sub$'New Code',1,4)=="FlAM"] = "Amelia Island"
meta_sub$State[substr(meta_sub$'New Code',1,4)=="FlFP"] = "Fort Pierce"
meta_sub$State[substr(meta_sub$'New Code',1,4)=="FlSB"] = "St Theresa Beach"

pc <- prcomp(t(dat))
scores <- data.frame(pc$x) #scores for all individuals through all PCs
scores$pop = meta_sub$State
agg <- aggregate(scores[,-ncol(scores)],by= list(cluster=scores$pop), mean)
colnames(agg)<-paste(colnames(agg),"center",sep="_")
colnames(agg)[1] ="pop"
scores <- merge(scores, agg, by = "pop", all = TRUE, sort = FALSE)
rownames(scores) = meta_sub$'New Code'
scores <- scores[, c("pop", "PC1", "PC2", "PC3", "PC4", "PC1_center", "PC2_center", "PC3_center", "PC4_center")]
cols <- brewer.pal(n = length(unique(scores$pop)), name = "Set1")
#(((pc$sdev)^2)/sum((pc$sdev)^2))*100

p2 <- ggplot(scores, aes(x=PC1, y=PC2, colour=pop))+
  geom_point(size=2, alpha = 0.5)+
    scale_color_manual(values = cols)+ 
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
    geom_segment(aes(x=PC1, y=PC2, xend=PC1_center, yend=PC2_center, colour = pop))+
  ggtitle("All Sampling Locations")+
  theme_classic()+
  labs(y= "PC2 (6.86%)", x= "PC1 (16.77%)", colour = "Sampled Populations")


#find sample IDs for outliers
tmp = scores[scores$PC1>10 & scores$PC1<18 & scores$PC2 > -10 & scores$PC2<0,]

pdf("output/PC.pdf")
print(p2 + annotate("text",x=tmp$PC1+1,y=tmp$PC2,label=rownames(tmp),adj=0))
dev.off()

#