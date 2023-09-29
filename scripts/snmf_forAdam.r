# snmf admix for adam
rm(list=ls())
out = read.delim("data/diop_80perc_genome.95perc.GT_noNAs.txt")

library(LEA)
library(readxl)
out2<- t(out)
inds <- rownames(out2)
loc <- colnames(out2)
meta <- read_xlsx("data/IndMeta_Final.xlsx")
meta <- meta[match(inds,meta$'New Code'),] # sort to match the matrix
meta$pops = substr(inds,1,4)
write.table(out2,file="data/dat.str",col.names = F, row.names = F,quote = F)
struct2geno("data/dat.str", ploidy=2, FORMAT = 1)

project = NULL
project = snmf("data/dat.str.geno", 
               K=2:10,
               entropy = TRUE,
               repetitions = 10,
               project = "new") 
project = load.snmfProject("data/dat.str.snmfProject")
pdf("output/snmf_cross-val.pdf")
plot(project, col = "blue", pch = 19, cex = 1.2)
dev.off()
best = which.min(cross.entropy(project, K = 6))
qmatrix <- Q(project,K=6,run=best)
pop_order = factor(meta$pops) # by pop
levels(pop_order)[levels(pop_order)=="SHLS"] <- "SCSH"
levels(pop_order)[levels(pop_order)=="SHHS"] <- "SCSH"
levels(pop_order)[levels(pop_order)=="CHLS"] <- "SCCH"
levels(pop_order)[levels(pop_order)=="CHHS"] <- "SCCH"
levels(pop_order)[levels(pop_order)=="BRLS"] <- "SCBR"
levels(pop_order)[levels(pop_order)=="BRHS"] <- "SCBR"
levels(pop_order)[levels(pop_order)=="WBHS"] <- "SCWB"
levels(pop_order)[levels(pop_order)=="BBHS"] <- "SCBB"
levels(pop_order)[levels(pop_order)=="MaWF"] <- "MA"
levels(pop_order)[levels(pop_order)=="MaBH"] <- "MA"
levels(pop_order)[levels(pop_order)=="MaDB"] <- "MA"
levels(pop_order)[levels(pop_order)=="VaWA"] <- "VA"
levels(pop_order)[levels(pop_order)=="NjRG"] <- "NJ"
levels(pop_order)[levels(pop_order)=="NcMV"] <- "NC"
levels(pop_order)[levels(pop_order)=="NcBE"] <- "NC"
levels(pop_order)[levels(pop_order)=="NcFR"] <- "NC"
pop_order = factor(pop_order,levels=levels(pop_order)[c(6,5,4,2,10,3,1,12,8,11,9,7)])
meta2 = meta[order(pop_order),]
meta2$pop_order = pop_order
write.csv(meta2,"output/SNMF_K=6_meta.csv",quote=F)
qmatrix_sort = qmatrix[order(pop_order),]
write.csv(qmatrix_sort,"output/SNMF_K=6_qmatrix.csv", row.names = meta2$'New Code')
barplot(t(qmatrix_sort),border=NA,space=0,xlab="Individuals",ylab="Admixture", xlim= c(0,235),
        legend=TRUE,legend.text= c("1","2","3","4","5","6"),args.legend = list(bty="n",x="right",ncol=1),
        col=c("black","purple","green","orange","red","lightgrey"))
