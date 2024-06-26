# snmf admix for adam
rm(list=ls())
out = read.delim("data/diop_233ind_3162snp.GT.txt",sep="\t")

library(LEA)
out2<- t(out)
inds <- rownames(out2)
loc <- colnames(out2)
meta = read.csv("data/233ind_meta.csv")
meta <- meta[match(inds,meta$New.Code),] # sort to match the matrix
meta$pops = substr(inds,1,4)

### run once
#write.table(out2,file="data/diop_233ind_3162snp.str",col.names = F, row.names = F,quote = F)
#struct2geno("data/diop_233ind_3162snp.str", ploidy=2, FORMAT = 1)
#project = NULL
#project = snmf("data/diop_233ind_3162snp.str.geno", 
#               K=2:10,
#               entropy = TRUE,
#               repetitions = 10,
#               project = "new") 
project = load.snmfProject("data/diop_233ind_3162snp.str.snmfProject")
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
pop_order2 = pop_order[order(pop_order)]
qmatrix_sort = qmatrix[order(pop_order),]

## make a pretty k=6
pdf("output/snmf_k=6.pdf",width=11,height=6)
barplot(t(qmatrix_sort),border=NA,space=0,xlab="Individuals",ylab="Admixture", xlim= c(0,235),ylim=c(-0.2,1.15),
        col=c("black","purple","green","orange","red","lightgrey"))
        #legend.text= c("1","2","3","4","5","6"),args.legend = list(bty="n",x="right",ncol=1))
endLine <- as.vector(tapply((1:nrow(qmatrix)),pop_order2,max))
segments(x0=endLine,y0=1,x1=endLine,y1=1.1,col="black",lwd=2)
meanPop <- as.vector(tapply((1:length(pop_order2)),pop_order2,mean))
text(levels(pop_order2),x=meanPop,y=1.07,cex=0.5,srt=90)
text(c("FlSB10","NcFR10","CHHS26"),x=(1:nrow(qmatrix))[match(c("FlSB10","NcFR10","CHHS26"),meta2$New.Code)],
        y=-0.07,cex=.3,srt=90)
dev.off()

