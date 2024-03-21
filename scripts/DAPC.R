#DAPC
library(adegenet)
rm(list=ls())
dat = read.delim("data/diop_233ind_3162snp.GT.txt")
meta = read.csv("data/233ind_meta.csv")
#meta <- read_xlsx("data/IndMeta_Final.xlsx")
#meta = meta[match(colnames(dat),meta$'New Code'),]
meta$pops = substr(meta$New.Code,1,4)
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

out_tr = t(dat)
out_tr = out_tr[order(pop_order),]
five <- find.clusters(out_tr, n.clust = 5,n.pca=100, n.da=5)#chose to keep 200 PCs, no 
dapc1 <- dapc(out_tr, five$grp,n.pca=150, n.da=5)
mycols <- c("black","purple","green","orange","red","blue")
pdf("output/dapc.pdf",width=8,height=5)
scatter(dapc1, posi.da = "bottomright", posi.leg = "bottomleft", bg ="white", pch = 20, cell=0, cstar = 0, col= mycols, solid = 0.4, cex=3, clab = 0, leg = TRUE, txt.leg = paste("Cluster", 1:5))
predict_clusters = predict(dapc1,out_tr)
fig = barplot(t(predict_clusters$posterior),col=mycols,border="NA",axisnames=F,ylim=c(0,1.1), space= c(0,0))
xbar = tapply(fig,sort(pop_order2),mean)
mtext(levels(pop_order2),at=xbar,cex=0.5,line=-1)
xmax = tapply(fig,sort(pop_order),max)+.6
segments(xmax,0,xmax,1,lwd=4,col="black")

print(predict_clusters$posterior[predict_clusters$posterior[,5]==1,])
dev.off()
