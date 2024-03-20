# cophylogeny
# make mtDNA and nucDNA NJ trees and output newick-formatted files
library(phytools)
library(readxl)
rm(list=ls())
tr_nuc = read.tree("data/diop_43ind_iqt.contree")
tr_mt = read.tree("data/diop_43ind_co1.contree")
obj <- cophylo(root(tr_mt,"MaWF01"),root(tr_nuc,"MaWF01"),rotate=T)
pdf("output/cophylogeny_ML.pdf")
plot(obj)
text(x=c(.25,-.25),y=c(1.03,1.03),c("nucDNA","mtDNA"))

par()
dev.off()
