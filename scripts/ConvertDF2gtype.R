### make gtype file

### LD and Ne
#library(strataG)
#library(ggplot2)
library(readxl)
rm(list=ls())
snp <- read.delim("data/diop_233ind_3162snp.GT.txt",sep="\t") 
#meta = read_xlsx("./Revised Diop Meta.xlsx")
#colnames(snp) = meta$'New Code'[match(colnames(snp),meta$plateID)]
#snp = snp[,!is.na(colnames(snp))]
#snp = snp[,!colnames(snp)%in%c("NA.","CHHS18","CHHS08")]

snp <- t(snp)# need col = loci; row = ind
### convert into 2 alleles per locus
nloci <- ncol(snp)
############################
### do this code once ######
snp2 <- c()
for (i in 1:nloci)
{
  tmp <- snp[,i]
  key <- data.frame(snp=c(0,1,2),a=c(1,1,2),b=c(1,2,2))
  a <- key$a[match(tmp,key$snp)]
  b <- key$b[match(tmp,key$snp)]
  snp2 <- cbind(snp2,a,b)
}
colnames(snp2) <- sort(c(paste(1:nloci,"a",sep=""),paste(1:nloci,"b",sep="")))
pop <- substr(rownames(snp),1,4)
dat <- data.frame(rownames(snp),pop,snp2)
colnames(dat) <- c("Ind","Pop",colnames(snp2))
gi <- df2gtypes(dat, ploidy = 2, id.col = 1, strata.col = 2, loc.col = 3)
saveRDS(gi,file="diop_90perc_genome.95perc.GT_noNAs_gtype") ### gtype-formatted data

pws = pairwiseTest(gi)
pairwiseMatrix(pws,"wcFst")
pairwiseSummary(pws,"wcFst")
