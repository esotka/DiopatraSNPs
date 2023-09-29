library(vcfR)
vcf = read.vcfR("data/diop_80perc_genome.95perc.recode.vcf.gz")
gt2 <- extract.gt(vcf,return.alleles=F)

gt2[gt2=="0/0"] = 0
gt2[gt2=="0/1"] = 1
gt2[gt2=="1/1"] = 2

table(gt2)

### replace NAs with the median genotype 
out <- c()
for (i in 1:nrow(gt2))
{
  tmp <- as.numeric(gt2[i,])
  tmp[is.na(tmp)] = median(tmp,na.rm=T)
  out <- rbind(tmp,out)
}
colnames(out) <- substr(colnames(gt2),1,3)
rownames(out) <- rownames(gt2)

library(readxl)
meta <- read_xlsx("data/IndMeta_Final.xlsx")
colnames(out) = meta$'New Code'[match(colnames(out),meta$plateID)]
out2 = out[,!is.na(colnames(out))] # blank
out2 = out2[,!colnames(out2)=="CHHS08"] # duplicate
out2 = out2[,!colnames(out2)=="CHHS18"] # duplicate
write.table(out2,"data/diop_80perc_genome.95perc.GT_noNAs.txt",sep="\t",quote=F)




