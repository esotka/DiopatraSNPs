library(vcfR)
rm(list=ls())
vcf = read.vcfR("old/diop_90perc_genome.95perc.recode.vcf.gz")
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
meta <- read_xlsx("old/Revised Diop Meta.xlsx")
colnames(out) = meta$'New Code'[match(colnames(out),meta$plateID)]
out2 = out[,!is.na(colnames(out))] # blank
#out2 = out2[,!colnames(out2)=="CHHS08"] # duplicate
#out2 = out2[,!colnames(out2)=="CHHS18"] # duplicate
write.table(out2,"data/diop_233ind_3162snp.GT.txt",sep="\t",quote=F)

meta2 = meta[match(colnames(out2),meta$'New Code'),c(1,5:10,13:15)]
write.csv(meta2,"data/233ind_meta.csv",row.names=F)

## modify vcf names

vcf2 = vcf
colnames(vcf2@gt)[-1] = substr(colnames(vcf2@gt)[-1],1,3)
vcf3 = vcf2[,c(1,match(meta2$plateID,colnames(vcf2@gt)))]
colnames(vcf3@gt) = c("FORMAT",meta2$'New Code'[match(colnames(vcf3@gt)[-1],meta2$plateID)])
write.vcf(vcf3,"data/diop_233ind_3162snp.vcf.gz")
