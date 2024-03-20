### phenotype link to loci
library(readxl)
library(lattice)
rm(list=ls())
meta = read_xlsx("data/Revised Diop Meta.xlsx")
dat <- read.delim('data/diop_90perc_genome.95perc.GT_noNAs.txt',header=T)
colnames(dat) = meta$'New Code'[match(colnames(dat),meta$plateID)]
dat = dat[,!is.na(colnames(dat))]
inds <- colnames(dat)

pheno = read_xlsx("data/SC Phenotypes.xlsx")
print(table(pheno$'New Code'%in%inds))
#FALSE  TRUE 
#   36   154
dat_sub = dat[,inds%in%pheno$'New Code']
cols = pheno$Phenotype[match(colnames(dat_sub),pheno$'New Code')]
stats = c()
for (i in 1:nrow(dat_sub))
{
    tmp = table(as.numeric(dat_sub[i,]),cols)
    p = chisq.test(tmp)$p.value
    df = chisq.test(tmp)$parameter[[1]]
    stats = rbind(stats,data.frame(snp=rownames(dat_sub)[i],df,p))
    }

print(histogram(~p | as.character(df),stats,breaks=100))

outlier = stats[stats$p<0.0001 & stats$df == 4,]
print(table(as.numeric(dat_sub[as.numeric(rownames(outlier))[1],]),cols))
print(table(as.numeric(dat_sub[as.numeric(rownames(outlier))[2],]),cols))
print(table(as.numeric(dat_sub[as.numeric(rownames(outlier))[3],]),cols))
