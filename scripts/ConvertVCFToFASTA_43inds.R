### convert vcf of genotype calls from bcftools view (CC, CT, etc..) into single letter IUPAC codes and then make a fasta file
rm(list=ls())
library(vcfR) # read.vcfR()
library(seqinr) # bma() Computing an IUPAC nucleotide symbol
library(ape)
library(readxl)
dat <- read.vcfR("data/diop_231ind_3162snp.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=T)
loci <- rownames(gt2)

#### use subset of samples that have been sequenced at mtDNA (Sotka et al 2023)
mt_meta = read.csv("data/231ind_meta.csv")
library(ape)
dna = read.dna("data/Diopatra_43ind_co1.fas","fasta")
gt2_sub = gt2[,match(labels(dna),colnames(gt2))]

### subset to loci to use

library(parallel)

## by locus.

out <- mclapply(1:nrow(gt2_sub), function(i)
  {
tmp1 <- paste(substr(gt2_sub[i,],1,1),substr(gt2_sub[i,],3,3),sep="")
tmp1 <- ifelse(nchar(tmp1)==1,"N",tmp1)
tmp2 <- mclapply(as.list(tmp1),FUN=function(j) bma(s2c(j)))
})

out2 <- as.data.frame(matrix(unlist(out),nrow=43)) # 43 individuals
out2.noNAs <- out2[,colSums(is.na(out2))==0] ### get rid of all NAs

write.table(x = out2.noNAs,"data/Diopatra_43inds.txt",quote=F,row.names = F,col.names = F,sep = "")

# now pull out individuals with 20% missing data
dna <- read.table("data/Diopatra_43inds.txt")

# now make a fasta file

ind = colnames(gt2_sub)
for (i in 1:length(ind))
{
  tmp <- read.table("data/Diopatra_43inds.txt",skip=i-1,nrows=1)
  if (i == 1)
  {
    write.table(x = paste(">",ind[i],"\n",tmp,"\n",sep=""),"data/Diopatra_43ind.fas",quote=F,row.names = F,col.names = F)
  }
  else
  {
    write.table(x = paste(">",ind[i],"\n",tmp,"\n",sep=""),"data/Diopatra_43ind.fas",quote=F,row.names = F,col.names = F,append = T)
  }
}

