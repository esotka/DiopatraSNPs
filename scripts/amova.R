### AMOVA on adults
### method 1 = pegas::vegan on mpgl numbers.
library(pegas)
library(readxl)
rm(list=ls())
meta = read.csv("data/233ind_meta.csv")
dat <- read.delim('data/diop_233ind_3162snp.GT.txt',header=T)
inds <- colnames(dat)

pop = substr(inds,1,4)
groups = rep("midAtl",length(inds)) # most pops
groups[pop=="FlSB"] = "GulfMex"
groups[pop=="FlFP"] = "SE_FL"
groups[pop%in%c("MaBH","MaDB","MaWF")] = "Mass"


st.d <- dist(t(dat))
print(m1 <- amova(st.d~groups/pop))
write.pegas.amova(m1,"output/amova.csv")
sig2 <- setNames(m1$varcomp$sigma2,rownames(m1$varcomp))
print(getPhi(sig2))

