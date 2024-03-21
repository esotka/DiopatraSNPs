## IBD - only mid-Atlantic (no southern FL, western FL, MA)

rm(list=ls())
iso <- read.csv("data/Isolation_by_distance_final.csv")
iso$pop1 = substr(iso$Pops,1,4)
iso$pop2 = substr(iso$Pops,8,11)
iso$midAtl_TF_1 = !(iso$pop1%in%c("FlFP","FlSB","MaBH","MaWF"))
iso$midAtl_TF_2 = !(iso$pop2%in%c("FlFP","FlSB","MaBH","MaWF"))
iso2 = iso[iso$midAtl_TF_1 & iso$midAtl_TF_2,]

fit <- lm(formula = Fst~ km, data=iso2)
f = ggplot(iso2, aes(x = km, y = Fst)) +
  geom_point() +
  stat_smooth(method = "lm")+
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
  #  label.y = 0.015)+
  theme_classic()+
  labs(x = "Distance (km)")
summary(fit)

pdf("output/IBD.pdf",width=5,height=5)
print(f)
dev.off()

### make matrix
fst_mat = matrix(data=NA,nrow=11,ncol=11,byrow = F)
colnames(fst_mat) = unique(iso2$pop1); rownames(fst_mat) = unique(iso2$pop1)
km_mat = matrix(data=NA,nrow=11,ncol=11,byrow = F)
colnames(km_mat) = unique(iso2$pop1); rownames(km_mat) = unique(iso2$pop1)

for(i in 1:11) # columns
{
    for(j in 1:11) ### rows
    {
        if (i == j | i > j) {fst_mat[i,j]==NA; km_mat[i,j]==NA}
        else {
            fst_mat[j,i] = iso2$Fst[iso2$pop1==colnames(fst_mat)[i] & iso2$pop2==rownames(fst_mat)[j]]
            km_mat[j,i] = iso2$km[iso2$pop1==colnames(km_mat)[i] & iso2$pop2==rownames(km_mat)[j]]
            }
    }}
#diag(fst_mat) = 0; diag(km_mat) = 0
#library(ape); mantel.test(fst_mat, km_mat)
library(ade4); print(mantel.rtest(as.dist(fst_mat),as.dist(km_mat),nrepet=10000))
