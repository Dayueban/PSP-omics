# for mediation analysis
# 参考：https://towardsdatascience.com/doing-and-reporting-your-first-mediation-analysis-in-r-2fe423b92171
rm(list = ls())
library(mediation)


setwd("E:\\Rcode\\11.29_RNA-Seq\\biodeep\\PSP\\meta_rna_cor")
rna=read.csv("contrary_in_HFD_PSP_all(2).csv",header = T)
which(colnames(rna)=="PSP5")
rna=rna[,1:18]
colnames(rna)[1]="name"

meta=read.csv("2Metab_代谢物鉴定定量列表ND-HFD-PSP.csv",header = T)
which(colnames(meta)=="PSP5") #6 22
meta=meta[,c(2,6:22)]


genus=read.csv("2Genus_otu_table_genus_ND_HFD_PSP_tidy.csv",header = T)
colnames(genus)[1]="name"

dat=rbind(rna,meta,genus)
which(dat$name=="Loxl3")  #7
which(dat$name=="Aminoadipic acid")  #161
which(dat$name=="Clostridium_sensu_stricto_1")  #577
dat=dat[c(7,161,577),]
dat_t=as.data.frame(t(dat))
colnames(dat_t)=dat_t[1,]
dat_t=dat_t[-1,]

df=as.data.frame(lapply(dat_t,as.numeric))
df1 <- as.data.frame(scale(df))



fit.totaleffect=lm(Clostridium_sensu_stricto_1~Loxl3,df1)
summary(fit.totaleffect)

# Step #2: The effect of the IV onto the mediator
fit.mediator=lm(Aminoadipic.acid~Loxl3,df1)
summary(fit.mediator)

# Step #3: The effect of the mediator on the dependent variable
fit.dv=lm(Clostridium_sensu_stricto_1~Loxl3+Aminoadipic.acid,df1)
summary(fit.dv)

# Step #4: Causal Mediation Analysis
results = mediate(fit.mediator, fit.dv, treat='Loxl3', 
                  mediator='Aminoadipic.acid', boot=T)
summary(results)





