rm(list = ls())
library(openxlsx)
library(ggplot2)
library(plyr)
library(psych) # 加载该包可使用corr.test函数
library(reshape2)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(plotly)
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

serum=read.xlsx("20230504-PSP生化指标数据整理_wy.xlsx",sheet = 3)
rownames(serum)=serum$Sample
serum=serum[,-1]


finaldat=cbind(dat_t,serum)
finaldat=as.data.frame(lapply(finaldat,as.numeric))


corr_df <- corr.test(finaldat,method = "spearman", adjust = "BH", alpha = .05)


cordat=corr_df[["r"]]
write.csv(cordat,file = "integrated_analysis_cor.csv")




link_list <- read.xlsx("00Integrated_analysis_link.xlsx",sheet = 1)
node_list <- read.xlsx("00Integrated_analysis_link.xlsx",sheet = 2)

link_list$IDsource <- match(link_list$source, node_list$node) - 1
link_list$IDtarget <- match(link_list$target, node_list$node) - 1

#plotly 包的桑基图


RColorBrewer::brewer.pal(20, "RdBu")

p <- plot_ly(
  type = 'sankey', orientation = 'h',
  textfont=20,
  node = list(
   color = node_list$color,
    pad = 10, thickness = 40,
    line = list(color = 'black', width = 1)
  ),
  
  link = list(
    source = link_list$IDsource, target = link_list$IDtarget,
    value = link_list$weight, color = link_list$color
    
  )
)
p <- p %>% layout(
  font = list(
    size = 25,color="black"
  ),
  xaxis = list(showgrid = F, zeroline = F),
  yaxis = list(showgrid = F, zeroline = F)
)
p




