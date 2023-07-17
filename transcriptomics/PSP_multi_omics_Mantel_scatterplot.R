rm(list = ls())
setwd("E:\\Rcode\\11.29_RNA-Seq\\biodeep\\PSP\\meta_rna_cor")
# 将基因和16S数据合并成一个表，然后同时和代谢数据进行ggcor分析
library(tidyverse)
library(ggcor)
#detach("package:linkET");
#library(linkET)
rna=read.csv("contrary_in_HFD_PSP_all(2).csv",header = T)
colnames(rna)
rna=rna[,1:18]
rna_t=as.data.frame(t(rna))
colnames(rna_t)=rna_t[1,]
rna_t=rna_t[-1,]
rna_t=as.data.frame(lapply(rna_t,as.numeric))

meta=read.csv("2Metab_代谢物鉴定定量列表ND-HFD-PSP.csv",header = T)
meta1=read.csv("差异代谢物定量列表ND-HFD.csv",header = T)
meta1=meta1[,c("name","mz","rt")]
meta2=read.csv("差异代谢物定量列表HFD-PSP.csv",header = T)
meta2=meta2[,c("name","mz","rt")]
metaall=rbind(meta1,meta2)
metaall=metaall[duplicated(metaall$name),]

for(i in 1:nrow(metaall)){
  datall_1=meta[which(meta$name==metaall$name[i]),]
  datall_2=meta[which(meta$name==metaall$name[i+1]),]
  meta_final=rbind(datall_1,datall_2,meta_final)
}
meta_final=meta_final[!duplicated(meta_final$name),]
meta_final=meta_final[,c(2,6:22)]
meta_final_t=as.data.frame(t(meta_final))
colnames(meta_final_t)=meta_final_t[1,]
meta_final_t=meta_final_t[-1,]
meta_final_t=as.data.frame(lapply(meta_final_t,as.numeric))


genus=read.csv("2Genus_otu_table_genus_ND_HFD_PSP_tidy.csv",header = T)
genus1=read.csv("Genus_HFD_PSP_P0.05.csv",header = T)
genus1=genus1[,c("X","HFD1")]
genus2=read.csv("Genus_ND_HFD_P0.05.csv",header = T)
genus2=genus2[,c("X","HFD1")]
genusall=rbind(genus1,genus2)
genusall=genusall[duplicated(genusall$X),]
for(i in 1:nrow(genusall)){
  datall_a=genus[which(genus$X==genusall$X[i]),]
  datall_b=genus[which(genus$X==genusall$X[i+1]),]
  genus_final=rbind(datall_a,datall_b,genus_final)
}
genus_final=genus_final[!duplicated(genus_final$X),]
genus_final_t=as.data.frame(t(genus_final))
colnames(genus_final_t)=genus_final_t[1,]
genus_final_t=genus_final_t[-1,]
genus_final_t=as.data.frame(lapply(genus_final_t,as.numeric))
#genus_final_t=genus_final_t[,c(43,2:42)]
##

metabo_16S=cbind(meta_final_t,genus_final_t)
rna_t1=as.data.frame(rna_t[,which(colnames(rna_t)=="Loxl3")])
colnames(rna_t1)="Loxl3"
rna_16S=cbind(rna_t1,genus_final_t)

meta_final_t1=as.data.frame(meta_final_t[,which(colnames(meta_final_t)=="Aminoadipic acid")])
colnames(meta_final_t1)="Aminoadipic.acid"
rna_meta=cbind(rna_t1,meta_final_t1)


mantel <- fortify_mantel(rna_meta,genus_final_t, mantel.fun = 'mantel.randtest',
                         spec.dist.method = 'bray', env.dist.method = 'euclidean',
                         spec.select = list(Loxl3 = 1,
                                            Aminoadipic.acid= 2))



mantel2 <- mutate(mantel,
                  r = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf),
                          labels = c('<0.25', '0.25-0.5', '>=0.5'), right = FALSE),
                  p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                labels = c('<0.001', '0.001-0.01', '0.01-0.05', '>=0.05'), right = FALSE))

corr <- correlate(genus_final_t, cor.test = TRUE) %>% 
  as_cor_tbl(type = "upper", show.diag = FALSE)
df <- combination_layout(mantel2, cor_tbl = corr)


options(ggcor.link.inherit.aes = FALSE)
p22 <- quickcor(corr) +
  geom_square() +
  geom_mark(data = get_data(type = 'upper', show.diag = FALSE), sep = '\n', size = 2, sig.thres = 0.05) +
  geom_link(aes(color = p.value, size = r), data = df) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  geom_start_point(fill = "#d50864", shape = 23, size = 4, data = df) +
  geom_end_point(fill = "dark blue", shape = 16, size = 2.5, data = df) +
  geom_start_label(aes(x = x - 0.5), hjust = 1, size = 8, data = df) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"),
                       limits = c(-1, 1),
                       breaks = seq(-1,1,0.5)) +
  #geom_end_label(aes(x = xend + 0.5), hjust = 0, size = 3.8, data = df) +
  expand_axis(x = -10) +
  #scale_fill_gradient2(midpoint = 0, low = 'blue', mid = 'white', high = 'red', space = 'Lab' ) +
  scale_color_manual(values = c('#4798b3', '#ff8c00', '#999999')) +
  #geom_diag_label(angle = 45) +
  remove_axis('x') +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2,keywidth = 1.5,keyheight = 1.5),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 4), 
                               order = 1,keywidth = 2,keyheight = 2),
         fill = guide_colorbar(title = "Pearson's r", order = 3,keywidth = 1.5,keyheight = 1.5)) +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=18)) +
  theme(legend.key = element_blank())

p22
genus_final_t1=as.data.frame(genus_final_t[,which(colnames(genus_final_t)=="Clostridium_sensu_stricto_1")])
colnames(meta_final_t1)="Clostridium_sensu_stricto_1"



tiff("wy_mantel_Genes_Meta_Genus.tiff",width=8000, height=6000,res=600,compression = 'lzw')
p22
dev.off()









library(ggplot2)
library(ggpubr)
library(ggpmisc)

genus_final_t1=as.data.frame(genus_final_t[,which(colnames(genus_final_t)=="Clostridium_sensu_stricto_1")])
colnames(genus_final_t1)="Clostridium_sensu_stricto_1"

class(Scatterplot_final$Clostridium_sensu_stricto_1)
Scatterplot_final=cbind(rna_meta,genus_final_t1)
Scatterplot_final=as.data.frame(lapply(Scatterplot_final,as.numeric))
tiff("spearman_Scatterplot_gene_meta_.tiff",width=4000, height=3000,res=600,compression = 'lzw')
ggplot(Scatterplot_final,aes(x=Scatterplot_final$Loxl3,y=Scatterplot_final$Aminoadipic.acid))+geom_point(color="#0073C2",size=3.5)+
  geom_smooth(method = "lm", color = "black", fill = "lightgray",se=F,lwd=0.8) +theme_classic()+
  stat_cor(method = "spearman",size=10)+
  xlab("Loxl3")+ylab("Aminoadipic.acid")+
  theme(axis.title.x =element_text(size=25), axis.title.y=element_text(size=25))+
  theme(axis.text.x = element_text(size = 15,color = "black"))+
  theme(axis.text.y = element_text(size = 15,color = "black"))
dev.off()

