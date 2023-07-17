rm(list = ls())
library(stringr)
library(qpcR)
library(readxl)
library(UpSetR)

#metabo <- read_xlsx(path = "Con_M1_M2_M3_venn.xlsx", sheet = 1)
#a <- read.csv("Table/venn_neg_NCHC2NCHE.csv", header = T, check.names = F)
#a <- read.csv("Table/venn_pos_NCHC2NCHE.csv", header = T, check.names = F)
# ND和PSP组间差异的
a <- read.csv("ND vs PSP/Serum_metabo_wilcoxn_test_BH_final_updown.csv", header = T, row.names = 1,check.names = F)
# ND和HFD组间差异的
b <- read.csv("ND vs HFD/1-差异分析/差异代谢物定量列表ND-HFD.csv", header = T, row.names = 1,check.names = F)
# HFD和PSP组间差异的
c <- read.csv("HFD vs PSP/1-差异分析/差异代谢物定量列表HFD-PSP.csv", header = T, row.names = 1,check.names = F)

# 下面这几个步骤的目的是为了cbind的时候短的那个能用NA补齐
aa <- rownames(a)
bb <- rownames(b)
cc <- rownames(c)
n <- max(length(aa),length(bb), length(cc))
length(aa) <- n
length(bb) <- n
length(cc) <- n

# 合并a,b,c数据框的行名
abc <- cbind(aa, bb, cc)
abc <- as.data.frame(abc)
names(abc) <- c("ND_PSP", "ND_HFD", "HFD_PSP")


# ND and PSP
a1 <- abc$ND_PSP
a1 <- a1[!is.na(a1)] # 去除空值
# ND and HFD
a2 <- abc$ND_HFD
a2 <- a2[!is.na(a2)] # 去除NA
# HFD and PSP
a3 <- abc$HFD_PSP
a3 <- a3[!is.na(a3)]

#' 查看NC vs HC间差异代谢物相对于NC vs HE有多少独有的
#NCHC <- setdiff(a1,a2)
#write.table(NCHC, "Table/Unique_in_NCHC.txt", quote = FALSE, row.names = FALSE,
#            col.names = FALSE)
#' 查看NC vs HE间差异代谢物相对于NC vs HC有多少独有的
#NCHE <- setdiff(a2, a1)
#write.table(NCHE, "Table/Unique_in_NCHE.txt", quote = FALSE, row.names = FALSE,
#            col.names = FALSE)

## second part: draw venn plot
# load library
library(stringr)
library(qpcR)
library(VennDiagram)

## Application on rap Lyrics
# protein with upregulated
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)
par(mar=c(2,1,2,1)) #sets the bottom, left, top and right margins

#求三个向量的交集
intersect_genus <- Reduce(intersect,  list(v1 = a1,
                                           v2 = a2,
                                           v3 = a3))
intersect_genus2 <- intersect(a2, a3) # 求a2和a3的交集

#[1] "Clostridium_sensu_stricto_1" "[Eubacterium]_siraeum_group" "Roseburia"                   "Erysipelotrichaceae_UCG-003"
#[5] "Blautia"                     "Sellimonas"                  "Bifidobacterium"             "Faecalitalea"               
#[9] "Staphylococcus"

venn.diagram(
  x = list(
    a1, a2, a3
  ),
  category.names = c("ND_PSP" , "ND_HFD", "HFD_PSP"),
  filename = "heatmap_plot/venn_ND_HFD_PSP_metabo.png",
  output = TRUE ,
  imagetype="png" ,
  height = 1200 , 
  width = 1200 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#246b93", "#f7aa5d","#2baeb5"),
  fill = c(alpha("#246b93",0.8), alpha("#f7aa5d",0.8), alpha("#2baeb5",0.8)),
  cex = 0.8,
  fontfamily = "sans",
  fontface = "bold",
  cat.cex = 0.6,
  #cat.default.pos = "inter",
  cat.pos = c(-30, 30, 180),
  margin = 0.1,
  cat.dist = c(0.055, 0.055, 0.055),
  cat.fontfamily = "sans",
  #cat.col = c("#7292AA","#90CBD3","#EAA58E"), # 标签的颜色
  cat.col = c("black", "black", "black"),
  cat.fontface = "bold",
  #rotation = 1 #绘制两个韦恩图的时候该参数不能用，否则报错
  rotation = 1
)

#绘制两个集合
venn.diagram(
  x = list(
   a2, a3
  ),
  category.names = c("ND_HFD", "HFD_PSP"),
  filename = "heatmap_plot/venn_ND_HFD_PSP2_metabo.png",
  output = TRUE ,
  imagetype="png" ,
  height = 1200 , 
  width = 1200 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#f7aa5d","#2baeb5"),
  fill = c(alpha("#f7aa5d",0.8), alpha("#2baeb5",0.8)),
  cex = 0.8,
  fontfamily = "sans",
  fontface = "bold",
  cat.cex = 0.6,
  #cat.default.pos = "inter",
  cat.pos = c(-90, -90),
  margin = 0.1,
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  #cat.col = c("#7292AA","#90CBD3","#EAA58E"), # 标签的颜色
  cat.col = c("black", "black"),
  cat.fontface = "bold",
  rotation.degree = 90,
  #rotation = 1 #绘制两个韦恩图的时候该参数不能用，否则报错
)


# 将有交集的n个genus做boxplot图
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(rstatix)
library(vegan)
library(picante)
library("RColorBrewer")
library(ggsignif)
library(agricolae)
library(dplyr)
library(multcompView)
library(cowplot)
b <- read.csv("ND vs HFD vs PSP/1-差异分析/差异代谢物定量列表.csv", 
              header = T, row.names = 1, check.names = F)
b1 <- b[, 6:ncol(b)]
b_t <- data.frame(t(b1), check.names = F)
b1 <- data.frame(b_t[, intersect_genus2])
names(b1) <- intersect_genus2

genus_label <- names(b1)
# 因为代谢物名称太长或者有特殊符号原因导致涉及计算的过程出错，将代谢物名称替换为简单的
names(b1) <- paste(rep("G", length(intersect_genus2)), 1:length(intersect_genus2), sep = "_")

Start = 1
Stop = ncol(b1)
gvec <- vector("list",length = length(Start:Stop))
group_list <- rep(c("ND", "HFD", "PSP"), c(6,6,5))
b1$group <- group_list
b1$group <- factor(b1$group,
                   levels = c("ND", "HFD", "PSP"))

colors_index <- c("#00BCD2", "#A60D0D","#006A44")
names(colors_index) <- c("ND", "HFD", "PSP")
colors_index <- as.data.frame(colors_index)
colors_index$group <- rownames(colors_index)

for(i in Start:Stop){
  m=names(b1)[i]
  index=b1[,c(m,"group")]
  model = aov(index[[m]] ~ group, data=index)
  # LSD test for stat label
  #out = LSD.test(model,"group", p.adj="BH") # alternative fdr
  out = LSD.test(model,"group", p.adj="none") # alternative fdr
  
  stat = out$groups
  stat$colors = colors_index
  index$stat=stat[as.character(index$group),]$groups
  index$colors=stat[as.character(index$group),]$colors
  max=max(index[,c(m)])
  min=min(index[,c(m)])
  x = index[,c("group",m)]
  y = x %>% group_by(group) %>% summarise_(Max=paste('max(',m,')',sep=""))
  y=as.data.frame(y)
  rownames(y)=y$group
  index$y=y[as.character(index$group),]$Max + (max-min)*0.05
  
  batch_neu <- ggplot(data = index, 
                      aes_string(x="group",y=m, fill="group"))+
    scale_fill_manual(values = c("#00BCD2", "#A60D0D","#006A44"))+
    geom_jitter(size = 1.2)+
    #geom_violin(size=1.5, alpha=.6)+
    geom_boxplot() +
    xlab("Treatment group")+
    ylab("Relitive Abundance")+
    ggtitle(genus_label[i])+
    theme_classic()+
    theme(legend.position = "none")+
    geom_text(data=index, aes(x=group, y=y, color="black", label= stat),size=4) +
    scale_y_continuous(trans = "log10")+
    theme(axis.title.y = element_text(size = 8,face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 7),
          axis.text.x = element_text(size = 8, color = "black", angle = 45,
                                     hjust = 0.95, vjust = 0.95),
          axis.text.y = element_text(size = 10, color = "black"))
  ggsave(batch_neu, file = paste0("Plots/plot_", genus_label[i], ".png"),
         units="in", width=4, height=3, dpi=300)
  gvec[[i-Start+1]] <- batch_neu
}

pp <- plot_grid(plotlist = gvec, nrow = 1) # R package cowplot
ggsave(pp, file = "heatmap_plot/metabo_intersect_ND_HFD_PSP_boxplot4.png", height = 3, width = 10, dpi = 300)






