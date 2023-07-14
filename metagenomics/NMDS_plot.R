rm(list = ls())
# options(stringsAsFactors = FALSE)
library(ade4)
library(MASS)
library(ggplot2)
library(ggpubr)
library(devtools)
library("FactoMineR")
library("factoextra")
library(ggbiplot)
library(vegan)
library(readxl)
##-------------------------------Work--------------------------##
# read table
OTU_df <- read_excel("H_otu_table.xlsx", sheet = 1)
OTU_df <- as.data.frame(OTU_df)

# set 1st column as rownames
rownames(OTU_df) <- OTU_df$Group
OTU_df$Group <- NULL

OTU_df <- sweep(OTU_df, 2, colSums(OTU_df),'/')*100
#OTU_df <- scale(OTU_df, )
#' Standard deviation standardization
stand_stand <- function(x){
  (x - mean(x,na.rm = T)) / sd(x, na.rm = T)
}
#' standardization of the alpha data frame
OTU_df2 <- data.frame(apply(OTU_df, 1, stand_stand))
#OTU_df2 <- as.data.frame(t(OTU_df))
#-------------------------- part one--PCA analysis -------------------------#
#计算bray_curtis距离
distance <- vegdist(OTU_df2, method = 'bray') 

#NMDS排序分析，k = 2预设两个排序轴
nmds <- metaMDS(distance, k = 2)

#获得应力值（stress）
stress <- nmds$stress

#将绘图数据转化为数据框
df <- as.data.frame(nmds$points)
#与分组数据合并
groups <- rep(c("H_FMT", "H_Con", "H_AB"), c(6,6,6))
df <- cbind(df, groups)


## 直接到这一步，直接读取已经计算好的NMDS值
df <- read.csv("../04.BetaDiversity/NMDS/NMDS_scores.csv", header = T, row.names = 1, check.names = F)

df$group <- "NA"
df$group[grep("^ND",rownames(df))] <- "ND"
df$group[grep("^HFD",rownames(df))] <- "HFD"
df$group[grep("^Cur",rownames(df))] <- "Cur"
df$group[grep("^PSP",rownames(df))] <- "PSP"
df$group[grep("^Atro",rownames(df))] <- "Atro"

# 保留ND、HFD和PSP组
df1 <- df[df$group %in% c("ND","HFD","PSP"), ]
# 分组标签按照所需要的顺序
df1$group <- factor(df1$group, levels = c("ND","HFD","PSP"))


#(5)作图2
library(ggExtra)
### 出图
pp <- ggplot(df1, aes(NMDS1, NMDS2))+ 
  geom_point(aes(color = group), size = 2.5, alpha = 0.8)+######设置处理group，点的大小
  scale_shape_manual(aes(values = group)) + 
  #scale_fill_manual(values = c("#00BCD2", "#A60D0D", "#006890","#006A44","#6F27B8"))+ ##底色设置
  scale_fill_manual(values = c("#00BCD2", "#A60D0D", "#006A44"))+ ##底色设置
  #scale_color_manual(values = c("#00BCD2", "#A60D0D", "#006890","#006A44","#6F27B8"))+ ##颜色设置
  scale_color_manual(values = c("#00BCD2", "#A60D0D","#006A44"))+ ##颜色设置
  stat_chull(geom = "polygon", aes(group = group, color = group, fill = group), alpha = 0.1) +
  #annotate("text", x = -0.85, y = 0.5, label = paste0("Pvalue < ", 0.05)) +
  #geom_polygon(data = border, aes(fill = group), alpha = 0.1, show.legend = FALSE) +
  theme_bw()+
  theme(legend.position = "right", legend.direction = "vertical")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color = "black",size =16),
        axis.title = element_text(color = "black", size = 18),
        legend.key.height=unit(0.8,"cm"),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(lineheight=0.8,face="bold",size=18),
        legend.title=element_blank())
pp

ggsave("Plots/NMDS_bc_20230425.png", plot = pp, width = 7, height = 6, dpi = 300)
