rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(scales)
library(patchwork)
library(cowplot)

a <- read.csv("taxa_abundance/Relative_group/otu_table_group.p.relative.csv", header = T, check.names = F)
# 计算相对丰度
#a1 <- sweep(a[, 2:ncol(a)], 2, colSums(a[, 2:ncol(a)]),'/')*100
#a1 <- cbind(a$Phylum, a1)
a1 <- a[,2:ncol(a)]*100
a1 <- cbind(a$Taxonomy, a1)
names(a1)[1] <- "Phylum"
# 将数据进行长宽矩阵转换

a_t <- melt(a1)

names(a_t) <- c("Phylum", "Group", "Abundance")
a_t$Group <- factor(a_t$Group, levels = c("ND", "HFD", "Cur", "PSP", "Atro"))
a_t$Phylum <- factor(a_t$Phylum, levels = c("Others", "Actinobacteriota","Chloroflexi","Acidobacteriota",
                                            "Desulfobacterota","Actinobacteria","Proteobacteria","unidentified_Bacteria",
                                            "Verrucomicrobiota","Bacteroidota","Firmicutes"))
# ND, HFD and PSP were showed first
a_t <- a_t[a_t$Group %in% c("ND","HFD","PSP"),]
# plot
p <- ggplot(a_t, aes(x = Group, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity",
           position = "stack") +
  theme_classic() + 
  theme(axis.text.y = element_text(size = 14,colour = "black"),
        axis.text.x = element_text(size = 12,angle = 60, hjust = 0.95, vjust = 0.95, colour = "black"),
        axis.title = element_text(size = 16, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 16, color = "black", face = "bold"),
        legend.key.height=unit(0.5,"cm"),
        legend.key.width=unit(0.5,"cm"),
        axis.ticks=element_line(size=1.2),
        panel.background = element_rect(colour = "black", size = 1)) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  labs(title = "Phylum") +
  #scale_y_continuous(labels = scales::percent, expand = c(0, 0)) + # 以百分数进行Y轴的展示
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("#b0b9b9","#ffe327","#2f4e87","#f69001","#f0eedf",
                               "#aed4e9","#f4a69a","#3ba889","#4593c3","#c5942e",
                               "#C696C4","#006b31","#fddd00",
                               "#d80b13")) 
p  
ggsave("Plots/Plylum_percent_barplot2.png", plot = p, width = 5, height = 6, dpi = 300)  


## for Genus
b <- read.csv("taxa_abundance/Relative_group/otu_table_group.g.relative.csv", header = T, check.names = F)
b1 <- b[,2:ncol(b)]*100
b1 <- cbind(b$Taxonomy, b1)
names(b1)[1] <- "Genus"
# 将数据进行长宽矩阵转换

b_t <- melt(b)

names(b_t) <- c("Genus", "Group", "Abundance")
b_t$Group <- factor(b_t$Group, levels = c("ND", "HFD", "Cur", "PSP", "Atro"))
b_t$Genus <- factor(b_t$Genus, levels = c("Others","Ligilactobacillus","Marvinbryantia","Roseburia","Colidextribacter",
                                            "Turicibacter","Clostridium_sensu_stricto_1","[Eubacterium]_xylanophilum_group",
                                            "[Eubacterium]_ruminantium_group","Limosilactobacillus",
                                            "Parabacteroides","Lactobacillus","Bacteroides","Lachnospiraceae_NK4A136_group",
                                            "Akkermansia",
                                            "Romboutsia"))
# ND, HFD and PSP were showed first
b_t <- b_t[b_t$Group %in% c("ND","HFD","PSP"),]
# plot
p1 <- ggplot(b_t, aes(x = Group, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity",
           position = "stack") +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 12,angle = 60, hjust = 0.95, vjust = 0.95, colour = "black"),
        axis.text.y = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 16, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 16, color = "black", face = "bold"),
        legend.key.height=unit(0.5,"cm"),
        legend.key.width=unit(0.5,"cm"),
        axis.ticks=element_line(size=1.2),
        panel.background = element_rect(colour = "black", size = 1)) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  labs(title = "Genus") + 
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) + # 以百分数进行Y轴的展示
  
  scale_fill_manual(values = c("#b0b9b9","#ffe327","#2f4e87","#f69001","#f0eedf",
                               "#aed4e9","#f4a69a","#3ba889","#4593c3","#f18e0c",
                               "#262a35","#c5942e","#C696C4","#006b31","#fddd00",
                               "#d80b13")) 
# 拼图，参考：https://mp.weixin.qq.com/s/qycXFTSt0VtHmm7z3NQcpg
p_total <- p + p1 +
  plot_annotation(tag_levels = "A")

ggsave("Plots/barplot_phylum_genus2.png", plot = p_total, width = 8, height = 6, dpi = 300)
