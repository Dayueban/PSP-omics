# represent the ratio of Firmicutes and Bacteroidetes
rm(list = ls())

# load packages
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(scales)
library(patchwork)
library(cowplot)
library(ggpubr)

# load dataset
# feces
a <- read.csv("taxa_abundance/Relative/Firmicutes_vs_Bacteroidetes_ratio.csv", header = T, check.names = F)

# saliva
#b <- read.csv("saliva_metag/Saliva_Firmicutes_vs_Bacteroidetes_ratio.csv", header = T, check.names = F)

# combine two dataset based on row
#ab <- rbind(a, b)

# add group information
ab <- a
ab <- ab %>%
  mutate(log2 = log2(ab$`Firmicutes/Bacteroidetes`))
ab$group <- NA
ab$group[grep("^ND", ab$ID)] <- "ND"
ab$group[grep("^HFD", ab$ID)] <- "HFD"
ab$group[grep("^Cur", ab$ID)] <- "Cur"
ab$group[grep("^PSP", ab$ID)] <- "PSP"
ab$group[grep("^Atro", ab$ID)] <- "Atro"

# ND, HFD and PSP groups were shown
ab <- ab[ab$group %in% c("ND", "HFD", "PSP"), ]

ab$group <- factor(ab$group, levels = c("ND","HFD","PSP"))

# add bodysite information
#ab$bodysite <- NA
#ab$bodysite[grep("*Fe", ab$ID)] <- "Feces"
#ab$bodysite[grep("*Sa", ab$ID)] <- "Saliva"


#my_comparisons_a <- list(c("ND","HFD"), c("ND", "Cur"), c("ND", "PSP"), c("ND", "Atro"),
#                         c("HFD","Cur"), c("HFD","PSP"), c("HFD","Atro"), c("Cur","PSP"),
#                         c("Cur","Atro"), c("PSP","Atro"))

my_comparisons_a <- list(c("ND","HFD"), c("ND", "PSP"), c("HFD","PSP"))
p1 <- ggplot(ab, aes(x = group,y = log2, fill = group)) +
  geom_violin(aes(fill = group)) +
  geom_boxplot(width=0.1, color = "white", size = 0.75) +
  #facet_grid(~bodysite) +
  
  #Change the order of items in the legend
  theme_classic() +
  scale_x_discrete(limits=c("ND", "HFD", "PSP")) +
  stat_compare_means(comparisons = my_comparisons_a,method = "t.test", label = "p")+
  #stat_compare_means(label.y = 3500)+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16),
        panel.background = element_rect(colour = "black", size = 1)) +
  ylab("Log2 Firmicutes/Bacteroidetes ratio") +
  xlab("") +
  labs(title = "Feces") +
  #strip.background = element_rect(fill=c("#C2C4C5")), # 分面的背景色设置
  #strip.text = element_text(size = 15,colour = "black", face = "bold")) +  # 分面的字体设置
  scale_color_manual(values=c("#00BCD2", "#A60D0D","#006A44")) +
  scale_fill_manual(values = c("#00BCD2", "#A60D0D","#006A44")) + 
  geom_jitter(position=position_jitter(0.2))
p1
ggsave("Plots/Firmucutes_vs_Bacteroidetes_ratio2.png", plot = p1, width = 4, height = 6)


