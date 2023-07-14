rm(list = ls())
options(stringsAsFactors = FALSE)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(rstatix)
library(vegan)
library(picante)
library("RColorBrewer")
display.brewer.all()

# stool
alpha_df <- read.csv("../03.AlphaDiversity/alpha_diversity_index.csv",header = T,sep = ',', row.names = 1)
alpha_all <- alpha_df[, c(2,4)]
#httr::set_config( httr::config( ssl_verifypeer = 0L ) )
#violin plot
alpha_all$group <- "NA"
#rownames(alpha_all) <- as.character(rownames(alpha_all))
alpha_all$group[grep("^ND",rownames(alpha_all))] <- "ND"
alpha_all$group[grep("^HFD",rownames(alpha_all))] <- "HFD"
alpha_all$group[grep("^Cur",rownames(alpha_all))] <- "Cur"
alpha_all$group[grep("^PSP",rownames(alpha_all))] <- "PSP"
alpha_all$group[grep("^Atro",rownames(alpha_all))] <- "Atro"

# 先只展示ND、HFD和PSP组
alpha_all <- alpha_all[alpha_all$group %in% c("ND","HFD","PSP"), ]
alpha_all$group <- factor(alpha_all$group, levels = c("ND","HFD","PSP"))

# do alpha analysis without group Severe
#alpha_all <- alpha_all[5:nrow(alpha_all), ]

#my_comparisons_a <- list(c("ND","HFD"), c("ND", "Cur"), c("ND", "PSP"), c("ND", "Atro"),
#                         c("HFD","Cur"), c("HFD","PSP"), c("HFD","Atro"), c("Cur","PSP"),
#                         c("Cur","Atro"), c("PSP","Atro"))

my_comparisons_a <- list(c("ND","HFD"), c("ND", "PSP"), c("HFD","PSP"))
#my_comparisons_a <- list(c("Health", "Mild"))
#stat_t <- t_test(data = alpha_all, Chao1~group)
#stat_t <- add_significance(stat_t, 'p')
#stat_t.test <-  add_xy_position(stat_t, x = 'group', dodge = 0.8)
p1 <- ggplot(alpha_all, aes(x=group,y=chao1,fill=group)) +
  geom_boxplot(position = position_dodge(1))
#Change the order of items in the legend
p1 <- p1 + theme_classic()
p1 <- p1 +  scale_x_discrete(limits=c("ND", "HFD", "PSP")) +
  stat_compare_means(comparisons = my_comparisons_a,method = "t.test", label = "p")+
  #stat_compare_means(label.y = 3500)+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.position = "none")+
  scale_color_manual(values=c("#00BCD2", "#A60D0D", "#006A44")) +
  scale_fill_manual(values = c("#00BCD2", "#A60D0D", "#006A44")) +
  geom_jitter(position=position_jitter(0.2))
p1




##Shannon index
p2 <- ggplot(alpha_all, aes(x=group,y=shannon,fill=group)) +
  geom_boxplot(position = position_dodge(1))
#Change the order of items in the legend
p2 <- p2 + theme_classic()
p2 <- p2 +  scale_x_discrete(limits=c("ND", "HFD", "PSP"))+
  stat_compare_means(comparisons = my_comparisons_a, method = "t.test",label = "p")+
  #stat_compare_means(label.y = 3500)+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.position = "none")
p2 <- p2 + scale_color_manual(values=c("#00BCD2", "#A60D0D","#006A44")) +
  scale_fill_manual(values= c("#00BCD2", "#A60D0D","#006A44"))
p2 <- p2 + geom_jitter(position=position_jitter(0.2))
p2
#Arrange on one page
p_total <- ggarrange(p1,p2,labels = "A")
ggsave("Plots/Feces_alpha_compare20230419.png", plot = p_total, width = 6, height = 6, dpi = 300)
