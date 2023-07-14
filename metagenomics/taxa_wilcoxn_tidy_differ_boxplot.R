###################################################################
###################################################################
###################################################################
# (3) wilcoxon t-test
# 以ND、HFD、PSP三组先进行分析
rm(list = ls())
library(reshape)
library(ggpubr)
library(ggsignif)
library(tidyverse)
library(cowplot)
a <- read.csv("taxa_abundance/Absolute/otu_table.g.absolute.csv",header = T, row.names = 1,
              check.names = F)
# (2)数据处理
sum1 <- apply(a[,1:ncol(a)], 1, sum)
sum2 <- apply(a[,1:ncol(a)], 1, function(x)length(x[x == 0]))
sumall <- sum(sum1)
num_sample <- ncol(a)-1

## 
a$ratio1 <- (sum1 / sumall) * 100
a$ratio1 <- paste(sprintf("%.2f", a$ratio1), "%", sep = "")
a$ratio2 <- (sum2 / num_sample)
a$ratio2 <- paste(sprintf("%.2f", a$ratio2), "%", sep = "")

## to exclude species which with abundance less than 0.01% and present 
## in at least 20% of the study cohort
keep <- intersect(rownames(a)[which(a$ratio1 > 0.05)],
                  rownames(a)[which(a$ratio2 < 0.8)])

## 保留keep变量里面的那些species
a2 <- a[rownames(a) %in% keep, ]
a2$ratio1 <- NULL
a2$ratio2 <- NULL

#a2 <- a2[, c(1:12,20:24)]
# calculating the relative abundance of bacteria at genus level
a2 <- sweep(a2, 2, colSums(a2), '/')*100
write.csv(a2, "taxa_abundance/Absolute/otu_table_genus_ND_HFD_PSP_tidy.csv")
nr <- nrow(a2)
Pvalue_ND_HFD <- c(rep(0,nr))
Pvalue_HFD_PSP <- c(rep(0,nr))
Pvalue_ND_PSP <- c(rep(0,nr))
Pvalue_ND_HFD_PSP <- c(rep(0,nr))
# add group information
group_list1 <- rep(c("ND", "HFD"), c(6,6))
group_list2 <- rep(c("HFD", "PSP"), c(6,5))
group_list3 <- rep(c("ND", "PSP"), c(6,5))
group_list4 <- rep(c("ND","HFD","PSP"), c(6,6,5))
#log2_FC <- c(rep(0,nr))
#FC <- c(rep(0,nr))
#group_list <- rep(c("Ctrl","RSV"), c(6,7))  
for (i in 1:nr){
  Pvalue_ND_HFD[i] <- wilcox.test(as.numeric(a2[i,c(grep("^ND", names(a2)),grep("^HFD",names(a2)))]) ~ 
                                    group_list1, data = a2)$p.value
  Pvalue_HFD_PSP[i] <- wilcox.test(as.numeric(a2[i,c(grep("^HFD", names(a2)),grep("^PSP",names(a2)))]) ~ 
                                     group_list2, data = a2)$p.value
  Pvalue_ND_PSP[i] <- wilcox.test(as.numeric(a2[i,c(grep("^ND", names(a2)),grep("^PSP",names(a2)))]) ~ 
                                    group_list3, data = a2)$p.value
  Pvalue_ND_HFD_PSP[i] <- kruskal.test(as.numeric(a2[i,c(grep("^ND", names(a2)),
                                                         grep("^HFD", names(a2)),
                                                         grep("^PSP", names(a2)))]) ~ group_list4)$p.value
  #Pvalue[i] <- wilcox.test(as.numeric(a3[i,]) ~ group_list, data = a3)$p.value
  #Pvalue[i] <- kruskal.test(as.numeric(a[i,]) ~ group_list, data = a)$p.value
  #FC[i] <- (mean(as.numeric(a3[i,1:6]))+0.001)/ (mean(as.numeric(a3[i,7:13]))+0.001)
  #log2_FC[i] <- log2((mean(as.numeric(a3[i,1:6]))+0.001)/ (mean(as.numeric(a3[i,7:13]))+0.001))
}

#Pvalue <- as.numeric(Pvalue)  #在操作的过程中发现秩和检验的包输出来不是列表的形式，这里变成列表
fdr.ND_HFD <- p.adjust(Pvalue_ND_HFD,method="BH",length(Pvalue_ND_HFD))    #p.adjust就是计算FDR的包，这个可要记得了
fdr.HFD_PSP <- p.adjust(Pvalue_HFD_PSP,method="BH",length(Pvalue_HFD_PSP))
fdr.ND_PSP <- p.adjust(Pvalue_ND_PSP,method="BH",length(Pvalue_ND_PSP))
fdr.ND_HFD_PSP <- p.adjust(Pvalue_ND_HFD_PSP,method="BH",length(Pvalue_ND_HFD_PSP))
#a3_res <- cbind(a3,Pvalue,fdr.w, FC,log2_FC)
a3 <- cbind(a2,Pvalue_ND_HFD,Pvalue_HFD_PSP, Pvalue_ND_PSP, Pvalue_ND_HFD_PSP,
            fdr.ND_HFD,fdr.HFD_PSP, fdr.ND_PSP,fdr.ND_HFD_PSP)
write.csv(a3, file="taxa_abundance/Absolute/ND_HFD_PSP_genus_tidy_wilcoxn_test_BH20230424-3.csv") # 校正后没有差异，计算完相对丰度后呢？

a3 <- read.csv("taxa_abundance/Absolute/ND_HFD_PSP_genus_tidy_wilcoxn_test_BH20230424-3.csv", header = T,
               row.names = 1, check.names = F)
### 作图
a4 <- a3[a3$fdr.ND_HFD_PSP < 0.05, ]
a5 <- a4[, 1:17]
a6 <- cbind(rownames(a5), a5)
names(a6)[1] <- "Genus"
a6_t <- melt(a6)

a6_t$group <- "NA"
#rownames(alpha_all) <- as.character(rownames(alpha_all))
a6_t$group[grep("^ND",a6_t$variable)] <- "ND"
a6_t$group[grep("^HFD",a6_t$variable)] <- "HFD"
a6_t$group[grep("^PSP",a6_t$variable)] <- "PSP"
a6_t$group <- factor(a6_t$group, levels = c("ND","HFD","PSP"))

my_comparisons_a <- list(c("ND","HFD"), c("ND", "PSP"), c("HFD","PSP"))

p_0 <- ggplot(a6_t, aes(x=Genus,y=value,fill=group)) +
  geom_boxplot(position = position_dodge(1))
# 参考：https://mp.weixin.qq.com/s/j_mGrI_nvH8wxNU1yH8R8A

## 使用ggpurb包 ------------------------
# 默认使用anova检验三组之间的差异：
# 批量t检验：
stat.test <- a6_t %>%
  group_by(Genus) %>%
  pairwise_t_test(
    value ~ group, paired = FALSE, 
    p.adjust.method = "none"
  )
# 添加显著性标记
stat.test <- stat.test %>% add_xy_position(x = "Genus")

p_3 <- p_0 + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  bracket.size = 0.3, # 粗细
  tip.length = 0.01  # 两边竖线的长度
) + ggtitle("Genus")
p_3

#Change the order of items in the legend
p1 <- p1 + theme_classic()
p1 <- p1 +stat_compare_means(comparisons = my_comparisons_a,method = "t.test", label = "p")+
  #stat_compare_means(label.y = 3500)+
  theme(axis.text.x=element_text(angle=60,color="black",vjust = 0.95,hjust = 0.95,size=14, face="bold"),
        axis.title = element_text(size = 16, color = "black"),
        legend.position = "none")+
  scale_color_manual(values=c("#00BCD2", "#A60D0D", "#006A44")) +
  scale_fill_manual(values = c("#00BCD2", "#A60D0D", "#006A44")) +
  geom_jitter(position=position_jitter(0.2))
p1 # 这个方法暂时行不通

## 按单个菌来进行作图
a5_t <- as.data.frame(t(a5))
a5_t <- a5_t[,c(1,3,14,16)] #只选择差异比较明显的
genus_label <- names(a5_t)
names(a5_t)[1:dim(a5_t)[2]] <- paste(rep("G", dim(a5_t)[2]), seq(1,dim(a5_t)[2],1), sep = "_")
a5_t$group <- "NA"
a5_t$group[grep("^ND",rownames(a5_t))] <- "ND"
a5_t$group[grep("^HFD",rownames(a5_t))] <- "HFD"
a5_t$group[grep("^PSP",rownames(a5_t))] <- "PSP"
a5_t$group <- factor(a5_t$group, levels = c("ND","HFD","PSP"))

my_comparisons_a <- list(c("ND","HFD"), c("ND", "PSP"), c("HFD","PSP"))
#colors_index <- c("#00BCD2", "#A60D0D","#006A44")
#names(colors_index) <- c("ND", "HFD", "PSP")
#colors_index <- as.data.frame(colors_index)
#colors_index$group <- rownames(colors_index)
library(ggpubr)
library(ggsignif)
Start = 1
Stop = ncol(a5_t)-1
gvec <- vector("list",length = length(Start:Stop))
for(i in Start:Stop){
  m=names(a5_t)[i]
  index=a5_t[,c(m,"group")]
  #index$group <- factor(index$group, levels = c("ND","HFD","PSP"))
  batch_neu <- ggplot(data = index, 
                      aes_string(x="group",y= m, fill="group"))+
    scale_fill_manual(values = c("#A60D0D", "#00BCD2","#006A44"))+
    scale_x_discrete(limits=c("ND", "HFD", "PSP"))+
    stat_compare_means(comparisons = my_comparisons_a,method = "wilcox.test", 
                       label = "p.signif",size=5,vjust = 0.5)+
    geom_jitter(size = 1.2)+
    #geom_violin(size=1.5, alpha=.6)+
    geom_boxplot() +
    xlab("Treatment group")+
    ylab("Relitive Abundance")+
    ggtitle(genus_label[i])+
    theme_classic()+
    theme(legend.position = "none")+
    #geom_text(data=index, aes(x=group, y=y, color="black", label= stat),size=4) +
    #scale_y_continuous(labels=function(y) format(y,scientific=FALSE))+
    theme(axis.title.y = element_text(size = 10,face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 10),
          axis.text.x = element_text(size = 8, color = "black", angle = 45,
                                     hjust = 0.95, vjust = 0.95),
          axis.text.y = element_text(size = 10, color = "black"))
  ggsave(batch_neu, file = paste0("Plots/plot_", genus_label[i], ".png"),
         units="in", width=4, height=4, dpi=300)
  gvec[[i-Start+1]] <- batch_neu
}

pp <- plot_grid(plotlist = gvec, nrow = 1) # R package cowplot
ggsave(pp, file = "Plots/genus_intersect_boxplot3.png", height = 3, width = 12, dpi = 300)


#*******************second part--drawing barplot***************#
df <- a3[a3$Pvalue_HFD_PSP < 0.05, ]
df1 <- df[, 1:17]
df2 <- df1[, c(grep("^HFD",names(df1)), grep("^PSP", names(df1)))]
write.csv(df2, "Genus_HFD_PSP_P0.05.csv")
df2_t <- as.data.frame(t(df2))
df2_t <- cbind(rownames(df2_t), df2_t)
names(df2_t)[1] <- "ID" 
# data2 <- read.csv("taxonomy/Species/species_differ_plot.csv",header = T, check.names = F)
library(reshape2)
## reshape dataframe, from long sheet to wild sheet
data2_rs <- melt(df2_t,id.vars = "ID",variable.name = "Taxa", value.name = "Relative_abundance")
# add a new cloumn as the group information
data2_rs$group <- "NA"
data2_rs$ID <- as.character(data2_rs$ID)
data2_rs$group[grep("^HFD",data2_rs$ID)] <- "HFD"
data2_rs$group[grep("^PSP",data2_rs$ID)] <- "PSP"
data2_rs$group <- factor(data2_rs$group, levels = c("HFD","PSP"))

# draw barplot a
library(ggplot2)
#par(mar=c(2,2,2,2))
p <- ggplot(data=data2_rs,aes(x=Taxa,y=Relative_abundance,fill=group))+geom_boxplot(notch = FALSE)
p <- p + theme_classic()
p <- p + theme(axis.text.x=element_text(angle=60,color="black",vjust = 0.95,hjust = 0.95,size=8, face="bold"),
               axis.text.y=element_text(size = 12, colour = "black"),
               axis.title.y=element_text(size = 14,face = "bold"))
p <- p + theme(plot.background=element_blank(),panel.background=element_blank())
p <- p + guides(fill=guide_legend(title=NULL))
p <- p + theme(legend.position = "top")
p <- p + theme(plot.title = element_text(size = 16))
p <- p + xlab("")
p <- p + labs(y = "Relative abundance (log10)")
p <- p + scale_fill_manual(values=c("#A60D0D", "#006A44"),
                           labels = c("HFD", "PSP"))
p <- p + theme(axis.line.x=element_line(colour = "black"),axis.line.y=element_line(colour = "black"),
               legend.key.height=unit(1,"cm"),
               legend.key.width=unit(1,"cm"),
               legend.text=element_text(lineheight=0.8,face="bold",size=14),
               legend.title=element_text(size=16))
#set the margin of plot
p <- p + theme(plot.margin = unit(c(1,0.5,1,0.5),"cm"))
p <- p + scale_y_continuous(trans = "log10")
#p <- p + scale_y_continuous(limits=c(miny*1.1,maxy*1.1))
#p <- p + geom_jitter() ##添加数据点
ggsave("Plots/Genus_diff_HFD_PSP20230424.png", device = "png", width = 5, height = 6, dpi = 300)


# ND与HFD
df <- a3[a3$Pvalue_ND_HFD < 0.05, ]
df1 <- df[, 1:17]
df2 <- df1[, c(grep("^ND",names(df1)), grep("^HFD", names(df1)))]
write.csv(df2, "Genus_ND_HFD_P0.05.csv")
df2_t <- as.data.frame(t(df2))
df2_t <- cbind(rownames(df2_t), df2_t)
names(df2_t)[1] <- "ID" 
# data2 <- read.csv("taxonomy/Species/species_differ_plot.csv",header = T, check.names = F)
library(reshape2)
## reshape dataframe, from long sheet to wild sheet
data2_rs <- melt(df2_t,id.vars = "ID",variable.name = "Taxa", value.name = "Relative_abundance")
# add a new cloumn as the group information
data2_rs$group <- "NA"
data2_rs$ID <- as.character(data2_rs$ID)
data2_rs$group[grep("^ND",data2_rs$ID)] <- "ND"
data2_rs$group[grep("^HFD",data2_rs$ID)] <- "HFD"
data2_rs$group <- factor(data2_rs$group, levels = c("ND", "HFD"))
# draw barplot a
library(ggplot2)
#par(mar=c(2,2,2,2))
p <- ggplot(data=data2_rs,aes(x=Taxa,y=Relative_abundance,fill=group))+geom_boxplot(notch = FALSE)
p <- p + theme_classic()
p <- p + theme(axis.text.x=element_text(angle=60,color="black",vjust = 0.95,hjust = 0.95,size=8, face="bold"),
               axis.text.y=element_text(size = 12, colour = "black"),
               axis.title.y=element_text(size = 14,face = "bold"))
p <- p + theme(plot.background=element_blank(),panel.background=element_blank())
p <- p + guides(fill=guide_legend(title=NULL))
p <- p + theme(legend.position = "top")
p <- p + theme(plot.title = element_text(size = 16))
p <- p + xlab("")
p <- p + labs(y = "Relative abundance (log10)")
p <- p + scale_fill_manual(values=c("#00BCD2","#A60D0D"),
                           labels = c("ND", "HFD"))
p <- p + theme(axis.line.x=element_line(colour = "black"),axis.line.y=element_line(colour = "black"),
               legend.key.height=unit(1,"cm"),
               legend.key.width=unit(1,"cm"),
               legend.text=element_text(lineheight=0.8,face="bold",size=14),
               legend.title=element_text(size=16))
#set the margin of plot
p <- p + theme(plot.margin = unit(c(1,0.5,1,0.5),"cm"))
p <- p + scale_y_continuous(trans = "log10")
#p <- p + scale_y_continuous(limits=c(miny*1.1,maxy*1.1))
#p <- p + geom_jitter(width = 0.2, height = 0.4) ##添加数据点
ggsave("Plots/Genus_diff_ND_HFD20230424.png", device = "png", width = 8, height = 6, dpi = 300)

# ND与PSP
df <- a3[a3$Pvalue_ND_PSP < 0.05, ]
df1 <- df[, 1:17]
df2 <- df1[, c(grep("^ND",names(df1)), grep("^PSP", names(df1)))]
write.csv(df2, "Genus_ND_PSP_P0.05.csv")
df2_t <- as.data.frame(t(df2))
df2_t <- cbind(rownames(df2_t), df2_t)
names(df2_t)[1] <- "ID" 
# data2 <- read.csv("taxonomy/Species/species_differ_plot.csv",header = T, check.names = F)
library(reshape2)
## reshape dataframe, from long sheet to wild sheet
data2_rs <- melt(df2_t,id.vars = "ID",variable.name = "Taxa", value.name = "Relative_abundance")
# add a new cloumn as the group information
data2_rs$group <- "NA"
data2_rs$ID <- as.character(data2_rs$ID)
data2_rs$group[grep("^ND",data2_rs$ID)] <- "ND"
data2_rs$group[grep("^PSP",data2_rs$ID)] <- "PSP"
data2_rs$group <- factor(data2_rs$group, levels = c("ND", "PSP"))
# draw barplot a
library(ggplot2)
#par(mar=c(2,2,2,2))
p <- ggplot(data=data2_rs,aes(x=Taxa,y=Relative_abundance,fill=group))+geom_boxplot(notch = FALSE)
p <- p + theme_classic()
p <- p + theme(axis.text.x=element_text(angle=60,color="black",vjust = 0.95,hjust = 0.95,size=8, face="bold"),
               axis.text.y=element_text(size = 12, colour = "black"),
               axis.title.y=element_text(size = 14,face = "bold"))
p <- p + theme(plot.background=element_blank(),panel.background=element_blank())
p <- p + guides(fill=guide_legend(title=NULL))
p <- p + theme(legend.position = "top")
p <- p + theme(plot.title = element_text(size = 16))
p <- p + xlab("")
p <- p + labs(y = "Relative abundance (log10)")
p <- p + scale_fill_manual(values=c("#00BCD2","#006A44"),
                           labels = c("ND", "PSP"))
p <- p + theme(axis.line.x=element_line(colour = "black"),axis.line.y=element_line(colour = "black"),
               legend.key.height=unit(1,"cm"),
               legend.key.width=unit(1,"cm"),
               legend.text=element_text(lineheight=0.8,face="bold",size=14),
               legend.title=element_text(size=16))
#set the margin of plot
p <- p + theme(plot.margin = unit(c(1,0.5,1,0.5),"cm"))
p <- p + scale_y_continuous(trans = "log10")
#p <- p + scale_y_continuous(limits=c(miny*1.1,maxy*1.1))
#p <- p + geom_jitter(width = 0.2, height = 0.4) ##添加数据点
ggsave("Plots/Genus_diff_ND_PSP20230424.png", device = "png", width = 8, height = 6, dpi = 300)

