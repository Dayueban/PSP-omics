# for PCA analsis with Lipid data
rm(list = ls())

# load packages
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
library(ggExtra)
library(ropls)

# load data
a <- read.csv("代谢物鉴定定量列表.csv", header = T, row.names = 1, check.names = F)
a1 <- a[, 5:ncol(a)]

a_t <- data.frame(t(a1))
a_t_log10 <- log10(a_t)

# PCA analysis for negative mode
a.bc <- vegdist(a_t_log10, na.rm = TRUE)  

a.bc.ma <- as.matrix(a.bc)
a.bc.df <- as.data.frame(a.bc.ma)
a.pca <- PCA(a.bc.df,graph = FALSE)

#' add the group information according to the rownames
groups <- rep("NA", ncol(a1))
groups[grep("^ND",names(a1))] <- "ND"
groups[grep("^HFD",names(a1))] <- "HFD"
groups[grep("^Cur",names(a1))] <- "Cur"
groups[grep("^PSP",names(a1))] <- "PSP"
groups[grep("^Atro",names(a1))] <- "Atro"



# 只取前面两个主成分
b <- as.data.frame(a.pca$ind$coord[,1:2])
b_PC <- a.pca$eig[,2][1:2]

# comp 1   comp 2 
# 45.20342 42.00901
b$group <- groups
# 分组标签按照所需要的顺序
b$group <- factor(b$group, levels = c("ND", "HFD", "Cur", "PSP", "Atro"))
names(b)[1:2] <- c("PCA1", "PCA2")
# 先只看ND, HFD和PSP三个组
b <- b[b$group %in% c("ND","HFD","PSP"), ]


# 另外一种作图方式
library(doBy)
se <- function(x) sd(x)/sqrt(length(x))
conf_95 <- function(x) t.test(x, conf.level = 0.95)$conf.int[2]
stat <- summaryBy(PCA1 + PCA2 ~ group, b, FUN = c(mean, sd, se, conf_95))

# (3)添加各组样本与其质心点的连线
# 各组样本与质心的连线
# 对于质心点，本示例直接使用各组的均值，已经在上一步中计算得到
plot_data <- merge(b, stat, by = 'group')

# (4) 作图1
p <- ggplot(plot_data) +
  geom_segment(aes(x = PCA1.mean, y = PCA2.mean, xend = PCA1, yend = PCA2, color = group), show.legend = FALSE) +
  geom_point(aes(x = PCA1, y = PCA2, color = group)) +
  #scale_color_manual(values = c("#00BCD2", "#A60D0D", "#006890","#006A44","#6F27B8")) +
  scale_color_manual(values = c("#00BCD2", "#A60D0D", "#006A44")) +
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'),
        panel.grid = element_blank(), legend.key = element_blank(),
        axis.text = element_text(color = 'black'), axis.ticks = element_line(color = 'black')) +
  xlab(paste("PCA1 ( ",round(b_PC[1], digits = 2),"%"," )", sep = "")) + ##坐标轴设置
  ylab(paste("PCA2 ( ",round(b_PC[2], digits = 2),"%"," )", sep = ""))

# 或者在质心点位置添加样本分组标签
p1 <- p +
  geom_point(data = stat, aes(x = PCA1.mean, y = PCA2.mean, color = group), shape = 22, fill = 'white', size = 10, show.legend = FALSE) +
  geom_text(data = stat, aes(x = PCA1.mean, y = PCA2.mean, color = group, label = group), size = 4, show.legend = FALSE) +
  stat_ellipse(aes(x = PCA1, y = PCA2, color = group),linetype = 2, level = 0.95, show.legend = FALSE) +
  #scale_color_manual(values = c("#00BCD2", "#A60D0D", "#006890","#006A44","#6F27B8")) +
  scale_color_manual(values = c("#00BCD2", "#A60D0D","#006A44")) +
  geom_vline(xintercept = 0, lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = 0, lty=4,col="black",lwd=0.5) +
  theme(legend.position = 'top',
        axis.title = element_text(size = 12, color = "black", face = "bold"),
        axis.text = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black"))

# 添加组间ANOSIM检验的结果，以表格展示
#anosim_all <- data.frame("Treatment"=c("Con Vs. AB","Con Vs. FMT","AB Vs. FMT","FMT Vs.DC"),
#                         "R"=c(0.892, 0.548, 0.012,0.26),
#                         "Pvalue"=c(0.015,0.006,0.314,0.07),stringsAsFactors = FALSE)

#anosim_all_p <- ggtexttable(anosim_all,rows = NULL, theme=ttheme("minimal", 
#                                                                 base_size = 12))

# 图形拼接
#pp <- p1 + annotation_custom(ggplotGrob(anosim_all_p),xmin = 3, ymin = 2.2, xmax = 5)

ggsave("metabo_PCA_plot20230423.png", plot = p1, width = 6, height = 6, dpi = 300)


# PLSDA analysis
plsda = opls(a_t_log10, groups, predI = 3)

# sample scores plot
sample.score = plsda@scoreMN %>% 
  as.data.frame() %>%
  mutate(group = groups)

# 先展示ND,HFD和PSP三组
sample.score <- sample.score[sample.score$group %in% c("ND","HFD","PSP"), ]
# 分组标签按照所需要的顺序
#sample.score$group <- factor(sample.score$group, levels = c("ND", "HFD", "Cur", "PSP", "Atro"))
sample.score$group <- factor(sample.score$group, levels = c("ND", "HFD", "PSP"))

p1 <- ggplot(sample.score, aes(p2, p3, color = group)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_point() +
  geom_point(aes(-10,-10), color = 'white') +
  labs(x = 'LDA2',y = 'LDA3') +
  stat_ellipse(level = 0.95, linetype = 'solid', 
               size = 1, show.legend = FALSE) +
  #scale_color_manual(values = c("#00BCD2", "#A60D0D", "#006890","#006A44","#6F27B8")) +
  scale_color_manual(values = c("#00BCD2", "#A60D0D","#006A44")) +
  theme_bw() +
  theme(legend.position = "top",
        legend.text = element_text(color = 'black',size = 12),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size = 15),
        axis.title = element_text(color = 'black',size = 15),
        axis.ticks = element_line(color = 'black'))
ggsave("metabo_PLSDA_plot20230423.png", plot = p1, width = 6, height = 6, dpi = 300)


library(scatterplot3d)
#install.packages("rgl")
png(filename = "metabo_LDA_3d_20230423.png", width = 6, height = 6, units = "in", res = 300)
par(mai = c(0.5, 0.5, 0.5, 0.5))
scatterplot3d(sample.score[,1:3], # 第1-3主成分
              # 颜色长度要和样本长度一样，且对应！
              color = rep(c("#00BCD2", "#A60D0D","#006A44"),
                          c(6,6,5)),
              pch = 15,
              lty.hide = 2,
              xlab = "LDA1", ylab = "LDA2", zlab = "LDA3",
              angle = 45
)
#par(usr=c(-4, 4, -4, 4, -4, 4))
legend("topleft",legend=c('ND','HFD', 'PSP'),
       fill=c("#00BCD2", "#A60D0D","#006A44"),box.col=NA)
dev.off()

