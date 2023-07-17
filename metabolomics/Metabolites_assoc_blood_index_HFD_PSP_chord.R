rm(list = ls())
#temp_nm2 <- "filtered"
#library (xlsx) ### Saving to spreadsheet
library (data.table) ### Fast read of large files into R
library (WGCNA) ### -	Clustering software. Previously reported work done using v1.34
library (flashClust) ### Clustering software
library (ppcor) ### Partial Spearman correlations, for confounder analysis. Previously reported work done using v1.0
library (gplots) ### Plotting
library (cowplot) ### Plotting; to arrange several plots on the same page
library (ggplot2) ### Plotting
library (plyr) ### Data transformations
library(stringr) ### Mainly Dataframe/Character manipulation
library(tidyverse)
library(circlize)

# load working datasheet
options (stringsAsFactors = FALSE)
#work_df <- read.csv("Ileum_OTU.csv",header = T,row.names = 1,check.names = F)
a <- read.csv("HFD vs PSP/1-差异分析/差异代谢物定量列表HFD-PSP.csv", header = T,
              row.names = 1, check.names = F)

dat <- as.data.frame(t(a[,7:ncol(a)]))
dat <- log10(dat + 0.0001)


# 差异的生化指标
temp_dir <- "F:\\其他人员\\朱亚玲\\PSP-姜黄素-生化指标测定\\"
b <- read.csv(paste(temp_dir,"/Blood_biochemical_indicators.csv",sep = ""), 
              header = T,row.names = 1, check.names = F)
b2 <- b[7:nrow(b), ]
# 对生化指标进行Z-score转化，消除量纲的影响
b_temp <- scale(b2)
b_temp <- as.data.frame(b_temp)

library(plyr)
library(psych)
library(reshape2)
corr_df <- corr.test(dat, b_temp, method = "spearman", adjust = "BH", alpha = .05)
corr_df_cor <- corr_df$r
corr_df_p <- corr_df$p.adj

#corr_df_cor <- corr_df_cor[abs(corr_df_cor) > 0.6 | corr_df_p < 0.05 ]

#' jiont the r and p matrix
asso_df <- rbind(corr_df_cor,corr_df_p)
#' transpose the data frame
asso_df <- as.data.frame(t(asso_df))
# rename the column variables
colnames(asso_df) <- paste(colnames(asso_df),rep(c("cor","P_adj"),c(dim(dat)[2],dim(dat)[2])),sep = "_")
#' re-order the asso_df according to the adjusted p value
#asso_df <- asso_df[order(asso_df$P_adj,decreasing = FALSE),]
# write.csv(asso_df, "M3C_metabo_transc_asso.csv",row.names = TRUE)

## use ggplot2
# Reset rownames
corr_df_cor <- data.frame(row=rownames(corr_df_cor),corr_df_cor,check.names = F) # create a column called "row" 
rownames(corr_df_cor) <- NULL
corr_df_p <- data.frame(row=rownames(corr_df_p),corr_df_p,check.names = F) # create a column called "row" 
rownames(corr_df_p) <- NULL
# Melt
nbar.m <- melt(corr_df_cor)
nbap.m <- melt(corr_df_p)
# Classify (you can classify differently for nbar and for nbap also)         
nbar.m$value2<-cut(nbar.m$value,breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),include.lowest=TRUE, label=c("(-0.75,-1)","(-0.5,-0.75)","(-0.25,-0.5)","(0,-0.25)","(0,0.25)","(0.25,0.5)","(0.5,0.75)","(0.75,1)")) # the label for the legend
nbap.m$value2<-cut(nbap.m$value,breaks=c(-Inf, 0.001, 0.01, 0.05),label=c("***", "**", "*")) 
nbar.m<-cbind.data.frame(nbar.m,nbap.m$value,nbap.m$value2) # adding the p value and its cut to the first dataset of R coefficients
names(nbar.m)[5]<-paste("valuep") # change the column names of the dataframe 
names(nbar.m)[6]<-paste("signif.")

# 另外一种作图方法
nbar.m <- nbar.m %>%
  mutate(value_abs = abs(value))
# 将物种名进行裁剪
nbar.m$row <- str_split(nbar.m$row, ";", simplify = TRUE)[,2]
# ab(value) > 0.6的才进行展示
##nbar.m0.6 <- nbar.m[which(nbar.m$value_abs > 0.6), ] # 相关系数大于0.6的
nbar.m <- nbar.m[!is.na(nbar.m$signif.), ] # 显著的
write.csv(nbar.m,"assoc_results/metabo_HFD_PSP_assoc_blood.csv")
#' part two: chord diagram plot
# 添加一列数值1，分别计算相同代谢物的个数，用于后续排序
df3 <- nbar.m
#write.csv(df3, file = "metabo_assoc_clinical_param.csv")
# remove column which the |correlation coefficients| less than 0.3
#df4 <- df3[-which(abs(df3$value) < 0.3) ,]
df4 <- df3[,1:3]
names(df4) <- c("Genus", "blood_index", "Coeff")
df4$value <- rep(1, nrow(df4))
# 根据列名taxo将相同细菌凑一起，然后根据value计算个数

data_total <- df4 %>%
  group_by(Genus) %>% 
  transmute(Total=sum(value))
# 根据列名blood_index将相同代谢物凑一起，然后根据value计算个数
data_total2 <- df4 %>%
  group_by(blood_index) %>% 
  transmute(Total2=sum(value))
# 将data_total中的Total列合并到df4数据框中
df4 <- cbind(df4, data_total$Total, data_total2$Total2)
names(df4)[c(ncol(df4)-1, ncol(df4))] <- c("Total", "Total2")

# change the first column by paste metabo names with numbers in column "Total"
df4$Genus <- paste(df4$Genus, "\t(",df4$Total, ")",sep = "")
df4$blood_index <- paste(df4$blood_index, "\t(", df4$Total2,")", sep = "")
# sort df4 based on the numbers of links
df4 <- df4[order(df4$Total2, decreasing = T),]
df4 <- df4[order(df4$Total, decreasing = T),]
# remove column 'value' and 'Total'
df4$value <- NULL
df4$Total <- NULL
df4$Total2 <- NULL

## set colours for segments
#df4 <- df4[-which(abs(df4$Corr_score) < 5), ] # only |Corr_score| > 3 were used for exhibition
circos.clear()
# 第一种图是将links根据正负相关显示的
#pdf(file = paste(temp_nm2,"_otus_metabo_circlize.pdf",sep = ""), height = 10, width = 10)
png(file = "assoc_results/metabo_blood_index_HFD_PSP_circlize202305010.png", height = 2400, width = 2400,
    res = 300)
#circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
#Put horizontally or vertically symmetric
# gap.after是为了将Otu和代谢物之间分一个gap
circos.par(start.degree = 0, gap.degree = 1, points.overflow.warning = FALSE,
           gap.after = c(rep(1.5, length(unique(df4$Genus))-1), 5, 
                         rep(1.5, length(unique(df4$blood_index))-1), 5),
           circle.margin = c(0.6,0.2,0.7,0.1)) #controls the margins on the left, right, bottom and top sides of the circle
#par(mar = rep(-2, 4))
#par(mar = c(-8,-8,0,0))

# 给代谢物上色，而Otu全部给予灰色
grid.col <- setNames(c(topo.colors(length(unique(df4$Genus))),rep("#BEBEBE", length(unique(df4$blood_index)))), 
                     c(unique(df4$Genus), unique(df4$blood_index)))
# now, plot the image with rotated labels
chordDiagram(df4,  
             preAllocateTracks = 1, 
             annotationTrack = "grid",
             annotationTrackHeight = c(0.03, 0.1),
             grid.col = grid.col,
             col = ifelse(df4$Coeff > 0, "#F76F72", "#A8D1DF"), # link为正负相关可以用两种颜色表示
             link.sort = TRUE, link.decreasing = FALSE)
#directional = 1, 
#direction.type = c("diffHeight", "arrows")
#link.arr.type = "big.arrow")

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.4), cex = 1)
}, bg.border = NA)
title(paste("Metabolites assocaited with" ,"blood indicators", sep = " "))
legend("topleft", pch = 15, bty = "n",col = c("#F76F72", "#A8D1DF"), title = "Association types",
       legend = c("Positive","Negative"), cex = 1.2)
dev.off()





