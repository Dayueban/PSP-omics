rm(list = ls());options(stringsAsFactors=FALSE)
library(openxlsx)
library(ggforce)
library(ggplot2)
library(dplyr)
NAME <- "HFD VS ND"
#setwd("F:/MJH/！科研相关/`yyds Dr.Zhu YL/paper/11_PSP Omics/PSP-RNA-seq/2volcano") 
setwd("F:/MJH/！科研相关/`yyds Dr.Zhu YL/paper/11_PSP Omics/PSP-RNA-seq/0DEseq_WY") 
exp <- read.csv("HFD_ND_deseq2.csv",header = T)
exp <- exp[,-6] #删去padj
exp=na.omit(exp) #删除NA值
exp <- exp[!duplicated(exp$ID),] #去除重复值
rownames(exp) <- exp$ID 
#修改列名，避免下面改代码
colnames(exp)[5] <- "pval"
exp$log2FoldChange=as.numeric(exp$log2FoldChange) #将log2FoldChange的数据类型转换为数值型（好像在这没必要？）
which(exp$log2FoldChange=="Inf")#这一步是为了查看是否有遗漏的缺失值
exp$log10pval=-log10(exp$pval)#创建新列-log10pval
results_frame <- exp

###---------------差异基因分析---FC1.0_P0.05###---------------
results_frame$up_down <- "NS"#先都设成不显著，再把显著的拿出来重新赋值。 #NS:Not Significant
results_frame[results_frame$log2FoldChange < -1.0&  
                results_frame$pval <0.05 &
                !is.na(results_frame$pval),]$up_down <- "DR" 

results_frame[results_frame$log2FoldChange > 1.0 &  
                results_frame$pval <0.05 &
                !is.na(results_frame$pval), ]$up_down <- "UR"

table(results_frame$up_down)



#火山图准备-#转换为因子指定绘图顺序；---
results_frame$up_down <- factor(results_frame$up_down,levels=c("UR","DR","NS")) #设置各组在图例中颜色设置的顺序 

###---------------volcon-ch###---------------
#将up和down的gene取出来，写成csv导出
result <- results_frame[results_frame$pval < 0.05 & abs(results_frame$log2FoldChange) > 1.0,]
result_order=result[order(result$pval,decreasing=F),]  ##筛选差异基因，并按照p值排序
nrow(result_order) # 264通过查看行数确定符合标准的差异基因数目
up <- result_order[result_order$pval < 0.05 & result_order$log2FoldChange>1.0,]
nrow(up) ## 117 查看up数目
down <- result_order[result_order$pval < 0.05 & result_order$log2FoldChange < -1.0,]
nrow(down) ## 147 down数目

#write.csv(up,paste0("000",NAME," fc1.0 p0.05_up 77 .csv"),row.names=F,quote=F)
#常规用csv导出即可，但是这次导出成csv出现了排序乱码，所以用xlsx
#write.xlsx(up,paste0(NAME,"_up_117_with_symbols.xlsx"))
#write.xlsx(down,paste0(NAME,"_down_147_with_symbols.xlsx"))
#write.xlsx(result_order,paste0(NAME," fc1.0 p0.05_ALLDIFF 264 .xlsx"))

###-------------------火山图绘制###---------------
library(ggplot2)
ggplot() + 
  geom_point(data=results_frame, #确定数据来源
             aes(x=log2FoldChange, #log2FoldChange为x轴
                 y=-log10(pval),  #pvalue 取-log10的值为纵坐标
                 colour=up_down, shape=up_down))  ####按照up_down分组来显示颜色,形状

range(results_frame$log2FoldChange)
range(results_frame$'log10pval')
g <-
  ggplot() + 
  geom_point(data=results_frame, 
             aes(x=log2FoldChange, 
                 y=-log10(pval), 
                 colour=up_down),size=2.5) +  #按照up_down分组来显示颜色
  scale_y_continuous(limits=c(0,9),breaks=c(0,2,4,6,8)) + #设置纵坐标极限,以及刻度的位置
  scale_x_continuous(limits=c(-7,7),breaks=c(-6,-4,-2,0,2,4,6)) +  
  scale_color_manual(breaks=c("UR","NS","DR"), 
                     values=c("#A60D0D","#aaaaaa","#00BCD2"),  #ND #00BCD2,HFD #A60D0D,Cur #006890,PSP #006A44, Atro #7027B8
                     #("#4A1985","#d8d8d8","#F8B606"),
                       #c("#CC0000","#aaaaaa","#2F5688"),  ##F8B606","#4A1985","#d8d8d8
                     #ND dark green ,#HFD red, #Cur #FF8000, #PSP #31b3fe, #Atro #e987fd
                    # 备选 ："orange","#bdbdbd","green","red","darkgreen"                   
                     #设置各组的颜色定义，这个颜色顺序设置与因子顺序一致，如果不设置为因子，则按组名字母顺序
                     #labels=c("Up Regulated","Down Regulated","NoSig"),
                     name="") +#设置图例的标题，这里冒号之间为空，就不显示标题
  xlab("log2(Fold Change)") +    #设置横轴标题
  ylab("-log10(pvalue)") + 
  ggtitle(NAME)+
  theme_bw() +   #通过这个来调到白底的背景，这与默认的不同
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background=element_rect(fill="transparent",colour="gray"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title=element_text(size=20,face="bold"))+
  theme(plot.title = element_text(size=20, face="bold"))+
  geom_hline(aes(yintercept=1.3),linetype=5,col="black")+  #颜色可自行修改 ： black ？？##设置横向辅助线
  geom_vline(aes(xintercept=c(-1,1)),linetype=5,col="black")+ #设置纵向辅助线
  theme(legend.position="none", #legend.title = element_blank(),
        # legend.position= c(0.8,0.9),#调节图例出现的位置。0~1之间，改变这个值，找到合适的位置
        # legend.direction="horizontal",   #图例的方向，"vertical"为水平排列。
        # legend.background=element_rect(fill="transparent",colour = NA),  #图例背景
        # legend.text=element_text(colour="black",size=30),#图例标记的字体大小,根据自己的需要调大调小
        axis.text.x = element_text(size=15,colour="black"),
        axis.text.y = element_text(size=15, angle=60),
        axis.title.x = element_text(size=16),  
        axis.title.y = element_text(size=16),
        plot.title = element_text(size=17,face="bold"))

g
#---筛选top10 up 和 down 基因，用于标记---#--
up_label=up[1:10,]
down_label=down[1:10,]
to_be_labeled <- rbind(up_label,down_label) #rbind 合并表格
nrow(to_be_labeled) 
#write.csv(to_be_labeled,file = paste(4,NAME,"dif_labeled.csv"),row.names=F)

#！！！注意-赶时间可以直接运行p6----
#突出差异表达基因#注意这个方法可以用于标注感兴趣的基因
#代码来源：https://mp.weixin.qq.com/s?__biz=MzIyOTY3MDA3MA==&mid=2247514615&idx=1&sn=d7801a40582c4ff97162fb6024df0fe9&chksm=e8bdd42edfca5d38e9240f5746679da5a92296215eb7ae6ed39a99d6618242cb2fb2a6c804dc&scene=126&sessionid=1681393763&subscene=227&key=44e84b829ca4f5d5a7cd4a368e1848eda1cd2f667b611e804e0743bc94d726d6af4985f27e6730e85f1a36383e6f554af977bd4173e694dbbebc76d20d9cde7d48ad8784d207e6061b2800a2046a5e01ba591a1485ad4445579592a3f4067dfd0d7e81b553ae246a679aa05b54a2598a48cc3dfd6712142d5a15734342c2cccf&ascene=0&uin=MjE1OTg5NzE2Mg%3D%3D&devicetype=Windows+10+x64&version=63090217&lang=zh_CN&countrycode=CN&exportkey=n_ChQIAhIQMe36NzJXWOhd%2BeljRlm29BLgAQIE97dBBAEAAAAAAPTGFDN%2BSkIAAAAOpnltbLcz9gKNyK89dVj04gENo6T9JvkhrFSIsDOvKg946V5iJCmWK1MIIVzjYHajv2DcDlxDOG6cQA76TSKoIggSIgDGgKGQXEqHrTDjvICLPqZuDuhNU5ypHY5OW4VxNhUz6TRuQbS1zEHYKiPAKmUJ%2FRuHAA9pRHxlPmjwwT%2BIAycH0oIQogdmlpnu39E2h5gHZTFojMG5fc7zlphVhpRpc6in0hBLpzz%2F1t8BUVaLpqBzQ7XjXO6IBVh8GXXTIQq7hL8uOrwc&acctmode=0&pass_ticket=MC5kbQUCHO8ucrA3Drwga%2Bh5HUMyN04RMxu7Ffu1iB6IO9mvMDfE7pUra5gCfD3%2FwiUwVOtaffdN3RaZkSXT6A%3D%3D&wx_header=1
gg <- g +
  geom_point(data = up_label,
             aes(x = log2FoldChange, y = -log10(pval)),
             color = '#006A44', size = 3.5, alpha = 1) +#ND #00BCD2,HFD #A60D0D,Cur #006890,PSP #006A44, Atro #7027B8
  geom_point(data = down_label,
             aes(x = log2FoldChange, y = -log10(pval)),
             color = '#A60D0D', size = 3.5, alpha = 1)
gg
#标注gene
#添加目标标签：
library(ggrepel)
library(patchwork)
ggg <- gg +
  geom_text_repel(data = up_label,
                  aes(x = log2FoldChange, y = -log10(pval), label = ID),
                  seed = 233,
                  size = 3.5,
                  color = '#006A44',
                  min.segment.length = 0,
                  force = 2,
                  force_pull = 2,
                  box.padding = 0.6,
                  max.overlaps = Inf
  ) + #上调相关标签(同一分组/配色)
  geom_text_repel(data = down_label,
                  aes(x = log2FoldChange, y = -log10(pval), label = ID),
                  seed = 233,
                  size = 3.5,
                  color = '#A60D0D',
                  min.segment.length = 0,
                  force = 2,
                  force_pull = 2,
                  box.padding = 0.6,
                  max.overlaps = Inf
  ) #下调相关标签
ggg
#选择将标签对齐排列：

p6 <- g +
  geom_text_repel(data = up_label,
                  aes(x = log2FoldChange, y = -log10(pval), label = ID),
                  seed = 233,
                  size = 4.5,
                  color = '#A60D0D',
                  min.segment.length = 0,
                  force = 2,
                  force_pull = 2,
                  box.padding = 0.1,
                  max.overlaps = Inf,
                  segment.linetype = 3,
                  segment.color = '#A60D0D',
                  segment.alpha = 1,
                  nudge_x = 3 - up_label$log2FoldChange,
                  direction = "y", #
                  hjust = 0 #0右对齐，1左对齐，0.5居中
  ) +
  geom_text_repel(data = down_label,
                  aes(x = log2FoldChange, y = -log10(pval), label = ID),
                  seed = 233,
                  size = 4.5,
                  color = '#00BCD2',
                  min.segment.length = 0,
                  force = 2,
                  force_pull = 2,
                  box.padding = 0.1,
                  max.overlaps = Inf,
                  segment.linetype = 3,
                  segment.color = '#00BCD2',
                  segment.alpha = 1,
                  nudge_x = -5 - down_label$log2FoldChange,
                  direction = "y",
                  hjust = 1
  )
p6
ggsave(paste0(NAME,"_vocalno_diff fc1 p0.05-3(ggsave).tiff"))
tiff(paste0(NAME,"_vocalno_diff fc1 p0.05.tiff"),width=5000, height=4500,res=600,compression = 'lzw')
p6
dev.off()










#---基因名称标记#--
library(ggrepel)
G <- g + geom_text_repel(data=to_be_labeled, 
                         aes(x=log2FoldChange, #log2FoldChange 列用作横坐标
                             y=-log10(pval)), 
                         label=to_be_labeled$ID,#标记
                         size=5)
tiff(paste0(NAME,"_vocalno_diff fc1 p0.05-2.tiff"),width=5000, height=4500,res=600,compression = 'lzw')
G
dev.off()
#到这就基本结束了----------



#top10 标记
up_label=up[1:5,]
down_label=down[1:5,]
to_be_labeled2 <- rbind(up_label,down_label) #rbind 合并表格
nrow(to_be_labeled2) 
#write.csv(to_be_labeled,file = paste(4,NAME,"dif_labeled.csv"),row.names=F)
p <- g+ geom_mark_circle(data = to_be_labeled2,
                          aes(x=log2FoldChange, #log2FoldChange 列用作横坐标
                              y=-log10(pval),fill = to_be_labeled2$ID, label = to_be_labeled2$ID),
                          expand = unit(2, "mm"),
                          show.legend = F,
                          label.fontsize = 15,
                          con.cap = 0)
p
ggsave(paste0(NAME,"_vocalno_diff fc1 p0.05-3.tiff"))





























#---基因名称标记#--
library(ggrepel)
G <- g + geom_text_repel(data=topgene, 
                         aes(x=log2FoldChange, #log2FoldChange 列用作横坐标
                             y=-log10(pvalue)), 
                         label=topgene$ID,#标记
                         size=4)
#tiff(paste0(NAME,"_vocalno_diff1.tiff"),width=5000, height=4500,res=600,compression = 'lzw')
G



#圈注感兴趣的基因：cebpb &ppard
interest1 <- topgene[12,]
interest2 <- topgene[20,]
interest <- rbind(interest1,interest2)
p <- G + geom_mark_circle(data = interest,
                            aes(x=log2FoldChange, #log2FoldChange 列用作横坐标
                                y=-log10(pvalue),fill = interest$ID, label = interest$ID),
                            expand = unit(3, "mm"),
                            show.legend = F,
                            label.fontsize = 10,
                            con.cap = 0)
p
ggsave(paste0(NAME,"_vocalno_diff_MJH.tiff"))
dev.off()
