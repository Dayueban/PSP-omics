# More rigorous screening of metabolites
rm(list = ls())
library(pheatmap)
a <- read.csv("HFD vs PSP/1-差异分析/差异代谢物定量列表.csv", header = T,
              row.names = 1, check.names = F)

b <- a[, c(1, 7:ncol(a))]
#write.csv(b2, "Model1 vs Control/Heatmap_differ_metabo.csv") # 筛选出来的代谢物可以进行富集分析

# scale before heatmap clustering
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
#fec_conj <- read.csv("MSMS/msms_Intensity2.csv")
#fec_conj_num <- read.csv("MSMS/msms_Intensity2_1.csv")
#type <- c(rep("glucuronide",13),rep("sulfate",20))
# 运行一次结果后发现M015个体和SV组聚成一个cluster，因此建议删除
#fec_conj <- fec_conj[,-16]
#fec_conj_num <- fec_conj_num[,-15]
sample <- c(rep("HFD",6),rep("PSP",5))
#type <- b$type
fecPhase <- b
#rownames(fecPhase) <- rownames(b2)
fecPhase$type <- NULL
fecPhase_t <- t(fecPhase)
#fecPhase$change <- NULL
fecPhase_norm <- apply(fecPhase_t, 2, cal_z_score)
#fecPhase_samplerow <- data.frame(sample)
#fecPhase_samplerow <- data.frame(type)
#fecPhase_typecol <- data.frame(type)
fecPhase_typecol <- data.frame(sample)
#row.names(fecPhase_samplerow) <- rownames(fecPhase)
row.names(fecPhase_typecol) <- names(fecPhase)
ann_colors <- list(sample = c(HFD = "#A60D0D", PSP = "#006A44"))
#ann_colors <- list(sample = c(Control = "#E16663", Model3 = "#EAA58E"))
#fecPhase_norm_t <- t(fecPhase_norm)
fecPhase_pheatmap <- pheatmap(fecPhase_norm,
                              color = colorRampPalette(c("navy", "#FEF9E7", "firebrick3"))(500),
                              #annotation_col = fecPhase_samplerow,
                              annotation_row  = fecPhase_typecol,
                              cutree_cols = 2,
                              cutree_rows = 2,
                              #cluster_cols = FALSE,
                              cluster_rows = FALSE,
                              border_color = "NA",
                              fontsize_row = 7,
                              fontsize_col = 8,
                              annotation_colors = ann_colors
)
save_pheatmap_png <- function(x, filename, width=2000, height=1400, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(fecPhase_pheatmap, "heatmap_plot/fecPhase_pheatmap_HFD_PSP2.png")
