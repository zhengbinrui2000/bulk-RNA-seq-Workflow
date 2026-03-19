#!/usr/bin/env Rscript

# =============================================================================
# 脚本名称: 样本QC - 相关性热图与PCA分析 (qc_and_pca.R)
# 功能描述:
#   该脚本是转录组数据分析中一个关键的质量控制(QC)步骤。它执行两个核心任务：
#   1. 计算所有样本之间的皮尔逊相关性，并绘制成热图。这有助于快速识别
#      样本间的相似性，检查生物学重复是否聚在一起，以及发现可能的异常样本。
#   2. 对样本进行主成分分析(PCA)，并绘制二维散点图。PCA是一种降维技术，
#      能展示样本在主要变异方向上的分布，是检查批次效应和样本聚类情况的
#      有力工具。
#
# 使用方法:
#   Rscript qc_and_pca.R <fpkm矩阵> <表型数据> <输出相关性热图.pdf> <输出PCA图.pdf>
# =============================================================================


# --- 1. 加载所需R包 ---
# 清除当前工作环境中的所有变量，确保脚本从一个干净的状态开始运行。
rm(list = ls())

# 加载所有必需的R包
library(ggplot2)     # 用于高级数据可视化 (主要用于绘制PCA图)。
library(ggrepel)     # 提供geom_text_repel，用于在图上添加不相互重叠的文本标签。
library(ggstatsplot) # 提供高级统计可视化函数，这里用于快速绘制精美的相关性热图。
library(grDevices)   # 用于pdf文字嵌入
library(showtext)    # 用于加载字体
font_add("Arial", "arial.ttf") # 加载Arial字体
showtext_auto() # 自动激活showtext，让它接管后续所有图形设备的文本渲染。
base_font_family <- "Arial" # 定义base_font_family变量为Arial

# --- 2. 解析命令行参数 ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("错误: 需要4个命令行参数。\n用法: Rscript qc_and_pca.R <fpkm矩阵> <表型数据> <输出相关性热图.pdf> <输出PCA图.pdf>", call. = FALSE)
}

# 按顺序为变量赋值，提高代码可读性
fpkm_matrix_file           <- args[1] # 输入：FPKM表达矩阵文件路径
phenotypic_data_file       <- args[2] # 输入：样本表型/分组文件路径
output_correlation_heatmap <- args[3] # 输出：相关性热图文件路径
output_pca_plot            <- args[4] # 输出：PCA图文件路径


# --- 3. 数据读取与预处理 ---
# 读取FPKM表达矩阵，`check.names = FALSE` 防止R自动修改可能包含'-'等特殊字符的样本名。
data <- read.table(fpkm_matrix_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 读取样本分组信息文件。
group <- read.table(phenotypic_data_file, header = TRUE, sep = "\t", check.names = FALSE)

# 过滤基因：保留在至少75%的样本中表达量大于0的基因。
keep <- rowSums(data > 0) >= floor(0.75 * ncol(data))
filter_fpkm_matrix <- data[keep, ]

# 对FPKM值进行log2(x+1)转换。
data_log2 <- log2(filter_fpkm_matrix + 1)

# 自动获取分组列名，约定使用表型文件的第二列作为分组依据。
if (ncol(group) < 2) {
  stop(paste("错误: 表型文件 '", phenotypic_data_file, "' 至少需要两列 (第一列样本名, 第二列分组信息)。"))
}
phenotype_col <- colnames(group)[2]



# --- 4. 生成并保存相关性矩阵热图 ---
p_corr <- ggcorrmat(
  data = data_log2,
  matrix.type = "upper",
  type = "parametric",
  colors = c("#2B6688", "#469393", "#84C2AE"),
  title = "Sample Correlation Heatmap"
) +
# 为ggcorrmat生成的图层添加主题设置，以修改字体
theme(
  plot.title = element_text(family = base_font_family, hjust = 0.5, size = 16, face = "bold"),
  axis.text.x = element_text(family = base_font_family, angle = 45, hjust = 1), # 旋转轴标签以防重叠
  axis.text.y = element_text(family = base_font_family)
)

# 使用 ggsave 保存为PDF
ggsave(output_correlation_heatmap, plot = p_corr, device = cairo_pdf, width = 8, height = 7, units = "in")
cat(paste("INFO: Correlation heatmap saved to '", output_correlation_heatmap, "'.\n", sep=""))


# --- 5. PCA 分析与可视化 ---
# 对数据进行主成分分析。
pca_res <- prcomp(t(data_log2), scale. = TRUE)

# 提取PCA的得分（scores）。
scores <- as.data.frame(pca_res$x)
scores$Sample <- rownames(scores)

# 数据合并
colnames(group)[1] <- "Sample"
scores <- merge(scores, group, by = "Sample", all.x = TRUE)

# 提取PC1和PC2的方差解释率。
explained_variance <- summary(pca_res)$importance[2, 1:2] * 100

# 使用ggplot2绘制PCA图。
pca_plot <- ggplot(scores, aes(x = PC1, y = PC2, color = .data[[phenotype_col]], shape = .data[[phenotype_col]])) +
  geom_point(size = 4, alpha = 0.7) +
  geom_text_repel(
    aes(label = Sample),
    size = 3,
    box.padding = 0.5,
    max.overlaps = Inf,
    show.legend = FALSE,
    family = base_font_family # 为ggrepel的标签也指定字体
  ) +
  scale_color_manual(values = c("#299D8F", "#E9C46A")) +
  scale_shape_manual(values = c(16, 15)) +
  labs(
    title = "PCA Analysis",
    subtitle = "PC1 vs PC2",
    x = paste0("PC1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PC2 (", round(explained_variance[2], 2), "%)"),
    color = "Group",
    shape = "Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    # 在这里为所有文本元素指定字体
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16, family = base_font_family),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30", family = base_font_family),
    axis.title = element_text(face = "bold", family = base_font_family),
    axis.text = element_text(family = base_font_family),
    legend.title = element_text(face = "bold", family = base_font_family),
    legend.text = element_text(family = base_font_family),
    legend.key = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

# 保存PCA图为PDF。
ggsave(output_pca_plot, plot = pca_plot, device = cairo_pdf, width = 8, height = 7, units = "in")
cat(paste("INFO: PCA plot saved to '", output_pca_plot, "'.\n", sep=""))

# 脚本结束时可以关闭showtext
showtext_auto(FALSE)
