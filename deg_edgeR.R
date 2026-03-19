#!/usr/bin/env Rscript

# =============================================================================
# 脚本名称: RNA-Seq 差异表达分析与可视化 (edgeR)
# 功能描述:
#   该脚本是一个完整的差异表达分析流程，专为在Nextflow等自动化管道中使用而设计。
#   1. 使用edgeR包对原始基因计数矩阵进行差异表达分析。
#   2. 根据用户定义的阈值（logFC, FDR）筛选并标记差异表达基因（DEGs）。
#   3. 生成并保存详细的差异表达结果表。
#   4. 绘制并保存信息丰富的火山图，高亮标记最显著的DEGs。
#   5. 对Top 20的DEGs绘制基因表达热图，直观展示其在不同样本组间的表达模式。
#
# 使用方法 (在Nextflow等流程中被自动调用):
#   Rscript this_script.R <counts> <pheno> <ctrl> <exp> <logFC> <FDR> <out_deg> <out_volcano.pdf> <out_heatmap.pdf>
# =============================================================================


# --- 1. 环境准备与包加载 ---
# 清除当前工作环境中的所有变量，确保脚本从一个干净的状态开始运行。
rm(list = ls())

# 加载所有必需的R包
library(edgeR)          # 核心包，用于进行差异表达分析。
library(ggplot2)        # 用于绘制火山图，是数据可视化的基石。
library(dplyr)          # 提供高效的数据处理和转换工具（如 mutate, filter, arrange）。
library(RColorBrewer)   # 提供高质量的调色板，用于热图。
library(ComplexHeatmap) # 提供功能强大且灵活的热图绘制工具。
library(ggrepel)        # 用于在ggplot图中添加不重叠的文本标签。
library(tidyverse)      # 一个包含ggplot2, dplyr等多个包的集合，提供一致的数据科学生态。
library(grid)           # R的基础绘图系统，ComplexHeatmap可能会用到。
library(grDevices)   # 用于pdf文字嵌入
library(showtext)    # 用于字体管理
font_add("Arial", "arial.ttf") # 加载Arial字体
showtext_auto() # 自动激活showtext，让它接管后续所有图形设备的文本渲染。
base_font_family <- "Arial" # 定义base_font_family变量为Arial

# --- 2. 解析命令行参数 ---
# 从命令行捕获所有传入的参数。这是脚本与外部（如Nextflow）交互的方式。
args <- commandArgs(trailingOnly = TRUE)

# --- 为每个参数赋予有意义的变量名 ---
raw_counts_file          <- args[1]  # 1. 输入：原始基因计数矩阵文件路径 (*.txt)
phenotypic_data_file     <- args[2]  # 2. 输入：样本分组信息文件路径 (*.tsv)
control_group_param      <- args[3]  # 3. 参数：对照组的名称 (例如 "Control")
experimental_group_param <- args[4]  # 4. 参数：实验组的名称 (例如 "Treatment")
logFC_threshold_param    <- as.numeric(args[5]) # 5. 参数：log2 Fold Change 筛选阈值
FDR_threshold_param      <- as.numeric(args[6]) # 6. 参数：FDR (校正后p值) 筛选阈值
deg_results_file         <- args[7]  # 7. 输出：差异表达基因结果文件路径 (*.csv)
volcano_plot_file        <- args[8]  # 8. 输出：火山图文件路径 (*.pdf)
heatmap_plot_file        <- args[9]  # 9. 输出：热图文件路径 (*.pdf)

# --- 3. 数据读取与预处理 (edgeR流程) ---
# 读取原始计数矩阵。row.names = 1 表示第一列是基因ID。
raw_counts <- read.table(raw_counts_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
# 读取分组信息文件。
group <- read.csv(phenotypic_data_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
# 将分组信息列转换为因子(factor)类型，这是R中进行统计建模的必要步骤。
group$information <- factor(group$information, levels = c(control_group_param, experimental_group_param))

# 创建 DGEList 对象，这是edgeR进行分析的基础数据结构。
d <- DGEList(counts = raw_counts, group = group$information)
# 使用`filterByExpr`函数过滤掉在所有样本中表达量都极低的基因。
# 这是一个自动化的过程，可以去除噪音，增加后续统计检验的效能。
keep <- filterByExpr(d)
d <- d[keep, , keep.lib.sizes = FALSE]
# 使用TMM方法计算归一化因子，以校正不同样本间文库大小和RNA组成的差异。
d <- calcNormFactors(d)


# --- 4. 差异表达分析 ---
# 创建设计矩阵。`~0 + ...` 表示无截距模型，直接估计每个组的均值，便于后续定义组间比较。
design <- model.matrix(~0 + group$information)
# 将设计矩阵的列名设置为实际的组名，提高可读性。
colnames(design) <- levels(group$information)
# 估计离散度(dispersion)。这是edgeR模型中的关键一步，用于量化基因表达的生物学变异。
d <- estimateDisp(d, design)

# 使用广义线性模型(GLM)进行拟合。
fit <- glmFit(d, design)
# 动态构建对比(contrast)字符串，例如 "Treatment-Control"。
contrast_string <- paste(experimental_group_param, "-", control_group_param)
# 使用`makeContrasts`函数创建对比矩阵，明确告诉edgeR要进行哪个比较。
contrast <- makeContrasts(contrasts = contrast_string, levels = design)
# 进行似然比检验(LRT)来确定每个基因在指定对比下的差异显著性。
lrt <- glmLRT(fit, contrast = contrast)
# 提取所有基因的检验结果，并转换为数据框。
DEG_edgeR <- as.data.frame(topTags(lrt, n = nrow(d)))


# --- 5. 结果整理与保存 ---
# 使用dplyr的`case_when`函数，根据用户定义的阈值，为每个基因添加一个`change`状态列。
DEG_edgeR <- DEG_edgeR %>%
  mutate(change = case_when(
    FDR < FDR_threshold_param & logFC > logFC_threshold_param ~ "up",   # 上调
    FDR < FDR_threshold_param & logFC < -logFC_threshold_param ~ "down", # 下调
    TRUE ~ "stable" # 不满足以上条件的为稳定
  ))

# 将完整的差异表达分析结果保存到CSV文件中。
write.csv(DEG_edgeR, deg_results_file)

# 为了方便后续绘图，将行名（基因ID）转换为一个名为 'gene' 的列。
new_df <- cbind(gene = rownames(DEG_edgeR), DEG_edgeR)


# --- 6. 绘制火山图 ---
# 过滤掉极端值，防止这些点扭曲火山图的整体视觉效果，使图形更具可读性。
df_filtered <- new_df %>%
  filter(logFC >= -12 & logFC <= 12) %>%
  filter(-log10(PValue) <= 50) %>%
  filter(!is.na(logFC), !is.na(PValue)) # 移除任何包含NA值的行

# 筛选出最显著的前5个上调基因，用于在图上进行标记。
top_upregulated <- df_filtered %>%
  filter(logFC > 0) %>%
  arrange(PValue) %>% # 按P值从小到大排序
  slice_head(n = 5)

# 筛选出最显著的前5个下调基因。
top_downregulated <- df_filtered %>%
  filter(logFC < 0) %>%
  arrange(PValue) %>%
  slice_head(n = 5)

# --- 使用ggplot2构建火山图 ---
p1 <- ggplot(data = df_filtered, aes(x = logFC, y = -log10(PValue))) +
  # 主体散点图：所有基因。
  geom_point(alpha = 0.6, aes(size = -log10(PValue), color = logFC)) +
  # 设置一个漂亮的颜色渐变。
  scale_color_gradientn(colours = c("#2166ac", "#4393c3", "#ffffbf", "#de77ae", "#c51b7d"), values = seq(0, 1, 0.2)) +

  # 高亮层1：突出显示Top 5上调基因。
  geom_point(data = top_upregulated, aes(size = -log10(PValue)), shape = 21, color = "#000000", fill = "#ff7f00") +
  geom_text_repel(data = top_upregulated, aes(label = gene), box.padding = 0.5, max.overlaps = 50, family = base_font_family) +

  # 高亮层2：突出显示Top 5下调基因。
  geom_point(data = top_downregulated, aes(size = -log10(PValue)), shape = 21, color = "#000000", fill = "#ff7f00") +
  geom_text_repel(data = top_downregulated, aes(label = gene), box.padding = 0.5, max.overlaps = 50, family = base_font_family) +

  # 调整图例和坐标轴。
  scale_size(range = c(1, 6)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)), limits = c(0, 50)) +
  xlim(c(-12, 12)) +

  # 添加阈值线：垂直线代表logFC阈值，水平线代表FDR阈值。
  geom_vline(xintercept = c(-logFC_threshold_param, logFC_threshold_param), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(FDR_threshold_param), lty = 4, col = "black", lwd = 0.8) +

  # 设置标签和主题。
  xlab('log2 fold change') + ylab('-log10 PValue') +
  theme_bw() +
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = base_font_family),
        legend.background = element_rect(color = '#808080', linetype = 1),
        legend.title = element_text(family = base_font_family),
        legend.text = element_text(family = base_font_family),
        axis.text = element_text(size = 12.5, color = "#000000", family = base_font_family),
        axis.title = element_text(size = 15, color = "#000000", family = base_font_family)) +
  coord_cartesian(clip = "off") # 防止图形元素被坐标轴边缘裁剪。

# 将火山图（ggplot对象）保存为PDF
ggsave(volcano_plot_file, plot = p1, device = cairo_pdf, width = 8, height = 12, units = "in")


# --- 7. 绘制差异基因热图 ---
# 筛选出Top 20的差异基因（不分上下调，按FDR和logFC绝对值排序）。
deg_top20 <- DEG_edgeR %>%
  filter(change != "stable") %>%
  arrange(FDR, desc(abs(logFC))) %>%
  head(20)

# 设置只有在存在差异基因时才绘制热图。
if (nrow(deg_top20) > 0) {
  # 从原始计数矩阵中提取这Top 20基因的表达数据。
  rawdata_heatmap_top20 <- raw_counts[rownames(raw_counts) %in% rownames(deg_top20), ]
  # 准备热图的列注释（即样本分组信息）。
  annotation_col <- data.frame(group = group$information)
  rownames(annotation_col) <- colnames(rawdata_heatmap_top20)

  # 动态设置注释颜色，确保颜色与传入的组名正确对应。
  annotation_colors <- list(group = setNames(c("#339DB5", "#C9352B"), c(control_group_param, experimental_group_param)))

  # 对表达数据进行Z-score标准化。`t(scale(t(...)))` 是一个关键技巧，用于对行(基因)进行标准化。
  # 这使得我们能够比较不同基因在样本间的相对表达模式，而不是它们的绝对表达量。
  rawdata_heatmap_scaled_top20 <- t(scale(t(as.matrix(rawdata_heatmap_top20)), center = TRUE, scale = TRUE))

  # 使用ComplexHeatmap绘制热图，并通过gpar()对象指定字体为Arial
  p2 <- Heatmap(rawdata_heatmap_scaled_top20,
                name = "Expression Level", # 图例名称
                col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), # 设置颜色梯度
                show_row_names = TRUE,      # 显示基因名
                row_names_gp = gpar(fontfamily = base_font_family, fontsize = 10), # 为行名(基因)设置字体
                show_column_names = FALSE,  # 不显示样本名（通常在注释条中显示）
                cluster_rows = TRUE,        # 对行（基因）进行聚类
                cluster_columns = TRUE,     # 对列（样本）进行聚类
                top_annotation = HeatmapAnnotation(
                    df = annotation_col,
                    col = annotation_colors,
                    # 为注释图例的文本设置字体
                    annotation_legend_param = list(
                        title_gp = gpar(fontfamily = base_font_family, fontsize = 12),
                        labels_gp = gpar(fontfamily = base_font_family, fontsize = 10)
                    )
                ),
                heatmap_legend_param = list(
                    title = "Expression",
                    at = c(-2, 0, 2),
                    labels = c("Low", "Medium", "High"),
                    # 为热图主图例的标题和标签设置字体
                    title_gp = gpar(fontfamily = base_font_family, fontsize = 12),
                    labels_gp = gpar(fontfamily = base_font_family, fontsize = 10)
                )
  )

  # 将热图保存为PDF。
  cairo_pdf(heatmap_plot_file, width = 8, height = 10)
  draw(p2) # ComplexHeatmap对象需要用`draw()`函数来绘制。
  dev.off()

} else {
  # 如果没有发现任何差异基因，则创建一个空的输出文件。
  # 这可以防止Nextflow等自动化流程因找不到预期的输出文件而报错。
  file.create(heatmap_plot_file)
}

# 脚本结束时，可选地关闭showtext的自动渲染功能
showtext_auto(FALSE)
