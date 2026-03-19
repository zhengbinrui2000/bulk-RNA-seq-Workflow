#!/usr/bin/env Rscript

# =============================================================================
# 脚本名称: FPKM表达谱分布小提琴图绘制 (fpkm_violin_plot.R)
# 功能描述:
#   - 读取FPKM格式的基因表达矩阵和样本分组信息。
#   - 对数据进行过滤，去除低表达基因，并进行log2(x+1)转换。
#   - 使用ggplot2绘制小提琴图，并叠加箱线图，以直观地比较不同组别
#     样本整体的基因表达水平分布。
#
# 使用方法:
#   Rscript fpkm_violin_plot.R <fpkm矩阵> <表型数据> <输出文件名.pdf> <对照组名> <实验组名>
# =============================================================================


# --- 1. 加载所需R包 ---
# 清除当前工作环境中的所有变量，确保脚本从一个干净的状态开始运行。
rm(list = ls())

# 加载所有必需的R包
library(ggplot2)   # 用于高级数据可视化，是绘图的核心。
library(tidyr)     # 用于数据整理，特别是将宽数据转换为长数据 (pivot_longer)。
library(dplyr)     # 用于数据处理和转换，提供了一套高效、易读的函数。
library(grDevices) # 用于pdf文字嵌入
library(showtext)  # 用于加载字体
font_add("Arial", "arial.ttf") # 假定arial.ttf在系统字体路径中
showtext_auto() # 自动激活showtext，让它接管后续所有图形设备的文本渲染。
base_font_family <- "Arial" # 定义base_font_family变量为Arial

# --- 2. 解析命令行参数 ---
# `commandArgs(trailingOnly = TRUE)` 获取在脚本名之后的所有参数。
args <- commandArgs(trailingOnly = TRUE)
# 检查参数数量是否正确，否则停止执行并给出使用提示。
if (length(args) != 5) {
  stop("错误: 需要5个命令行参数。\n用法: Rscript fpkm_violin_plot.R <fpkm矩阵> <表型数据> <输出文件名.pdf> <对照组名> <实验组名>", call. = FALSE)
}

# 按顺序为变量赋值
fpkm_matrix_file        <- args[1] # 输入：FPKM表达矩阵文件路径
phenotypic_data_file    <- args[2] # 输入：样本表型/分组文件路径
output_file             <- args[3] # 输出：图片文件路径
control_group_name      <- args[4] # 输入：对照组名称 (用于设定绘图顺序)
experimental_group_name <- args[5] # 输入：实验组名称 (用于设定绘图顺序)


# --- 3. 数据读取与预处理 ---
# 读取FPKM表达矩阵。
data <- read.table(fpkm_matrix_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 读取样本分组信息文件。
group <- read.table(phenotypic_data_file, header = TRUE, sep = "\t", check.names = FALSE)

# 约定：表型文件第一列为样本名，第二列为分组信息。为保证后续代码的健壮性，这里进行重命名。
if (ncol(group) < 2) {
  stop(paste("错误: 表型文件 '", phenotypic_data_file, "' 至少需要两列 (第一列样本名, 第二列分组信息)。"))
}
colnames(group)[1] <- "sample"
colnames(group)[2] <- "group"
phenotype_col <- "group" # 将分组列的名称存储在变量中，方便后续引用
cat(paste("INFO: 使用表型文件第一列作为样本名，第二列作为分组信息。\n", sep=""))


# --- 4. 数据处理与转换 ---
# 过滤基因：保留在至少75%的样本中表达量大于0的基因。
# 这是一个常用的过滤步骤，可以去除在大部分样本中都未表达的基因，减少噪音。
# `data > 0` 会生成一个布尔矩阵，`rowSums` 对每行求和（即计算每个基因为TRUE的样本数）。
keep <- rowSums(data > 0) >= floor(0.75 * ncol(data))
filtered_fpkm_matrix <- data[keep, ]

# 对FPKM值进行log2(x+1)转换。
# `+1`是为了避免对0取对数（log(0)为负无穷）。
# log转换可以使数据分布更接近正态分布，减小极端值的影响，使可视化效果更好。
data_log2 <- log2(filtered_fpkm_matrix + 1)

# 转置数据框，使行为样本，列为基因，然后转换为data.frame。
data_log2_t <- as.data.frame(t(data_log2))
# 将行名（即样本名）添加为一个新的列，以便后续与分组信息合并。
data_log2_t$sample <- rownames(data_log2_t)

# 使用 `merge` 函数，通过共同的 'sample' 列，将表达数据和分组信息合并在一起。
merged_data <- merge(data_log2_t, group, by = "sample")

# 将宽数据转换为长数据。这是`ggplot2`绘图的推荐数据格式。
# 原始数据是“宽”的：每个基因一列。转换后是“长”的：所有基因名在一列，所有表达值在另一列。
long_data <- merged_data %>%
pivot_longer(
    cols = -c(sample, !!phenotype_col), # 选择除 'sample' 和分组列之外的所有列（即所有基因列）进行转换。
    names_to = "Gene",                  # 新生成的、包含所有基因名的列，其列名称为 "Gene"。
    values_to = "Expression"            # 新生成的、包含所有表达值的列，其列名称为 "Expression"。
  )

# 确保分组列是因子（factor）类型，并使用`levels`参数设置好组别顺序。
# 这可以控制x轴上的组别显示顺序，而不是按默认的字母顺序。
long_data[[phenotype_col]] <- factor(long_data[[phenotype_col]], levels = c(control_group_name, experimental_group_name))
# 过滤掉分组信息为NA的行，防止绘图时出错或出现不必要的"NA"组。
long_data <- long_data %>% filter(!is.na(.data[[phenotype_col]]))


# --- 5. 结果可视化 ---
P <- ggplot(long_data, aes(x = .data[[phenotype_col]], y = Expression, fill = .data[[phenotype_col]])) +
  # 绘制小提琴图，它结合了箱线图和密度图的特点，可以很好地展示数据的分布形状和密度。
  # `trim = FALSE`: 不修剪小提琴的尾部，完整展示数据范围。
  # `alpha = 0.7`: 设置透明度，使图形更好看。
  geom_violin(trim = FALSE, alpha = 0.7, color = "black") +

  # 手动设置填充颜色，`#F8766D`是红色系，`#7CAE00`是绿色系，是ggplot2的默认色。
  scale_fill_manual(values = c("#F8766D", "#7CAE00")) +

  # 在小提琴图内部叠加一个箱线图，以清晰地显示中位数、四分位数等统计信息。
  # `width = 0.1`: 设置箱线图的宽度，使其比小提琴窄。
  # `outlier.shape = NA`: 不显示异常值点，因为小提琴图本身已经展示了数据的完整分布。
  geom_boxplot(width = 0.1, position = position_dodge(0.9), alpha = 0.3, outlier.shape = NA) +

  # 设置坐标轴和图表标题。
  xlab("Group") +
  ylab(expression(paste("FPKM (log"["2"], "+1)"))) + # 使用 expression 获得更美观的log2格式
  ggtitle("Gene Expression Distribution Between Groups") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", color = "grey20", family = base_font_family), # X轴刻度文字
    axis.text.y = element_text(size = 12, color = "grey20", family = base_font_family), # Y轴刻度文字
    axis.title.x = element_text(size = 14, face = "bold", color = "black", family = base_font_family), # X轴标题
    axis.title.y = element_text(size = 14, face = "bold", color = "black", family = base_font_family), # Y轴标题
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, family = base_font_family), # 图表主标题居中、加粗
    legend.text = element_text(family = base_font_family), # 如果图例有文字，也一并修改
    legend.position = "right", # 图例位置
    legend.title = element_blank() # 不显示图例标题
  )

# --- 6. 保存图片 ---
ggsave(output_file, plot = P, device = cairo_pdf, width = 8, height = 8, units = "in")

# 绘图完成后可以关闭 showtext
showtext_auto(FALSE)

cat(paste("INFO: Plotting complete. PDF saved to '", output_file, "'\n", sep=""))
