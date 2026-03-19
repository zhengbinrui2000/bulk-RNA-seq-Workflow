#!/usr/bin/env Rscript

# =============================================================================
# 脚本名称: KEGG 通路富集分析与可视化 (为适配Nextflow明确输出而修改)
# 功能描述:
#   该脚本对差异表达基因进行KEGG通路富集分析，并使用GOplot包创建圈图。
#   1. 根据物种名加载注释库并匹配KEGG物种代码。
#   2. 读取差异基因列表，转换基因为ENTREZID。
#   3. 使用clusterProfiler进行KEGG富集分析。
#   4. 将富集结果和基因的logFC信息整合，生成并保存KEGG圈图到指定文件。
#   5. 如果没有富集结果，则创建空的占位符文件。
#
# 使用方法 (在Nextflow中被自动调用):
#   Rscript <this_script.R> <deg_file> <species> <results_out.csv> <plot_out.pdf>
# =============================================================================

# --- 1. 环境准备与参数解析 ---
rm(list = ls())

library(tidyverse)        # 用于数据处理和可视化 (包含了dplyr, tidyr, ggplot2等)
library(AnnotationHub)    # 用于获取基因注释数据库
library(AnnotationDbi)    # 用于处理注释数据
library(clusterProfiler)  # 用于富集分析的核心包
library(GOplot)           # 用于绘制圈图
library(grDevices)        # 用于pdf文字嵌入
library(showtext)         # 用于字体管理
font_add("Arial", "arial.ttf") # 加载Arial字体
showtext_auto() # 自动激活showtext，让它接管后续所有图形设备的文本渲染。
par(family = "Arial") # 为R基础绘图系统设置字体。par()会影响当前图形设备后续的所有绘图操作。

# ===================================================================
# Section 0: 命令与环境设置
# ===================================================================

# --- 命令行参数解析 ---
args <- commandArgs(trailingOnly = TRUE)

deg_file        <- args[1] # 输入1: 差异基因文件
species         <- args[2] # 输入2: 物种名
results_csv_out <- args[3] # 输入3: 结果表格的输出路径 (.csv)
plot_pdf_out    <- args[4] # 输入4: 圈图的输出路径 (.pdf)

# ===================================================================
# Section 2: 数据加载与物种选择
# ===================================================================
# --- 根据物种名选择注释数据库 (OrgDb) ---
if (species == "Homo_sapiens") {
    library(org.Hs.eg.db)
    Orgdb <- org.Hs.eg.db
} else if (species == "Mus_musculus") {
    library(org.Mm.eg.db)
    Orgdb <- org.Mm.eg.db
} else if (species == "Arabidopsis_thaliana") {
    library(org.At.tair.db)
    Orgdb <- org.At.tair.db
} else if (species == "Anopheles") {
    library(org.Ag.eg.db)
    Orgdb <- org.Ag.eg.db
} else if (species == "Bos_taurus") {
    library(org.Bt.eg.db)
    Orgdb <- org.Bt.eg.db
} else if (species == "Caenorhabditis_elegans") {
    library(org.Ce.eg.db)
    Orgdb <- org.Ce.eg.db
} else if (species == "Canis_familiaris") {
    library(org.Cf.eg.db)
    Orgdb <- org.Cf.eg.db
} else if (species == "Drosophila_melanogaster") {
    library(org.Dm.eg.db)
    Orgdb <- org.Dm.eg.db
} else if (species == "Danio_rerio") {
    library(org.Dr.eg.db)
    Orgdb <- org.Dr.eg.db
} else if (species == "Escherichia_coli_K12") {
    library(org.EcK12.eg.db)
    Orgdb <- org.EcK12.eg.db
} else if (species == "Escherichia_coli_Sakai") {
    library(org.EcSakai.eg.db)
    Orgdb <- org.EcSakai.eg.db
} else if (species == "Gallus_gallus") {
    library(org.Gg.eg.db)
    Orgdb <- org.Gg.eg.db
} else if (species == "Macaca_mulatta") {
    library(org.Mmu.eg.db)
    Orgdb <- org.Mmu.eg.db
} else if (species == "Pan_troglodytes") {
    library(org.Pt.eg.db)
    Orgdb <- org.Pt.eg.db
} else if (species == "Rattus_norvegicus") {
    library(org.Rn.eg.db)
    Orgdb <- org.Rn.eg.db
} else if (species == "Saccharomyces_cerevisiae") {
    library(org.Sc.sgd.db)
    Orgdb <- org.Sc.sgd.db
} else if (species == "Sus_scrofa") {
    library(org.Ss.eg.db)
    Orgdb <- org.Ss.eg.db
} else if (species == "Xenopus_laeviss") {
    library(org.Xl.eg.db)
    Orgdb <- org.Xl.eg.db
} else if (species == "Oryza_sativa") {
    Orgdb <- loadDb("/opt/annotationhub/org.Oryza_sativa_Japonica_Group.eg.sqlite")
} else if (species == "Zea_mays") {
    Orgdb <- loadDb("/opt/annotationhub/org.Zea_mays.eg.sqlite")
} else {
    stop("未知的物种或缺少对应的注释包!")
}

# --- 根据物种名选择KEGG物种代码 ---
kegg_organism_map <- c(
    "Homo_sapiens" = "hsa", "Mus_musculus" = "mmu", "Rattus_norvegicus" = "rno",
    "Drosophila_melanogaster" = "dme", "Caenorhabditis_elegans" = "cel",
    "Saccharomyces_cerevisiae" = "sce", "Arabidopsis_thaliana" = "ath",
    "Danio_rerio" = "dre", "Bos_taurus" = "bta", "Gallus_gallus" = "gga",
    "Xenopus_laeviss" = "xla", "Anopheles" = "aga", "Escherichia_coli_K12" = "eco",
    "Pan_troglodytes" = "ptr", "Canis_familiaris" = "cfa", "Sus_scrofa" = "ssc",
    "Oryza_sativa" = "osa", "Zea_mays" = "zma"
)
organism <- kegg_organism_map[species]
if (is.na(organism)) {
    stop("未知的物种或缺少对应的KEGG物种代码映射!")
}

# --- 加载并筛选差异表达基因(DEG) ---
deg_filter <- function(df) {
  df[df$change != "stable", ]
}
DEG_edgeR <- read.csv(deg_file, row.names = 1)
DEG_edgeR <- rownames_to_column(DEG_edgeR, var = "ID")
DEGs_filtered <- deg_filter(DEG_edgeR)


# ===================================================================
# Section 3: ID转换与数据整合
# ===================================================================
Converted_IDs <- bitr(DEGs_filtered$ID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = Orgdb)
DEGs_for_plot <- DEGs_filtered %>%
  left_join(Converted_IDs, by = c("ID" = "SYMBOL")) %>%
  na.omit() %>%
  dplyr::select(-ID) %>%
  dplyr::rename(ID = ENTREZID) %>%
  dplyr::select(ID, logFC, change, PValue, FDR)


# ===================================================================
# Section 4: KEGG 富集分析
# ===================================================================
kegg_enrich_result <- enrichKEGG(
  gene         = DEGs_for_plot$ID,
  organism     = organism,
  keyType      = "kegg",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
ekp_df <- as.data.frame(kegg_enrich_result)


# ===================================================================
# Section 5: 结果处理与可视化 (GOplot)
# ===================================================================
if (nrow(ekp_df) > 0) {

  print("成功找到KEGG富集结果，正在进行可视化...")

  # --- 为`GOplot`准备富集结果数据 ---
  ekp_for_plot <- ekp_df %>%
    dplyr::rename(
      term     = "Description",
      adj_pval = "p.adjust",
      genes    = "geneID"
    ) %>%
    mutate(genes = gsub("/", ",", genes))

  circ_kegg <- circle_dat(ekp_for_plot, DEGs_for_plot)

  # --- 使用从参数接收的文件路径保存富集结果表格 ---
  write.csv(ekp_for_plot, file = results_csv_out, row.names = FALSE)

  # --- 绘制并保存圈图 (PDF格式) ---
  circ_kegg$term <- sapply(circ_kegg$term, function(x) paste(strwrap(x, width = 40), collapse = "\n"))
  nsub_value <- min(10, nrow(ekp_for_plot))

  # 使用cairo_pdf设备保存为PDF，以确保字体正确嵌入
  cairo_pdf(filename = plot_pdf_out, width = 12, height = 9)
  GOCircle(
    circ_kegg,
    nsub       = nsub_value,
    title      = "KEGG Enrichment Circle Plot",
    label.size = 5,
    zsc.col    = c('#E4B112', '#CCDDAE', '#8CBA54'),
    lfc.col    = c('#8CBA54', '#E4B112')
  )
  dev.off()

  print(paste("KEGG富集分析可视化完成，文件已保存。"))

} else {
  print("未找到任何显著富集的KEGG通路。正在创建空的占位符文件。")

  # 创建空的CSV结果文件
  write.csv(data.frame(), file = results_csv_out)
  # 创建空的PDF图表文件
  file.create(plot_pdf_out)
}

# 脚本结束时，可选地关闭showtext的自动渲染功能
showtext_auto(FALSE)
