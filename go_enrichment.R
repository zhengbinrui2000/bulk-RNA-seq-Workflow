#!/usr/bin/env Rscript

# =============================================================================
# 脚本名称: GO 富集分析与可视化 (为适配Nextflow明确输出而修改)
# 功能描述:
#   该脚本承接差异表达分析的结果，进行GO功能富集分析。
#   1. 根据物种名称加载对应的Bioconductor注释数据库。
#   2. 读取差异基因列表，转换基因为ENTREZID。
#   3. 对GO的三个本体（BP, CC, MF）进行富集分析。
#   4. 将每个GO本体的富集结果保存到由Nextflow指定的CSV文件中。
#   5. 如果某个本体没有富集结果，则创建一个空的CSV文件作为占位符。
#   6. 合并富集条目并使用ggplot2绘制条形图，保存到Nextflow指定的文件中。
#   7. 如果没有任何富集结果，则创建一个空的PDF文件作为占位符。
#
# 使用方法 (在Nextflow中被自动调用):
#   Rscript this_script.R <deg_file> <species> <bp_out.csv> <cc_out.csv> <mf_out.csv> <plot_out.pdf>
# =============================================================================


# --- 1. 环境准备与参数解析 ---
rm(list = ls())

library(tidyverse)        # 提供一套完整的数据科学工具，核心是dplyr和ggplot2。
library(AnnotationHub)    # 用于访问和管理Bioconductor的注释资源（如此处用于加载本地数据库）。
library(AnnotationDbi)    # 提供处理Bioconductor注释数据库的核心功能（如select, mapIds）。
library(clusterProfiler)  # GO和KEGG富集分析的核心包。
library(grDevices)   # 用于pdf文字嵌入
library(showtext)    # 用于字体管理
font_add("Arial", "arial.ttf") # 加载Arial字体
showtext_auto() # 自动激活showtext，让它接管后续所有图形设备的文本渲染。
base_font_family <- "Arial" # 定义base_font_family变量为Arial


# --- 命令行参数解析 ---
args <- commandArgs(trailingOnly = TRUE)
deg_file     <- args[1] # 输入1: 差异表达基因结果文件
species      <- args[2] # 输入2: 物种名称
bp_csv_out   <- args[3] # 输入3: BP 结果的输出路径 (.csv)
cc_csv_out   <- args[4] # 输入4: CC 结果的输出路径 (.csv)
mf_csv_out   <- args[5] # 输入5: MF 结果的输出路径 (.csv)
plot_pdf_out <- args[6] # 输入6: 最终图表的输出路径 (.pdf)


# --- 根据物种名称，动态加载对应的注释包 ---
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
    stop("未知的或不支持的物种!")
}

# --- 2. 数据加载与预处理 ---
deg_filter <- function(df) {
  df[df$change != "stable", ]
}
DEG_edgeR <- read.csv(deg_file, row.names = 1)
DEG_edgeR <- rownames_to_column(DEG_edgeR, var = "ID")
DEGs <- deg_filter(DEG_edgeR)
Converted_ID <- bitr(DEGs$ID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = Orgdb)


# --- 3. GO富集分析 (BP, CC, MF) ---
ego_BP <- tryCatch({ enrichGO( gene = Converted_ID$ENTREZID, OrgDb = Orgdb, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE) }, error = function(e) { message("GO BP 分析失败: ", e$message); return(NULL) })
ego_CC <- tryCatch({ enrichGO( gene = Converted_ID$ENTREZID, OrgDb = Orgdb, keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE) }, error = function(e) { message("GO CC 分析失败: ", e$message); return(NULL) })
ego_MF <- tryCatch({ enrichGO( gene = Converted_ID$ENTREZID, OrgDb = Orgdb, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE) }, error = function(e) { message("GO MF 分析失败: ", e$message); return(NULL) })

# 对每个分类，如果结果存在，则写入结果；否则，写入一个空文件以防止nextflow报错。
if (!is.null(ego_BP) && nrow(as.data.frame(ego_BP)) > 0) {
    write.csv(as.data.frame(ego_BP), file = bp_csv_out, row.names = FALSE)
} else {
    write.csv(data.frame(), file = bp_csv_out) # 创建一个带表头的空文件
}

if (!is.null(ego_CC) && nrow(as.data.frame(ego_CC)) > 0) {
    write.csv(as.data.frame(ego_CC), file = cc_csv_out, row.names = FALSE)
} else {
    write.csv(data.frame(), file = cc_csv_out)
}

if (!is.null(ego_MF) && nrow(as.data.frame(ego_MF)) > 0) {
    write.csv(as.data.frame(ego_MF), file = mf_csv_out, row.names = FALSE)
} else {
    write.csv(data.frame(), file = mf_csv_out)
}


# --- 3. 结果整合与处理，为绘图做准备 ---
results_list <- list()
display_number <- c(10, 10, 10)

if (!is.null(ego_BP) && nrow(ego_BP) > 0) {
    ego_result_BP <- as.data.frame(ego_BP)[1:min(display_number[1], nrow(ego_BP)), ]
    ego_result_BP$type <- "biological process"
    results_list$BP <- ego_result_BP
}
if (!is.null(ego_CC) && nrow(ego_CC) > 0) {
    ego_result_CC <- as.data.frame(ego_CC)[1:min(display_number[2], nrow(ego_CC)), ]
    ego_result_CC$type <- "cellular component"
    results_list$CC <- ego_result_CC
}
if (!is.null(ego_MF) && nrow(ego_MF) > 0) {
    ego_result_MF <- as.data.frame(ego_MF)[1:min(display_number[3], nrow(ego_MF)), ]
    ego_result_MF$type <- "molecular function"
    results_list$MF <- ego_result_MF
}


# --- 4. 结果可视化与保存 ---
if (length(results_list) > 0) {
    go_enrich_df_full <- do.call(rbind, results_list)
    go_enrich_df <- data.frame(
        ID = go_enrich_df_full$ID,
        Description = go_enrich_df_full$Description,
        GeneNumber = go_enrich_df_full$Count,
        type = factor(go_enrich_df_full$type, levels = c("biological process", "cellular component", "molecular function"))
    )
    for(i in 1:nrow(go_enrich_df)){
        description_splite <- strsplit(go_enrich_df$Description[i], split = " ")
        word_count <- min(5, length(description_splite[[1]]))
        description_collapse <- paste(description_splite[[1]][1:word_count], collapse = " ")
        go_enrich_df$Description[i] <- description_collapse
        go_enrich_df$Description <- gsub(pattern = "NA", "", go_enrich_df$Description)
    }
    go_enrich_df$type_order <- factor(go_enrich_df$Description, levels = rev(go_enrich_df$Description))
    COLS <- c("biological process" = "#66C3A5", "cellular component" = "#8DA1CB", "molecular function" = "#FD8D62")

    p <- ggplot(data = go_enrich_df, aes(x = type_order, y = GeneNumber, fill = type)) +
        geom_bar(stat = "identity", width = 0.8) +
        scale_fill_manual(values = COLS, name = "Ontology", drop = TRUE) +
        coord_flip() +
        xlab("GO term") +
        ylab("Gene_Number") +
        labs(title = "The Most Enriched GO Terms") +
        theme_bw() +
        # 在 theme() 中为所有文本元素指定字体
        theme(
            axis.text.y = element_text(size = 10, face = "bold", family = base_font_family),
            axis.text.x = element_text(family = base_font_family),
            axis.title = element_text(family = base_font_family),
            plot.title = element_text(hjust = 0.5, face = "bold", family = base_font_family),
            legend.title = element_text(face = "bold", family = base_font_family),
            legend.text = element_text(family = base_font_family)
        )

    # 使用ggsave保存为PDF
    ggsave(plot_pdf_out, plot = p, device = cairo_pdf, width = 10, height = 8, units = "in")

} else {
    # 如果没有结果，则在指定的绘图输出路径创建一个空文件。
    print("没有发现显著富集的GO条目。正在创建空的占位符图表文件。")
    file.create(plot_pdf_out)
}

# 脚本结束时，可选地关闭showtext的自动渲染功能
showtext_auto(FALSE)
