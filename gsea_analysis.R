# 清空工作环境并加载所需包
rm(list = ls())
library(AnnotationDbi)
library(clusterProfiler)
library(tidyverse)
library(enrichplot)
library(grDevices)     # 用于pdf文字嵌入
library(showtext)      # 用于字体管理

# 加载并设置全局字体
font_add("Arial", "arial.ttf") # 确保容器或系统中有这个字体文件
showtext_auto() # 自动激活showtext，让它接管后续所有图形设备的文本渲染

# 命令行参数解析
args <- commandArgs(trailingOnly = TRUE)
deg_file          <- args[1] # 输入1: 差异基因文件
species           <- args[2] # 输入2: 物种名称
go_gsea_out       <- args[3] # 输入3: 基于GO数据库的GSEA结果(.csv)
go_gsea_pdf_out   <- args[4] # 输入4: 基于GO数据库的GSEA可视化结果(.pdf)
kegg_gsea_out     <- args[5] # 输入5: 基于KEGG数据库的GSEA结果(.csv)
kegg_gsea_pdf_out <- args[6] # 输入6: 基于KEGG数据库的GSEA可视化结果(.pdf)

# 根据物种选择适当的注释包
if (species == "Homo_sapiens") {
    library(org.Hs.eg.db)  # 加载人类注释包
    Orgdb <- org.Hs.eg.db
} else if (species == "Mus_musculus") {
    library(org.Mm.eg.db)  # 加载小鼠注释包
    Orgdb <- org.Mm.eg.db
} else if (species == "Arabidopsis_thaliana") {
    library(org.At.tair.db)  # 加载拟南芥注释包
    Orgdb <- org.At.tair.db
} else if (species == "Anopheles") {
    library(org.Ag.eg.db)  # 加载按蚊注释包
    Orgdb <- org.Ag.eg.db
} else if (species == "Bos_taurus") {
    library(org.Bt.eg.db)  # 加载牛注释包
    Orgdb <- org.Bt.eg.db
} else if (species == "Caenorhabditis_elegans") {
    library(org.Ce.eg.db)  # 加载秀丽隐杆线虫注释包
    Orgdb <- org.Ce.eg.db
} else if (species == "Canis_familiaris") {
    library(org.Cf.eg.db)  # 加载犬注释包
    Orgdb <- org.Cf.eg.db
} else if (species == "Drosophila_melanogaster") {
    library(org.Dm.eg.db)  # 加载苍蝇注释包
    Orgdb <- org.Dm.eg.db
} else if (species == "Danio_rerio") {
    library(org.Dr.eg.db)  # 加载斑马鱼注释包
    Orgdb <- org.Dr.eg.db
} else if (species == "Escherichia_coli_K12") {
    library(org.EcK12.eg.db)  # 加载大肠杆菌K12注释包
    Orgdb <- org.EcK12.eg.db
} else if (species == "Escherichia_coli_Sakai") {
    library(org.EcSakai.eg.db)  # 加载大肠杆菌Sakai注释包
    Orgdb <- org.EcSakai.eg.db
} else if (species == "Gallus_gallus") {
    library(org.Gg.eg.db)  # 加载鸡注释包
    Orgdb <- org.Gg.eg.db
} else if (species == "Macaca_mulatta") {
    library(org.Mmu.eg.db)  # 加载恒河猴注释包
    Orgdb <- org.Mmu.eg.db
} else if (species == "Pan_troglodytes") {
    library(org.Pt.eg.db)  # 加载黑猩猩注释包
    Orgdb <- org.Pt.eg.db
} else if (species == "Rattus_norvegicus") {
    library(org.Rn.eg.db)  # 加载大鼠注释包
    Orgdb <- org.Rn.eg.db
} else if (species == "Saccharomyces_cerevisiae") {
    library(org.Sc.sgd.db)  # 加载酵母注释包
    Orgdb <- org.Sc.sgd.db
} else if (species == "Sus_scrofa") {
    library(org.Ss.eg.db)  # 加载猪注释包
    Orgdb <- org.Ss.eg.db
} else if (species == "Xenopus_laeviss") {
    library(org.Xl.eg.db)  # 加载爪蟾注释包
    Orgdb <- org.Xl.eg.db
} else if (species == "Oryza_sativa") {
    Orgdb <- loadDb("/opt/annotationhub/org.Oryza_sativa_Japonica_Group.eg.sqlite")  # 水稻注释包（SQLite格式）
} else if (species == "Zea_mays") {
    Orgdb <- loadDb("/opt/annotationhub/org.Zea_mays.eg.sqlite")  # 玉米注释包（SQLite格式）
} else {
    stop("Unknown species!")
}

# 根据物种名选择对应的KEGG ID
if (species == "Homo_sapiens") {
  organism <- "hsa"
} else if (species == "Mus_musculus") {
  organism <- "mmu"
} else if (species == "Arabidopsis_thaliana") {
  organism <- "ath"
} else if (species == "Anopheles") {
  organism <- "ag"
} else if (species == "Bos_taurus") {
  organism <- "bta"
} else if (species == "Caenorhabditis_elegans") {
  organism <- "cel"
} else if (species == "Canis_familiaris") {
  organism <- "cfa"
} else if (species == "Drosophila_melanogaster") {
  organism <- "dme"
} else if (species == "Danio_rerio") {
  organism <- "dre"
} else if (species == "Escherichia_coli_K12") {
  organism <- "eco"
} else if (species == "Escherichia_coli_Sakai") {
  organism <- "ecs"
} else if (species == "Gallus_gallus") {
  organism <- "gga"
} else if (species == "Macaca_mulatta") {
  organism <- "mmu"
} else if (species == "Pan_troglodytes") {
  organism <- "ptr"
} else if (species == "Rattus_norvegicus") {
  organism <- "rno"
} else if (species == "Saccharomyces_cerevisiae") {
  organism <- "sce"
} else if (species == "Sus_scrofa") {
  organism <- "ssc"
} else if (species == "Xenopus_laeviss") {
  organism <- "xla"
} else if (species == "Oryza_sativa") {
  organism <- "osa"
} else if (species == "Zea_mays") {
  organism <- "zma"
} else {
  stop("Unknown species or missing KEGG ID mapping!")
}

deg_filter <- function(df) {
  df[df$change != "stable", ]
}

DEG_edgeR <- read.csv(deg_file, row.names = 1)
DEG_edgeR <- rownames_to_column(DEG_edgeR, var = "ID")
DEGs <- deg_filter(DEG_edgeR)

# 'select()' returned 1:1 mapping between keys and columns
Converted_ID <- bitr(DEGs$ID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = Orgdb)

merged_df <- merge(DEGs, Converted_ID, by.x = "ID", by.y = "SYMBOL")

# --- 数据准备和格式转换  ---
# 1. 创建只包含 ENTREZID 和 logFC 的数据框
new_df <- merged_df[, c("ENTREZID", "logFC")]

# 2. 创建 GSEA 所需的已排序、已命名的向量
geneList <- new_df$logFC
names(geneList) <- new_df$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)
# 去除重复的ENTREZID,保留第一个
geneList <- geneList[!duplicated(names(geneList))]

# --- 运行 GO GSEA 分析 ---
go <- gseGO(geneList = geneList,
            OrgDb = Orgdb,
            ont = 'BP',
            pvalueCutoff = 0.05,
            minGSSize = 10,
            maxGSSize = 500,
            pAdjustMethod = 'BH',
            keyType = "ENTREZID"
)

# --- 保存GO GSEA结果 ---
# 检查是否有富集结果
if (!is.null(go) && nrow(go@result) > 0) {
    # 保存表格
    write.csv(go@result, file = go_gsea_out, row.names = FALSE)

    # 开启PDF设备，并指定字体
    pdf(go_gsea_pdf_out, height = 5, width = 6)
    # 绘制最显著的一条通路图
    print(gseaplot2(go, geneSetID = 1, title = go$Description[1]))
    # 关闭PDF设备
    dev.off()
} else {
    # 如果没有结果，创建一个空的CSV和一个提示性的PDF
    write.csv(data.frame(), file = go_gsea_out)
    pdf(go_gsea_pdf_out)
    plot.new()
    text(0.5, 0.5, "No significant GO enrichment found.")
    dev.off()
}


# --- 运行 KEGG GSEA 分析 ---
# **已修正**: 将硬编码的 "zma" 替换为从上方逻辑中获取的 'organism' 变量
kk <- gseKEGG(geneList,
              organism = organism, 
              pvalueCutoff = 0.05,
              pAdjustMethod = 'BH',
              minGSSize = 10,
              maxGSSize = 500,
              keyType = "kegg" # 通常gseKEGG使用ENTREZID, keyType为kegg
)

# --- 保存KEGG GSEA结果 ---
# 检查是否有富集结果
if (!is.null(kk) && nrow(kk@result) > 0) {
    # 保存表格
    write.csv(kk@result, file = kegg_gsea_out, row.names = FALSE)
    
    # 开启PDF设备，并指定字体
    pdf(kegg_gsea_pdf_out, height = 5, width = 6)
    # 绘制最显著的一条通路图
    print(gseaplot2(kk, geneSetID = 1, title = kk$Description[1]))
    # 关闭PDF设备
    dev.off()
} else {
    # 如果没有结果，创建一个空的CSV和PDF
    write.csv(data.frame(), file = kegg_gsea_out)
    pdf(kegg_gsea_pdf_out)
    plot.new()
    text(0.5, 0.5, "No significant KEGG enrichment found.")
    dev.off()
}
