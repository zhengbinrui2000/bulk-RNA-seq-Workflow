import yaml
import os
import pandas as pd
from snakemake.io import glob_wildcards

# 加载配置文件
with open("config.yaml", "r") as config_file:
    config = yaml.safe_load(config_file)

# 定义路径（替换 {PROJECT_DIR} 变量）
SPECIES = config["SPECIES"]
PROJECT_DIR = config["PROJECT_DIR"]
GENOTYPIC_DATA = config["GENOTYPIC_DATA"].format(PROJECT_DIR=PROJECT_DIR)
PHENOTYPIC_DATA = config["PHENOTYPIC_DATA"].format(PROJECT_DIR=PROJECT_DIR)
LOGS_DIR = config["LOGS_DIR"].format(PROJECT_DIR=PROJECT_DIR)
UPSTREAM_ANALYSIS = config["UPSTREAM_ANALYSIS"].format(PROJECT_DIR=PROJECT_DIR)
DOWNSTREAM_ANALYSIS = config["DOWNSTREAM_ANALYSIS"].format(PROJECT_DIR=PROJECT_DIR)
ANNOTATION_FILE = config["ANNOTATION_FILE"].format(PROJECT_DIR=PROJECT_DIR)
GENOME_FILE = config["GENOME_FILE"].format(PROJECT_DIR=PROJECT_DIR)
SINGULARITY_IMAGE = config["SINGULARITY_IMAGE"].format(PROJECT_DIR=PROJECT_DIR)

# 读取表型数据并获取样本ID列表
phenotypic_data = pd.read_csv(f"{PHENOTYPIC_DATA}", index_col=0)
SAMPLE_IDS = phenotypic_data.index.tolist()

# 打印调试信息
print(f"SAMPLE_IDS: {SAMPLE_IDS}")

# 主规则
rule all:
    input:
        # 质量控制结果
        expand(f"{UPSTREAM_ANALYSIS}/fastqc/original/{{sample}}.R1_fastqc.html", sample=SAMPLE_IDS),
        expand(f"{UPSTREAM_ANALYSIS}/fastqc/filtered/{{sample}}.R1.filtered_fastqc.html", sample=SAMPLE_IDS),
        expand(f"{UPSTREAM_ANALYSIS}/cleaned/{{sample}}.R1.filtered.fastq.gz", sample=SAMPLE_IDS),
        f"{UPSTREAM_ANALYSIS}/multiqc_result/fastp/multiqc_report.html",

        # 索引和注释
        expand(f"{PROJECT_DIR}/references/hisat2_index/genome_index.{{suffix}}", suffix=["1.ht2", "2.ht2", "3.ht2", "4.ht2", "5.ht2", "6.ht2", "7.ht2", "8.ht2"]),

        # 比对结果
        expand(f"{UPSTREAM_ANALYSIS}/alignment/bam_file/{{sample}}.sorted.bam", sample=SAMPLE_IDS),
        expand(f"{UPSTREAM_ANALYSIS}/alignment/stats/{{sample}}.stats.txt", sample=SAMPLE_IDS),
        f"{UPSTREAM_ANALYSIS}/multiqc_result/alignment/multiqc_report.html",

        # 定量结果
        expand(f"{UPSTREAM_ANALYSIS}/quantification/stringtie/{{sample}}.gtf", sample=SAMPLE_IDS),
        expand(f"{UPSTREAM_ANALYSIS}/quantification/stringtie/{{sample}}.abundance.tsv", sample=SAMPLE_IDS),
        f"{PROJECT_DIR}/references/fixed_output.gtf",
        f"{UPSTREAM_ANALYSIS}/quantification/featurecounts/gene_count_matrix.txt",
        f"{UPSTREAM_ANALYSIS}/quantification/featurecounts/transcript_count_matrix.txt",
        f"{UPSTREAM_ANALYSIS}/quantification/stringtie/normalized_data/fpkm_matrix.csv",
        f"{UPSTREAM_ANALYSIS}/quantification/stringtie/normalized_data/tpm_matrix.csv",

        # 下游分析结果
        f"{DOWNSTREAM_ANALYSIS}/TGE_Analysis/gene_expression_violin.pdf",
        f"{DOWNSTREAM_ANALYSIS}/GEC_Analysis/correlation_heatmap.pdf",
        f"{DOWNSTREAM_ANALYSIS}/GEC_Analysis/PCA_Plot.pdf",
        f"{DOWNSTREAM_ANALYSIS}/DEG_Analysis_edgeR/DEG_edgeR.csv",
        f"{DOWNSTREAM_ANALYSIS}/DEG_Analysis_edgeR/edgeR_volcano.pdf",
        f"{DOWNSTREAM_ANALYSIS}/DEG_Analysis_edgeR/edgeR_heatmap.pdf",
        expand(f"{DOWNSTREAM_ANALYSIS}/GO_Enrichment/enrichment_result_{{ont}}.csv", ont=["BP", "CC", "MF"]),
        expand(f"{DOWNSTREAM_ANALYSIS}/GO_Enrichment/GOCircle_Plot_{{ont}}.pdf", ont=["BP", "CC", "MF"]),
        f"{DOWNSTREAM_ANALYSIS}/GO_Enrichment/extracted_data.csv",
        f"{DOWNSTREAM_ANALYSIS}/GO_Enrichment/DEGs_with_ENTREZID.csv",
        f"{DOWNSTREAM_ANALYSIS}/KEGG_Enrichment/Kegg_enrichment_result.txt",
        f"{DOWNSTREAM_ANALYSIS}/KEGG_Enrichment/KeggCircle_Plot.pdf",
        f"{DOWNSTREAM_ANALYSIS}/GSEA_Enrichment/GSEA_GO_Result.csv",
        f"{DOWNSTREAM_ANALYSIS}/GSEA_Enrichment/GSEA_GO_CirclePlot.pdf",
        f"{DOWNSTREAM_ANALYSIS}/GSEA_Enrichment/GSEA_KEGG_Result.csv",
        f"{DOWNSTREAM_ANALYSIS}/GSEA_Enrichment/GSEA_KEGG_CirclePlot.pdf"

# 规则1：FastQC 原始质量控制
rule fastqc_original:
    input:
        R1 = f"{GENOTYPIC_DATA}/{{sample}}.R1.fastq.gz",
        R2 = f"{GENOTYPIC_DATA}/{{sample}}.R2.fastq.gz"
    output:
        R1_fastqc = f"{UPSTREAM_ANALYSIS}/fastqc/original/{{sample}}.R1_fastqc.html",
        R2_fastqc = f"{UPSTREAM_ANALYSIS}/fastqc/original/{{sample}}.R2_fastqc.html"
    log:
        f"{LOGS_DIR}/upstream/fastqc_original/{{sample}}.log"
    params:
        output_dir = f"{UPSTREAM_ANALYSIS}/fastqc/original"
    shell:
        """
        apptainer exec {SINGULARITY_IMAGE} fastqc {input.R1} {input.R2} -o {params.output_dir} &>> {log}
        """

# 规则2：使用 Fastp 进行数据清洗
rule fastp:
    input:
        R1 = f"{GENOTYPIC_DATA}/{{sample}}.R1.fastq.gz",
        R2 = f"{GENOTYPIC_DATA}/{{sample}}.R2.fastq.gz"
    output:
        R1_filtered = f"{UPSTREAM_ANALYSIS}/cleaned/{{sample}}.R1.filtered.fastq.gz",
        R2_filtered = f"{UPSTREAM_ANALYSIS}/cleaned/{{sample}}.R2.filtered.fastq.gz",
        html_report = f"{UPSTREAM_ANALYSIS}/fastp/{{sample}}.fastp.html",
        json_report = f"{UPSTREAM_ANALYSIS}/fastp/{{sample}}.fastp.json"
    log:
        f"{LOGS_DIR}/upstream/fastp/{{sample}}.log"
    shell:
        """
        apptainer exec {SINGULARITY_IMAGE} fastp -i {input.R1} -I {input.R2} \
            -o {output.R1_filtered} -O {output.R2_filtered} \
            --html {output.html_report} --json {output.json_report} \
            --length_required 35 --qualified_quality_phred 20 &>> {log}
        """

# 规则3：FastQC 过滤后质量控制
rule fastqc_filtered:
    input:
        R1 = f"{UPSTREAM_ANALYSIS}/cleaned/{{sample}}.R1.filtered.fastq.gz",
        R2 = f"{UPSTREAM_ANALYSIS}/cleaned/{{sample}}.R2.filtered.fastq.gz"
    output:
        R1_fastqc = f"{UPSTREAM_ANALYSIS}/fastqc/filtered/{{sample}}.R1.filtered_fastqc.html",
        R2_fastqc = f"{UPSTREAM_ANALYSIS}/fastqc/filtered/{{sample}}.R2.filtered_fastqc.html"
    log:
        f"{LOGS_DIR}/upstream/fastqc_filtered/{{sample}}.log"
    params:
        output_dir = f"{UPSTREAM_ANALYSIS}/fastqc/filtered"
    shell:
        """
        apptainer exec {SINGULARITY_IMAGE} fastqc {input.R1} {input.R2} -o {params.output_dir} &>> {log}
        """

# 规则4：MultiQC 汇总报告
rule multiqc:
    input:
        original_fastqc = expand(f"{UPSTREAM_ANALYSIS}/fastqc/original/{{sample}}.R1_fastqc.html", sample=SAMPLE_IDS),
        filtered_fastqc = expand(f"{UPSTREAM_ANALYSIS}/fastqc/filtered/{{sample}}.R1.filtered_fastqc.html", sample=SAMPLE_IDS),
        fastp_reports = expand(f"{UPSTREAM_ANALYSIS}/fastp/{{sample}}.fastp.html", sample=SAMPLE_IDS)
    output:
        fastqc_original = f"{UPSTREAM_ANALYSIS}/multiqc_result/fastqc_original/multiqc_report.html",
        fastqc_filtered = f"{UPSTREAM_ANALYSIS}/multiqc_result/fastqc_filtered/multiqc_report.html",
        fastp_summary = f"{UPSTREAM_ANALYSIS}/multiqc_result/fastp/multiqc_report.html"
    log:
        f"{LOGS_DIR}/upstream/multiqc/multiqc.log"
    params:
        original_dir = f"{UPSTREAM_ANALYSIS}/fastqc/original",
        filtered_dir = f"{UPSTREAM_ANALYSIS}/fastqc/filtered",
        fastp_dir = f"{UPSTREAM_ANALYSIS}/fastp",
        multiqc_original_dir = f"{UPSTREAM_ANALYSIS}/multiqc_result/fastqc_original",
        multiqc_filtered_dir = f"{UPSTREAM_ANALYSIS}/multiqc_result/fastqc_filtered",
        multiqc_fastp_dir = f"{UPSTREAM_ANALYSIS}/multiqc_result/fastp"
    shell:
        """
        # MultiQC for original FastQC reports
        apptainer exec {SINGULARITY_IMAGE} multiqc {params.original_dir} \
            -o {params.multiqc_original_dir} &>> {log}

        # MultiQC for filtered FastQC reports
        apptainer exec {SINGULARITY_IMAGE} multiqc {params.filtered_dir} \
            -o {params.multiqc_filtered_dir} &>> {log}

        # MultiQC for Fastp reports
        apptainer exec {SINGULARITY_IMAGE} multiqc {params.fastp_dir} \
            -o {params.multiqc_fastp_dir} &>> {log}
        """

# 规则5：Hisat2生成索引文件
rule hisat2_index:
    input:
        genome = GENOME_FILE
    output:
        expand(f"{PROJECT_DIR}/references/hisat2_index/genome_index.{{suffix}}", suffix=["1.ht2", "2.ht2", "3.ht2", "4.ht2", "5.ht2", "6.ht2", "7.ht2", "8.ht2"])
    log:
        f"{LOGS_DIR}/upstream/hisat2/index.log"
    params:
        index_prefix = f"{PROJECT_DIR}/references/hisat2_index/genome_index"
    shell:
        """
        apptainer exec {SINGULARITY_IMAGE} hisat2-build {input.genome} {params.index_prefix} &>> {log}
        """

# 规则6：Hisat2进行对比
rule hisat2_alignment:
    input:
        R1 = f"{UPSTREAM_ANALYSIS}/cleaned/{{sample}}.R1.filtered.fastq.gz",
        R2 = f"{UPSTREAM_ANALYSIS}/cleaned/{{sample}}.R2.filtered.fastq.gz",
        index = expand(f"{PROJECT_DIR}/references/hisat2_index/genome_index.{{suffix}}", suffix=["1.ht2", "2.ht2", "3.ht2", "4.ht2", "5.ht2", "6.ht2", "7.ht2", "8.ht2"])
    output:
        sorted_bam = f"{UPSTREAM_ANALYSIS}/alignment/bam_file/{{sample}}.sorted.bam",
        stats = f"{UPSTREAM_ANALYSIS}/alignment/stats/{{sample}}.stats.txt"
    log:
        f"{LOGS_DIR}/upstream/hisat2_alignment/{{sample}}.log"
    params:
        index_prefix = f"{PROJECT_DIR}/references/hisat2_index/genome_index"
    shell:
        """
        # HISAT2 比对并直接生成排序后的 BAM 文件,无需预先生成sam文件
        apptainer exec {SINGULARITY_IMAGE} hisat2 -x {params.index_prefix} \
            -1 {input.R1} -2 {input.R2} 2>> {log} | \
        apptainer exec {SINGULARITY_IMAGE} samtools sort -o {output.sorted_bam} -
        
        # 生成统计信息
        apptainer exec {SINGULARITY_IMAGE} samtools flagstat {output.sorted_bam} > {output.stats}
        """
        
# 规则7：MultiQC 汇总比对报告
rule multiqc_alignment:
    input:
        # 收集所有样本的比对统计文件
        alignment_stats = expand(f"{UPSTREAM_ANALYSIS}/alignment/stats/{{sample}}.stats.txt", sample=SAMPLE_IDS)
    output:
        # MultiQC 汇总报告的输出文件
        alignment_multiqc = f"{UPSTREAM_ANALYSIS}/multiqc_result/alignment/multiqc_report.html"
    log:
        # 日志文件路径
        f"{LOGS_DIR}/upstream/multiqc_alignment/multiqc_alignment.log"
    params:
        # 输入目录（包含比对统计文件）
        alignment_dir = f"{UPSTREAM_ANALYSIS}/alignment/stats",
        # MultiQC 输出目录
        multiqc_output_dir = f"{UPSTREAM_ANALYSIS}/multiqc_result/alignment"
    shell:
        """
        # 创建 MultiQC 输出目录（如果不存在）
        mkdir -p {params.multiqc_output_dir}
        
        # 运行 MultiQC 来汇总比对统计文件
        apptainer exec {SINGULARITY_IMAGE} multiqc {params.alignment_dir} \
            -o {params.multiqc_output_dir} &>> {log}
        """

# 规则8：使用 AGAT 工具修复基因注释文件
rule agat_convert_gxf:
    input:
        gtf_file = ANNOTATION_FILE  # 原始的基因注释文件
    output:
        fixed_gtf = f"{PROJECT_DIR}/references/fixed_output.gtf"  # 修复后的基因注释文件
    log:
        f"{LOGS_DIR}/upstream/agat/agat_convert_gxf.log"
    shell:
        """
        apptainer exec {SINGULARITY_IMAGE} agat_convert_sp_gxf2gxf.pl -g {input.gtf_file} -o {output.fixed_gtf} &>> {log}
        """

# 规则9：使用 FeatureCounts 生成 gene_count_matrix 和transcript_count_matrix并处理输出格式
rule featurecounts:
    input:
        bam = expand(f"{UPSTREAM_ANALYSIS}/alignment/bam_file/{{sample}}.sorted.bam", sample=SAMPLE_IDS),  # 批量 BAM 文件
        annotation_gtf = f"{PROJECT_DIR}/references/fixed_output.gtf"  # 使用修复后的 GTF 文件
    output:
        raw_counts = f"{UPSTREAM_ANALYSIS}/quantification/featurecounts/gene_count_matrix.txt",  # 所有样本的基因表达矩阵
        transcript_counts = f"{UPSTREAM_ANALYSIS}/quantification/featurecounts/transcript_count_matrix.txt"  # 所有样本的转录本表达矩阵
    log:
        f"{LOGS_DIR}/upstream/featurecounts/featurecounts.log"  # 日志文件
    shell:
        """
        # 计算基因的表达矩阵，批量处理 BAM 文件
        apptainer exec {SINGULARITY_IMAGE} featureCounts -a {input.annotation_gtf} -t exon -g gene_id -p -o {output.raw_counts} {input.bam} &>> {log}
        
        # 计算转录本的表达矩阵，批量处理 BAM 文件
        apptainer exec {SINGULARITY_IMAGE} featureCounts -a {input.annotation_gtf} -t transcript -g transcript_id -p -o {output.transcript_counts} {input.bam} &>> {log}

        # 删除以 "#" 开头的行，并将路径和后缀去除
        # 对于 gene_count_matrix.txt
        grep -v "^#" {output.raw_counts} | \
        sed 's|{UPSTREAM_ANALYSIS}/alignment/bam_file/||g' | \
        sed 's/.sorted.bam//g' | \
        cut -f 1,7- > {output.raw_counts}.tmp && mv {output.raw_counts}.tmp {output.raw_counts}

        # 对于 transcript_count_matrix.txt
        grep -v "^#" {output.transcript_counts} | \
        sed 's|{UPSTREAM_ANALYSIS}/alignment/bam_file/||g' | \
        sed 's/.sorted.bam//g' | \
        cut -f 1,7- > {output.transcript_counts}.tmp && mv {output.transcript_counts}.tmp {output.transcript_counts}
        """

# 规则10：StringTie进行定量获得每个样本的.TSV文件
rule stringtie:
    input:
        sorted_bam = f"{UPSTREAM_ANALYSIS}/alignment/bam_file/{{sample}}.sorted.bam",
        annotation_gtf = f"{PROJECT_DIR}/references/fixed_output.gtf"  # 使用修复后的 GTF 文件
    output:
        gtf = f"{UPSTREAM_ANALYSIS}/quantification/stringtie/{{sample}}.gtf",
        abundance = f"{UPSTREAM_ANALYSIS}/quantification/stringtie/{{sample}}.abundance.tsv"
    log:
        f"{LOGS_DIR}/upstream/quantification/stringtie/{{sample}}.log"
    shell:
        """
        apptainer exec {SINGULARITY_IMAGE} stringtie {input.sorted_bam} \
            -G {input.annotation_gtf} \
            -e -o {output.gtf} \
            -A {output.abundance} &>> {log}
        """

# 规则11：基于.TSV文件归纳FPKM和TPM表达矩阵
rule extract_fpkm_tpm:
    input:
        abundance_files = expand(f"{UPSTREAM_ANALYSIS}/quantification/stringtie/{{sample}}.abundance.tsv", sample=SAMPLE_IDS)  # 确保输入是 .abundance.tsv 文件
    output:
        fpkm_matrix = f"{UPSTREAM_ANALYSIS}/quantification/stringtie/normalized_data/fpkm_matrix.csv",
        tpm_matrix = f"{UPSTREAM_ANALYSIS}/quantification/stringtie/normalized_data/tpm_matrix.csv",
        abundance_file_list = f"{UPSTREAM_ANALYSIS}/quantification/stringtie/abundance_file_list.txt"
    log:
        f"{LOGS_DIR}/upstream/quantification/extract_fpkm_tpm/extract_fpkm_tpm.log"
    shell:
        """
        # 确保输出目录存在
        mkdir -p $(dirname {output.abundance_file_list})
        
        # 生成 abundance 文件列表
        for file in {input.abundance_files}; do
            sample_id=$(basename "$file" .abundance.tsv)
            echo -e "$sample_id $file" >> {output.abundance_file_list}
        done

        # 运行提取 FPKM 和 TPM 的脚本
        apptainer exec --bind /home/data/t210528/scripts:/mnt/scripts {SINGULARITY_IMAGE} python /mnt/scripts/extract_fpkm_tpm.py3 \
            --abundance_files {input.abundance_files} \
            --fpkm_matrix {output.fpkm_matrix} \
            --tpm_matrix {output.tpm_matrix} \
            > {log} 2>&1
        """
 
# 规则12：组间基因总表达量分析
rule TGE_Analysis:
    input:
        fpkm_matrix = f"{UPSTREAM_ANALYSIS}/quantification/stringtie/normalized_data/fpkm_matrix.csv",
        phenotypic_data = f"{PHENOTYPIC_DATA}"
    output:
        violin_plot = f"{DOWNSTREAM_ANALYSIS}/TGE_Analysis/gene_expression_violin.pdf"
    log:
        f"{LOGS_DIR}/downstream/TGE_Analysis/violin_plot.log"
    shell:
        """
        # 运行 R 脚本生成小提琴图
        apptainer exec {SINGULARITY_IMAGE} Rscript /opt/scripts/TGE_Analysis.R \
            {input.fpkm_matrix} {input.phenotypic_data} {output.violin_plot} &>> {log}
        """

# 规则13：样本间相关性分析
rule GEC_Analysis:
    input:
        fpkm_matrix = f"{UPSTREAM_ANALYSIS}/quantification/stringtie/normalized_data/fpkm_matrix.csv",
        phenotypic_data = f"{PHENOTYPIC_DATA}"
    output:
        correlation_heatmap = f"{DOWNSTREAM_ANALYSIS}/GEC_Analysis/correlation_heatmap.pdf",
        output_pca_plot = f"{DOWNSTREAM_ANALYSIS}/GEC_Analysis/PCA_Plot.pdf"
    log:
        f"{LOGS_DIR}/downstream/GEC_Analysis/analysis.log"
    shell:
        """
        apptainer exec {SINGULARITY_IMAGE} Rscript /opt/scripts/GEC_Analysis.R \
            {input.fpkm_matrix} {input.phenotypic_data} \
            {output.correlation_heatmap} {output.output_pca_plot} &>> {log}
        """
        
# 规则14：基于edgeR包的差异分析
rule DEG_Analysis:
    input:
        raw_counts = f"{UPSTREAM_ANALYSIS}/quantification/featurecounts/gene_count_matrix.txt",
        phenotypic_data = f"{PHENOTYPIC_DATA}"
    output:
        deg_results = f"{DOWNSTREAM_ANALYSIS}/DEG_Analysis_edgeR/DEG_edgeR.csv",
        volcano_plot = f"{DOWNSTREAM_ANALYSIS}/DEG_Analysis_edgeR/edgeR_volcano.pdf",
        heatmap_plot = f"{DOWNSTREAM_ANALYSIS}/DEG_Analysis_edgeR/edgeR_heatmap.pdf"
    log:
        f"{LOGS_DIR}/downstream/DEG_Analysis_edgeR/analysis.log"
    shell:
        """    
        apptainer exec {SINGULARITY_IMAGE} Rscript /opt/scripts/DEG_Analysis_edgeR.R \
            {input.raw_counts} {input.phenotypic_data} \
            {output.deg_results} {output.volcano_plot} {output.heatmap_plot} &>> {log}
        """
        
# 规则15：GO 富集分析
rule GO_Enrichment:
    input:
        deg_edgeR = f"{DOWNSTREAM_ANALYSIS}/DEG_Analysis_edgeR/DEG_edgeR.csv"
    output:
        enrichment_results = expand(f"{DOWNSTREAM_ANALYSIS}/GO_Enrichment/enrichment_result_{{ont}}.csv", ont=["BP", "CC", "MF"]),
        circle_plots = expand(f"{DOWNSTREAM_ANALYSIS}/GO_Enrichment/GOCircle_Plot_{{ont}}.pdf", ont=["BP", "CC", "MF"]),
        extracted_data = f"{DOWNSTREAM_ANALYSIS}/GO_Enrichment/extracted_data.csv",
        entrezid_data = f"{DOWNSTREAM_ANALYSIS}/GO_Enrichment/DEGs_with_ENTREZID.csv"
    log:
        f"{LOGS_DIR}/downstream/GO_Enrichment/GO_Enrichment.log"
    params:
        output_dir = f"{DOWNSTREAM_ANALYSIS}/GO_Enrichment",
        species = SPECIES  # 物种名称传递给 R 脚本
    shell:
        """
        apptainer exec {SINGULARITY_IMAGE} Rscript /opt/scripts/GO_Enrichment.R \
            {input.deg_edgeR} {params.output_dir} {params.species} \
            {output.extracted_data} {output.entrezid_data} &>> {log}
        """
        
# 规则16：KEGG富集分析
rule KEGG_Enrichment:
    input:
        extracted_data = f"{DOWNSTREAM_ANALYSIS}/GO_Enrichment/extracted_data.csv",
        degs_file = f"{DOWNSTREAM_ANALYSIS}/GO_Enrichment/DEGs_with_ENTREZID.csv"
    output:
        enrichment_results = f"{DOWNSTREAM_ANALYSIS}/KEGG_Enrichment/Kegg_enrichment_result.txt",
        circle_plot = f"{DOWNSTREAM_ANALYSIS}/KEGG_Enrichment/KeggCircle_Plot.pdf"
    log:
        f"{LOGS_DIR}/downstream/KEGG_Enrichment/KEGG_Enrichment.log"
    params:
        output_dir = f"{DOWNSTREAM_ANALYSIS}/KEGG_Enrichment",
        species = SPECIES  # 物种名称传递给 R 脚本
    shell:
        """
        apptainer exec {SINGULARITY_IMAGE} Rscript /opt/scripts/KEGG_Enrichment.R \
            {input.extracted_data} {input.degs_file} {params.output_dir} {params.species} &>> {log}
        """

# 规则17：基于GO和KEGG数据库的GSEA分析
rule GSEA_Enrichment:
    input:
        deg_file = f"{DOWNSTREAM_ANALYSIS}/DEG_Analysis_edgeR/DEG_edgeR.csv",
        extracted_data = f"{DOWNSTREAM_ANALYSIS}/GO_Enrichment/extracted_data.csv",
        entrezid_data = f"{DOWNSTREAM_ANALYSIS}/GO_Enrichment/DEGs_with_ENTREZID.csv"
    output:
        go_result         = f"{DOWNSTREAM_ANALYSIS}/GSEA_Enrichment/GSEA_GO_Result.csv",
        go_circleplot     = f"{DOWNSTREAM_ANALYSIS}/GSEA_Enrichment/GSEA_GO_CirclePlot.pdf",
        kegg_result       = f"{DOWNSTREAM_ANALYSIS}/GSEA_Enrichment/GSEA_KEGG_Result.csv",
        kegg_circleplot   = f"{DOWNSTREAM_ANALYSIS}/GSEA_Enrichment/GSEA_KEGG_CirclePlot.pdf"
    log:
        f"{LOGS_DIR}/downstream/GSEA_Enrichment/GSEA_Enrichment.log"
    params:
        output_dir = f"{DOWNSTREAM_ANALYSIS}/GSEA_Enrichment",
        species    = SPECIES # 物种名称传递给 R 脚本
    shell:
        """
        apptainer exec {SINGULARITY_IMAGE} Rscript /opt/scripts/GSEA_Enrichment.R \
            {input.extracted_data} {input.entrezid_data} {params.species} {params.output_dir} &>> {log}
        """