#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// 原始数据通道
Channel.fromFilePairs("${params.genotypic_dir}/*_{1,2}.fastq.gz", size: 2)
    .set { raw_reads }

// 流程定义 (Processes)

// --- 质控 ---
process FASTQC_RAW {
    publishDir "${params.output_dir}/fastqc_raw", mode: 'copy'
    tag "${sample}"
    
    input:
    tuple val(sample), path(reads)
    
    output:
    path "*.{html,zip}", emit: raw_reports
    
    script: 
    """
    fastqc -t ${task.cpus} ${reads.join(' ')}
    """
}

// --- 过滤 ---
process FASTP_FILTER {
    publishDir "${params.output_dir}/fastp", mode: 'copy'
    tag "${sample}"
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}_filtered_{1,2}.fastq.gz"), emit: reads
    path "*.{html,json}", emit: fastp_reports
    
    script:
    def (read1, read2) = reads
    """
    fastp \\
        -i ${read1} \\
        -I ${read2} \\
        -o ${sample}_filtered_1.fastq.gz \\
        -O ${sample}_filtered_2.fastq.gz \\
        --html ${sample}_filtered_report.html \\
        --json ${sample}_filtered_report.json \\
        --thread ${task.cpus}
    """
}

// --- 二次质控 ---
process FASTQC_CLEAN {
    publishDir "${params.output_dir}/fastqc_clean", mode: 'copy'
    tag "${sample}"
    
    input:
    tuple val(sample), path(filtered_reads)
    
    output:
    path "*.{html,zip}", emit: clean_reports
    
    script:
    """
    fastqc -t ${task.cpus} ${filtered_reads.join(' ')}
    """
}

// --- 使用star进行索引构建 ---
process STAR_GENOME_GENERATE {
    publishDir "${params.star_index}", mode: 'copy'
    tag "Build_STAR_Index"

    input:
    path reference
    path annotation

    output:
    path "star_index_dir", emit: index

    script:
    """
    STAR --runThreadN ${task.cpus} \\
         --runMode genomeGenerate \\
         --genomeDir star_index_dir \\
         --genomeFastaFiles ${reference} \\
         --sjdbGTFfile ${annotation} \\
         --sjdbOverhang 149
    """
}

// --- 使用star进行比对 ---
process STAR_ALIGN {
    publishDir "${params.output_dir}/star_aligned/${sample}", mode: 'copy'
    tag "${sample}"

    input:
    tuple val(sample), path(reads)
    path star_index

    output:
    tuple val(sample), path("${sample}.Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(sample), path("${sample}.Log.final.out"), emit: log

    script:
    def (read1, read2) = reads
    """
    STAR --runThreadN ${task.cpus} \\
         --genomeDir ${star_index} \\
         --readFilesIn ${read1} ${read2} \\
         --readFilesCommand zcat \\
         --outFileNamePrefix "${sample}." \\
         --outSAMtype BAM SortedByCoordinate \\
         --outSAMunmapped Within \\
         --outSAMattributes Standard
    """
}

// --- 使用featurecounts进行定量 ---
process FEATURECOUNTS {
    publishDir "${params.output_dir}/featurecounts", mode: 'copy'
    tag "featureCounts_Quantification"

    input:
    path bams
    path annotation

    output:
    path "counts.txt", emit: counts_matrix
    path "counts.txt.summary", emit: summary

    script:
    """
    featureCounts \\
        -T ${task.cpus} \\
        -p \\
        -g gene_id \\
        -a ${annotation} \\
        -o counts.txt \\
        ${bams.join(' ')}
    """
}

// --- 对gene count进行格式处理并根据表型文件排序
process CLEAN_AND_REORDER_COUNTS {
    publishDir "${params.output_dir}/featurecounts", mode: 'copy'
    tag "清理并重排计数矩阵"

    input:
    path counts_file
    path phenotype_file

    output:
    path "cleaned_counts.tsv", emit: cleaned_counts

    script:
    """
    python3 /home/data/t210528/scripts/clean_and_reorder_counts.py \\
        --counts ${counts_file} \\
        --phenotype ${phenotype_file} \\
        --output cleaned_counts.tsv
    """
}

// --- 基于gene count获取FPKM/TPM ---
process QUANTIFICATION {
    publishDir "${params.output_dir}/quantification", mode: 'copy'
    tag "TPM_FPKM_Quantification"

    input:
    path counts_file
    path annotation

    output:
    path "fpkm_counts.tsv", emit: fpkm_counts
    path "tpm_counts.tsv" , emit: tpm_counts

    script:
    """
    python3 /home/data/t210528/scripts/calculate_fpkm_tpm.py \\
        --gtf ${annotation} \\
        --counts ${counts_file} \\
        --out_fpkm fpkm_counts.tsv \\
        --out_tpm tpm_counts.tsv
    """
}

// --- 汇总报告流程 ---
process MULTIQC_ALL {
    publishDir "${params.output_dir}", mode: 'copy'
    label 'multiqc'

    input:
    path raw_reports
    path fastp_reports
    path clean_reports
    path star_logs
    path featurecounts_summary

    output:
    path "multiqc_reports", emit: report

    script:
    """
    multiqc . -f -o multiqc_reports
    """
}

// --- 基于FPKM分析组间总表达量差异 ---
process FPKM_VIOLIN_PLOT {
    publishDir "${params.output_dir}/visualization", mode: 'copy'
    tag "ViolinPlot"

    input:
    path fpkm_matrix
    path phenotype_file
    
    output:
    path "Gene_Expression_Violin_Plot.pdf", emit: violin_plot

    script:
    """
    Rscript /home/data/t210528/scripts/fpkm_violin_plot.R \\
        ${fpkm_matrix} \\
        ${phenotype_file} \\
        Gene_Expression_Violin_Plot.pdf \\
        "${params.control_group}" \\
        "${params.experimental_group}"
    """
}

// --- 样本间相关性与PCA分析 ---
process QC_AND_PCA {
    publishDir "${params.output_dir}/visualization", mode: 'copy'
    tag "QC_and_PCA"

    input:
    path fpkm_matrix
    path phenotype_file

    output:
    path "correlation_heatmap.pdf", emit: correlation_heatmap
    path "pca_plot.pdf",            emit: pca_plot

    script:
    """
    Rscript /home/data/t210528/scripts/qc_and_pca.R \\
        ${fpkm_matrix} \\
        ${phenotype_file} \\
        correlation_heatmap.pdf \\
        pca_plot.pdf
    """
}

// --- 基于edgeR包的差异分析 ---
process DIFFERENTIAL_EXPRESSION {
    publishDir "${params.output_dir}/differential_expression", mode: 'copy'
    tag "DEG - ${params.experimental_group}_vs_${params.control_group}"

    input:
    path counts_matrix
    path phenotype_file

    output:
    path "DEG_results.csv",       emit: deg_results
    path "volcano_plot.pdf",      emit: volcano_plot
    path "heatmap_top20_DEG.pdf", emit: heatmap

    script:
    """
    Rscript /home/data/t210528/scripts/deg_edgeR.R \\
        ${counts_matrix} \\
        ${phenotype_file} \\
        "${params.control_group}" \\
        "${params.experimental_group}" \\
        ${params.logFC} \\
        ${params.FDR} \\
        DEG_results.csv \\
        volcano_plot.pdf \\
        heatmap_top20_DEG.pdf
    """
}

// --- GO富集分析 ---
process GO_ENRICHMENT {
    publishDir "${params.output_dir}/go_enrichment", mode: 'copy'
    tag "GO Enrichment"

    input:
    path deg_results
    val species

    output:
    path "GO_enrichment_bar_plot.pdf", emit: go_plot
    path "GO_results_BP.csv"         , emit: go_results_bp
    path "GO_results_CC.csv"         , emit: go_results_cc
    path "GO_results_MF.csv"         , emit: go_results_mf

    script:
    """
    Rscript /home/data/t210528/scripts/go_enrichment.R \\
        ${deg_results} \\
        ${species} \\
        GO_results_BP.csv \\
        GO_results_CC.csv \\
        GO_results_MF.csv \\
        GO_enrichment_bar_plot.pdf
    """
}

// --- KEGG富集分析 ---
process KEGG_ENRICHMENT {
    publishDir "${params.output_dir}/kegg_enrichment", mode: 'copy'
    tag "KEGG Enrichment"

    input:
    path deg_results
    val species

    output:
    path "KEGG_enrichment_results.csv", emit: kegg_results
    path "KEGG_circle_plot.pdf"       , emit: kegg_plot

    script:
    """
    Rscript /home/data/t210528/scripts/kegg_enrichment.R \\
        ${deg_results} \\
        ${species} \\
        KEGG_enrichment_results.csv \\
        KEGG_circle_plot.pdf
    """
}

// --- GSEA 分析 ---
process GSEA_ANALYSIS {
    publishDir "${params.output_dir}/gsea", mode: 'copy'
    tag "GSEA - GO & KEGG"

    input:
    path deg_results
    val species

    output:
    path "gsea_go_results.csv", emit: go_csv
    path "gsea_go_plot.pdf",    emit: go_pdf
    path "gsea_kegg_results.csv", emit: kegg_csv
    path "gsea_kegg_plot.pdf",    emit: kegg_pdf

    script:
    """
    Rscript /home/data/t210528/scripts/gsea_analysis.R \\
        ${deg_results} \\
        ${species} \\
        gsea_go_results.csv \\
        gsea_go_plot.pdf \\
        gsea_kegg_results.csv \\
        gsea_kegg_plot.pdf
    """
}


// ===================================================================================
// ==                                 工作流编排 (Workflow)                                 ==
// ===================================================================================
workflow {
    // 1. 构建 STAR 索引文件 (运行一次)
    STAR_GENOME_GENERATE(
        params.reference,
        params.annotation
    )

    // 2. 原始数据质控
    FASTQC_RAW(raw_reads)

    // 3. 数据过滤
    FASTP_FILTER(raw_reads)

    // 4. 过滤后数据质控
    FASTQC_CLEAN(FASTP_FILTER.out.reads)

    // 5. 使用 STAR 进行比对
    STAR_ALIGN(
        FASTP_FILTER.out.reads,
        STAR_GENOME_GENERATE.out.index
    )

    // 6. 使用 featureCounts 进行定量
    FEATURECOUNTS(
        STAR_ALIGN.out.bam.map { it[1] }.collect(),
        params.annotation
    )

    // 7. 对gene count矩阵进行格式处理和重排序
    CLEAN_AND_REORDER_COUNTS(
        FEATURECOUNTS.out.counts_matrix,
        file(params.phenotype_dir)
    )

    // 8. 计算 FPKM 和 TPM
    QUANTIFICATION(
        CLEAN_AND_REORDER_COUNTS.out.cleaned_counts,
        params.annotation
    )

    // 9. 收集所有报告并统一调用 MultiQC
    MULTIQC_ALL(
        FASTQC_RAW.out.raw_reports.collect(),
        FASTP_FILTER.out.fastp_reports.collect(),
        FASTQC_CLEAN.out.clean_reports.collect(),
        STAR_ALIGN.out.log.map { it[1] }.collect(),
        FEATURECOUNTS.out.summary.collect()
    )
    
    // 10. 组间总表达量差异分析
    FPKM_VIOLIN_PLOT(
        QUANTIFICATION.out.fpkm_counts,
        file(params.phenotype_dir)
    )

    // 11. 样本间相关性与PCA分析
    QC_AND_PCA(
        QUANTIFICATION.out.fpkm_counts,
        file(params.phenotype_dir)
    )
    
    // 12. 使用edgeR包进行差异分析
    DIFFERENTIAL_EXPRESSION(
        CLEAN_AND_REORDER_COUNTS.out.cleaned_counts,
        file(params.phenotype_dir)
    )

    // 13. 进行GO富集分析
    GO_ENRICHMENT(
        DIFFERENTIAL_EXPRESSION.out.deg_results,
        params.species
    )

    // 14. 进行KEGG富集分析
    KEGG_ENRICHMENT(
        DIFFERENTIAL_EXPRESSION.out.deg_results,
        params.species
    )

    // 15. 进行GSEA分析
    GSEA_ANALYSIS(
        DIFFERENTIAL_EXPRESSION.out.deg_results,
        params.species
    )
}
