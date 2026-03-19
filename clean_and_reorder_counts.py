#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
一个用于清理并重排序featureCounts输出文件的命令行脚本。

该脚本会处理由featureCounts生成的原始计数文件，通过移除路径和文件扩展名来清理样本列名，
然后根据一个给定的表型文件（phenotype file）中指定的顺序对样本列进行重新排序。

此脚本设计为通过命令行运行，需要三个参数：
输入计数文件的路径、表型文件的路径以及输出文件的路径。

使用示例:
    python3 your_script_name.py --counts raw_counts.txt --phenotype samplesheet.tsv --output cleaned_counts.tsv
"""

import pandas as pd
import argparse
import os

def clean_and_reorder_counts(counts_file, phenotype_file, output_file):
    """
    清理featureCounts的列名，并根据表型文件对列进行重排序。

    此函数读取一个featureCounts的原始计数矩阵和一个表型文件。它会从计数文件
    头部的BAM文件路径中提取出干净的样本名，然后根据表型文件中提供的样本顺序
    来重新排列计数矩阵的列。

    参数 (Args):
        counts_file (str): featureCounts输出的原始计数文件路径。
        phenotype_file (str): 包含期望样本顺序的表型文件路径（第一列应为样本名）。
        output_file (str): 清理并重排后的计数矩阵将被保存到的文件路径。

    异常 (Raises):
        ValueError: 如果表型文件中列出的任何样本在清理后的计数文件中都找不到，
                    则会引发此异常。
    """
    # --- 步骤 1: 读取表型文件，获取期望的样本顺序 ---
    # 表型文件应为制表符分隔（tab-separated）。
    phenotype_df = pd.read_csv(phenotype_file, sep='\t')
    # 假设第一列包含样本名，并移除可能存在的前后空格。
    sample_order = phenotype_df.iloc[:, 0].str.strip().tolist()
    print(f"从表型文件获取的期望样本顺序: {sample_order}")

    # --- 步骤 2: 读取featureCounts的原始计数文件 ---
    # featureCounts的输出通常包含以'#'开头的注释行。
    # `comment='#'` 参数可以确保在读取时跳过这些行。
    counts_df = pd.read_csv(counts_file, sep='\t', comment='#')

    # --- 步骤 3: 清理列名 ---
    # 原始列名是BAM文件的完整路径，例如：
    # /path/to/your/bams/Control1.Aligned.sortedByCoord.out.bam
    # 我们的目标是把它简化为样本名，例如：'Control1'。

    # 创建一个字典，用于映射旧列名到新的、清理后的列名。
    rename_dict = {}
    # featureCounts输出的前6列是基因的元数据（Geneid, Chr等）。
    # 样本计数数据从第7列开始。
    for col in counts_df.columns[6:]:
        # 从完整路径中获取基本文件名。
        # 例如: "Control1.Aligned.sortedByCoord.out.bam"
        base_name = os.path.basename(col)
        # 移除特定的后缀以分离出干净的样本名。
        # 例如: "Control1"
        # 注意: 如果你的文件名格式不同，你可能需要修改这里被替换的字符串。
        sample_name = base_name.replace('.Aligned.sortedByCoord.out.bam', '')
        rename_dict[col] = sample_name

    # 在DataFrame上执行重命名操作。
    counts_df.rename(columns=rename_dict, inplace=True)
    print(f"清理后的列名: {counts_df.columns.tolist()}")

    # --- 步骤 4: 根据表型文件对列进行重排序 ---
    # 最终的列顺序应始终以 'Geneid' 列开始，然后是表型文件中定义的样本列。
    final_columns = ['Geneid'] + sample_order

    # 验证所有来自表型文件的样本是否存在于计数数据中。
    # 这可以防止因文件不匹配而导致的错误。
    missing_cols = [col for col in final_columns if col not in counts_df.columns]
    if missing_cols:
        # 如果有任何样本缺失，就引发一个带有描述性信息的错误。
        raise ValueError(f"错误：表型文件中的以下样本在计数文件中不存在: {missing_cols}")

    # 选择并重排序DataFrame的列。
    reordered_df = counts_df[final_columns]

    # --- 步骤 5: 将处理完成的最终结果写入新的TSV文件 ---
    # `index=False` 防止pandas将DataFrame的行索引写入到输出文件中。
    reordered_df.to_csv(output_file, sep='\t', index=False)
    print(f"成功创建已清理和重排序的文件: {output_file}")


if __name__ == '__main__':
    # 这个代码块设置了命令行界面。
    # 它只在脚本被直接执行时才会运行。
    parser = argparse.ArgumentParser(
        description="清理并根据表型文件重排序featureCounts的输出。"
    )

    # 定义命令行参数。
    parser.add_argument(
        '--counts',
        type=str,
        required=True,
        help='featureCounts原始计数文件的路径 (例如: counts.txt)。'
    )
    parser.add_argument(
        '--phenotype',
        type=str,
        required=True,
        help='包含样本顺序的表型文件路径 (例如: a.tsv)。'
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='输出的清理后计数文件的路径 (例如: cleaned_counts.tsv)。'
    )
    
    # 解析用户提供的参数。
    args = parser.parse_args()
    
    # 使用解析后的参数调用主函数。
    clean_and_reorder_counts(args.counts, args.phenotype, args.output)
