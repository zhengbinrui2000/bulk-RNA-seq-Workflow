#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
一个用于从featureCounts计数矩阵和GTF注释文件计算FPKM和TPM的脚本。

该脚本执行以下操作：
1.  解析GTF文件，计算每个基因所有外显子（exon）的总长度。
2.  读取一个清理过的基因计数矩阵（例如，由前一个脚本生成）。
3.  使用基因长度和计数值来计算每个基因在每个样本中的FPKM（每千碱基每百万读数片段数）。
4.  使用基因长度和计数值来计算每个基因在每个样本中的TPM（每百万转录本数）。
5.  将计算出的FPKM和TPM矩阵分别保存到指定的输出文件中。

使用示例:
    python3 your_script_name.py --gtf gencode.v38.annotation.gtf --counts cleaned_counts.tsv --out_fpkm fpkm_matrix.tsv --out_tpm tpm_matrix.tsv
"""

import pandas as pd
from collections import defaultdict
import argparse

# --- 1. 设置命令行参数解析 ---
# 这部分代码定义了脚本运行时需要从用户那里接收哪些输入。
parser = argparse.ArgumentParser(description="从featureCounts矩阵和GTF文件计算FPKM和TPM。")
parser.add_argument('--gtf', type=str, required=True, help='GTF注释文件的路径。')
parser.add_argument('--counts', type=str, required=True, help='来自featureCounts的清理后计数文件的路径。')
parser.add_argument('--out_fpkm', type=str, required=True, help='FPKM矩阵的输出路径。')
parser.add_argument('--out_tpm', type=str, required=True, help='TPM矩阵的输出路径。')

args = parser.parse_args()

# --- 2. 从 GTF 文件中提取每个基因的外显子总长度 ---
# FPKM 和 TPM 的计算需要基因的有效长度，通常使用基因所有外显子的长度之和。
print(f"信息: 正在从 {args.gtf} 读取基因长度...")
# `defaultdict(int)` 创建一个特殊的字典，如果访问一个不存在的键，它会自动将其值初始化为0。
# 这在累加长度时非常方便，无需检查键是否存在。
gene_lengths = defaultdict(int)
with open(args.gtf, 'r') as f:
    for line in f:
        # GTF文件中的注释行以'#'开头，应跳过。
        if line.startswith('#'):
            continue

        # GTF是制表符分隔的格式，将其分割成字段列表。
        fields = line.strip().split('\t')

        # 我们只关心'exon'（外显子）类型的特征，因为这是表达并被测序的部分。
        if fields[2] == 'exon':
            gene_id = ''
            # 第9列（属性列）包含基因ID，但格式复杂，需要解析。
            # 格式通常是 'key1 "value1"; key2 "value2"; ...'
            for item in fields[8].split(';'):
                item = item.strip() # 去除前后空格
                if item.startswith('gene_id'):
                    # 找到以'gene_id'开头的条目后，提取引号内的ID。
                    gene_id = item.split('"')[1]
                    break # 找到gene_id后就不需要再查找了

            # 如果成功提取到gene_id，则计算该外显子的长度并累加到对应基因上。
            if gene_id:
                # GTF坐标是1-based且包含两端，所以长度是 (终点 - 起点 + 1)。
                exon_length = abs(int(fields[4]) - int(fields[3])) + 1
                gene_lengths[gene_id] += exon_length

# --- 3. 加载经过清理的 count 矩阵 ---
print(f"信息: 正在从 {args.counts} 读取计数矩阵...")
# `index_col=0` 指定文件中的第一列（通常是Geneid）作为DataFrame的索引。
counts_df = pd.read_csv(args.counts, sep='\t', index_col=0)

# 为了确保计算正确，基因长度的顺序必须与count矩阵的基因顺序完全一致。
# 下面的代码根据counts_df的索引顺序，从gene_lengths字典中创建一个新的pandas Series。
# `.get(gene, 0)`确保如果某个基因在GTF中没找到，其长度为0，而不是报错。
ordered_gene_lengths = pd.Series([gene_lengths.get(gene, 0) for gene in counts_df.index], index=counts_df.index)

# 过滤掉长度为0的基因，这可能是因为它们在GTF文件中没有被注释，或者没有外显子。
# 这样做可以避免后续计算中出现“除以零”的错误。
original_gene_count = len(counts_df)
ordered_gene_lengths = ordered_gene_lengths[ordered_gene_lengths > 0]
# 同时，也要从原始的counts_df中移除这些基因，以保持数据一致。
counts_df = counts_df.loc[ordered_gene_lengths.index]
print(f"信息: 移除了 {original_gene_count - len(counts_df)} 个长度为零的基因。")

# --- 4. 计算 FPKM ---
# FPKM公式: (原始读数 * 10^9) / (基因长度 * 该样本的总读数)
print("信息: 正在计算 FPKM...")
# 步骤 a: 计算 RPK (Reads Per Kilobase)，即每千碱基的读数。
# 这是公式中的 (原始读数 / 基因长度) 部分，其中基因长度以kb为单位。
# `axis=0` 表示按行操作，即用每一行的计数值除以该行对应的基因长度。
rpk = counts_df.div(ordered_gene_lengths / 1000, axis=0)

# 步骤 b: 计算 "per million" 缩放因子。
# 这是公式中的 (该样本的总读数 / 1,000,000) 部分。
# `counts_df.sum()` 默认按列求和，得到每个样本的总读数。
scaling_factor_fpkm = counts_df.sum() / 1_000_000

# 步骤 c: 计算最终的 FPKM。
# 用RPK值除以缩放因子。
# `axis=1` 表示按列操作，即用每一列的RPK值除以该列（样本）对应的缩放因子。
fpkm = rpk.div(scaling_factor_fpkm, axis=1)

# --- 5. 计算 TPM ---
# TPM公式: (RPK / 所有基因RPK之和) * 1,000,000
print("信息: 正在计算 TPM...")
# 步骤 a: 计算RPK (上一步已经计算过了，这里直接使用`rpk`变量)。

# 步骤 b: 计算每个样本所有RPK值的总和。
rpk_sum_per_sample = rpk.sum() # 默认按列求和

# 步骤 c: 计算TPM的缩放因子。
# TPM的逻辑是先将RPK归一化，使其总和为1，然后再乘以一百万。
# `rpk_sum_per_sample / 1_000_000` 就是这个缩放因子。
scaling_factor_tpm = rpk_sum_per_sample / 1_000_000

# 步骤 d: 计算最终的 TPM。
# 用RPK值除以TPM缩放因子。
tpm = rpk.div(scaling_factor_tpm, axis=1)

# --- 6. 保存结果 ---
print(f"信息: 正在保存 FPKM 矩阵到 {args.out_fpkm}...")
# 将结果保存为制表符分隔的文件，`sep='\t'`。
fpkm.to_csv(args.out_fpkm, sep='\t')

print(f"信息: 正在保存 TPM 矩阵到 {args.out_tpm}...")
tpm.to_csv(args.out_tpm, sep='\t')

print("信息: 全部完成。")
