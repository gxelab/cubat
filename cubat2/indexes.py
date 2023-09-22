import numpy as np
import pandas as pd
from Bio.SeqUtils import seq1
from Bio.Data import CodonTable
from collections import Counter
from codon_table import get_codon_table_by_index
from sequences_loader import sequences_to_codonframe


def calculate_rscu(codonframe, codon_table_index=1):
    """
    计算序列的相对密码子使用频率 (RSCU)。

    参数:
    - codonframe (pd.DataFrame): 包含每个序列的ID以及每个密码子的数量的DataFrame。
    - codon_table_index (int): 选择的密码子表的序号，默认为1。

    返回:
    - rscuSheet (pd.DataFrame): 包含 RSCU 值的DataFrame。
    """

    # 获取选择的密码子表
    codon_table = get_codon_table_by_index(codon_table_index)

    # 剔除终止密码子列
    stop_codons = codon_table.stop_codons
    rscuSheet = codonframe.drop(columns=stop_codons, errors='ignore')

    # 创建 rscuSheet 用于存储 RSCU 值
    rscuSheet = rscuSheet.copy()

    # 计算 RSCU
    for col in rscuSheet.columns[:]:
        codon_count = rscuSheet[col]
        print(col)
        print(codon_count)
        target_aa = codon_table.forward_table.get(col)
        print(target_aa)
        if target_aa:
            codons_for_aa = [codon for codon, amino_acid in codon_table.forward_table.items() if amino_acid == target_aa]
            aa_count = sum(rscuSheet[codon] for codon in codons_for_aa)
            aa_multiplicity = len(codons_for_aa)
            rscu = aa_multiplicity * codon_count / aa_count
            rscuSheet[col] = rscu

    return rscuSheet


def calculate_relative_adapt(file_path, codon_table_index=1):
    """
    计算相对适应度。

    参数:
    - file_path (str): 包含序列的文件路径。
    - codon_table_index (int): 选择的密码子表的序号，默认为1。

    返回:
    - relAdapt (pd.Series): 包含每个密码子的相对适应度的 Series。
    """

    # 从文件中读取序列并计算密码子数量
    codonframe = sequences_to_codonframe(file_path, codon_table=codon_table_index)

    # 将密码子数量纵向加和
    codonframe_sum = codonframe.iloc[:, :].sum()
    print(codonframe_sum)
    # 获取选择的密码子表
    codon_table = get_codon_table_by_index(codon_table_index)

     # 创建字典以存储每个氨基酸的密码子数量总和
    amino_acid_codon_counts = {}
    for codon, amino_acid in codon_table.forward_table.items():
        if amino_acid not in amino_acid_codon_counts:
            amino_acid_codon_counts[amino_acid] = 0
        amino_acid_codon_counts[amino_acid] += codonframe_sum[codon]
    
    # 计算相对适应度
    relAdapt = pd.Series()
    for codon in codon_table.forward_table:
        amino_acid = codon_table.forward_table[codon]
        relAdapt[codon] = codonframe_sum[codon] / amino_acid_codon_counts[amino_acid]

    return relAdapt



# 示例用法：
if __name__ == "__main__":
    file_path = "../test_data/Drosophila_melanogaster.fna"
    codonframe = sequences_to_codonframe(file_path)
    # RSCU测试代码
    if False:
        rscuSheet = calculate_rscu(codonframe, codon_table_index=1)
        print(rscuSheet)
        rscuSheet.to_csv("test.csv")

    # RSCU测试代码
    if False:
        relAdapt = calculate_relative_adapt(file_path, codon_table_index=1)
        print(relAdapt)
        relAdapt.to_csv("test.csv")
    