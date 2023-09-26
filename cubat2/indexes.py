import time
import numpy as np
import pandas as pd
from Bio.SeqUtils import seq1
from Bio.Data import CodonTable
from collections import Counter
from sequences_loader import sequences_to_codonframe
from codon_table import get_codon_table_by_index, codon_table_completion

def calculate_rscu(codonframe, codon_table_index=1):
    """
    计算序列的相对密码子使用频率 (RSCU)。

    参数:
    - codonframe (pd.DataFrame): 包含每个序列的ID以及每个密码子的数量的DataFrame。
    - codon_table_index (int): 选择的密码子表的序号, 默认为1。

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
    - codon_table_index (int): 选择的密码子表的序号, 默认为1。

    返回:
    - relAdapt (pd.Series): 包含每个密码子的相对适应度的 Series。
    """

    # 从文件中读取序列并计算密码子数量
    codonframe = sequences_to_codonframe(file_path, codon_table=codon_table_index)

    # 将密码子数量纵向加和
    codonframe_sum = codonframe.iloc[:, :].sum()
    # codonframe_sum.to_csv("test.csv")
    # print(codonframe_sum)

    # 创建字典以存储每个氨基酸的密码子数量总和
    amino_acid_codon_counts = {}
    for codon, amino_acid in codon_table_completion(codon_table_index).items():
        if amino_acid not in amino_acid_codon_counts:
            amino_acid_codon_counts[amino_acid] = 0
        amino_acid_codon_counts[amino_acid] += codonframe_sum[codon]
    
    # 计算相对适应度
    relAdapt = pd.Series()
    for codon in codon_table_completion(codon_table_index):
        amino_acid = codon_table_completion(codon_table_index)[codon]
        relAdapt[codon] = codonframe_sum[codon] / amino_acid_codon_counts[amino_acid]

    return relAdapt


def relative_adapt_loader(csv_file, codon_table_index=1):
    """
    根据输入的 CSV 文件内容判断是密码子数量数据还是相对适应度数据，然后计算相对适应度并返回结果。

    参数:
    - csv_file (str): 输入的 CSV 文件路径。

    返回:
    - relAdapt (pd.Series): 包含每个密码子的相对适应度的 Series。
    """

    # 从 CSV 文件中读取数据
    codonframe_sum = pd.read_csv(csv_file, index_col=0).iloc[:, 0]
    
    # 判断数据类型
    # TODO 判断这里可以加一个QA
    # print(codonframe_sum)
    if codonframe_sum.apply(lambda x: isinstance(x, int) and x >= 0).all():
    # 创建字典以存储每个氨基酸的密码子数量总和
        amino_acid_codon_counts = {}
        for codon, amino_acid in codon_table_completion(codon_table_index).items():
            if amino_acid not in amino_acid_codon_counts:
                amino_acid_codon_counts[amino_acid] = 0
            amino_acid_codon_counts[amino_acid] += codonframe_sum[codon]
    
        # 计算相对适应度
        relAdapt = pd.Series()
        for codon in codon_table_completion(codon_table_index):
            amino_acid = codon_table_completion(codon_table_index)[codon]
            relAdapt[codon] = codonframe_sum[codon] / amino_acid_codon_counts[amino_acid]
    else:
        relAdapt = pd.Series(codonframe_sum)
    
    return relAdapt


def calculate_CAI(sequences, relAdapt, codon_table_index=1):
    """
    计算密码子适应性指数 (CAI)。

    参数:
    - sequences (pd.DataFrame): 包含密码子数量数据的DataFrame。
    - relAdapt (pd.DataFrame, pd.Series or str): 可以是包含相对适应度数据的DataFrame、Series或CSV文件的路径。
    - codon_table_index (int): 选择的密码子表的序号, 默认为1。

    返回:
    - cai (pd.DataFrame): 包含 CAI 和 CAI2 值的DataFrame。
    """

    # 根据输入的 relAdapt 类型调用相应的函数获取相对适应度数据
    if isinstance(relAdapt, pd.DataFrame) or isinstance(relAdapt, pd.Series):
        relAdapt = relAdapt
    elif isinstance(relAdapt, str):
        if relAdapt.endswith('.csv'):
            # 如果输入是CSV文件路径，则调用 relative_adapt_loader 函数
            relAdapt = relative_adapt_loader(relAdapt, codon_table_index=codon_table_index)
        else:
            # 如果输入是文件路径但不是CSV文件，则调用 calculate_relative_adapt 函数
            relAdapt = calculate_relative_adapt(relAdapt, codon_table_index=codon_table_index)

    # 获取选择的密码子表
    codon_table = get_codon_table_by_index(codon_table_index)

    # 计算 CAI
    cai = pd.DataFrame(index=sequences.index, columns=["CAI", "CAI2"])
    cai.fillna(0, inplace=True)

    # 使用向量化操作计算 CAI
    for col in sequences.columns:
        codon_counts = sequences[col]
        target_aa = codon_table.forward_table.get(col)

        if target_aa:
            # 使用向量化操作计算 CAI 和 CAI2
            cai["CAI2"] += codon_counts * relAdapt[col]
            cai["CAI"] += codon_counts * np.log(relAdapt[col])

    # 使用向量化操作计算 CAI2 和 CAI
    cai["CAI2"] /= sequences.sum(axis=1)
    cai["CAI"] = np.exp(cai["CAI"] / sequences.sum(axis=1))

    return cai



# 示例用法：
if __name__ == "__main__":
    start_time = time.time()
    file_path = "../test_data/Drosophila_melanogaster.fna"
    codonframe = sequences_to_codonframe(file_path)
    # RSCU测试代码
    if False:
        rscuSheet = calculate_rscu(codonframe, codon_table_index=1)
        print(rscuSheet)
        rscuSheet.to_csv("test.csv")

    # calculate_relative_adapt测试代码
    if False:
        relAdapt = calculate_relative_adapt(file_path, codon_table_index=1)
        print(relAdapt)
        relAdapt.to_csv("test.csv")

    # relative_adapt_loader测试代码
    if False:
        csv_file = "test.csv"  # 替换为你的 CSV 文件路径
        relAdapt = relative_adapt_loader(csv_file)
        print(relAdapt)
    
    # CAI测试代码
    if True:
        relAdapt = calculate_relative_adapt(file_path, codon_table_index=1)
        cai = calculate_CAI(codonframe, relAdapt)
        print(cai)

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"函数执行时间：{execution_time}秒")