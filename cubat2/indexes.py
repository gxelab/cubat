import time
import numpy as np
import pandas as pd
from collections import defaultdict
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
    根据输入的 CSV 文件内容判断是密码子数量数据还是相对适应度数据,然后计算相对适应度并返回结果。

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
            # 如果输入是CSV文件路径,则调用 relative_adapt_loader 函数
            relAdapt = relative_adapt_loader(relAdapt, codon_table_index=codon_table_index)
        else:
            # 如果输入是文件路径但不是CSV文件,则调用 calculate_relative_adapt 函数
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



def calculate_enc(codonframe, codon_table_index=1):
    """
    计算基因的enc值。

    参数:
    - codonframe (pd.DataFrame): 包含密码子频率信息的DataFrame,由sequences_to_codonframe函数生成。
    - codon_table_index (int): 选择的密码子表的序号,默认为1。

    返回:
    - enc (pd.DataFrame): 包含每个基因的enc值的DataFrame,以"enc"列命名。
    """

    # 获取选择的密码子表
    codon_table = get_codon_table_by_index(codon_table_index)

    # 获取同义密码子的列表
    codon_list = list(codon_table.forward_table.keys())

    # 创建一个字典用于存储每个同义密码子类别的信息
    subfam_dict = defaultdict(list)
    for codon in codon_list:
        subfam = codon_table.forward_table[codon]
        subfam_dict[subfam].append(codon)


    # 计算每个同义密码子类别内的 F_CF 和 n_j
    fcf_dict = {}
    nj_dict = {}
    for subfam, codons in subfam_dict.items():
        subfam_frame = codonframe[codons]
        n = subfam_frame.sum(axis=1)
        mx = subfam_frame.max(axis=1)
        ncol_mx = subfam_frame.shape[1]
        fcf = (mx + 1) / (n + ncol_mx)
        nj = mx + 1
        fcf_dict[subfam] = fcf
        nj_dict[subfam] = nj

    # 计算 enc
    enc = pd.DataFrame(index=codonframe.index, columns=["enc"], dtype=float)
    enc.fillna(0, inplace=True)

    ss = len(subfam_dict)  # 同义密码子类别的数量

    N_single = 1 / ss  # 只有一个同义密码子类别的基因的贡献

    if ss >= 2:
        N_double = 2 * (ss - 1) / (ss * (ss - 2))  # 有两个同义密码子类别的基因的贡献
    else:
        N_double = 0

    if ss >= 3:
        N_triple = 6 / (ss * (ss - 1) * (ss - 2))  # 有三个同义密码子类别的基因的贡献
    else:
        N_triple = 0

    if ss >= 4:
        N_quad = 24 / (ss * (ss - 1) * (ss - 2) * (ss - 3))  # 有四个同义密码子类别的基因的贡献
    else:
        N_quad = 0

    for gene_id in codonframe.index:
        fcf_gene = [fcf_dict[subfam][gene_id] for subfam in subfam_dict.keys()]
        nj_gene = [nj_dict[subfam][gene_id] for subfam in subfam_dict.keys()]

        enc.loc[gene_id, "enc"] = N_single + N_double * sum(fcf_gene) + N_triple * sum(nj_gene) + N_quad * sum(nj_gene)

    return enc


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
    if False:
        relAdapt = calculate_relative_adapt(file_path, codon_table_index=1)
        cai = calculate_CAI(codonframe, relAdapt)
        print(cai)
    
    # ENC测试代码
    if True:
        enc = calculate_enc(codonframe)
        print(enc)


    end_time = time.time()
    execution_time = end_time - start_time
    print(f"执行时间：{execution_time}秒")