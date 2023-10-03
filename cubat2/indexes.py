import time
import numpy as np
import pandas as pd
from itertools import product
from collections import defaultdict
from sequences_loader import sequences_to_codonframe
from codon_table import get_codon_table_by_index, codon_table_completion, codon_table_subfamily, codon_pair_table

def calculate_RSCU(codonframe, codon_table_index=1):
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
    - codon_table_index (int): 选择的密码子表的序号, 默认为1。

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


def calculate_ENC(codonframe, codon_table_index=1):
    """
    计算基因的ENC(Effective Number of Codons)值。

    参数:
    - codonframe (pd.DataFrame): 包含每个序列的ID以及每个密码子的数量的DataFrame。
    - codon_table_index (int): 选择的密码子表的序号, 默认为1。

    返回:
    - enc (pd.DataFrame): 包含每个基因的ENC值的DataFrame。
    """
    
    # 获取选择的密码子表
    codon_table = get_codon_table_by_index(codon_table_index)
    subfamily_table = codon_table_subfamily(codon_table_index)
    
    # 剔除终止密码子列
    stop_codons = codon_table.stop_codons
    codonframe = codonframe.drop(columns=stop_codons, errors='ignore')

    # 计算同义密码子的亚家族（subfamily）
    subfamily_dict = {}
    for codon, subfamily in zip(subfamily_table['codon'], subfamily_table['subfamily']):
        subfamily_dict[codon] = subfamily
    

   # 根据亚家族分组密码子
    codon_list = {}
    for index, row in subfamily_table.iterrows():
        subfamily = row['subfamily']
        codon = row['codon']
        if subfamily not in codon_list:
            codon_list[subfamily] = []
        codon_list[subfamily].append(codon)

    # 计算 ENC
    f_cf = pd.DataFrame(index=codonframe.index)
    n_cf = pd.DataFrame(index=codonframe.index)

    for subfamily, codons in codon_list.items():
        mx = codonframe[codons].values
        n = np.sum(mx, axis=1)
        p = (mx + 1) / (n.reshape(-1, 1) + len(codons))
        f_cf[subfamily] = np.sum(np.power(p, 2), axis=1)
        n_cf[subfamily] = n + len(codons)


    subfamily_group = [len(codons) for codons in codon_list.values()]
    N_single = np.sum(np.array(subfamily_group) == 1)
    N_double = np.sum(np.array(subfamily_group) == 2) * np.sum(n_cf.iloc[:, np.array(subfamily_group) == 2], axis=1) / np.sum(n_cf.iloc[:, np.array(subfamily_group) == 2] * f_cf.iloc[:, np.array(subfamily_group) == 2], axis=1)
    N_triple = np.sum(np.array(subfamily_group) == 3) * np.sum(n_cf.iloc[:, np.array(subfamily_group) == 3], axis=1) / np.sum(n_cf.iloc[:, np.array(subfamily_group) == 3] * f_cf.iloc[:, np.array(subfamily_group) == 3], axis=1)
    N_quad = np.sum(np.array(subfamily_group) == 4) * np.sum(n_cf.iloc[:, np.array(subfamily_group) == 4], axis=1) / np.sum(n_cf.iloc[:, np.array(subfamily_group) == 4] * f_cf.iloc[:, np.array(subfamily_group) == 4], axis=1)
    Nc = N_single + N_double + N_triple + N_quad
    enc = pd.DataFrame(Nc, columns=['enc'], index=codonframe.index)
    return enc


def calculate_tRNA_relative_adapt(tRNA_GCN, codon_table_index=1,
                  penalty_coefficient={'WC': 0, 'GU': 0.6295, 'UG': 0.7861, 
                                       'AI': 0.9075, 'CI': 0.4659, 'UI': 0}):
    """
    计算基因的TAI(Effective Number of Codons)值。

    参数:
    - tRNA_GCN (str): 输入的 CSV 文件路径。
    - codon_table_index (int): 选择的密码子表的序号, 默认为1。
    - penalty_coefficient (dic): 各类配对的惩罚系数(原文是selective constraint) TODO 自定义惩罚系数表

    返回:
    - tRNA_rel_adapt_table (pd.DataFrame): 包含每个codon的tRNA_rel_adapt_norm。
    """
    # 从文件中读取 tRNA_GCN 数据
    tRNA_GCN = pd.read_csv(tRNA_GCN)
    pair_table = codon_pair_table(codon_table_index)
    
    # 初始化 tRNA_GCN、penalty_coefficient 和 tRNA_rel_adapt 列
    pair_table['tRNA_GCN'] = 0
    # penalty_coefficient 和 tRNA_rel_adapt 列初始化为浮点型
    pair_table['penalty_coefficient'] = 0.0
    pair_table['tRNA_rel_adapt'] = 0.0
    # 只保留有拷贝数的tRNA
    # pair_table = pair_table[pair_table['anti_codon'].isin(list(tRNA_GCN['antiCodon']))].reset_index(drop=True)

    for index, row in pair_table.iterrows():
        anti_codon = row['anti_codon']
        tRNA_GCN_value = tRNA_GCN[tRNA_GCN['antiCodon'] == anti_codon]['copyCount'].values
        if len(tRNA_GCN_value) == 1:
            pair_table.at[index, 'tRNA_GCN'] = int(tRNA_GCN_value[0])
        elif len(tRNA_GCN_value) > 2:
            raise ValueError(f"Error in tRNA_GCN data: Expected 1 value for anti-codon {anti_codon}, "
                            f"but found {len(tRNA_GCN_value)} values.")
        
        pair_table.at[index, 'penalty_coefficient'] = penalty_coefficient[row['pair_type']]
    
        # 计算 tRNA_rel_adapt 并存储在 tRNA_rel_adapt 列中
        pair_table.at[index, 'tRNA_rel_adapt'] = (1 - pair_table.at[index, 'penalty_coefficient']) * pair_table.at[index, 'tRNA_GCN']
    
    tRNA_rel_adapt_table = pair_table.groupby('codon')['tRNA_rel_adapt'].sum().reset_index()
    # 按 'codon' 列汇总并重置索引
    # tRNA_rel_adapt_table = tRNA_rel_adapt_table.groupby('codon')['tRNA_rel_adapt'].sum().reset_index()

    tRNA_rel_adapt = tRNA_rel_adapt_table['tRNA_rel_adapt']
    tRNA_rel_adapt_table['tRNA_rel_adapt_norm'] = tRNA_rel_adapt / tRNA_rel_adapt.max()

    stop_codons_len = len(get_codon_table_by_index(codon_table_index).stop_codons)
    number_of_0 = sum(tRNA_rel_adapt==0)# 利用了Ture=1的特性，计算有几个值为0
    wi_mean = sum(tRNA_rel_adapt_table['tRNA_rel_adapt_norm']) / (64 - stop_codons_len - number_of_0)
    tRNA_rel_adapt_table['tRNA_rel_adapt_norm'] = tRNA_rel_adapt_table['tRNA_rel_adapt_norm'].apply(lambda x: (wi_mean if x == 0 else x))
    return tRNA_rel_adapt_table


def calculate_TAI(codonframe, tRNA_relative_adapt):
    """
    计算基因的TAI(tRNA adaptation index)值。

    参数:
    - codonframe (pd.DataFrame): 包含每个序列的ID以及每个密码子的数量的DataFrame。
    - tRNA_relative_adapt (pd.DataFrame): 包含每个codon tRNA_rel_adapt_norm的DataFrame。

    返回:
    - tai(pd.DataFrame): 包含每个基因的ENC值的DataFrame。
    """
    codonframe = codonframe[tRNA_relative_adapt['codon']] # 忽略终止密码子
    tai_values = np.exp(np.dot(codonframe, np.log(tRNA_relative_adapt['tRNA_rel_adapt_norm'])) / np.sum(codonframe, axis=1))
    
    tai = pd.DataFrame(tai_values, columns=['tai'], index=codonframe.index)
    return tai
    
# 示例用法：
if __name__ == "__main__":
    start_time = time.time()
    file_path = "../test_data/Drosophila_melanogaster.fna"
    codonframe = sequences_to_codonframe(file_path)
    # RSCU测试代码
    if False:
        rscuSheet = calculate_RSCU(codonframe, codon_table_index=1)
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
    if False:
        enc = calculate_ENC(codonframe)
        print(enc)

    # calculate_tRNA_relative_adapt测试代码
    if False:
        tRNA_relative_adapt = calculate_tRNA_relative_adapt(tRNA_GCN='../test_data/Drosophila_melanogaster_trna_counts.csv')
        print(tRNA_relative_adapt)
    
    # calculate_tRNA_relative_adap测试代码
    if False:
        tRNA_relative_adapt = calculate_tRNA_relative_adapt(tRNA_GCN='../test_data/Drosophila_melanogaster_trna_counts.csv')
        # print(tRNA_relative_adapt)
        # tRNA_relative_adapt.to_csv("test.csv")
        tai = calculate_TAI(codonframe, tRNA_relative_adapt)
        print(tai)

    
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"执行时间：{execution_time}秒")

