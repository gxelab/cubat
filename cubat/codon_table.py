import pandas as pd
from Bio.Seq import Seq
from Bio.Data import CodonTable


def get_codon_table_by_index(codon_table_index=1):
    """
    根据密码子表的序号获取密码子表以及每个密码子表的起始密码子和终止密码子。

    参数:
    codon_table_index (str): 密码子表的序号。

    返回:
    codon_table (Bio.Data.CodonTable): 密码子表对象。
    start_codon (str): 起始密码子。
    stop_codons (list): 终止密码子列表。
    """
    # 获取密码子表的所有可能表格
    all_codon_tables = CodonTable.unambiguous_dna_by_id

    # 确保给定的序号在有效范围内
    if codon_table_index < 1 or codon_table_index > len(all_codon_tables):
        raise ValueError("Invalid codon table index")

    # 获取指定序号的密码子表
    codon_table = all_codon_tables[int(codon_table_index)]

    # # 获取起始密码子和终止密码子
    # start_codon = codon_table.start_codons
    # stop_codons = codon_table.stop_codons
    # TODO 自定义密码子表
    return codon_table


def codon_table_completion(codon_table_index=1):
    """
    由于biopyhton的密码子表不包含终止密码子, 所以需要一个函数来返回完整密码子表。

    参数:
    - codon_table_index (str): 密码子表的编号，默认为 '1'。

    返回:
    - complete_codon_table (pd.DataFrame): 完整的密码子表。
    """
    codon_table = CodonTable.unambiguous_dna_by_id[int(codon_table_index)].forward_table
    codons = ["TTT", "TTC", "TTA", "TTG",
              "TCT", "TCC", "TCA", "TCG",
              "TAT", "TAC", "TAA", "TAG",
              "TGT", "TGC", "TGA", "TGG",
              "CTT", "CTC", "CTA", "CTG",
              "CCT", "CCC", "CCA", "CCG",
              "CAT", "CAC", "CAA", "CAG",
              "CGT", "CGC", "CGA", "CGG",
              "ATT", "ATC", "ATA", "ATG",
              "ACT", "ACC", "ACA", "ACG",
              "AAT", "AAC", "AAA", "AAG",
              "AGT", "AGC", "AGA", "AGG",
              "GTT", "GTC", "GTA", "GTG",
              "GCT", "GCC", "GCA", "GCG",
              "GAT", "GAC", "GAA", "GAG",
              "GGT", "GGC", "GGA", "GGG"]
    
    complete_codon_table = codon_table.copy()
    for codon in codons:
        if codon not in complete_codon_table:
            complete_codon_table[codon] = "*"

    return complete_codon_table


def codon_table_subfamily(codon_table_index= 1):
    """
    根据密码子表编号获取密码子表的子家族信息。

    参数:
    - codon_table_index (str): 密码子表的编号，默认为 '1'。

    返回:
    - codon_table_subfamily (pd.DataFrame): 包含密码子子家族信息的DataFrame, 包括 aa_code, amino_acid, codon, 和 subfam 列。
    """
    codon_table = CodonTable.unambiguous_dna_by_id[int(codon_table_index)]
    codon_table_dict = {
        'codon': [],
        'amino_acid': [],
        'subfamily': []
    }

    for codon, amino_acid in codon_table.forward_table.items():
        subfamily = f"{amino_acid}_{codon[:2]}"
        codon_table_dict['codon'].append(codon)
        codon_table_dict['amino_acid'].append(amino_acid)
        codon_table_dict['subfamily'].append(subfamily)

    codon_table_subfamily = pd.DataFrame(codon_table_dict)
    return codon_table_subfamily


def codon_pair_table(codon_table_index= 1):
    """
    根据密码子表编号获取密码子表的配对信息。

    参数:
    - codon_table_index (str): 密码子表的编号，默认为 '1'。

    返回:
    - codon_pair_table (pd.DataFrame): 包含密码子配对信息的DataFrame, 包含codon, anti_codon, type(配对类型)。
    """
    codon_table = CodonTable.unambiguous_dna_by_id[int(codon_table_index)] # 这里使用的是RNA的密码子表
    codon_table_dict = {
        'strict_base': [], # 前两位碱基严格配对
        'wobble_base': [], # 摆动配对的碱基
        'codon': [],
        'anti_codon': [],
        'pair_type': [],
        'codon_amino_acid': [],
        'anti_amino_acid': []
    }

    for codon, amino_acid in codon_table.forward_table.items():
        # 前两位碱基严格配对，后一位碱基为摆动配对
        strict_base, wobble_base = codon[0:2], codon[2]
        codon_table_dict['strict_base'].append(strict_base)
        codon_table_dict['wobble_base'].append(wobble_base)

        codon_table_dict['codon'].append(codon)
        anti_codon = str(Seq(codon).reverse_complement()) # 反向互补序列
        codon_table_dict['anti_codon'].append(anti_codon)
        
        codon_table_dict['pair_type'].append('WC') # (standard)'
        codon_table_dict['codon_amino_acid'].append(amino_acid)
        codon_table_dict['anti_amino_acid'].append(amino_acid)

    standard_table = pd.DataFrame(codon_table_dict)

    wobble_table = pd.DataFrame(columns=standard_table.columns.tolist())
    for index, row in standard_table.iterrows():
        if row['wobble_base'] in ['T', 'G']:
            # 这里是GU或者UG配对的anti codon
            anti_codon = ('G' if row['wobble_base'] == 'T' else 'T') + row['anti_codon'][1:3]
            wobble_name = 'UG' if row['wobble_base'] == 'T' else 'GU'# (wobble)'
        elif row['wobble_base'] in ['A', 'C', 'T']:
            # 这里是AI, CI, 或者UI配对的anti codon
            anti_codon = 'A' + row['anti_codon'][1:3]
            
            wobble_name = row['wobble_base'] + 'I'# (wobble)'
        # 如果是终止密码子那么就跳过     
        if str(Seq(anti_codon).reverse_complement()) in codon_table.stop_codons:
            continue

        # 相关信息按行添加进wobble_table
        wobble_table.loc[index] = [row['strict_base'], 
                                row['wobble_base'],
                                row['codon'],
                                anti_codon,
                                wobble_name,
                                row['codon_amino_acid'],
                                codon_table.forward_table[str(Seq(anti_codon).reverse_complement())]
                                ]

    # 删除codon和anti codon氨基酸不一样的情况
    wobble_table = wobble_table[wobble_table['codon_amino_acid'] == wobble_table['anti_amino_acid']]

    # 合并dataframe并且进行整理
    codon_pair_table = pd.concat([standard_table, wobble_table], axis=0)
    codon_pair_table = codon_pair_table.drop(columns=['strict_base', 'wobble_base', 'anti_amino_acid'])
    codon_pair_table.rename(columns={'codon_amino_acid': 'amino_acid'}, inplace=True)
    codon_pair_table = codon_pair_table.reset_index(drop=True)

    return codon_pair_table


# 测试代码：
if __name__ == "__main__":
    codon_table_index = 1  # 例如，选择标准密码子表
    codon_table = codon_pair_table(codon_table_index)
    print(codon_table)
    