import pandas as pd
from Bio.Data import CodonTable


def get_codon_table_by_index(index=1):
    """
    根据密码子表的序号获取密码子表以及每个密码子表的起始密码子和终止密码子。

    参数:
    index (int): 密码子表的序号。

    返回:
    codon_table (Bio.Data.CodonTable): 密码子表对象。
    start_codon (str): 起始密码子。
    stop_codons (list): 终止密码子列表。
    """
    # 获取密码子表的所有可能表格
    all_codon_tables = CodonTable.unambiguous_dna_by_id

    # 确保给定的序号在有效范围内
    if index < 1 or index > len(all_codon_tables):
        raise ValueError("Invalid codon table index")

    # 获取指定序号的密码子表
    codon_table = all_codon_tables[index]

    # # 获取起始密码子和终止密码子
    # start_codon = codon_table.start_codons
    # stop_codons = codon_table.stop_codons
    # TODO 自定义密码子表
    return codon_table


def codon_table_completion(index=1):
    codon_table = CodonTable.unambiguous_dna_by_id[index].forward_table
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


def codon_table_subfamily(index= 1):
    """
    根据密码子表编号获取密码子表的子家族信息。

    参数:
    - index (str): 密码子表的编号，默认为 '1'。

    返回:
    - codon_table_subfamily (pd.DataFrame): 包含密码子子家族信息的DataFrame, 包括 aa_code, amino_acid, codon, 和 subfam 列。
    """
    codon_table = CodonTable.unambiguous_dna_by_id[int(index)]
    codon_table_dict = {
        'codon': [],
        'amino_acid': [],
        'subfamily': []
    }

    for codon, amino_acid in codon_table.forward_table.items():
        aa_code = CodonTable.ambiguous_dna_by_id[int(index)].forward_table[codon]
        subfamily = f"{amino_acid}_{codon[:2]}"
        codon_table_dict['codon'].append(codon)
        codon_table_dict['amino_acid'].append(amino_acid)
        codon_table_dict['subfamily'].append(subfamily)

    codon_table_subfamily = pd.DataFrame(codon_table_dict)
    return codon_table_subfamily


# 测试代码：
if __name__ == "__main__":
    index = 1  # 例如，选择标准密码子表
    codon_table = codon_table_subfamily(index)
    print(codon_table)
    