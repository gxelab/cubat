import re
from collections import Counter

import numpy as np
import pandas as pd
from Bio import SeqIO
from codon_table import get_codon_table_by_index


def read_sequences(file_path):
    """
    从序列文件中读取所有序列并返回一个列表。

    参数:
    file_path (str): 序列文件的路径。

    返回:
    seq_records (list): 所有读取的序列对象的列表。
    """

    try:
        seq_records = list(SeqIO.parse(file_path, "fasta"))
        return seq_records
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {file_path}")


def qa(seq_records, codon_table=1, min_len=6, check_len=True, check_start=True, check_stop=True, check_istop=True):
    """
    读取序列的同时对序列进行质量鉴定。

    参数:
    seq_records (list): 包含序列记录的列表。
    get_codon_table_by_index (index): 密码子表对象。
    min_len (int): 最小序列长度,默认为6。
    check_len (bool): 是否检查序列长度是否为3的倍数,默认为True。
    check_start (bool): 是否检查起始密码子,默认为True。
    check_stop (bool): 是否检查终止密码子,默认为True。
    check_istop (bool): 是否检查内部终止密码子,默认为True。

    返回:
    problem_seqs (list): 经过质量鉴定和处理后的序列记录列表。
    warnings (list): 警告信息列表。
    """

    codon_table = get_codon_table_by_index(codon_table)
    warnings = []
    warnings_types = []
    problem_seqs = []

    for record in seq_records:
        seq = record.seq

        # 检查序列长度是否满足最小长度要求
        if len(seq) < min_len:
            warnings.append(f"Sequence '{record.id}': Length is less than the minimum required length ({min_len}).")
            warnings_types.append("Shorter than minimum length")
            problem_seqs.append(record)

        # 检查序列长度是否为3的倍数（检查CDS长度）
        if check_len and len(seq) % 3 != 0:
            warnings.append(f"Sequence '{record.id}': The length is not a multiple of 3.")
            warnings_types.append("Not a multiple of 3")
            problem_seqs.append(record)

        # 检查起始密码子
        if check_start and seq[:3].upper() not in codon_table.start_codons:
            warnings.append(f"Sequence '{record.id}': The start codon does not match the selected codon table.")
            warnings_types.append("Start codon mismatch")
            problem_seqs.append(record)

        # 检查终止密码子
        if check_stop:
            valid_stop_codons = codon_table.stop_codons
            if not any(seq[-3:].upper() == stop_codon for stop_codon in valid_stop_codons):
                warnings.append(f"Sequence '{record.id}': The stop codon does not match the selected codon table.")
                warnings_types.append("Stop codon mismatch")
                problem_seqs.append(record)

        
        # 检查内部终止密码子
        if check_istop:
            valid_stop_codons = codon_table.stop_codons
            seq_str = str(seq)
            positions = []

            for stop_codon in valid_stop_codons:
                stop_len = len(stop_codon)
                for i in range(0, len(seq_str) - stop_len + 1, 3):
                    codon = seq_str[i:i+stop_len]
                    if codon == stop_codon:
                        positions.append(i)

        if len(positions) > 1:
            warnings.append(f"Sequence '{record.id}' Contains internal stop codons.")
            warnings_types.append("Internal stop codons")
            problem_seqs.append(record)


    # # 打印警告信息（黄色文本）
    # for warning in warnings:
    #     print("\033[93mWarning:", warning, "\033[0m")

    return problem_seqs, warnings, warnings_types



def count_codon(sequence):
    '''
    避免重复不断复制和重建DataFrame,
    收集每个密码子计数字典,返回,然后一次性创建一个DataFrame.

    参数:
    sequence (seq): 某个序列。

    返回:
    codon_dict (dic): 密码子计数字典
    '''
    codon_num = {
    'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0,
    'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0,
    'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
    'TGT': 0, 'TGC': 0, 'TGA': 0, 'TGG': 0,
    'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0,
    'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
    'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0,
    'CGT': 0, 'CGC': 0, 'CGA': 0, 'CGG': 0,
    'ATT': 0, 'ATC': 0, 'ATA': 0, 'ATG': 0,
    'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0,
    'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0,
    'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
    'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0,
    'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0,
    'GAT': 0, 'GAC': 0, 'GAA': 0, 'GAG': 0,
    'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0
}
    # 使用Counter计算密码子数量,它内部优化了计数操作。
    # 使用正则表达式拆分序列:免去循环,较为高效。
    codon_counts = Counter(re.findall(r'.{3}', sequence))
    codon_num.update(codon_counts)

    codon_dict = dict(codon_num)
    return codon_dict

def sequences_to_codonframe(seq_records, codon_table=1, exclude_qa_problems=False, display_warnings=False):
    """
    将DNA序列记录列表转换为密码子计数的DataFrame。
    该函数接受DNA序列记录列表和一个可选的排除QA问题序列的标志,然后执行以下操作：
    1. 初始化密码子列表和空的DataFrame。
    2. 遍历输入的DNA序列记录列表。
    3. 对于每个记录,如果启用了排除QA问题序列并且该记录在QA阶段被鉴定为问题序列,则跳过。
    4. 否则,使用count_codon函数计算DNA序列的密码子数量,并将其存储为字典。
    5. 将ID添加到字典中。
    6. 将字典添加到数据列表中。
    7. 最后,将数据列表转换为DataFrame,其中包含每个序列的ID以及每个密码子的数量。

    参数:
    - seq_records (list): 包含DNA序列记录的列表,每个记录包括ID和序列。
    - exclude_qa_problems (bool): 是否排除在QA阶段被鉴定为问题序列,默认为False。

    返回:
    - codonframe (pd.DataFrame): 包含每个序列的ID以及每个密码子的数量的DataFrame。
    """

    # 这里只是用于获取密码子列表，免去在程序里再打一遍。
    codon_list = get_codon_table_by_index(codon_table ).forward_table.keys()
    codonframe = pd.DataFrame(columns=["ID"] + list(codon_list))
    warning_counts = {}  # 用于统计不同类型警告的数量

    data = []
    for record in seq_records:
        seq_id = record.id
        seq = str(record.seq)

        codon_counts = count_codon(seq)
        codon_counts["ID"] = seq_id

        _, warnings, warnings_types = qa([record], codon_table)  # 调用 qa 函数并获取警告
        
        if warnings and display_warnings:
            for warning in warnings:
                print("\033[93mWarning:", warning, "\033[0m")

        # 分类统计警告类型
        if warnings_types:
            warning_key = '; '.join(warnings_types)
            if warning_key in warning_counts:
                warning_counts[warning_key] += 1
            else:
                warning_counts[warning_key] = 1
        # 加入codonframe
        if not exclude_qa_problems or not warnings:
            data.append(codon_counts)

    # 打印不同类型的警告及其数量
    print("\033[93mWarning:\033[0m\n")
    for warning, count in warning_counts.items():
        print(f"'{warning}': {count} times")
    print("\n\nsequences has been loaded")
    codonframe = pd.DataFrame(data)
    codonframe = codonframe.set_index("ID")
    # TODO 加一个类似的dataframe看每条序列的问题
    return codonframe

# '''

# 用法示例：

# 选择文件路径
file_path = "../test_data/Drosophila_melanogaster.fna"
# file_path = "../test_data/test.fasta"

# qa
seq_records = read_sequences(file_path)
# problem_seqs, warnings = qa(seq_records, 1)


# # 遍历并打印质量鉴定后的序列名称
# for record in problem_seqs:
#     print("问题序列名称:", record.id)
# print(warnings)

# codonframe计算
# 将序列转换为codonframe（默认情况下不排除QA问题序列）
codonframe = sequences_to_codonframe(seq_records)
codonframe.to_csv("test.csv")
# # 打印codonframe
print(codonframe)
# '''
# print(get_codon_table_by_index(1))
