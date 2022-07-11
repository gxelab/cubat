import re
import time
import numpy as np
import pandas as pd
from Bio import SeqIO
# from Bio.Seq import Seq
from Bio.Data import CodonTable


def generate_info_dataframe(parse_data):
    id_list = []
    codon_frequency_list = []
    frequency_check = []
    start_codon_check = []
    for record in parse_data:
        id_list.append(record.id)
        codon_frequency_list.append(len(record.seq))
        start_codon_check.append(str(record.seq[0:3]))
        if len(record.seq) % 3 == 0:
            frequency_check.append("True")
        else:
            frequency_check.append("False")

    info_dataframe = pd.DataFrame({"id": id_list,
                                   "codon_frequency": codon_frequency_list,
                                   "Multiple of three?": frequency_check,
                                   "start codon": start_codon_check
                                   })
    return info_dataframe


def count_codon(seq_list):
    codon_num = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0, 'TAT': 0, 'TAC': 0,
                 'TAA': 0, 'TAG': 0, 'TGT': 0,
                 'TGC': 0, 'TGA': 0, 'TGG': 0, 'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0,
                 'CCG': 0,
                 'CAT': 0,
                 'CAC': 0, 'CAA': 0, 'CAG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0, 'CGG': 0, 'ATT': 0, 'ATC': 0, 'ATA': 0,
                 'ATG': 0,
                 'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0, 'AGT': 0, 'AGC': 0,
                 'AGA': 0,
                 'AGG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0, 'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0, 'GAT': 0,
                 'GAC': 0,
                 'GAA': 0, 'GAG': 0, 'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

    depart = len(seq_list) / 8
    times = int(depart)
    rest = int((depart - times) * 8)
    for i in range(times):
        codon_num[seq_list.pop()] += 1
        codon_num[seq_list.pop()] += 1
        codon_num[seq_list.pop()] += 1
        codon_num[seq_list.pop()] += 1
        codon_num[seq_list.pop()] += 1
        codon_num[seq_list.pop()] += 1
        codon_num[seq_list.pop()] += 1
    for i in range(rest):
        codon_num[seq_list.pop()] += 1

    codon_array = np.array(list(codon_num.values()))

    return codon_array


def generate_codon_dataframe(parse_data):
    parse_data = list(parse_data)
    times = len(parse_data)
    result_array = np.zeros((times, 64), dtype=int)
    extracted_name = []
    for t in range(0, times):
        extracted_name.append(parse_data[t].name)
        seq_str = str(parse_data[t].seq)
        seq_list = re.findall('.{3}', seq_str)
        result_array[t] = count_codon(seq_list)
    codon_dataframe = pd.DataFrame(result_array,
                                   columns=['TTT', 'TTC', 'TTA', 'TTG',
                                            'TCT', 'TCC', 'TCA', 'TCG',
                                            'TAT', 'TAC', 'TAA', 'TAG',
                                            'TGT', 'TGC', 'TGA', 'TGG',
                                            'CTT', 'CTC', 'CTA', 'CTG',
                                            'CCT', 'CCC', 'CCA', 'CCG',
                                            'CAT', 'CAC', 'CAA', 'CAG',
                                            'CGT', 'CGC', 'CGA', 'CGG',
                                            'ATT', 'ATC', 'ATA', 'ATG',
                                            'ACT', 'ACC', 'ACA', 'ACG',
                                            'AAT', 'AAC', 'AAA', 'AAG',
                                            'AGT', 'AGC', 'AGA', 'AGG',
                                            'GTT', 'GTC', 'GTA', 'GTG',
                                            'GCT', 'GCC', 'GCA', 'GCG',
                                            'GAT', 'GAC', 'GAA', 'GAG',
                                            'GGT', 'GGC', 'GGA', 'GGG'],
                                   index=extracted_name)

    return codon_dataframe


def rscu(codon_dataframe):
    stop_codon = []
    rscu_dataframe = pd.DataFrame(columns=codon_dataframe.columns, index=codon_dataframe.index)
    for codon in codon_dataframe.columns:
        try:
            for seq in codon_dataframe.index:
                # /
                amino_acid_value = CodonTable.standard_dna_table.forward_table[codon]
                rscu_amino_acid = \
                    [key for key, value in CodonTable.standard_dna_table.forward_table.items() if
                     value == amino_acid_value]
                codon_frequency = 0
                for rscu_codon in rscu_amino_acid:
                    codon_frequency += codon_dataframe.at[seq, rscu_codon]
                    if codon_frequency / len(rscu_amino_acid):
                        rscu_dataframe.at[seq, codon] = \
                            codon_dataframe.at[seq, codon] / (codon_frequency / len(rscu_amino_acid))

        except KeyError:
            stop_codon.append(codon)
    return rscu_dataframe


start = time.time()

Sars2 = SeqIO.parse("Test_Data/Sars_cov_2.ASM985889v3.cds.fasta", "fasta")
Sars2_codon_dataframe = generate_codon_dataframe(Sars2)
print(rscu(Sars2_codon_dataframe))


# print(Sars2_info_dataframe)
# print(Sars2_codon_dataframe)
print(time.time() - start)
