import pandas as pd
from Bio import SeqIO


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


def generate_codon_dataframe(parse_data):
    codon_dataframe = pd.DataFrame(
        columns=["TTT", "TTC", "TTA", "TTG",
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
                 "GGT", "GGC", "GGA", "GGG"])
    for record in parse_data:
        codon_dataframe.loc[record.id] = [0, 0, 0, 0,
                                          0, 0, 0, 0,
                                          0, 0, 0, 0,
                                          0, 0, 0, 0,
                                          0, 0, 0, 0,
                                          0, 0, 0, 0,
                                          0, 0, 0, 0,
                                          0, 0, 0, 0,
                                          0, 0, 0, 0,
                                          0, 0, 0, 0,
                                          0, 0, 0, 0,
                                          0, 0, 0, 0,
                                          0, 0, 0, 0,
                                          0, 0, 0, 0,
                                          0, 0, 0, 0,
                                          0, 0, 0, 0]
        for codon in [str(record.seq)[i:i + 3] for i in range(0, len(record.seq), 3)]:
            codon_dataframe.loc[record.id, codon == codon_dataframe.columns] += 1
    return codon_dataframe


Sars2 = SeqIO.parse("Sars_cov_2.ASM985889v3.cds.fasta", "fasta")
print(generate_info_dataframe(Sars2))
