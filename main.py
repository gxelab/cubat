import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.Data import CodonTable


class Cubat:
    i = 0
    correspondence = {}
    sequences = {}

    def __init__(self, data_path):
        self.data_path = data_path
        Cubat.get_data(self.data_path)

    @staticmethod
    def get_data(data_path):
        # Extract sequences and merge into one line string
        data = open(data_path, 'rb').read().decode('utf-8')  # Open file
        sequence = re.findall(r']((?:.|\n)*?)>', data)
        names = re.findall(r'\[((?:.|\n)*?)]', data)

        for dna_seq in sequence:
            Cubat.i += 1
            locals()['info' + str(Cubat.i)] = names[Cubat.i]
            locals()['seq' + str(Cubat.i)] = dna_seq
            locals()['seq' + str(Cubat.i)] = [seq.strip() for seq in locals()['seq' + str(Cubat.i)]]  # Remove '\n'
            locals()['seq' + str(Cubat.i)] = ''.join(locals()['seq' + str(Cubat.i)])  # Merge into one line string
            Cubat.correspondence.update({locals()['info' + str(Cubat.i)]: locals()['seq' + str(Cubat.i)]})
            Cubat.sequences.update({Cubat.i: locals()['seq' + str(Cubat.i)]})
        return Cubat.sequences

    @staticmethod
    def cut(obj, sec):
        # Function for cutting sequence
        return [obj[i:i + sec] for i in range(0, len(obj), sec)]

    @staticmethod
    def locate_start_codon(selected_start_codon, mrna_seq):
        start_codon = selected_start_codon + '.*'
        jump_to_starting_point = re.search(start_codon, str(mrna_seq))  # Jump to start point
        starting_point = jump_to_starting_point.group(0)
        return starting_point

    @staticmethod
    def locate_stop_codon(seq):
        if seq.count('UGA'):
            return seq.index('UGA')
        elif seq.count('UAA'):
            return seq.index('UAA')
        elif seq.count('UAG'):
            return seq.index('UAG')

    @staticmethod
    def process_mrna_seq(cut_seq, stop_codon_location):
        return cut_seq[0:stop_codon_location]

    @staticmethod
    def count_codon(processed_mrna_seq):
        # Convert mRNA_seq into dataframe of pandas
        codon_series = pd.value_counts(processed_mrna_seq)
        codon_dic = {'codon': codon_series.index, 'quantity': codon_series.values}
        codon_dataframe = pd.DataFrame(codon_dic)
        return codon_dataframe

    @staticmethod
    def generate_amino_dataframe(codon_dataframe):
        codon_table = CodonTable.unambiguous_rna_by_id[1].forward_table
        amino_acid = codon_dataframe['codon'].map(codon_table)
        amino_dic = {'amino_acid': amino_acid.values}
        amino_dataframe = pd.DataFrame(amino_dic)
        return amino_dataframe

    @staticmethod
    def generate_dataframe(data):
        dna_seq = Seq(data)  # Extract a sequence
        mrna_seq = dna_seq.transcribe()  # Transcribed into mRNA
        located_seq = Cubat.locate_start_codon('AUG', mrna_seq)
        cut_seq = Cubat.cut(located_seq, 3)  # Cutting RNA sequence
        stop_codon_location = Cubat.locate_stop_codon(cut_seq)
        processed_mrna_seq = Cubat.process_mrna_seq(cut_seq, stop_codon_location)
        codon_dataframe = Cubat.count_codon(processed_mrna_seq)
        amino_dataframe = Cubat.generate_amino_dataframe(codon_dataframe)
        merged_dataframe = codon_dataframe.merge(amino_dataframe, left_index=True, right_index=True)
        cols = merged_dataframe.columns[[0, 2, 1]]
        dataframe = merged_dataframe[cols]
        return dataframe

    @staticmethod
    def generate_dataframes():
        key = 0
        for key in Cubat.sequences:
            print(Cubat.generate_dataframe(Cubat.sequences[key]))
        return 'Complete, ' + str(key) + ' sequences in total.'


sars_cov_2 = Cubat('Test_Data/Sars_cov_2.ASM985889v3.cds.fasta')
print(sars_cov_2.generate_dataframes())

# sars_cov_2_seq1 = sars_cov_2.generate_dataframe(sars_cov_2.sequences[1])
# plt.rcParams['figure.figsize'] = (12, 14)
# codons = sars_cov_2_seq1['codon'].values.tolist()
# quantity = sars_cov_2_seq1['quantity'].values.tolist()
# sars_cov_2_seq1_pivot_table = pd.DataFrame([])
# for codon in codons:
#     sars_cov_2_seq1_pivot_table = pd.concat([sars_cov_2_seq1_pivot_table,
#                                              sars_cov_2_seq1.loc[sars_cov_2_seq1['codon'] == codon]])
# sars_cov_2_seq1_pivot_table = sars_cov_2_seq1_pivot_table.pivot_table(index='codon', values='quantity',
#                                                                       columns='amino_acid').fillna(0)
# sns.heatmap(sars_cov_2_seq1_pivot_table, cmap=plt.cm.Reds, linewidths=0.01)
# 
# plt.bar(codons, quantity, color='pink')
# plt.xticks(fontsize=7)
# plt.show()
