import re
# import os
import math
import pandas as pd
import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.Data import CodonTable


class Cubat:
    i = 0
    correspondence = {}
    sequences = {}
    sequence_number = {}
    total_sequences = ''

    def __init__(self, data_path, codon_table=CodonTable.unambiguous_rna_by_id[1]):
        self.data_path = data_path
        self.filename = re.search(r'/.*', data_path).group()[1:]
        # print(self.filename)
        self.codon_table_forward_table = codon_table.forward_table
        self.codon_table = codon_table
        Cubat.get_data(self, self.data_path)

    def get_data(self, data_path):
        # Extract sequences and merge into one line string
        data = open(data_path, 'rb').read().decode('utf-8')  # Open file
        sequence = re.findall(r']((?:.|\n)*?)>', data)
        names = re.findall(r'>((?:.|\n)*?)]', data)
        for dna_seq in sequence:
            locals()['info' + str(Cubat.i)] = names[Cubat.i] + ']'
            locals()['seq' + str(Cubat.i)] = dna_seq
            locals()['seq' + str(Cubat.i)] = [seq.strip() for seq in locals()['seq' + str(Cubat.i)]]  # Remove '\n'
            locals()['seq' + str(Cubat.i)] = ''.join(locals()['seq' + str(Cubat.i)])  # Merge into one line string
            # Cubat.total_sequences.update = ''
            Cubat.correspondence.update({locals()['info' + str(Cubat.i)]: locals()['seq' + str(Cubat.i)]})
            Cubat.sequence_number.update({locals()['info' + str(Cubat.i)]: str(Cubat.i)})
            Cubat.sequences.update({Cubat.i: locals()['seq' + str(Cubat.i)]})
            Cubat.total_sequences = Cubat.total_sequences + ''.join(locals()['seq' + str(Cubat.i)])
            Cubat.i += 1
        Cubat.correspondence.update({(self.filename + '(total sequences)'): Cubat.total_sequences})
        Cubat.sequence_number.update({(self.filename + '(total sequences)'): 'total'})
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

    # @staticmethod
    # def locate_stop_codon(seq, stop_codon):
    #     # return seq.index(stop_codon)
    #     if seq.count('UGA1'):
    #         print(seq.index('UGA1'))
    #         return seq.index('UGA1')
    #     elif seq.count('UAA1'):
    #         print(seq.index('UAA'))
    #         return seq.index('UAA')
    #     elif seq.count('UAG'):
    #         print(seq.index('UAG'))
    #         return seq.index('UAG')

    # @staticmethod
    # def process_mrna_seq(cut_seq, stop_codon_location):
    #     return cut_seq[0:stop_codon_location]

    @staticmethod
    def count_codon(processed_mrna_seq):
        # Convert mRNA_seq into dataframe of pandas
        codon_series = pd.value_counts(processed_mrna_seq)
        codon_dic = {'codon': codon_series.index, 'Obsi': codon_series.values}
        # Obsi = Observed number of occurrences of codon 'i'
        codon_dataframe = pd.DataFrame(codon_dic)
        return codon_dataframe

    def generate_amino_dataframe(self, codon_dataframe):
        # codon_table = CodonTable.unambiguous_rna_by_id[1].forward_table
        amino_acid = codon_dataframe['codon'].map(self.codon_table_forward_table)
        amino_dic = {'amino_acid': amino_acid.values}
        amino_dataframe = pd.DataFrame(amino_dic)
        return amino_dataframe

    def generate_dataframe(self, data_info):
        data = Cubat.correspondence[data_info]
        dna_seq = Seq(data)  # Extract a sequence
        mrna_seq = dna_seq.transcribe()  # Transcribed into mRNA
        if len(dna_seq) % 3 != 0:  # Determine whether the sequence is a multiple of 3.
            print(data_info)
            print('\033[33m%s\033[0m' % 'Warning: the sequence is not a multiple of 3.')
        cut_seq = Cubat.cut(str(mrna_seq), 3)  # Cutting RNA sequence
        codon_dataframe = Cubat.count_codon(cut_seq)
        # Determine whether the sequence starts with the start codon in the codon table.
        if cut_seq[0] not in self.codon_table.start_codons:
            print(data_info)
            print('\033[33m%s\033[0m' % 'Warning: The first codon is not the start codon in the codon table.')

            # Otherwise, a function for locating the start codon will be called.
            start_codon = input('Please enter the start codon you need.')
            try:
                located_seq = Cubat.locate_start_codon(start_codon, mrna_seq)
            except AttributeError:
                print('\033[31m%s\033[0m' % 'Error: the start codon was not found.')
            else:
                cut_seq = Cubat.cut(located_seq, 3)  # Cutting RNA sequence

        # Determine whether the sequence ends with the stop codon in the codon table.
        if cut_seq[-1] not in self.codon_table.stop_codons:
            print(data_info)
            print('\033[33m%s\033[0m' % 'Warning: last codon is not stop codon.')

            # Otherwise, a function for locating the stop codon will be called.
            stop_codon = input('Please enter the stop codon you need.')
            try:
                stop_codon_location = cut_seq.index(stop_codon) + 1
                processed_mrna_seq = cut_seq[0:stop_codon_location]
            except ValueError:
                print('\033[31m%s\033[0m' % 'Error: the stop codon was not found.')
            else:
                codon_dataframe = Cubat.count_codon(processed_mrna_seq)

        amino_dataframe = Cubat.generate_amino_dataframe(self, codon_dataframe)
        merged_dataframe = codon_dataframe.merge(amino_dataframe, left_index=True, right_index=True)
        cols = merged_dataframe.columns[[0, 2, 1]]
        dataframe = merged_dataframe[cols]
        return dataframe

    def generate_dataframe_total(self):
        return Cubat.generate_dataframe(self, (self.filename + '(total sequences)'))

    def generate_dataframes(self, output_path):
        for key in Cubat.correspondence:
            file_path = output_path + Cubat.sequence_number[key] + '.csv'
            open(file_path, 'w')
            Cubat.generate_rscu_dataframe(self, Cubat.generate_dataframe(self, key), key) \
                .to_csv(file_path, sep=',', index=False, header=True)
            # Cubat.generate_dataframe(self, key).to_csv(file_path, sep=',', index=False, header=False)
            # print(key, '\n', Cubat.generate_dataframe(self, key))  # Cubat.sequences[key]))
        return 'Completed'  # 'Completed, ' + str(key) + ' sequences in total.'

    def info_rscu(self, rscu_codon, sequence_info):
        rscu_dataframe = Cubat.generate_dataframe(self, sequence_info)
        return Cubat.rscu(self, rscu_codon, rscu_dataframe)

    def rscu(self, rscu_codon, rscu_dataframe):
        try:
            rscu_amino_acid = list(rscu_dataframe.loc[rscu_dataframe['codon'] == rscu_codon, 'amino_acid'])[0]
            the_amino_dataframe = rscu_dataframe[(rscu_dataframe['amino_acid'] == rscu_amino_acid)]
            amino_sum = the_amino_dataframe['Obsi'].sum()
            # amino_sum_dataframe = the_amino_dataframe.groupby(by=['amino_acid"])["Obsi"].sum()
            if isinstance(rscu_amino_acid, str):
                if rscu_amino_acid == '*':
                    num_encoding_the_codons = len(self.codon_table.stop_codons)

                else:
                    num_encoding_the_codons = \
                        len([k for k, v in self.codon_table_forward_table.items() if v == rscu_amino_acid])
                # The number of all codons encoding this amino acid
                average_amino_acid_encoding = amino_sum / num_encoding_the_codons
                # The average number of uses of all codons encoding the amino ac
                the_codon_usage = the_amino_dataframe.loc[rscu_dataframe["codon"] == rscu_codon, "Obsi"].iloc[0]
                # Query the codon usage
                rscu_value = the_codon_usage / average_amino_acid_encoding
            else:
                rscu_value = 1.0

        except IndexError:
            # print('\033[31m%s\033[0m' % 'Error: the codon was not found.')
            return 'NaN'
        else:
            return rscu_value

    def generate_rscu_dataframe(self, codon_dataframe, sequence_info):
        rscu_list = []
        codon_list = codon_dataframe['codon'].tolist()
        for rscu_codon in codon_list:
            rscu_list.append(Cubat.info_rscu(self, rscu_codon, sequence_info))
        rscu_dataframe = codon_dataframe.merge(pd.DataFrame({'RSCU': rscu_list}), left_index=True, right_index=True)
        return rscu_dataframe

    @staticmethod
    def generate_pivot_table(dataframe, values):
        codons = dataframe['codon'].values.tolist()
        # values_list = dataframe[values].values.tolist()
        pivot_table = pd.DataFrame([])
        for codon in codons:
            pivot_table = pd.concat([pivot_table, dataframe.loc[dataframe['codon'] == codon]])
        pivot_table = pivot_table.pivot_table(index='codon', values=values, columns='amino_acid').fillna(0)
        print(type(pivot_table))
        return pivot_table

    # @staticmethod
    # def cai(seq_dataframe, high_expression_gene_table):
    #     return '---'


sars_cov_2 = Cubat('Test_Data/Sars_cov_2.ASM985889v3.cds.fasta')
# sars_cov_2_total = sars_cov_2.generate_dataframe_total()
sars_cov_2_seq1 = sars_cov_2.generate_dataframe('ENSSAST00005000002.1 cds primary_assembly:ASM985889v3:MN908947.3:266:21555:1 gene:ENSSASG00005000002.1 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:ORF1ab description:ORF1a polyprotein;ORF1ab polyprotein [Source:NCBI gene (formerly Entrezgene);Acc:43740578]')
sars_cov_2_seq1 = sars_cov_2_seq1.fillna('*')

Homo_sapiens_table = pd.read_excel('Reference_species/Homo_sapiens.xlsx')
Homo_sapiens_table = Homo_sapiens_table.fillna('*')
table_index = 0
for table_codon in Homo_sapiens_table['codon']:
    dna_codon = Seq(table_codon)  # Extract a sequence
    rna_codon = dna_codon.transcribe()  # Transcribed into mRNA
    # Homo_sapiens_table.loc[Homo_sapiens_table['codon'][table_index]] = str(table_codon)
    Homo_sapiens_table.loc[Homo_sapiens_table['codon'] == table_codon, 'codon'] = str(rna_codon)
    table_index += 1

wij_dataframe = pd.DataFrame(columns=('codon', 'wij'))
for wij_codon in sars_cov_2_seq1['codon']:
    rac_amino_acid = list(sars_cov_2_seq1[(sars_cov_2_seq1['codon'] == wij_codon)]['amino_acid'])[0]
    large_codon = Homo_sapiens_table.loc[Homo_sapiens_table['amino_acid'] == rac_amino_acid].sort_values(by='Obsi', ascending=False).reset_index(drop=True)['codon'][0]
    wij = int(Homo_sapiens_table[Homo_sapiens_table.codon == wij_codon]['Obsi']) / int(Homo_sapiens_table[Homo_sapiens_table.codon == large_codon]['Obsi'])
    wij_dataframe = wij_dataframe.append(pd.DataFrame({'codon': [wij_codon], 'wij': [wij]}), ignore_index=True)
cai_numerator = 0
cai_denominator = 0
for fij_codon in sars_cov_2_seq1['codon']:
    cai_numerator = cai_numerator + int(sars_cov_2_seq1[sars_cov_2_seq1.codon == fij_codon]['Obsi']) * np.log(float(wij_dataframe[wij_dataframe.codon == fij_codon]['wij']))
    cai_denominator = cai_denominator + int(sars_cov_2_seq1[sars_cov_2_seq1.codon == fij_codon]['Obsi'])
print(cai_numerator)
print(cai_denominator)
cai = math.exp(cai_numerator/cai_denominator)
print(cai)
