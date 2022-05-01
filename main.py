import re
# import os
# import math
import pandas as pd
import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio import SeqIO
# import Bio.Data.CodonTable


class Cubat:
    name = []
    sequences = []
    total_sequences = ''
    part_name = []
    length=0
    times=0
    genecode = 0

    def __init__(self, data_path, gene_code=1):
        self.data_path = data_path
        fasta_list = list(SeqIO.parse(data_path, "fasta"))
        self.length=len(fasta_list)
        self.genecode=gene_code
        for i in fasta_list:
            self.name.append(i.name)
            self.sequences.append(i.seq)
        codon_table = CodonTable.unambiguous_dna_by_id[gene_code]
        self.codon_table_forward_table = codon_table.forward_table
        self.codon_table = codon_table
        self.part_name=[]
        Cubat.get_data(self, self.data_path)

#     def get_data(self, data_path):
#         # Extract sequences and merge into one line string
#         data = open(data_path, 'rb').read().decode('utf-8')  # Open file
#         sequence = re.findall(r']((?:.|\n)*?)>', data)
#         names = re.findall(r'>((?:.|\n)*?)]', data)
#         for dna_seq in sequence:
#             locals()['info' + str(Cubat.i)] = names[Cubat.i] + ']'
#             locals()['seq' + str(Cubat.i)] = dna_seq
#             locals()['seq' + str(Cubat.i)] = [seq.strip() for seq in locals()['seq' + str(Cubat.i)]]  # Remove '\n'
#             locals()['seq' + str(Cubat.i)] = ''.join(locals()['seq' + str(Cubat.i)])  # Merge into one line string
#             # Cubat.total_sequences.update = ''
#             Cubat.correspondence.update({locals()['info' + str(Cubat.i)]: locals()['seq' + str(Cubat.i)]})
#             Cubat.sequence_number.update({locals()['info' + str(Cubat.i)]: str(Cubat.i)})
#             Cubat.sequences.update({Cubat.i: locals()['seq' + str(Cubat.i)]})
#             Cubat.total_sequences = Cubat.total_sequences + ''.join(locals()['seq' + str(Cubat.i)])
#             Cubat.i += 1
#         Cubat.correspondence.update({(self.filename + '(total sequences)'): Cubat.total_sequences})
#         Cubat.sequence_number.update({(self.filename + '(total sequences)'): 'total'})
#         return Cubat.sequences

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

#     @staticmethod
#     def count_codon(processed_mrna_seq):
#         # Convert mRNA_seq into dataframe of pandas
#         codon_series = pd.value_counts(processed_mrna_seq)
#         codon_dic = {'codon': codon_series.index, 'Obsi': codon_series.values}
#         # Obsi = Observed number of occurrences of codon 'i'
#         codon_dataframe = pd.DataFrame(codon_dic)
#         return codon_dataframe

    @staticmethod
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

        while len(seq_list) > 0:
            codon_num[seq_list.pop(0)] += 1
        codon_array = np.array(list(codon_num.values()))
        return codon_array


    def generate_amino_dataframe(self, codon_dataframe):
        # codon_table = CodonTable.unambiguous_rna_by_id[1].forward_table
        amino_acid = codon_dataframe['codon'].map(self.codon_table_forward_table)
        amino_dic = {'amino_acid': amino_acid.values}
        amino_dataframe = pd.DataFrame(amino_dic)
        return amino_dataframe

#     def generate_dataframe(self, data_info):
#         data = Cubat.correspondence[data_info]
#         dna_seq = Seq(data)  # Extract a sequence
#         mrna_seq = dna_seq.transcribe()  # Transcribed into mRNA
#         if len(dna_seq) % 3 != 0:  # Determine whether the sequence is a multiple of 3.
#             print(data_info)
#             print('\033[33m%s\033[0m' % 'Warning: the sequence is not a multiple of 3.')
#         cut_seq = Cubat.cut(str(mrna_seq), 3)  # Cutting RNA sequence
#         codon_dataframe = Cubat.count_codon(cut_seq)
#         # Determine whether the sequence starts with the start codon in the codon table.
#         if cut_seq[0] not in self.codon_table.start_codons:
#             print(data_info)
#             print('\033[33m%s\033[0m' % 'Warning: The first codon is not the start codon in the codon table.')

#             # Otherwise, a function for locating the start codon will be called.
#             start_codon = input('Please enter the start codon you need.')
#             try:
#                 located_seq = Cubat.locate_start_codon(start_codon, mrna_seq)
#             except AttributeError:
#                 print('\033[31m%s\033[0m' % 'Error: the start codon was not found.')
#             else:
#                 cut_seq = Cubat.cut(located_seq, 3)  # Cutting RNA sequence

#         # Determine whether the sequence ends with the stop codon in the codon table.
#         if cut_seq[-1] not in self.codon_table.stop_codons:
#             print(data_info)
#             print('\033[33m%s\033[0m' % 'Warning: last codon is not stop codon.')

#             # Otherwise, a function for locating the stop codon will be called.
#             stop_codon = input('Please enter the stop codon you need.')
#             try:
#                 stop_codon_location = cut_seq.index(stop_codon) + 1
#                 processed_mrna_seq = cut_seq[0:stop_codon_location]
#             except ValueError:
#                 print('\033[31m%s\033[0m' % 'Error: the stop codon was not found.')
#             else:
#                 codon_dataframe = Cubat.count_codon(processed_mrna_seq)

#         amino_dataframe = Cubat.generate_amino_dataframe(self, codon_dataframe)
#         merged_dataframe = codon_dataframe.merge(amino_dataframe, left_index=True, right_index=True)
#         cols = merged_dataframe.columns[[0, 2, 1]]
#         dataframe = merged_dataframe[cols]
#         return dataframe.fillna('*')

    def get_sequences(self, data_info=None, method="all"):
        data_input = []


        if method == "part":
            if type(data_info) == int:
                data_input.append(self.sequences[data_info])
                self.part_name.append(self.name[data_info])

            if type(data_info) == str:
                try:
                    index=Cubat.name.index(data_info)
                except ValueError:
                    print("ERROR.Its that we cant find the name you signed in" )
                data_input.append(self.sequences[index])
                self.part_name.append(data_info)

            if type(data_info) == list:
                for i in data_info:
                    if type(i) == str:
                        try:
                            index = Cubat.name.index(i)
                        except ValueError:
                            print("ERROR.Its that the list may have a name that couldnt be found")
                        data_input.append(Cubat.sequences[index])
                        self.part_name.append(i)
                    if type(i) == int:
                        data_input.append(Cubat.sequences[i])
                        self.part_name.append(self.name[i])

        elif method == 'all':
            data_input = list(self.sequences)
            self.part_name = list(self.name)
        Cubat.times = len(data_input)
        self.times = len(data_input)

        return data_input


    def process_data(self,data_input, QC=False):
        lack_start = np.zeros( self.times, dtype=np.bool_)
        rearr = np.zeros((self.times, 64), dtype=float)
        for t in range(0, self.times):
            seq_str=str(data_input[t])
            seq_list = re.findall('.{3}', seq_str)
            if QC and seq_list[0] not in self.codon_table.start_codons:
                lack_start[t] = True
            rearr[t] = self.count_codon(seq_list)
        result = (rearr,)
        if QC:
            not_mut3 = (np.sum(rearr, axis=1) % 3 != 0)
            inner_end = (np.sum(rearr[:, np.array(code_end)], axis=1) != 0)
            result += (not_mut3,) + (inner_end,) + (lack_start,)
        return result

    def generate_frequency_dataframe(self, count_arr):
        codon_dataframe = pd.DataFrame(count_arr[0],
                                       columns=['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC',
                                                'TAA', 'TAG', 'TGT', 'TGC', 'TGA', 'TGG',
                                                'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC',
                                                'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG',
                                                'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC',
                                                'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG',
                                                'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC',
                                                'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG'], index=self.part_name)
        return codon_dataframe

    def generate_dataframe(self, data_input, quality_control=True,enc=True):
        result_dataframe=pd.DataFrame()

        if quality_control:
            count_arr, lack_start, not_mut3, inner_stop = Cubat.process_data(self,data_input, QC=True)
            concat=np.array([lack_start,not_mut3,inner_stop]).T
            QC_dataframe=pd.DataFrame(concat,columns=['lack_start_codon','not_multipler_of_3','inner_stop_codon'],index=self.part_name)
            result_dataframe=pd.concat([result_dataframe,QC_dataframe],axis=1)


        else:
            count_arr = Cubat.process_data(self,data_input,QC=False)

        if enc:
            enc_arr = Cubat.enc(self, count_arr,1).T
            enc_dataframe=pd.DataFrame(enc_arr,columns=['enc'],index=self.part_name)
            result_dataframe = pd.concat([result_dataframe, enc_dataframe],axis=1)

        return result_dataframe


    def enc(self,frequency_array, genecode=1):
        if type(frequency_array)==tuple:
            frequency_array=frequency_array[0]
        process_arr = np.delete(frequency_array, code_del[genecode], axis=1)
        FcF_array, an_array = np.empty((21, self.times), dtype=float), np.empty((21, self.times), dtype=float)
        fly = len(all_24familys)
        num1xa = process_arr
        for i in range(0, fly):
            F_array = num1xa[:, codon24[i]:codon24[i + 1]]
            an = np.sum(F_array, axis=1)
            an_array[i] = an
            FcF_array[i] = np.sum(np.square(F_array + 1), axis=1) / np.square(an + all_24familys[i])

        family2, family3, family4 = (all_24familys == 2), (all_24familys == 3), (all_24familys == 4)
        family2sum = code_K[genecode][0] * np.sum(an_array[family2], axis=0) / np.sum(
            np.multiply(FcF_array[family2], an_array[family2]), axis=0)
        family3_prevent0 = np.sum(np.multiply(FcF_array[family3], an_array[family3]), axis=0)
        family3_prevent0[family3_prevent0 == 0] = 1
        family3sum = code_K[genecode][1] * np.sum(an_array[family3], axis=0) / family3_prevent0
        family4sum = code_K[genecode][2] * np.sum(an_array[family4], axis=0) / np.sum(
            np.multiply(FcF_array[family4], an_array[family4]), axis=0)

        enc_array = code_ns[genecode] + family2sum + family3sum + family4sum

        return enc_array


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

    @staticmethod
    def cai(seq_dataframe, high_expression_gene_table=pd.read_excel('Reference_species/Homo_sapiens.xlsx')):
        high_expression_gene_table = high_expression_gene_table.fillna('*')
        table_index = 0
        for table_codon in high_expression_gene_table['codon']:
            dna_codon = Seq(table_codon)  # Extract a sequence
            rna_codon = dna_codon.transcribe()  # Transcribed into mRNA
            # high_expression_gene_table.loc[high_expression_gene_table['codon'][table_index]] = str(table_codon)
            high_expression_gene_table.loc[high_expression_gene_table['codon'] == table_codon, 'codon'] = str(rna_codon)
            table_index += 1

        wij_dataframe = pd.DataFrame(columns=('codon', 'wij'))
        for wij_codon in seq_dataframe['codon']:

            wij_amino_acid = list(seq_dataframe[(seq_dataframe['codon'] == wij_codon)]['amino_acid'])[0]
            large_codon = high_expression_gene_table.loc[
                high_expression_gene_table['amino_acid'] == wij_amino_acid
            ].sort_values(by='Obsi', ascending=False).reset_index(drop=True)['codon'][0]

            wij = int(high_expression_gene_table[high_expression_gene_table.codon == wij_codon]['Obsi']) / int(
                high_expression_gene_table[high_expression_gene_table.codon == large_codon]['Obsi'])

            wij_dataframe = wij_dataframe.append(pd.DataFrame({'codon': [wij_codon], 'wij': [wij]}), ignore_index=True)
        cai_numerator = 0
        cai2_numerator = 0
        cai_denominator = 0
        for fij_codon in seq_dataframe['codon']:
            cai_numerator = cai_numerator + int(seq_dataframe[seq_dataframe.codon == fij_codon]['Obsi']) * np.log(
                float(wij_dataframe[wij_dataframe.codon == fij_codon]['wij']))
            cai2_numerator = cai2_numerator + int(seq_dataframe[seq_dataframe.codon == fij_codon]['Obsi']) * float(
                wij_dataframe[wij_dataframe.codon == fij_codon]['wij'])
            cai_denominator = cai_denominator + int(seq_dataframe[seq_dataframe.codon == fij_codon]['Obsi'])
        cai = np.exp(cai_numerator / cai_denominator)
        cai2 = cai2_numerator / cai_denominator
        return cai, cai2

    @staticmethod
    def csc(seq_info, mrna_half_life_info):
        data_total = seq_info.generate_dataframe_total()

        # Generate a dictionary of sequence information and mrna half-life.
        csc_index = 0
        csc_correspondence = {}
        for info in list(seq_info.correspondence.keys())[0:5]:
            csc_correspondence[info] = list(mrna_half_life_info['Total Half-life'])[csc_index]
            csc_index = csc_index + 1

        # Initialize and allocate memory space.
        for csc_codon in data_total['codon']:
            locals()['csc_' + str(csc_codon)] = []
            locals()['mrna_hf_' + str(csc_codon)] = []

        # Store the frequency of codons in the list.
        for info in list(seq_info.correspondence.keys())[0:5]:
            for csc_codon in seq_info.generate_dataframe(info)['codon']:
                locals()['csc_' + str(csc_codon)].append(
                    int(seq_info.generate_dataframe(info).loc[
                            seq_info.generate_dataframe(info).codon == csc_codon]['Obsi']))
                locals()['mrna_hf_' + str(csc_codon)].append(csc_correspondence[info])

        csc_dataframe = pd.DataFrame(columns=('codon', 'csc'))
        for csc_codon in data_total['codon']:
            locals()['codon_Obsi_list_' + str(csc_codon)] = np.array(locals()['csc_' + str(csc_codon)])
            locals()['mrna_hf_list_' + str(csc_codon)] = np.array(locals()['mrna_hf_' + str(csc_codon)])
            if len(locals()['mrna_hf_list_' + str(csc_codon)]) > 1 \
                    and len(locals()['codon_Obsi_list_' + str(csc_codon)]) > 1:
                csc = np.corrcoef(locals()['codon_Obsi_list_' + str(csc_codon)],
                                  locals()['mrna_hf_list_' + str(csc_codon)])[0, 1]
                csc_dataframe = csc_dataframe.append(pd.DataFrame({'codon': [csc_codon],
                                                                   'csc': [csc]}), ignore_index=True)
        return csc_dataframe


brewer_yeast = Cubat('Test_Data/CSC_test/csc_test_total5.fna')
mrna_hl = pd.read_excel('Test_Data/CSC_test/mrna_life_time.xlsx')  # Half-life of mrna
print(brewer_yeast.csc(brewer_yeast, mrna_hl))
