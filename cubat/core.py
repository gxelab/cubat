import re
import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import os
import ntpath
import statsmodels.api as sm


# import matplotlib.pyplot as plt


class GenecodeCompute:
    stop_codon_location = None
    number_of_stopcodon = None
    single_codon_location = None
    number_of_single = None
    delete_location = None
    number_of_delete = None
    codon_family24 = None
    all_24families = None
    number_of_eachfamily = None
    sorting = None
    codon_familyami = None
    all_amifamilies = None

    @staticmethod
    def get_genecode(genecode):
        codon_table = CodonTable.unambiguous_dna_by_id[genecode].forward_table
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
        complete_codon_table = codon_table
        for codon in codons:
            try:
                codon_table[codon]
            except KeyError:
                complete_codon_table.update({codon: "*"})
        codon_series = pd.Series(complete_codon_table)[codons]
        complete_codon_table = dict(codon_series)
        amino_acid_arr = np.array(list(complete_codon_table.values()))
        return amino_acid_arr

    @staticmethod
    def compute_genecode(ami_arr):
        # get stop codons
        stop_codon_location__ = np.nonzero(ami_arr == "*")[0]
        number_of_stopcodon__ = len(stop_codon_location__)

        # get single aminoacid
        mask_forward = np.empty(64, dtype=np.bool_)
        mask_backward = np.empty(64, dtype=np.bool_)
        mask_forward[:1] = True
        mask_forward[1:] = (ami_arr[1:] != ami_arr[:-1])
        mask_backward[-1:] = True
        mask_backward[:-1] = ami_arr[:-1] != ami_arr[1:]
        locate_forward = np.nonzero(mask_forward)[0]
        locate_backward = np.nonzero(mask_backward)[0]
        single_codon_location__ = locate_forward[locate_forward == locate_backward]
        for stop_codon in stop_codon_location__:
            single_codon_location__ = single_codon_location__[single_codon_location__ != stop_codon]

        # get codons need to be deleted when computing(stop and single)
        delete_location__ = np.concatenate((single_codon_location__, stop_codon_location__))
        number_of_delete__ = len(delete_location__)
        number_of_single__ = len(single_codon_location__)

        # get 2,4 codon family
        # note that it's a must to remove single aminoacid and stop codon first
        mask_syn = np.delete(mask_forward, delete_location__)
        codon_family24__ = np.concatenate(np.nonzero(mask_syn) + ([mask_syn.size],))
        all_24families__ = np.diff(codon_family24__)
        family_type, number_of_eachfamily__ = np.unique(all_24families__, return_counts=True)

        # get amino acid codon family
        # note that sorting and deleting are a must before computing
        ami_delete = np.delete(ami_arr, delete_location__)
        perm__ = np.argsort(ami_delete)
        ami_sort = ami_delete[perm__]
        mask_ami = np.empty(64 - number_of_delete__, dtype=np.bool_)
        mask_ami[:1] = True
        mask_ami[1:] = (ami_sort[1:] != ami_sort[:-1])
        codon_familyami__ = np.concatenate(np.nonzero(mask_ami) + ([mask_ami.size],))
        all_amifamilies__ = np.diff(codon_familyami__)
        return (stop_codon_location__, number_of_stopcodon__, single_codon_location__, number_of_single__,
                delete_location__, number_of_delete__, codon_family24__, all_24families__, number_of_eachfamily__,
                perm__, codon_familyami__, all_amifamilies__)

    def __init__(self, genecode):
        amino_acid_array = GenecodeCompute.get_genecode(genecode)
        data_tuple = GenecodeCompute.compute_genecode(amino_acid_array)
        self.stop_codon_location = data_tuple[0]
        self.number_of_stopcodon = data_tuple[1]
        self.single_codon_location = data_tuple[2]
        self.number_of_single = data_tuple[3]
        self.delete_location = data_tuple[4]
        self.number_of_delete = data_tuple[5]
        self.codon_family24 = data_tuple[6]
        self.all_24families = data_tuple[7]
        self.number_of_eachfamily = data_tuple[8]
        self.sorting = data_tuple[9]
        self.codon_familyami = data_tuple[10]
        self.all_amifamilies = data_tuple[11]

        return


def optimized_codon(processed_array, enc_array, genecode_data):
    # call genecode database
    codon_familyami = genecode_data.codon_familyami
    all_amifamilies = genecode_data.all_amifamilies

    # do the correlation analysis
    fop_opt = np.zeros(18, dtype=int)
    zvalue_array = np.zeros(59, dtype=float)
    cir = 0
    for i in range(0, len(codon_familyami) - 1):
        one_ami = processed_array[:, codon_familyami[i]:codon_familyami[i + 1]]
        ami_count = np.sum(one_ami, axis=1)
        ami_count[ami_count == 0] = 1
        probability = one_ami.T / ami_count.T
        binomial_array = np.zeros(6, dtype=float)
        for m in range(0, all_amifamilies[i]):
            endog = probability[m]
            exog = enc_array
            exog = sm.add_constant(exog, prepend=True)
            res = sm.GLM(endog=endog, var_weights=ami_count, exog=exog, missing=True,
                         family=sm.families.Binomial()).fit()
            binomial_array[m] = res.tvalues[1]
            zvalue_array[cir] = res.tvalues[1]
            cir += 1
        locate = np.argmin(binomial_array)
        locate += codon_familyami[i]
        fop_opt[i] = locate
    cbi_opt = np.nonzero(zvalue_array < (-np.sqrt(processed_array.shape[0]) / 3))[0]
    return_tuple = (fop_opt, cbi_opt)
    return return_tuple


def codon_table_completion(genecode):
    codon_table = CodonTable.unambiguous_dna_by_id[genecode].forward_table
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
    complete_codon_table = codon_table
    for codon in codons:
        try:
            codon_table[codon]
        except KeyError:
            complete_codon_table.update({codon: "*"})
    codon_series = pd.Series(complete_codon_table)[codons]
    complete_codon_table = dict(codon_series)

    return complete_codon_table


def info_preprocess(input_path, file_format):
    if file_format not in ['csv', 'save']:
        id_list = []
        sequence_length_list = []
        frequency_check = []
        start_codon_check = []
        seq_list = []
        description = []
        parse_data = SeqIO.parse(input_path, file_format)
        data = list(parse_data)

        # todo 这可能是主要的限速步骤，有机会看看能不能不用for loop。
        for record in data:
            id_list.append(record.name)
            sequence_length_list.append(len(record.seq))
            start_codon_check.append(str(record.seq[0:3]))
            seq_list.append(str(record.seq))
            description.append(record.description)

            if len(record.seq) % 3 == 0:
                frequency_check.append("True")
            else:
                frequency_check.append("False")
    else:
        saved_dataframe = pd.read_csv(input_path, index_col=0)
        seq_list = saved_dataframe[saved_dataframe.columns[4:68]].to_numpy()
        id_list = list(saved_dataframe.index)
        description = saved_dataframe[saved_dataframe.columns[3]].tolist()
        sequence_length_list = saved_dataframe[saved_dataframe.columns[2]].tolist()
        frequency_check = saved_dataframe[saved_dataframe.columns[1]].tolist()
        start_codon_check = saved_dataframe[saved_dataframe.columns[0]].tolist()

    info_tuple = (id_list, seq_list, description, sequence_length_list, frequency_check, start_codon_check)

    return info_tuple


def count_codon(sequence):
    # the sequence should be a list.
    # create a dictionary for counting
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

    # count the frequency of each codon in a sequence,
    # which will be divided into groups of 8 for higher running velocity.
    counter = Counter(sequence)
    codon_num.update(counter)

    codon_array = np.array(list(codon_num.values()))
    return codon_array


def frequency_count(seq_data):
    if type(seq_data) == tuple:
        seq_data = list(seq_data)
    times = len(seq_data)
    result_array = np.zeros((times, 64), dtype=int)
    for t in range(0, times):
        seq_str = seq_data[t]
        seq_list = re.findall('.{3}', seq_str)
        result_array[t] = count_codon(seq_list)
    return result_array


def rscu_compute(processed_array, genecode_data):
    # call genecode_data
    codons = ['TTT', 'TTC', 'TTA', 'TTG',
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
              'GGT', 'GGC', 'GGA', 'GGG']
    stop_location = genecode_data.stop_codon_location
    single_location = genecode_data.single_codon_location
    delete_location = genecode_data.delete_location
    codon_familyami = genecode_data.codon_familyami
    all_amifamilies = genecode_data.all_amifamilies
    sorting = genecode_data.sorting

    # process codon for generating dataframe
    # processed_array = np.delete(frequency_array, delete_location, axis=1)[:, sorting]
    processed_codon = np.delete(np.array(codons), delete_location)[sorting]
    single_codon = np.array(codons)[single_location]
    stop_codon = np.array(codons)[stop_location]

    # do the computing
    result_array = np.zeros(processed_array.shape)
    length = processed_array.shape[0]
    for i in range(0, len(codon_familyami) - 1):
        one_ami = processed_array[:, codon_familyami[i]:codon_familyami[i + 1]]
        ami_count = np.sum(one_ami, axis=1)
        ami_count[ami_count == 0] = 1
        rscuij = all_amifamilies[i] * one_ami.T / ami_count
        result_array[:, codon_familyami[i]:codon_familyami[i + 1]] = rscuij.T
    result_dataframe = pd.DataFrame(result_array, columns=processed_codon)
    for single in single_codon:
        result_dataframe[single] = np.ones(length)
    nan_array = np.zeros(length)
    nan_array[:] = np.nan
    for stop in stop_codon:
        result_dataframe[stop] = nan_array

    # sorting the dataframe
    rscu_dataframe = result_dataframe.loc[:, codons]

    return rscu_dataframe


def enc_compute(frequency_array, genecode_data):
    # call the genecode database
    families = genecode_data.all_24families
    codon_location = genecode_data.codon_family24
    amount_of_eachfamily = genecode_data.number_of_eachfamily
    amount_of_single = genecode_data.number_of_single
    delete_location = genecode_data.delete_location

    # previous work for enc computing, including getting the amount of sequences,
    # processing frequency array and create new array.
    number_of_sequence = frequency_array.shape[0]
    process_arr = np.delete(frequency_array, delete_location, axis=1)
    FcF_array, an_array = np.empty((21, number_of_sequence), dtype=float), np.empty((21, number_of_sequence),
                                                                                    dtype=float)

    # do the computing
    number_of_families = len(families)
    for i in range(0, number_of_families):
        F_array = process_arr[:, codon_location[i]:codon_location[i + 1]]
        an = np.sum(F_array, axis=1)
        an_array[i] = an
        FcF_array[i] = np.sum(np.square(F_array + 1), axis=1) / np.square(an + families[i])

    family2, family3, family4 = (families == 2), (families == 3), (families == 4)
    family2sum = amount_of_eachfamily[0] * np.sum(an_array[family2], axis=0) / np.sum(
        np.multiply(FcF_array[family2], an_array[family2]), axis=0)
    family3_prevent0 = np.sum(np.multiply(FcF_array[family3], an_array[family3]), axis=0)
    family3_prevent0[family3_prevent0 == 0] = 1
    family3sum = amount_of_eachfamily[1] * np.sum(an_array[family3], axis=0) / family3_prevent0
    family4sum = amount_of_eachfamily[2] * np.sum(an_array[family4], axis=0) / np.sum(
        np.multiply(FcF_array[family4], an_array[family4]), axis=0)
    enc_array = amount_of_single + family2sum + family3sum + family4sum
    return enc_array


def ite_compute(frequency_array, ref_heg, ref_bg):
    # ref_heg means frequency of high expressed genes in the reference. ref_bg means that of other genes.
    # compute weight for each codon, which will be put in an array.

    Wte = np.array(np.zeros((29, 2), dtype=float))
    for i in range(0, 58, 2):
        Wte[int(i / 2), 0] = ref_heg[i] / ref_bg[i]
        Wte[int(i / 2), 1] = ref_heg[i + 1] / ref_bg[i + 1]
        Wte[int(i / 2)] *= 1 / max([ref_heg[i] / ref_bg[i], ref_heg[i + 1] / ref_bg[i + 1]])
    Wte = Wte.flatten()

    # compute ite

    num = np.delete(frequency_array, [10, 11, 14, 15, 34, 35], axis=1)
    mut = np.multiply(np.log(Wte), num)
    row_sum = np.sum(mut, axis=1)
    ITE = np.exp(row_sum / np.sum(num, axis=1))
    ite_array = ITE

    return ite_array


def chi_squared(processed_array, genecode_data):
    codon_familyami = genecode_data.codon_familyami
    all_amifamilies = genecode_data.all_amifamilies

    times = processed_array.shape[0]
    all_ami = np.zeros((18, times), dtype=float)
    for i in range(0, len(all_amifamilies)):
        codon_count = processed_array[:, codon_familyami[i]:codon_familyami[i + 1]]
        ami_count = np.sum(codon_count, axis=1)
        ami_count[ami_count == 0] = 1
        every_codon = np.square((codon_count.T / ami_count) - (1 / all_amifamilies[i]))
        all_codon = every_codon * ami_count * all_amifamilies[i]
        all_ami[i] = np.sum(all_codon, axis=0)
    sum_arr = np.sum(all_ami, axis=0)
    chi2_array = sum_arr / np.sum(processed_array, axis=1)

    return chi2_array


def cbi_compute(processed_array, genecode_data, opt_codon_cbi):
    # process optimized codons. remind that opt_codon_fop could be in two format,
    # a 0-1 array judging whether a codon is an optimized codon or not or an integer array ranged
    # from 0~58 containing the location of optimized codon in the order of
    # amino acid from A to Y("GCA to TAT" less stop codons, W and M.).

    # call the genecode database
    codon_familyami = genecode_data.codon_familyami
    all_amifamilies = genecode_data.all_amifamilies
    delete_location = genecode_data.delete_location
    sorting = genecode_data.sorting

    # process the parameter
    opt_array_cbi = np.array(opt_codon_cbi)
    if any(opt_array_cbi > 1):
        processed_opt_cbi = opt_codon_cbi
    else:
        opt_array_cbi = np.delete(opt_array_cbi, delete_location)[sorting]
        processed_opt_cbi = np.nonzero(opt_array_cbi)[0]

    # do the computing
    all_codon = np.sum(processed_array, axis=1)
    opt_codon = np.sum(processed_array[:, processed_opt_cbi], axis=1)
    rand_number = np.zeros(processed_array.shape[0], dtype=float)
    cbi_list = list(processed_opt_cbi)
    for i in range(0, len(codon_familyami) - 1):
        opt_ami = 0
        one_ami = processed_array[:, codon_familyami[i]:codon_familyami[i + 1]]
        ami_count = np.sum(one_ami, axis=1)
        if len(cbi_list) != 0:
            while cbi_list[0] < codon_familyami[i + 1]:
                cbi_list.pop(0)
                opt_ami += 1
                if len(cbi_list) == 0:
                    break

        weight = opt_ami / all_amifamilies[i]
        rand_number += ami_count * weight

    cbi_array = (opt_codon - rand_number) / (all_codon - rand_number)

    return cbi_array


def fop_compute(processed_array, opt_codon_fop):
    # process optimized codons. remind that opt_codon_fop could be in two format,
    # a 0-1 array judging whether a codon is an optimized codon or not or an integer array ranged from
    # 0~58 containing the location of optimized codon in the order of
    # amino acid from A to Y("GCA to TAT" less stop codons, W and M.).
    opt_array_fop = np.array(opt_codon_fop)
    if any(opt_array_fop > 1):
        processed_opt_fop = opt_codon_fop
    else:
        processed_opt_fop = np.nonzero(opt_array_fop)[0]

    # do the computing
    opt_frequency_array = processed_array[:, processed_opt_fop]
    fop_array = np.sum(opt_frequency_array, axis=1) / np.sum(processed_array, axis=1)

    return fop_array


def tai_compute(frequency_array, tRNA_GCN, tai_s=None):
    # note that in the tRNA_GCN, the value is still in the order of "TTT" to "GGG", but means a tRNA or an anti-codon.
    # for instance, "TTT" means tRNA or anti-codon is "TTT" instead of a codon.

    # compute the weight for each codon.

    if type(tai_s) == pd.DataFrame:
        # todo 注意如果给矩阵
        SI_C = tai_s.loc['SI_C']
        SI_A = tai_s.loc['SI_A']
        SG_U = tai_s.loc['SG_U']
        SU_G = tai_s.loc['SU_G']
    else:
        SI_C = 0.28
        SI_A = 0.9999
        SG_U = 0.41
        SU_G = 0.68

    tGCN_by_codon = []
    tRNA = ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TAG', 'TAA', 'TGT', 'TGC', 'TGA',
            'TGG', 'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC',
            'CGA', 'CGG', 'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT',
            'AGC', 'AGA', 'AGG', 'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG',
            'GGT', 'GGC', 'GGA', 'GGG']
    for i in tRNA:
        tGCN_by_codon.append(tRNA_GCN[tRNA.index(Seq(i).reverse_complement())])
    Ws = np.array(np.zeros((1, 64), dtype=float))
    for i in range(0, 61, 4):
        Ws[0, i] = 1 * tGCN_by_codon[i] + (1 - SG_U) * tGCN_by_codon[i + 1]
        Ws[0, i + 1] = 1 * tGCN_by_codon[i + 1] + (1 - SI_C) * tGCN_by_codon[i]
        Ws[0, i + 2] = 1 * tGCN_by_codon[i + 2] + (1 - SI_A) * tGCN_by_codon[i]
        Ws[0, i + 3] = 1 * tGCN_by_codon[i + 3] + (1 - SU_G) * tGCN_by_codon[i + 2]
    Ws = np.delete(Ws, [10, 11, 14, 35], axis=1)
    Ws[Ws == 0] = np.mean(Ws != 0)
    Ws = Ws / Ws.max()

    # do the computing

    num = np.delete(frequency_array, [10, 11, 14, 35], axis=1)
    mut = np.multiply(np.log(Ws), num)
    tai_sum = np.sum(mut, axis=1)
    tAI = np.exp(tai_sum / np.sum(num, axis=1))
    tai_array = tAI

    return tai_array


def generate_cai_ref(processed_array):
    cai_ref = 0
    for i in processed_array:
        cai_ref = cai_ref + i
    return cai_ref


def cai_compute(processed_array, reference, genecode_data):
    delete_location = genecode_data.delete_location
    sorting = genecode_data.sorting
    codon_familyami = genecode_data.codon_familyami

    ref = np.delete(reference, delete_location)[sorting]
    W_array = np.zeros(ref.shape)
    for i in range(0, len(codon_familyami) - 1):
        one_ami = ref[codon_familyami[i]:codon_familyami[i + 1]]
        W = one_ami / one_ami.max()
        W_array[codon_familyami[i]:codon_familyami[i + 1]] = W
    all_codon = np.sum(processed_array, axis=1)
    multi_array = np.sum(np.multiply(processed_array, np.log(W_array)), axis=1)
    cai_array = np.exp(multi_array / all_codon)
    cai2_array = np.sum(np.multiply(processed_array, W_array), axis=1) / all_codon
    return_tuple = (cai_array, cai2_array)
    return return_tuple


def csc_compute(codon_dataframe, mrna_hl):
    # Half-life of mrna
    csc_dataframe = pd.DataFrame(columns=codon_dataframe.columns, index=["csc"])
    for codon in csc_dataframe:
        csc_dataframe[codon] = float(
            np.corrcoef(mrna_hl, codon_dataframe[codon].values.T)[[0], [1]])
    return csc_dataframe


class Analyze:
    enc_result = None
    ite_result = None
    X2_result = None
    fop_result = None
    cbi_result = None
    tai_result = None
    cai_result = None
    cai2_result = None
    rscu_result = None
    csc_result = None
    inputpath = ''
    output = ''
    prefix = ''
    codons = ['TTT', 'TTC', 'TTA', 'TTG',
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
              'GGT', 'GGC', 'GGA', 'GGG']

    def __init__(self, input_path, genecode=1, file_format='fasta', save=False,
                 quality_control=True, enc=False, ite=False, ite_ref='', X2=False,
                 fop=False, fop_opt='', cbi=False, cbi_opt='', tai=False, tai_gcn='', tai_s='', cai=False, cai_ref='',
                 rscu=False,
                 csc=False, csc_ref='', output=None, prefix=''):

        # read genecode and get related parameters
        genecode_data = GenecodeCompute(genecode=genecode)
        delete_location = genecode_data.delete_location
        sorting = genecode_data.sorting
        # enc_flag = False

        print('processing input file...')
        # use Seq.IO to read the file and do info preprocess.
        information = info_preprocess(input_path, file_format)
        Analyze.inputpath = input_path

        # count the frequency
        print('counting codon frequency...')
        # save_dataframe = pd.read_csv('Test_Data\Sars_cov_2.ASM985889v3.cds.fasta.save.csv', index_col=0)
        # frequency_array = save_dataframe[save_dataframe.columns[2:66]].to_numpy()
        # print(frequency_array[0], frequency_array1[0])
        # print(type(frequency_array[0]), type(frequency_array1[0]))

        if file_format not in ['csv', 'save']:
            frequency_array = frequency_count(information[1])
        else:
            frequency_array = information[1]
        frequency_dataframe = pd.DataFrame(frequency_array,
                                           columns=Analyze.codons,
                                           index=np.array(information[0]))
        # inner codon control
        inner_stop_codons = np.sum(frequency_array[:, genecode_data.stop_codon_location], axis=1)
        inner_stop_codons = (inner_stop_codons > 1)
        # compute indexes one by one
        print('computing indexes for gene...')
        index_dataframe = pd.DataFrame(index=np.array(information[0]))
        index_dataframe['sequences_length'] = information[3]
        index_dataframe['description'] = information[2]

        if quality_control:
            index_dataframe['multi_of_3'] = information[4]
            index_dataframe['start_codon'] = information[5]
            index_dataframe['inner_stopcodon'] = inner_stop_codons

        if prefix:
            prefix_name = prefix
        else:
            prefix_name = ntpath.basename(input_path)

        if output:
            out = os.path.join(output, prefix_name)
            self.output = output
        else:
            padir = os.path.dirname(input_path)
            out = os.path.join(padir, prefix_name)
            self.output = padir

        if save and file_format not in ['csv', 'save']:
            print("saving...")
            save_dataframe = frequency_dataframe
            save_dataframe.insert(loc=0, column='description', value=information[2])
            save_dataframe.insert(loc=0, column='sequence_length', value=information[3])
            save_dataframe.insert(loc=0, column='frequency_check', value=information[4])
            save_dataframe.insert(loc=0, column='start_codon_check', value=information[5])
            save_dataframe.insert(loc=len(save_dataframe.columns), column='seq_list', value=information[1])
            save_dataframe.to_csv(out + ".save.csv", sep=",", index=True)

        if enc:
            enc_array = enc_compute(frequency_array, genecode_data)
            index_dataframe['enc'] = enc_array
            self.enc_result = enc_array

        if ite:
            try:
                ite_ref = pd.read_csv(ite_ref, index_col=0, header=0).loc[Analyze.codons]
            except:
                raise FileNotFoundError('the reference for ite is unreachable.')
            # todo 自己寻找reference。

            ite_array = ite_compute(frequency_array, np.array(ite_ref.iloc[:, 0].values),
                                    np.array(ite_ref.iloc[:, 1].values))
            # todo 注意报错
            index_dataframe['ite'] = ite_array
            self.ite_result = ite_array

        if tai:
            try:
                tai_gcn = pd.read_csv(tai_gcn, index_col=0, header=0).loc[Analyze.codons]
            except:
                # todo 注意s的修改
                raise FileNotFoundError('the gene copy number for tai is unreachable.')
            try:
                tai_s = pd.read_csv(tai_s, index_col=0, header=0)
            except:
                raise FileNotFoundError('the selective constraint for tai is unreachable.')

            tai_array = tai_compute(frequency_array, np.array(tai_gcn.iloc[:, 0].values), tai_s)
            index_dataframe['tai'] = tai_array
            self.tai_result = tai_array

        # process frequency array to compute indexes require codon:amino_acid congruent relationship.
        processed_array = np.delete(frequency_array, delete_location, axis=1)[:, sorting]

        if X2:
            X2_array = chi_squared(processed_array, genecode_data)
            index_dataframe['X2'] = X2_array
            self.X2_result = X2_array

        if fop:
            try:
                fop_opt = pd.read_csv(fop_opt, index_col=0, header=0).loc[Analyze.codons]
            except:
                raise FileNotFoundError('the reference for fop is unreachable.')
            # todo 自己寻找reference。

            fop_array = fop_compute(frequency_array, np.array(fop_opt.iloc[:, 0].values))
            index_dataframe['fop'] = fop_array
            self.fop_result = fop_array

        if cbi:
            try:
                cbi_opt = pd.read_csv(cbi_opt, index_col=0, header=0).loc[Analyze.codons]
            except:
                raise FileNotFoundError('the reference for cbi is unreachable.')
            # todo 自己寻找reference。注意排序的问题。两次排序，先排成ACGT，再按照氨基酸排。

            cbi_array = cbi_compute(frequency_array, opt_codon_cbi=np.array(cbi_opt.iloc[:, 0].values),
                                    genecode_data=genecode_data)
            index_dataframe['cbi'] = cbi_array
            self.cbi_result = cbi_array

        if cai:
            try:
                if ".csv" in cai_ref:
                    cai_ref = pd.read_csv(cai_ref, index_col=0, header=0).loc[Analyze.codons]
                # else:
                #     cai_ref = generate_cai_ref()
            except:
                raise FileNotFoundError('the reference for cai is unreachable.')

            cai_array, cai2_array = cai_compute(processed_array, genecode_data=genecode_data,
                                                reference=np.array(cai_ref.iloc[:, 0].values))
            # reference_dataframe=pd.read_excel('Homo_sapiens.xlsx')
            index_dataframe['cai'] = cai_array
            index_dataframe['cai2'] = cai2_array

            self.cai_result = cai_array
            self.cai2_result = cai2_array

        # print(index_dataframe)

        print('counting indexes for codon...')
        if rscu:
            rscu_dataframe = rscu_compute(processed_array, genecode_data)
            rscu_dataframe.index = frequency_dataframe.index
            self.rscu_result = rscu_dataframe
            rscu_dataframe.to_csv(out + ".rscu.csv", sep=",", index=True)

        if csc:
            try:
                csc_ref = pd.read_csv(csc_ref, index_col=0, header=0)
            except:
                raise FileNotFoundError('the reference for csc is unreachable.')
            # todo 注意文件格式
            csc_dataframe = csc_compute(frequency_dataframe, mrna_hl=np.array(csc_ref.iloc[:, 0].values))
            self.csc_result = csc_dataframe
            csc_dataframe.to_csv(out + ".csc.csv", sep=",", index=True)

        # output file
        print('writing output files...')

        # when the out path existing already, a figure will be accumulated to the end of outpath,
        # which aims to avoid repetition for names if a single prefix correspond to a directory.
        while os.path.exists(out + ".index.csv"):
            try:
                serial = int(out[-1]) + 1
                out = out[:-1] + str(serial)
            except:
                out = out + str(1)
        self.prefix = ntpath.basename(out)

        index_dataframe.to_csv(out + ".index.csv", sep=",", index=True)
        # restore original working path
        print('done')
        return

    # def plot(self, rscu_barplot=False, output='', prefix=''):
    #     print('start plotting...')
    #     # decide output path and prefix
    #     if prefix:
    #         try:
    #             n = int(self.prefix[-1])
    #             prefix_name = prefix + str(n)
    #         except:
    #             prefix_name = prefix
    #     else:
    #         prefix_name = self.prefix
    #
    #     if output:
    #         out = os.path.join(output, prefix_name)
    #     else:
    #         out = os.path.join(self.output, prefix_name)

    # do the plotting
    # if rscu_barplot:
    #     if type(self.rscu_result) != pd.DataFrame:
    #         raise ValueError('rscu was not computed.')
    #     plt.bar(range(len(analyze.codons)), self.rscu_result.iloc[0], tick_label=analyze.codons)
    #     plt.title('''rscu_barplot''', fontsize=20)
    #     plt.savefig(out + 'rscu.jpg', dpi=500, bbox_inches='tight')
    #
    # print('done')
    # return


# if FAS type fasta.

if __name__ == '__main__':
    test = Analyze('Test_Data/Sars_cov_2.ASM985889v3.cds.fasta', file_format='fasta',
                   genecode=1, enc=True, ite=True,
                   ite_ref='example/ite_ref.csv', tai=True, tai_gcn='example/tai_gcn.csv', tai_s='example/tai_s.csv',
                   cai=True, save=True,
                   cai_ref='example/cai_ref.csv', cbi=True, cbi_opt='example/cbi_opt.csv',
                   fop=True, fop_opt='example/fop_opt.csv', X2=True,
                   rscu=True, output='Test_Data')
    # lll.plot(rscu_barplot=True)

# .save.csv', file_format='csv',
