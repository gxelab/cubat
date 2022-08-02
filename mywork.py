import re
import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from select_parameters import select_parameters
from optimized_codon import optimized_codon
from genecode_data import genecode_compute
from Bio import SeqIO
from Bio.Seq import Seq


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


def info_preprocess(parse_data):
    data = list(parse_data)

    id_list = []
    sequence_length_list = []
    frequency_check = []
    start_codon_check = []
    seq_list = []

    for record in data:
        id_list.append(record.name)
        sequence_length_list.append(len(record.seq))
        start_codon_check.append(str(record.seq[0:3]))
        seq_list.append(str(record.seq))
        if len(record.seq) % 3 == 0:
            frequency_check.append("True")
        else:
            frequency_check.append("False")

    info_tuple = (id_list, seq_list, sequence_length_list, frequency_check, start_codon_check)

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

    # count the frequency of each codon in a sequence, which will be divided into groups of 8 for higher running velocity.

    depart = len(sequence) / 8
    times = int(depart)
    rest = int((depart - times) * 8)
    for i in range(times):
        codon_num[sequence.pop()] += 1
        codon_num[sequence.pop()] += 1
        codon_num[sequence.pop()] += 1
        codon_num[sequence.pop()] += 1
        codon_num[sequence.pop()] += 1
        codon_num[sequence.pop()] += 1
        codon_num[sequence.pop()] += 1
        codon_num[sequence.pop()] += 1
    for i in range(rest):
        codon_num[sequence.pop()] += 1

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


def rscu_compute(codon_dataframe, genecode):
    stop_codon = []
    rscu_dataframe = pd.DataFrame(columns=codon_dataframe.columns, index=codon_dataframe.index)
    for codon in codon_dataframe.columns:
        try:
            amino_acid_value = CodonTable.unambiguous_dna_by_id[genecode].forward_table[codon]
            rscu_codons = \
                [key for key, value in CodonTable.unambiguous_dna_by_id[genecode].forward_table.items() if
                 value == amino_acid_value]
            rscu_dataframe[codon] = codon_dataframe[codon] / (
                    codon_dataframe[rscu_codons].sum(axis=1) / len(rscu_codons))

        except KeyError:
            stop_codon.append(codon)
    return rscu_dataframe


def enc_compute(frequency_array, genecode_data):
    # call the genecode database

    familys = genecode_data.all_24families
    codon_location = genecode_data.codon_family24
    amount_of_eachfamily = genecode_data.number_of_eachfamily
    amount_of_single = genecode_data.number_of_single
    delete_location = genecode_data.delete_location

    # previous work for enc computing, including getting the amount of sequences, processing frequency array and create new array.

    number_of_sequence = frequency_array.shape[0]
    process_arr = np.delete(frequency_array, delete_location, axis=1)
    FcF_array, an_array = np.empty((21, number_of_sequence), dtype=float), np.empty((21, number_of_sequence),
                                                                                    dtype=float)

    # do the computing

    number_of_families = len(familys)
    # num1xa = process_arr
    for i in range(0, number_of_families):
        F_array = process_arr[:, codon_location[i]:codon_location[i + 1]]
        an = np.sum(F_array, axis=1)
        an_array[i] = an
        FcF_array[i] = np.sum(np.square(F_array + 1), axis=1) / np.square(an + familys[i])

    family2, family3, family4 = (familys == 2), (familys == 3), (familys == 4)
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
    # process optimized codons. remind that opt_codon_fop could be in two format, a 0-1 array judging whether a codon is an
    # optimized codon or not or an integer array ranged from 0~58 containing the location of optimized codon in the order of
    # amino acid from A to Y("GCA to TAT" less stop codons, W and M.).

    # call the genceode database
    codon_familyami = genecode_data.codon_familyami
    all_amifamilies = genecode_data.all_amifamilies
    delect_location = genecode_data.delete_location
    sorting=genecode_data.sorting

    # process the parameter
    opt_array_cbi = np.array(opt_codon_cbi)
    if any(opt_array_cbi > 1):
        processed_opt_cbi = opt_codon_cbi
    else:
        opt_array_cbi=np.delete(opt_array_cbi,delect_location)[sorting]
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
    # process optimized codons. remind that opt_codon_fop could be in two format, a 0-1 array judging whether a codon is an
    # optimized codon or not or an integer array ranged from 0~58 containing the location of optimized codon in the order of
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


def tai_compute(frequency_array, tRNA_GCN, SI_C=0.28, SI_A=0.9999, SG_U=0.41, SU_G=0.68):
    # note that in the tRNA_GCN, the value is still in the order of "TTT" to "GGG", but means a tRNA or an anti-codon.
    # for instance, "TTT" means tRNA or anti-codon is "TTT" instead of a codon.

    # compute the weight for each codon.

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
    sum = np.sum(mut, axis=1)
    tAI = np.exp(sum / np.sum(num, axis=1))
    tai_array = tAI

    return tai_array


def generate_wij_dataframe(codon_dataframe, genecode):
    wij_dict = {}
    for codon in codon_dataframe.columns:
        amino_acid_value = codon_table_completion(genecode)[codon]
        wij_codons = \
            [key for key, value in codon_table_completion(genecode).items() if
             value == amino_acid_value]
        wij = \
            codon_dataframe[wij_codons].sum()[codon] / codon_dataframe[wij_codons].sum().sort_values(ascending=False)[0]
        wij_dict.update({codon: wij})
    return pd.DataFrame(pd.Series(wij_dict))


def cai_compute(codon_dataframe, reference):
    reset_reference_dataframe = reference.loc[codon_dataframe.columns]
    reset_reference_dataframe.columns = ["wij"]

    codon_dataframe.loc["All sequences"] = codon_dataframe.apply(lambda x: x.sum())
    cai_up_dataframe = pd.DataFrame(columns=codon_dataframe.columns, index=codon_dataframe.index)
    cai2_up_dataframe = pd.DataFrame(columns=codon_dataframe.columns, index=codon_dataframe.index)

    for codon in codon_dataframe.columns:
        cai_up_dataframe[codon] = codon_dataframe[codon] * np.log(reset_reference_dataframe.loc[codon, "wij"])
        cai2_up_dataframe[codon] = codon_dataframe[codon] * reset_reference_dataframe.loc[codon, "wij"]

    cai_up_dataframe["Total"] = cai_up_dataframe.apply(lambda x: x.sum(), axis=1)
    cai2_up_dataframe["Total"] = cai2_up_dataframe.apply(lambda x: x.sum(), axis=1)
    codon_dataframe["Total frequency"] = codon_dataframe.apply(lambda x: x.sum(), axis=1)
    cai_value = cai_up_dataframe["Total"] / codon_dataframe["Total frequency"]
    cai2_value = cai2_up_dataframe["Total"] / codon_dataframe["Total frequency"]

    for key in cai_value.keys():
        cai_value[key] = np.exp(cai_value[key])

    return (cai_value, cai2_value)


def csc_compute(codon_dataframe, mrna_hl_location):
    mrna_hl = pd.read_excel(mrna_hl_location)  # Half-life of mrna
    array_mrna = mrna_hl[["Total Half-life"]].values.T
    csc_dataframe = pd.DataFrame(columns=codon_dataframe.columns, index=["csc"])
    for codon in csc_dataframe:
        csc_dataframe[codon] = float(
            np.corrcoef(array_mrna, codon_dataframe[codon].values.T)[[0], [1]])
    return csc_dataframe


def codon_bias_analyze(file_path, genecode=1, parameters=select_parameters(skip_select=True), file_format='fasta',
                       quality_control=True, rscu=True, enc=True, ite=True, X2=True,
                       fop=True, cbi=True, tai=True, cai=True, csc=False,):

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
    enc_flag = False

    # read genecode and get related parameters

    genecode_data = genecode_compute(genecode=genecode)
    delete_location = genecode_data.delete_location
    sorting = genecode_data.sorting

    # call other parameters for index computing
    try:
        # this is just a test, to judge whether what was given is class('select_parameter').
        test1 = parameters.cai_reference_frequency
    except:
        try:
            parameters = select_parameters(skip_select=True, parameter_file=parameters)
        except ValueError:
            print('the parameter file needs to be an excel file or you may input a dataframe.')

    # input file

    file = SeqIO.parse(file_path, file_format)
    information = info_preprocess(file)

    # count the frequency

    frequency_array = frequency_count(information[1])
    frequency_dataframe = pd.DataFrame(frequency_array,
                                       columns=codons,
                                       index=np.array(information[0]))

    print(frequency_dataframe)

    # compute indexes one by one

    index_dataframe = pd.DataFrame(index=np.array(information[0]))
    index_dataframe['sequnece_length'] = information[2]

    if quality_control:
        index_dataframe['mult_of_3'] = information[3]
        index_dataframe['start_codon'] = information[4]

    if enc:
        enc_array = enc_compute(frequency_array, genecode_data)
        index_dataframe['enc'] = enc_array
        enc_flag = True

    if ite:
        ite_array = ite_compute(frequency_array, parameters.ite_heg, parameters.ite_bg)
        index_dataframe['ite'] = ite_array

    if tai:
        tai_array = tai_compute(frequency_array, parameters.tai_tGCN)
        index_dataframe['tai'] = tai_array

    # process frequency array to compute indexes require codon:amino_acid congruent relationship.

    processed_array = np.delete(frequency_array, delete_location, axis=1)[:, sorting]

    if X2:
        X2_array = chi_squared(processed_array, genecode_data)
        index_dataframe['X2'] = X2_array

    if fop and cbi:

        # judge whether correlation analysis were used

        if parameters.cbi_opt_codon and parameters.fop_opt_codon:
            cbi_opt_codon = parameters.cbi_opt_codon
            fop_opt_codon = parameters.fop_opt_codon

        else:
            if enc_flag:
                fop_optimized, cbi_optimized = optimized_codon(processed_array,
                                                               enc_array=enc_array,
                                                               genecode_data=genecode_data)
            else:
                fop_optimized, cbi_optimized = optimized_codon(processed_array,
                                                               enc_array=enc_compute(frequency_array, genecode_data),
                                                               genecode_data=genecode_data)
            if parameters.cbi_opt_codon:
                fop_opt_codon = fop_optimized
            elif parameters.fop_opt_codon:
                cbi_opt_codon = cbi_optimized
            else:
                fop_opt_codon = fop_optimized
                cbi_opt_codon = cbi_optimized

        fop_array = fop_compute(processed_array, opt_codon_fop=fop_opt_codon)
        index_dataframe["fop"] = fop_array

        cbi_array = cbi_compute(processed_array, genecode_data=genecode_data, opt_codon_cbi=cbi_opt_codon)
        index_dataframe['cbi'] = cbi_array

    elif fop:
        if parameters.fop_opt_codon:
            fop_opt_codon = parameters.fop_opt_codon
        else:
            if enc_flag:
                fop_optimized, cbi_optimized = optimized_codon(processed_array,
                                                               enc_array=enc_array,
                                                               genecode_data=genecode_data)
            else:
                fop_optimized, cbi_optimized = optimized_codon(processed_array,
                                                           enc_array=enc_compute(frequency_array, genecode_data),
                                                           genecode_data=genecode_data)
            fop_opt_codon = fop_optimized
        fop_array = fop_compute(processed_array, opt_codon_fop=fop_opt_codon)
        index_dataframe["fop"] = fop_array

    elif cbi:
        if parameters.fop_opt_codon:
            cbi_opt_codon = parameters.cbi_opt_codon
        else:
            if enc_flag:
                fop_optimized, cbi_optimized = optimized_codon(processed_array,
                                                               enc_array=enc_array,
                                                               genecode_data=genecode_data)
            else:
                fop_optimized, cbi_optimized = optimized_codon(processed_array,
                                                               enc_array=enc_compute(frequency_array, genecode_data),
                                                               genecode_data=genecode_data)
            cbi_opt_codon = cbi_optimized

        cbi_array = cbi_compute(processed_array, genecode_data=genecode_data, opt_codon_cbi=cbi_opt_codon)
        index_dataframe['cbi'] = cbi_array

    if cai:
        reference_dataframe = pd.DataFrame(columns=codons)
        reference_dataframe.loc[len(reference_dataframe)] = np.array(parameters.cai_reference_frequency)
        reference_weight = generate_wij_dataframe(reference_dataframe, genecode)
        # reference_dataframe=pd.read_excel('Homo_sapiens.xlsx')
        cai_dataframe, cai2_dataframe = cai_compute(frequency_dataframe, reference_weight)
        index_dataframe['cai'] = cai_dataframe
        index_dataframe['cai2'] = cai2_dataframe

    print(index_dataframe)

    if rscu:
        rscu_dataframe = rscu_compute(frequency_dataframe, genecode)
        if cai:
            rscu_dataframe = rscu_dataframe.drop(columns='Total frequency')
        print(rscu_dataframe)

    if csc:
        csc_dataframe = csc_compute(frequency_dataframe, mrna_hl_location=None)
        print(csc_dataframe)

    return


# lll = codon_bias_analyze('C:/Users/YuanYe/Desktop/Homo_sapiens.GRCh38.cds2.FAS',genecode=1, quality_control=True, rscu=False,
#                          fop=True, cbi=False)
