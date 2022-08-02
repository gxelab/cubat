import numpy as np
from Bio.Data import CodonTable
import pandas as pd


class genecode_compute():
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

        # get codons need to be deleted when computing(stop and single)
        get_delete_location = np.concatenate((single_codon_location__, stop_codon_location__))
        delete_location__ = np.unique(get_delete_location)
        number_of_delete__ = len(delete_location__)
        number_of_single__ = number_of_delete__ - number_of_stopcodon__

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

        return (
        stop_codon_location__, number_of_stopcodon__, single_codon_location__, number_of_single__, delete_location__,
        number_of_delete__
        , codon_family24__, all_24families__, number_of_eachfamily__, perm__, codon_familyami__, all_amifamilies__)

    def __init__(self, genecode):
        amino_acid_array = genecode_compute.get_genecode(genecode)
        data_tuple = genecode_compute.compute_genecode(amino_acid_array)
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
