import numpy as np

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
    def _compute_genecode(amino_acid_sequence,
                          TTT=None,
                          TTC=None, TTA=None, TTG=None, TCT=None, TCC=None, TCA=None, TCG=None, TAT=None, TAC=None,
                          TAA=None,
                          TAG=None, TGT=None, TGC=None, TGA=None, TGG=None, CTT=None, CTC=None, CTA=None, CTG=None,
                          CCT=None,
                          CCC=None, CCA=None, CCG=None, CAT=None, CAC=None, CAA=None, CAG=None, CGT=None, CGC=None,
                          CGA=None,
                          CGG=None, ATT=None, ATC=None, ATA=None, ATG=None, ACT=None, ACC=None, ACA=None, ACG=None,
                          AAT=None,
                          AAC=None, AAA=None, AAG=None, AGT=None, AGC=None, AGA=None, AGG=None, GTT=None, GTC=None,
                          GTA=None,
                          GTG=None, GCT=None, GCC=None, GCA=None, GCG=None, GAT=None, GAC=None, GAA=None, GAG=None,
                          GGT=None,
                          GGC=None, GGA=None, GGG=None):
        amino_acid_arr = np.array(list(amino_acid_sequence))
        order = np.array(
            [TTT, TTC, TTA, TTG, TCT, TCC, TCA, TCG, TAT, TAC, TAA, TAG, TGT, TGC, TGA, TGG, CTT, CTC, CTA, CTG, CCT,
             CCC,
             CCA, CCG, CAT, CAC, CAA, CAG, CGT, CGC, CGA, CGG, ATT, ATC, ATA, ATG, ACT, ACC, ACA, ACG, AAT, AAC, AAA,
             AAG,
             AGT, AGC, AGA, AGG, GTT, GTC, GTA, GTG, GCT, GCC, GCA, GCG, GAT, GAC, GAA, GAG, GGT, GGC, GGA, GGG])
        changed_location = np.nonzero(order)
        if len(changed_location) > 0:
            amino_acid_arr[changed_location] = order[order != None]

        return amino_acid_arr

    @staticmethod
    def _get_genecode(ami_arr):
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

       
        return (stop_codon_location__,number_of_stopcodon__,single_codon_location__,number_of_single__,delete_location__,number_of_delete__
                ,codon_family24__,all_24families__,number_of_eachfamily__,perm__,codon_familyami__,all_amifamilies__)

    def __init__(self,amino_acid_order):
        amino_acid_array=genecode_compute._compute_genecode(amino_acid_order)
        data_tuple=genecode_compute._get_genecode(amino_acid_array)
        self.stop_codon_location=data_tuple[0]
        self.number_of_stopcodon=data_tuple[1]
        self.single_codon_location=data_tuple[2]
        self.number_of_single=data_tuple[3]
        self.delete_location=data_tuple[4]
        self.number_of_delete=data_tuple[5]
        self.codon_family24=data_tuple[6]
        self.all_24families=data_tuple[7]
        self.number_of_eachfamily=data_tuple[8]
        self.sorting=data_tuple[9]
        self.codon_familyami=data_tuple[10]
        self.all_amifamilies=data_tuple[11]

        return
