import numpy as np


stop_codon_location={}
number_of_stopcodon= {}
single_codon_location= {}
number_of_single= {}
delete_location= {}
number_of_delete={}
codon_family24={}
all_24familys={}
number_of_eachfamily={}
sorting={}
codon_familyami= {}
all_amifamilys={}

def _compute_genecode(amino_acid_sequence='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', TTT=None,
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
        [TTT, TTC, TTA, TTG, TCT, TCC, TCA, TCG, TAT, TAC, TAA, TAG, TGT, TGC, TGA, TGG, CTT, CTC, CTA, CTG, CCT, CCC,
         CCA, CCG, CAT, CAC, CAA, CAG, CGT, CGC, CGA, CGG, ATT, ATC, ATA, ATG, ACT, ACC, ACA, ACG, AAT, AAC, AAA, AAG,
         AGT, AGC, AGA, AGG, GTT, GTC, GTA, GTG, GCT, GCC, GCA, GCG, GAT, GAC, GAA, GAG, GGT, GGC, GGA, GGG])
    changed_location = np.nonzero(order)
    if len(changed_location) > 0:
        amino_acid_arr[changed_location] = order[order != None]

    return amino_acid_arr


def _get_genecode(ami_arr,indexes):

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
    all_24familys__ = np.diff(codon_family24__)
    family_type, number_of_eachfamily__ = np.unique(all_24familys__, return_counts=True)

    # get amino acid codon family
    # note that sorting and deleting are a must before computing
    ami_delete = np.delete(ami_arr, delete_location__)
    perm__ = np.argsort(ami_delete)
    ami_sort = ami_delete[perm__]
    mask_ami = np.empty(64 - number_of_delete__, dtype=np.bool_)
    mask_ami[:1] = True
    mask_ami[1:] = (ami_sort[1:] != ami_sort[:-1])
    codon_familyami__ = np.concatenate(np.nonzero(mask_ami) + ([mask_ami.size],))
    all_amifamilys__ = np.diff(codon_familyami__)

    # row = [stop_codon_location__, number_of_stopcodon__, single_codon_location__, number_of_single__,
    #        delete_location__, number_of_delete__,
    #        codon_family24__, number_of_eachfamily__, perm__, codon_familyami__]
    # dict={'single_codon_location':stop_codon_location__}
    # with open('genecode_database.csv', 'a') as gd:
    #     write = csv.writer(gd,lineterminator = '\n')
    #     write.writerow(row)

    stop_codon_location[indexes]=stop_codon_location__
    number_of_stopcodon[indexes]=number_of_stopcodon__
    single_codon_location[indexes]=single_codon_location__
    number_of_single[indexes]=number_of_single__
    delete_location[indexes]=delete_location__
    number_of_delete[indexes]=number_of_delete__
    codon_family24[indexes]=codon_family24__
    all_24familys[indexes]=all_24familys__
    number_of_eachfamily[indexes]=number_of_eachfamily__
    sorting[indexes]=perm__
    codon_familyami[indexes]=codon_familyami__
    all_amifamilys[indexes]=all_amifamilys__

    return


def new_genecode(id_or_name,amino_acid_sequence='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', TTT=None,
                 TTC=None, TTA=None, TTG=None, TCT=None, TCC=None, TCA=None, TCG=None, TAT=None, TAC=None, TAA=None,
                 TAG=None, TGT=None, TGC=None, TGA=None, TGG=None, CTT=None, CTC=None, CTA=None, CTG=None, CCT=None,
                 CCC=None, CCA=None, CCG=None, CAT=None, CAC=None, CAA=None, CAG=None, CGT=None, CGC=None, CGA=None,
                 CGG=None, ATT=None, ATC=None, ATA=None, ATG=None, ACT=None, ACC=None, ACA=None, ACG=None, AAT=None,
                 AAC=None, AAA=None, AAG=None, AGT=None, AGC=None, AGA=None, AGG=None, GTT=None, GTC=None, GTA=None,
                 GTG=None, GCT=None, GCC=None, GCA=None, GCG=None, GAT=None, GAC=None, GAA=None, GAG=None, GGT=None,
                 GGC=None, GGA=None, GGG=None):
    return _get_genecode(_compute_genecode(amino_acid_sequence, TTT,
                                           TTC, TTA, TTG, TCT, TCC, TCA, TCG, TAT, TAC, TAA,
                                           TAG, TGT, TGC, TGA, TGG, CTT, CTC, CTA, CTG, CCT,
                                           CCC, CCA, CCG, CAT, CAC, CAA, CAG, CGT, CGC, CGA,
                                           CGG, ATT, ATC, ATA, ATG, ACT, ACC, ACA, ACG, AAT,
                                           AAC, AAA, AAG, AGT, AGC, AGA, AGG, GTT, GTC, GTA,
                                          GTG, GCT, GCC, GCA, GCG, GAT, GAC, GAA, GAG, GGT,
                                           GGC, GGA, GGG),indexes=id_or_name)


new_genecode(1,'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(2,'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG')
new_genecode(3,'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(4,'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(5,'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG')
new_genecode(6,'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(9,'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG')
new_genecode(10,'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG')
new_genecode(11,'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(12,'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(13,'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG')
new_genecode(14,'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG')
new_genecode(16,'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(21,'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG')
new_genecode(22,'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(23,'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(24,'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG')
new_genecode(25,'FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(26,'FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(27,'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(28,'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(29,'FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(30,'FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(31,'FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
new_genecode(33,'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG')


