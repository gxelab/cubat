import pandas as pd


reference_frequnecy_cai = {
    'Homo_sapiens_EMBOSS': [544814, 655524, 234253, 406803, 477158, 567853, 383463, 143585, 385678, 495846, 20594,
                            17361, 327792, 403140, 36263, 417781, 413535, 632996, 224870, 1295575, 558133, 645805,
                            538723, 226405, 340239, 483801, 382585, 1109029, 148414, 345088, 200253, 376554, 502403,
                            677714, 229688, 710465, 414756, 615448, 476341, 197169, 532711, 615538, 764918, 1035170,
                            384496, 626669, 373118, 375881, 349016, 470297, 223912, 919492, 595935, 912370, 510981,
                            242826, 702386, 820851, 922616, 1291190, 345134, 727678, 524128, 531888]}
tGCN_tai = {
    'Homo_sapiens_GRCh38_GtRNAdb': [12, 8, 0, 6, 6, 9, 1, 6, 5, 5, 3, 4, 6, 8, 4, 7, 15, 8, 0, 13, 5, 5, 7, 4, 20, 13,
                                    7, 9, 5, 4, 4, 4, 0, 0, 1,
                                    0, 0, 0, 0, 7, 15, 10, 0, 9, 9, 26, 9, 9, 25, 13, 13, 9, 8, 14, 29, 0, 3, 0, 10, 0,
                                    0, 0, 0, 0]}
reference_heg_ite = {
    'Homo_sapiens_EMBOSS': [544814.0, 655524.0, 234253.0, 406803.0, 477158.0, 567853.0, 383463.0, 143585.0, 385678.0,
                            495846.0,
                            327792.0, 403140.0, 413535.0, 632996.0, 224870.0, 1295575.0, 558133.0, 645805.0, 538723.0,
                            226405.0,
                            340239.0, 483801.0, 382585.0, 1109029.0, 148414.0, 345088.0, 200253.0, 376554.0, 502403.0,
                            677714.0,
                            414756.0, 615448.0, 476341.0, 197169.0, 532711.0, 615538.0, 764918.0, 1035170.0, 384496.0,
                            626669.0,
                            373118.0, 375881.0, 349016.0, 470297.0, 223912.0, 919492.0, 595935.0, 912370.0, 510981.0,
                            242826.0,
                            702386.0, 820851.0, 922616.0, 1291190.0, 345134.0, 727678.0, 524128.0, 531888.0, 17361.0,
                            20594.0,
                            36263.0]}
reference_bg_ite = {
    'Homo_sapiens_EMBOSS': [687589.0, 512749.0, 345215.0, 295841.0, 459081.0, 342347.0, 414992.0, 355639.0, 504963.0,
                            376561.0,
                            418699.0, 312233.0, 749622.0, 559009.0, 677630.0, 580714.0, 575017.0, 428803.0, 519794.0,
                            445452.0,
                            472034.0, 352006.0, 803248.0, 688366.0, 312558.0, 233081.0, 282540.0, 242131.0, 532065.0,
                            396773.0,
                            497528.0, 371017.0, 449746.0, 385423.0, 657750.0, 490499.0, 969364.0, 830724.0, 579225.0,
                            431940.0,
                            403343.0, 345656.0, 573163.0, 427420.0, 518118.0, 444016.0, 660594.0, 492619.0, 597152.0,
                            511746.0,
                            872554.0, 650683.0, 1192155.0, 1021651.0, 621672.0, 463594.0, 561968.0, 481594.0, 23436.0,
                            27347.0,
                            23436.0]}
optimized_codon_fop = {'correlation_analysis': []}
optimized_codon_cbi = {'correlation_analysis': []}


class select_parameters():
    cai_reference_frequency = reference_frequnecy_cai['Homo_sapiens_EMBOSS']
    tai_tGCN = tGCN_tai['Homo_sapiens_GRCh38_GtRNAdb']
    ite_heg = reference_heg_ite['Homo_sapiens_EMBOSS']
    ite_bg = reference_bg_ite['Homo_sapiens_EMBOSS']
    fop_opt_codon = optimized_codon_fop['correlation_analysis']
    cbi_opt_codon = optimized_codon_cbi['correlation_analysis']

    @staticmethod
    def parameter_select_cai():
        # select reference set for cai computing
        cai_dict = {1: 'Homo_sapiens_EMBOSS'}
        while True:
            try:
                # select reference set
                cai_choice = int(input("select reference set for cai computing"
                                       "\nYour choice is: "
                                       "\n(1) 'Homo_sapiens_reference_set_fromEMBOSS[DEFAULT]'\n"))
            except ValueError:
                print("Only one reference set could be chosen at a time.")
                continue
            if cai_choice not in cai_dict.keys():
                raise ValueError('the number you choose is out of range')

            select_parameters.cai_reference_frequency = reference_frequnecy_cai[cai_dict[cai_choice]]
            break
        print(f"reference set for cai computing was selected\n")
        pass
        # result +=(cai_reference_frequency,)
        # result +=(cai_reference_order,)

    @staticmethod
    def parameter_select_fop():
        # select optimized codons for fop computing
        fop_dict = {0: 'correlation_analysis'}
        while True:
            try:
                fop_choice = int(input("select optimized codons for fop computing"
                                       "\nYour choice is: "
                                       "\n(0) Correlation_analysis_for_optimized_codon[DEFAULT]"
                                       "\n"))
            except ValueError:
                print("Only one reference could be chosen at a time.")
                continue
            select_parameters.fop_opt_codon = optimized_codon_cbi[fop_dict[fop_choice]]
            break
        print(f"optimized codon for fop computing was selected\n")

    @staticmethod
    def parameter_select_tai():
        # select reference GCN of tRNAs for tai computing
        tai_dict = {1: 'Homo_sapiens_GRCh38_GtRNAdb'}
        while True:
            try:
                tai_choice = int(input("select reference GCN of tRNAs for tai computing"
                                       "\nYour choice is: "
                                       "\n(1) Homo_sapiens_GRCh38_fromGtRNAdb[DEFAULT]"
                                       "\n"))
            except ValueError:
                print("Only one reference GCN could be chosen at a time.")
                continue
            select_parameters.tai_tGCN = tai_dict[tai_choice]
            break
        print("the GCN of tRNA for tai computing was selected\n")

    @staticmethod
    def parameter_select_ite():
        # select reference set for ite computing
        ite_dict = {1: 'Homo_sapiens_EMBOSS'}
        while True:
            try:
                ite_choice = int(input("select reference set for ite computing"
                                       "\nYour choice is: "
                                       "\n(1) 'Homo_sapiens_reference_fromEMBOSS[DEFAULT]'\n"))
            except ValueError:
                print("Only one reference set could be chosen at a time.")
                continue
            select_parameters.ite_heg = reference_heg_ite[ite_dict[ite_choice]]
            select_parameters.ite_bg = reference_bg_ite[ite_dict[ite_choice]]
            break
        print(f"reference set for ite computing was selected\n")

    @staticmethod
    def parameter_select_cbi():
        # select optimized codons for fop computing
        cbi_dict = {0: 'correlation_analysis'}
        while True:
            try:
                cbi_choice = int(input("select optimized codons for cbi computing"
                                       "\nYour choice is: "
                                       "\n(0) Correlation_analysis_for_optimized_codon[DEFAULT]"
                                       "\n"))
            except ValueError:
                print("Only one reference could be chosen at a time.")
                continue
            select_parameters.cbi_opt_codon = optimized_codon_cbi[cbi_dict[cbi_choice]]
            break
        print(f"optimized codon for cbi computing was selected\n")

    @staticmethod
    def process_parameter_file(parameter_file):
        # the input file welcomes an excel file or a dataframe,
        # namely if you have a parameter file in form of csv, make it a dataframe first.

        try:
            parameter_dataframe = pd.read_excel(parameter_file, index_col=0)
        except:
            if type(parameter_file) == pd.DataFrame:
                parameter_dataframe = parameter_file
            else:
                raise AttributeError('the input parameter file should be an excel file(.xlsx) or a dataframe.')

        # sorting the dataframe.
        parameter_dataframe = parameter_dataframe.loc[["TTT", "TTC", "TTA", "TTG",
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
                                                       "GGT", "GGC", "GGA", "GGG"]]

        # change parameters
        if 'cai_ref' in parameter_dataframe.columns:
            select_parameters.cai_reference_frequency = parameter_dataframe['cai_ref'].values.tolist()

        if 'tai_tGCN' in parameter_dataframe.columns:
            select_parameters.tai_tGCN = parameter_dataframe['tai_tGCN'].values.tolist()

        if 'ite_heg' in parameter_dataframe.columns:
            select_parameters.ite_heg = parameter_dataframe['ite_heg'].values.tolist()

        if 'ite_bg' in parameter_dataframe.columns:
            select_parameters.ite_bg = parameter_dataframe['ite_bg'].values.tolist()

        if 'cbi_opt' in parameter_dataframe.columns:
            # todo 注意导入的文件与实际计算时的可能不同。使得cbi和fop在计算时能容纳不同的opt。
            select_parameters.cbi_opt_codon = parameter_dataframe['cbi_opt'].values.tolist()

        if 'fop_opt' in parameter_dataframe.columns:
            select_parameters.fop_opt_codon = parameter_dataframe['fop_opt'].values.tolist()
        return

    fundict = {1: parameter_select_cai, 2: parameter_select_tai, 3: parameter_select_ite, 4: parameter_select_fop,
               5: parameter_select_cbi}

    def __init__(self, parameter_file=None, cai_ref=None, tai_tGCN=None, ite_heg=None, ite_bg=None, fop_opt=None,
                 cbi_opt=None, skip_select=False):

        additional_request = [parameter_file, cai_ref, tai_tGCN, ite_heg, ite_bg, fop_opt, cbi_opt]

        # if skip_select is True, you will select parameters among the prestored database.
        if not skip_select:
            while True:
                try:
                    indexes = list(
                        map(int, input("Press the serial number of index of which parameter you want to change: "
                                       "\n(1) cai "
                                       "\n(2) tai "
                                       "\n(3) ite "
                                       "\n(4) fop "
                                       "\n(5) cbi "
                                       "\n").split()))
                except ValueError:
                    print("Please input number")
                    continue
                for index in indexes:
                    if 0 < index < 6:
                        select_parameters.fundict[index]()
                    else:
                        print("The number is not between 1 and 5!")
                break

        # if your own parameter data is input, parameters will be replace by the data. the priority is
        # single_parameter_list(use cai_ref=) > parameter_file > parameter_selected. when you try to change a single
        # parameter, for example using cai_ref, make it equal to a list or an array is both OK but not a dataframe.
        if additional_request:
            if parameter_file:
                select_parameters.process_parameter_file(parameter_file=parameter_file)

            if cai_ref:
                select_parameters.cai_reference_frequency = list(cai_ref)

            if tai_tGCN:
                select_parameters.tai_tGCN = list(tai_tGCN)

            if ite_heg:
                select_parameters.ite_heg = list(ite_heg)

            if ite_bg:
                select_parameters.ite_bg = list(ite_bg)

            if cbi_opt:
                select_parameters.cbi_opt_codon = list(cbi_opt)

            if fop_opt:
                select_parameters.fop_opt_codon = list(fop_opt)

        return

    @staticmethod
    def generate_parameter_file(sorting='TCAG'):
        # use different sorting method, the result xlsx file will have sample name of different sorting.
        # for example, "TCAG" for "TTT" to "GGG", "amino_acid" for "GCA" to "TGA".
        if sorting == 'TCAG':
            parameter_dataframe = pd.DataFrame(index=["TTT", "TTC", "TTA", "TTG",
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
        elif sorting == 'amino_acid':
            parameter_dataframe = pd.DataFrame(
                index=['GCA', 'GCC', 'GCG', 'GCT', 'TGC', 'TGT', 'GAC', 'GAT', 'GAA', 'GAG', 'TTC', 'TTT',
                       'GGA', 'GGC',
                       'GGG', 'GGT', 'CAC', 'CAT', 'ATA', 'ATC', 'ATT', 'AAA', 'AAG', 'CTA', 'CTC', 'CTG',
                       'CTT', 'TTA',
                       'TTG', 'ATG', 'AAC', 'AAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CAA', 'CAG', 'AGA', 'AGG',
                       'CGA', 'CGC',
                       'CGG', 'CGT', 'AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT', 'ACA', 'ACC', 'ACG', 'ACT',
                       'GTA', 'GTC',
                       'GTG', 'GTT', 'TGG', 'TAC', 'TAT', 'TAA', 'TAG', 'TGA'])
        else:
            raise ValueError("Counldn't find the sorting method you use.")

        parameter_dataframe['cai_ref'] = select_parameters.cai_reference_frequency
        parameter_dataframe['tai_tGCN'] = select_parameters.tai_tGCN
        ite_dataframe = pd.DataFrame(index=["TTT", "TTC", "TTA", "TTG",
                                            "TCT", "TCC", "TCA", "TCG",
                                            "TAT", "TAC",
                                            "TGT", "TGC", "TGG",
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
        ite_dataframe['ite_heg'] = select_parameters.ite_heg
        ite_dataframe['ite_bg'] = select_parameters.ite_bg
        parameter_dataframe = pd.concat([parameter_dataframe, ite_dataframe], axis=1)

        if select_parameters.cbi_opt_codon:
            parameter_dataframe['cbi_opt'] = select_parameters.cbi_opt_codon
        if select_parameters.fop_opt_codon:
            parameter_dataframe['fop_opt'] = select_parameters.fop_opt_codon

        parameter_dataframe.to_excel('parameter_file.xlsx')

        return


# c = select_parameters( parameter_file='parameter_file.xlsx', cai_ref=[1], skip_select=True)
# c.process_parameter_file()
# print(c.cai_reference_frequency)
# print(c.cbi_opt_codon)
# print(c.ite_heg)
