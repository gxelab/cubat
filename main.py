import re

import pandas as pd
from Bio.Seq import Seq
import numpy as np
import csv


def get_data(data_path, line_begin, line_end):
    # Extract sequences and merge into one line string
    import linecache
    lines = linecache.getlines(data_path)[line_begin:line_end]  # Extract sequences
    lines = [line.strip() for line in lines]  # Remove '\n'
    lines = ''.join(lines)  # Merge into one line string
    return lines


def cut(obj, sec):
    # Function for cutting sequence
    return [obj[i:i + sec] for i in range(0, len(obj), sec)]


def locate_start_codon(selected_start_codon, mRNA_seq):
    start_codon = selected_start_codon + '.*'
    jump_to_starting_point = re.search(start_codon, str(mRNA_seq))  # Jump to start point
    starting_point = jump_to_starting_point.group(0)
    return starting_point


def locate_stop_codon(cut_seq):
    if cut_seq.count('UGA'):
        return cut_seq.index('UGA')
    elif cut_seq.count('UAA'):
        return cut_seq.index('UAA')
    elif cut_seq.count('UAG'):
        return cut_seq.index('UAG')


def process_mRNA_seq(location_of_stop_codon):
    return cut_seq[0:location_of_stop_codon]


Sars_cov_2 = get_data('Test_Data/Sars_cov_2.ASM985889v3.cds.fasta', 1, 356)  # Extract a sequence
DNA_seq = Seq(Sars_cov_2)
mRNA_seq = DNA_seq.transcribe()  # Transcribed into mRNA
located_seq = locate_start_codon('AUG', mRNA_seq)
cut_seq = cut(located_seq, 3)  # Cutting RNA sequence
location_of_stop_codon = locate_stop_codon(cut_seq)
processed_mRNA_seq = process_mRNA_seq(location_of_stop_codon)
result = pd.value_counts(processed_mRNA_seq)
print(result)
`