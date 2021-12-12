import re
# import biopython
from Bio.Seq import Seq
import numpy as np
import csv
import linecache


def cut(obj, sec):
    # Function for cutting sequence
    return [obj[i:i + sec] for i in range(0, len(obj), sec)]


testData = linecache.getline('Test_SRA/SRR17184641.fastq', 2)  # Extract a sequence from SRR17184641
DNA_seq = Seq(testData)
mRNA_seq = DNA_seq.reverse_complement().transcribe()  # Transcribed into mRNA
startCodon = 'AUG.*'
jumpStartingPoint = re.search(startCodon, str(mRNA_seq))  # Jump to start point
StartingPoint = jumpStartingPoint.group(0)
result = cut(str(StartingPoint), 3)  # Cutting RNA sequence
print(result)
