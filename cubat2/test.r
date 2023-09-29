library(Biostrings)
library(cubar)
cds = readDNAStringSet("../test_data/Drosophila_melanogaster.fna")
cf = count_codons(cds)

