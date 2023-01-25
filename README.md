 ![](/docs/mstile-150x150.png)

### Introduction

CUBAT is a python program to help analyze codon usage bias(CUB). You can run it as a script or the command line.

### Installation

For install:

```
pip install cubat
```

For update:

```
pip install cubat --upgrade
```

You can easily calculate multiple codon bias indexes with CUBAT.​

The following is a comparison with other software that can calculate CUB indexes.

|                                       | CUBAT | codonW | DAMBE | EMBOSS |
| ------------------------------------- | ----- | ------ | ----- | ------ |
| RSCU(relative synonymous codon usage) | √     | √      | √     |        |
| Nc:(effective number of codons)       |       | √      |       | √      |
| Nc(effective number of codons,SYX13)  | √     |        | √     |        |
| CAI(codon adaptation index)           | √     | √      | √     | √      |
| CAI2(Xuhua Xia,2007)                  | √     |        | √     |        |
| CBI(codon bias index)                 | √     | √      | √     |        |
| Fop(frequency of optimal codons)      | √     | √      | √     |        |
| TAI(tRNA adaptation index)            | √     |        |       |        |
| CSC(codon stabilization coefficient)  | √     |        |       |        |
| the scaled χ2                         | √     |        |       |        |
| Amino acid usage                      | √     | √      | √     |        |
| Codon table replaceability            | √     |        | √     |        |
| cross-platform                        | √     |        | √     | √      |
