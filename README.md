
# CUBAT <img src='docs/mstile-150x150.png' align="right" height="139" />

Codon Usage Bias Analysis Toolkit (CUBAT) is a powerful python package for codon usage bias (CUB) analysis. CUBAT can be imported as a module for custom analyses and has a command line interface as well.

### Installation

You can install `cubat` with `pip` from command line:

```bash
pip install cubat
```

To update, run the following command:

```bash
pip install cubat --upgrade
```

### Main features
CUBAT can calculate various indices related to codon usage. Here is a summary of indices supported by CUBAT in comparison with other popular software.

|  Index                                 | CUBAT | codonW | DAMBE | EMBOSS |
| -------------------------------------  | ----- | ------ | ----- | ------ |
| RSCU (relative synonymous codon usage) | √     | √      | √     |        |
| Nc (effective number of codons)        |       | √      |       | √      |
| Nc (effective number of codons,SYX13)  | √     |        | √     |        |
| CAI (codon adaptation index)           | √     | √      | √     | √      |
| CAI2 (Xuhua Xia,2007)                  | √     |        | √     |        |
| CBI (codon bias index)                 | √     | √      | √     |        |
| Fop (frequency of optimal codons)      | √     | √      | √     |        |
| TAI (tRNA adaptation index)            | √     |        |       |        |
| CSC (codon stabilization coefficient)  | √     |        |       |        |
| Scaled χ<sup>2<sup>                    | √     |        |       |        |
| **Other Features**                     |
| Amino acid usage                       | √     | √      | √     |        |
| Custom genetic code                    | √     |        | √     |        |
| Cross-platform                         | √     |        | √     | √      |
