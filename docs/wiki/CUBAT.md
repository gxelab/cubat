---
title: CUBAT
permalink: /
layout: wiki
redirect_from:
 - wiki/
 - wiki/Website
 - wiki/Website_HTML
 - wiki/Main_Page
 - wiki/Main
 - wiki/CUBAT
---

### Introduction

CUBAT is a python program to help analyze codon usage bias(CUB). You can run it as a script or the command line.

### Installation

For install:

```bash
pip install cubat
```

For update:

```bash
pip install cubat --upgrade
```

### Usage
You can easily calculate multiple codon bias indexes with CUBAT.​

For example: You want calculate CAI(with reference of human), ENC and RSCU of SARS-COV-2.
All you need to do is provide a fasta file of SARS-COV-2 and human reference for cai(You can use the built-in data).

Run this:
```bash
cubat analyze --cr example/cai_ref.csv  -erc  Test_Data/Sars_cov_2.fasta
```
You will get two csv sheets with your desired results there.

### And more...
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
