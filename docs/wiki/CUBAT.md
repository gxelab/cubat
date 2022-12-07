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

You can calculate multiple codon bias indexes with CUBAT.​

* [RSCU(relative synonymous codon usage)](https://gxelab.github.io/CUBAT/wiki/indexes/RSCU)
* ENC(effective number of codons, also called the Nc)
* CAI(codon adaptation index)
* CBI(codon bias index)
* Fop (fraction of optimal codons)
* tAI (tRNA adaption index)
* CSC(codon stabilization coefficient)
* the scaled χ2

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
