---
title: CAI
permalink: indexes/CAI
layout: wiki
---

### Introduction

The efficiency of translating mRNA to protein depends partially on the coding strategy of an mRNA. It is reflected in codon usage bias, which is often measured by codon-specific and gene-specific indices. A representative of the first class is the relatively synonymous codon usage or RSCU (Sharp et al. 1986). A representative of the second class is the codon adaptation index, or CAI (Sharp and Li, 1987).

### How to calculate?

#### Old CAI

$$
w_{i j}=\frac{f_{i j . r e f}}{\operatorname{Maxf}_{i, v e f}}
$$

$$
CAI=\exp \left(\frac{\sum_{i=1}^{m} \sum_{j=1}^{n_{i}}\left[f_{i j} \ln \left(w_{i j}\right)\right]}{\sum_{i=1}^{m} \sum_{j=1}^{n_{i}} f_{i j}}\right)
$$

#### New CAI

$$
CAI2= \frac{\sum_{i=1}^{m} \sum_{j=1}^{n_{i}}\left[f_{i j} w_{i j}\right]}{\sum_{i=1}^{m} \sum_{j=1}^{n_{i}} f_{i j}}
$$

### Meaning of value

CAI refers to the adaptation coefficient of all codons encoding the protein for a certain generative if the gene uses the optimal codon. The CAI value is within the range of 0–1. A more considerable value indicates more robust adaptability. CAI value is widely used in the evaluation of gene expression levels.

> Sharp et al. 1986
> 
> Sharp and Li, 1987
> 
> An Improved Implementation of Codon Adaptation Index

‍
