---
title: ENC
permalink: indexes/ENC
layout: wiki
---

### Introduction

The codon usage table of a gene can be subdivided according to the number of synonymous codons belonging to each aa. Thus for a gene using the 'universal' code, there is two aa with only one codon choice, 9 with two, 1 with three, 5 with four, and 3 with six. These represent five SF types, designated SF types 1, 2, 3, 4, and 6, according to their number of synonymous codons. The measure of SCU bias developed here will involve combining contributions to overall bias from each of the five SF types. However, first, we must consider the more straightforward problem of measuring the contribution era single aa to the overall SCU bias of a gene.

An aa with four synonymous codons is analogous to a locus with four alleles. Equal codon usage would be equivalent to minimum homozygosity. Homozygosity (F) can be calculated from the squared allele (codon) frequencies:

### How to calculate?

#### Old Nc

$$
\hat{F}=\frac{n \sum_{i=1}^{k} p_{i}^{2}-1}{n-1}
$$

* n is the total usage of aa, $n = n_1 + .. + n_k$
* p is the usage frequency of the synonymous codons, $p_1 = \frac{n_1}{n}$ $p_2 = \frac{n_2}{n}$...

$$
\hat{N}_{c}=2+\left(9 / \bar{\hat{F}}_{2}\right)+\left(1 / \hat{F}_{3}\right)+\left(5 /\bar{\hat{F}}_{4}\right)+\left(3 / \bar{\hat{F}}_{6}\right)
$$

> The ‘effective number of codons’ used in a gene

#### New Nc

$$
F_{C F}=\sum_{i=1}^{m} p_{i}^{2}
$$

$$
N_{c}=N_{\mathrm{s}}+\frac{K_{2} \sum_{j}^{K_{2}} n_{j}}{\sum_{j=1}^{K_{2}}\left(n_{j} F_{C F . j}\right)}+\frac{K_{3} \sum_{j}^{K_{3}} n_{j}}{\sum_{j=1}^{K_{3}}\left(n_{j} F_{C F . j}\right)}+\frac{K_{4} \sum_{j}^{K_{4}} n_{j}}{\sum_{j=1}^{K_{4}}\left(n_{j} F_{C F . j}\right)}
$$

> An Improved Implementation of Effective Number of Codons (Nc)

‍
