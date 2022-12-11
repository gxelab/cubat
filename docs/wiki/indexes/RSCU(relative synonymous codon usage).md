---
title: RSCU
permalink: /CUBAT/wiki/indexes/RSCU
layout: wiki
redirect_from:
 - /CUBAT/indexes/RSCU
---

# Introduction

It is used to describe the relative usage of each codon in species (or part of species).

Its value ranges from 0 to 6.

Methionine, tryptophan, and stop codons are usually omitted.

# How to calculate?

$$
\operatorname{RSCU}_{i j}=\frac{x_{i j}}{\frac{1}{n_{i}} \sum_{j=1}^{n_{i}} x_{i j}}
$$

$$x_{i j}$$ is the number of occurrences of the jth codon for the ith amino acid.

$${\frac{1}{n_{i}} \sum_{j=1}^{n_{i}} x_{i j}}$$ is the average codon use of a codon family.

# Meaning of value

* 0 refers to the absence of this codon.
* 1 refer to this codon is no bias.
* 6 refers to this codon as the only codon used in a six-codon family.​

> 1. Codon usage in yeast: cluster analysis clearly differentiates highly and lowly expressed genes: page 3
> 2. *Synonymous but not the same: the causes and consequences of codon bias*: page 2

‍
