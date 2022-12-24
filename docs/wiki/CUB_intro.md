---
title: Codon Usage Bias Introduction
permalink: /CUB_intro
layout: wiki
---

### What is codon usage bias？

•The protein sequences are composed of only 20 standard amino acids, yet there are 64 different codons.

• 61 codons encoding for amino acids plus three stop codons.

So there exists the "degeneracy" of the codons, and codon usage in organisms is not average or random, which leads to codon usage bias(CUB). There are many reasons for the existence of CUB, such as the following.

1. Origin of the genetic code system.

2. Mutation, selection, and random drift.

3. GC content, RNA structure, protein structure, etc.

There are many theories explaining the origin of codons. In general, codons undergo adaptive selection to achieve a required level of optimization that meets protein synthesis requirements throughout the life of an organism.

### What does CUB Effects?

1. Affect the translation elongation speed and efficiency.

2. Affect Initiation(Ramp).

3. Affect translation fidelity. Regulates co-translation protein folding and protein function.

4. Determining mRNA level.

### Applications

1. Heterologous gene expression (Synthetic biology).

2. Codon optimization.

3. Development of attenuated viral vaccines.

4. Gene therapy.

### Codon bias indices

To study CUB, many scientists have proposed indexes for evaluating it. You can easily calculate the following indexes in CUBAT.

- [RSCU(relative synonymous codon usage)](https://gxelab.github.io/CUBAT/indexes/RSCU)
- [ENC(effective number of codons, also called the Nc)](https://gxelab.github.io/CUBAT/indexes/ENC)
- [CAI(codon adaptation index)](https://gxelab.github.io/CUBAT/indexes/CAI)
- [CBI(codon bias index)](https://gxelab.github.io/CUBAT/indexes/CBI)
- [Fop (fraction of optimal codons)](https://gxelab.github.io/CUBAT/indexes/FOP)
- [tAI (tRNA adaption index)](https://gxelab.github.io/CUBAT/indexes/TAI)
- [CSC(codon stabilization coefficient)](https://gxelab.github.io/CUBAT/indexes/CSC)
- [the scaled χ2](https://gxelab.github.io/CUBAT/indexes/X2)
