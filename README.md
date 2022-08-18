codonOptimizer.py
===============

Introduction
------------

codonOptimizer.py is a module with functions to calculate relative synonymous codon usage (RSCU) values, 
and to use RSCU values to create codon optimized versions of provided gene or coding sequences. 

Example use cases are provided below and can be run using a cloned version of this repository as your working directory. 

Requirements
------------

codonOptimizer.py relies on Biopython.
Installation requires setuptools.
The example below requires numpy, pandas, and matplotlib.

Installation
------------

Install using pip:
`pip install -i https://test.pypi.org/simple/ codonOptimizer==0.1.0`

Or clone this repo and import the module manually

Example
-------

### Comparing codon usage patterns in *Saccharomyces cerevisiae* and *Aureobassidium pullulans*
Collaborators were unable to express GFP-tagged proteins in the fungus *Aureobassidium pullulans*, and were concered that codon preferences in *A. pullulans* could be causing poor expression. The following code compares codon usage between *S. cerevisiae*, the organism for which this GFP was codon optimized previously, and *A. pullulans*. This example relies on the very popular librarires numpy, pandas, and matplotlib. 
```
import codonOptimizer as co

co.clean_CDS("README_files/ScCDS.fna", output_file="README_files/ApCDS.fna")
co.get_pref_codons("README_files/ApCDS.fna", output_file="README_files/Ap_pref_codons.fasta")

co.get_pref_codons("README_files/ScCDS.fna", output_file="README_files/Sc_pref_codons.fasta")

import pandas as pd

Ap_pref = pd.read_csv("README_files/Ap_pref_codons.fasta", names=["Codon", "CAI"], header=None)
Sc_pref = pd.read_csv("README_files/Sc_pref_codons.fasta", names=["Codon", "CAI"], header=None)

import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(5, 15))
im = ax.imshow(np.vstack((Ap_pref["CAI"],Sc_pref["CAI"])).T,cmap="Reds")
fig.colorbar(im)
ax.set_yticks(np.arange(len(Ap_pref["Codon"])))
ax.set_xticks(np.arange(2))
ax.set_yticklabels(labels=Ap_pref["Codon"])
ax.set_xticklabels(labels=("A. pullulans", "S. cerevisiae"))
ax.set_aspect('auto')
```
![](README_files/example.png)

### Codon optimizing green fluorescent protein for expression in *Aureobassidium pullulans*
Codon preferences appeared quite different between these two organisms, so here we convert every codon in the GFP sequence into the *A. pullulans* preferred codon for that amino acid. 
```
import codonOptimizer as co

co.clean_CDS("README_files/CDS.fna", output_file="README_files/ApCDS.fna")
pref_codon_dict = co.get_pref_codons("README_files/ApCDS.fna")

co.codon_optimize("README_files/GFP.fa", pref_codon_dict, output_file="README_files/ApOpt_GFP.fa")
```
