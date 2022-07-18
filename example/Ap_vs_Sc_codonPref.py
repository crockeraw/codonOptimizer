import codonOptimizer as co

Ap_records = co.clean_CDS("source-data/CDS.fna", output_file="source-data/ApCDS.fna")
co.get_pref_codons("source-data/ApCDS.fna", output_file="Ap_pref_codons.fasta")

co.get_pref_codons("source-data/Sc_CDS.fna", output_file="Sc_pref_codons.fasta")

import pandas as pd

Ap_pref = pd.read_csv("Ap_pref_codons.fasta", names=["Codon", "CAI"], header=None)
Sc_pref = pd.read_csv("Sc_pref_codons.fasta", names=["Codon", "CAI"], header=None)

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