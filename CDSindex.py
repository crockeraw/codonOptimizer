from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
import os

records = []
for record in SeqIO.parse("source-data/CDS.fna", "fasta"):
    if len(record.seq)%3 != 0:
        pass
    elif "N" in str(record.seq):
        pass
    else:
        records.append(record)

os.system("mkdir derived-data")

SeqIO.write(records, "derived-data/CDS_div_by3.fna", "fasta")

index = CodonAdaptationIndex()
index.generate_index("derived-data/CDS_div_by3.fna")

index.index