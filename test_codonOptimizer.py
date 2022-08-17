import Bio
from Bio import SeqIO
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
import os
import numpy as np

def clean_CDS(input_file, output_file=False):
    """
    Accepts fasta file containing coding sequences.\n
    Returns or writes file containing only coding sequences
    that are divisible by 3, do not contain ambiguous nucleotides, and begin with ATG. 
    """
    records = []
    n_dropped = 0 
    for record in SeqIO.parse(input_file, "fasta"):
        if len(record.seq)%3 != 0:
            n_dropped += 1 
            pass
        elif "N" in str(record.seq):
            n_dropped += 1 
            pass
        elif str(record.seq[0:3]) != "ATG":
            n_dropped += 1 
            pass
        else:
            records.append(record)
    print(f"{n_dropped} records were removed. {len(records)} were maintained.")
    if output_file:
        SeqIO.write(records, output_file, "fasta")
        print(f"cleaned CDS file saved as {output_file}")
    else: 
        return records

def get_pref_codons(input_file, output_file=False):
    """
    Accepts coding sequences file (should be cleaned first to remove unusual sequences). \n
    Returns or saves file with codon adaptation index (CAI) scores.
    """
    index = CodonAdaptationIndex()
    index.generate_index(input_file)

    # Values of 1 is preferred codon. Otherwise, relative to preffered.
    preferred_codons_dict = index.index

    if output_file:
        import csv
        with open(output_file, 'w') as f:
            for key in preferred_codons_dict.keys():
                f.write("%s,%s\n"%(key,preferred_codons_dict[key]))
        print(f"preferred codons file saved as {output_file}")
    else:
        return preferred_codons_dict

def codon_optimize(input_file, preferred_codons_dict, output_file=False, gene_id=None, gene_description=None):
    index = CodonAdaptationIndex()
    index.set_cai_index(preferred_codons_dict)
    
    select_rec = SeqIO.read(f"{input_file}", "fasta")
    seq = str(select_rec.seq)
    print(f"Initial CAI: {index.cai_for_gene(seq)}")

    locs = np.where(np.array(list(preferred_codons_dict.values())) == 1.0)[0]
    pref = {}
    for i in locs:
        key = list(preferred_codons_dict.keys())[int(i)]
        val = Bio.Seq.translate(key)
        pref[val] = key

    optimized= ''
    for i in range(0,len(seq),3):
        try:
            codon = seq[i:i+3]
            aa = Bio.Seq.translate(codon)
            pref_codon = pref[aa]
            optimized += pref_codon
        except:
            print(f"Ignoring incomplete codon beginning at nucleotide {i}")
    
    print(f"Codon Optimized CAI: {index.cai_for_gene(optimized)}")
    opt_record = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(optimized))
    if output_file:
        opt_record.id = str(gene_id)
        opt_record.description = str(gene_description)
        if output_file == True:
            SeqIO.write(f"codon-optimized_{output_file}", "fasta")
        else:
            SeqIO.write(output_file, "fasta")
    else:
        return opt_record
