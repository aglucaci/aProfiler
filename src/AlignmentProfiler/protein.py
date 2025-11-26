# alignmentprofiler/protein.py

from Bio.Seq import Seq
from Bio.Data import IUPACData
from collections import Counter
import math

ESSENTIAL_AA = {'H', 'I', 'L', 'K', 'M', 'F', 'T', 'W', 'V'}

def translate_alignment(alignment):
    proteins = []
    for record in alignment:
        seq = str(record.seq).replace("-", "")  # remove gaps
        coding_dna = Seq(seq)
        try:
            protein = str(coding_dna.translate(to_stop=True))
        except Exception:
            protein = ""
        proteins.append(protein)
    return proteins

def amino_acid_frequencies(protein_seqs):
    total_counts = Counter()
    total = 0
    for protein in protein_seqs:
        total_counts.update(protein)
        total += len(protein)
    freqs = {aa: count / total for aa, count in total_counts.items() if total > 0}
    return freqs

def essential_aa_stats(protein_seqs):
    essential = 0
    nonessential = 0
    for seq in protein_seqs:
        for aa in seq:
            if aa in ESSENTIAL_AA:
                essential += 1
            elif aa in IUPACData.protein_letters:  # valid amino acids
                nonessential += 1
    return {
        "essential_aa": essential,
        "nonessential_aa": nonessential,
        "ratio": essential / nonessential if nonessential > 0 else None
    }

def protein_diversity(protein_seqs):
    unique = set(protein_seqs)
    entropy = -sum((protein_seqs.count(p)/len(protein_seqs)) * math.log2(protein_seqs.count(p)/len(protein_seqs)) for p in unique)
    return {
        "num_unique_proteins": len(unique),
        "shannon_entropy": entropy
    }

