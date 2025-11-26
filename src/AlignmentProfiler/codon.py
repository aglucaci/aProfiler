# alignmentprofiler/codon.py

from collections import Counter, defaultdict
from Bio.Seq import Seq
from Bio.Data import CodonTable

def get_translated_codons(alignment):
    codon_matrix = []
    for record in alignment:
        seq = str(record.seq).replace('-', '')  # ungapped for translation
        codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3) if len(seq[i:i+3]) == 3]
        codon_matrix.append(codons)
    return codon_matrix

def codon_diversity(alignment):
    codon_matrix = get_translated_codons(alignment)
    site_diversity = []
    for i in range(min(map(len, codon_matrix))):  # loop over codon columns
        site_codons = [row[i] for row in codon_matrix]
        unique_codons = set(site_codons)
        site_diversity.append(len(unique_codons))
    return {
        "mean_codon_diversity": sum(site_diversity) / len(site_diversity),
        "max_codon_diversity": max(site_diversity),
        "sitewise_codon_diversity": site_diversity
    }

def rscu(alignment):
    codon_usage = Counter()
    amino_usage = defaultdict(Counter)
    table = CodonTable.unambiguous_dna_by_name["Standard"]
    codon_matrix = get_translated_codons(alignment)

    for codons in codon_matrix:
        for codon in codons:
            if codon in table.forward_table:
                aa = table.forward_table[codon]
                codon_usage[codon] += 1
                amino_usage[aa][codon] += 1

    rscu_table = {}
    for aa, codons in amino_usage.items():
        total = sum(codons.values())
        n = len(codons)
        for codon, count in codons.items():
            expected = total / n if n else 0
            rscu_table[codon] = round(count / expected, 3) if expected else 0

    return rscu_table

