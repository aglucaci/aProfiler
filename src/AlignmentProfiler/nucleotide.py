# alignmentprofiler/nucleotide.py

from collections import Counter

VALID_BASES = ['A', 'T', 'G', 'C']
AMBIGUOUS_BASES = ['N', 'X', 'Y', 'R', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V']

def nucleotide_frequencies(alignment):
    counts = Counter()
    for record in alignment:
        counts.update(str(record.seq).upper())

    total = sum(counts[base] for base in VALID_BASES)
    freqs = {base: counts[base] / total if total else 0 for base in VALID_BASES}
    freqs['total_bases'] = total
    return freqs

def ambiguous_character_count(alignment):
    counts = Counter()
    for record in alignment:
        seq = str(record.seq).upper()
        for base in AMBIGUOUS_BASES:
            counts[base] += seq.count(base)
    return dict(counts)

def gap_character_analysis(alignment):
    gap_total = sum(str(record.seq).count('-') for record in alignment)
    ungapped_total = sum(len(str(record.seq).replace('-', '')) for record in alignment)
    total = gap_total + ungapped_total
    return {
        "gap_fraction": gap_total / total if total else 0,
        "total_gaps": gap_total,
        "total_ungapped": ungapped_total
    }

