# alignmentprofiler/stats.py

import numpy as np
from collections import Counter

def basic_stats(alignment):
    return {
        "num_sequences": len(alignment),
        "alignment_length": alignment.get_alignment_length(),
        "num_gapped_sites": sum(1 for i in range(alignment.get_alignment_length()) if '-' in alignment[:, i]),
        "num_conserved_sites": sum(1 for i in range(alignment.get_alignment_length()) if len(set(alignment[:, i])) == 1),
        "gap_fraction_total": sum(seq.seq.count('-') for seq in alignment) / (len(alignment) * alignment.get_alignment_length())
    }

def column_stats(alignment):
    aln_len = alignment.get_alignment_length()
    gap_counts = []
    entropy = []

    for i in range(aln_len):
        column = alignment[:, i]
        freqs = Counter(column)
        total = sum(freqs.values())
        probs = [count / total for count in freqs.values() if count > 0]
        ent = -sum(p * np.log2(p) for p in probs)
        entropy.append(ent)
        gap_counts.append(freqs.get('-', 0) / total)

    return {
        "mean_entropy": np.mean(entropy),
        "max_entropy": np.max(entropy),
        "mean_gap_fraction": np.mean(gap_counts),
        "max_gap_fraction": np.max(gap_counts)
    }

def sequence_stats(alignment):
    stats = []
    aln_len = alignment.get_alignment_length()
    for record in alignment:
        seq = str(record.seq)
        gaps = seq.count('-')
        a_count = seq.upper().count('A')
        t_count = seq.upper().count('T')
        g_count = seq.upper().count('G')
        c_count = seq.upper().count('C')
        valid = len(seq.replace('-', ''))

        stats.append({
            "id": record.id,
            "length": len(seq),
            "gap_fraction": gaps / len(seq),
            "gc_content": (g_count + c_count) / valid if valid > 0 else 0,
            "at_content": (a_count + t_count) / valid if valid > 0 else 0,
            "a_count": a_count,
            "t_count": t_count,
            "g_count": g_count,
            "c_count": c_count,
            "num_gaps": gaps,
            "ambiguous_bases": sum(seq.upper().count(x) for x in ['N', 'X'])
        })
    return stats
