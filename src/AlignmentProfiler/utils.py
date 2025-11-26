# alignmentprofiler/utils.py

def is_valid_codon(codon):
    return len(codon) == 3 and all(base in "ATGCatgc" for base in codon)

def normalize_list(values):
    total = sum(values)
    return [v / total if total else 0 for v in values]

def sliding_window(seq, window_size, step=1):
    for i in range(0, len(seq) - window_size + 1, step):
        yield seq[i:i + window_size]

