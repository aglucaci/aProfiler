import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import Counter
from math import log2

def plot_entropy(alignment, save_path=None):
    aln_len = alignment.get_alignment_length()
    entropy = []

    for i in range(aln_len):
        column = alignment[:, i]
        freqs = Counter(column)
        total = sum(freqs.values())
        probs = [count / total for count in freqs.values() if count > 0]
        col_entropy = -sum(p * log2(p) for p in probs)
        entropy.append(col_entropy)

    plt.figure(figsize=(12, 4))
    plt.plot(entropy, lw=1)
    plt.title("Shannon Entropy per Column")
    plt.xlabel("Position")
    plt.ylabel("Entropy")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

def plot_gap_heatmap(alignment, save_path=None):
    gap_matrix = []

    for record in alignment:
        row = [1 if char == '-' else 0 for char in str(record.seq)]
        gap_matrix.append(row)

    plt.figure(figsize=(12, len(alignment) * 0.3))
    sns.heatmap(gap_matrix, cbar=False)
    plt.title("Gap Heatmap (1 = gap)")
    plt.xlabel("Position")
    plt.ylabel("Sequences")
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

def plot_codon_diversity(diversity_summary, save_path=None):
    sitewise_div = diversity_summary.get("sitewise_codon_diversity", [])
    if not sitewise_div:
        print("[!] No sitewise codon diversity data found.")
        return

    plt.figure(figsize=(12, 4))
    plt.plot(sitewise_div, lw=1)
    plt.title("Codon Diversity Per Codon Site")
    plt.xlabel("Codon Position")
    plt.ylabel("Unique Codons")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
        print(f"[+] Codon diversity plot saved to {save_path}")
    else:
        plt.show()

def plot_protein_diversity(diversity_summary, save_path=None):
    entropy = diversity_summary.get("shannon_entropy")
    num_unique = diversity_summary.get("num_unique_proteins")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(["Unique Proteins", "Shannon Entropy"], [num_unique, entropy])
    ax.set_ylabel("Value")
    ax.set_title("Protein Diversity Metrics")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        print(f"[+] Protein diversity plot saved to {save_path}")
    else:
        plt.show()
