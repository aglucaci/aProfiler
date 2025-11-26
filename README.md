
<p align="center">
  <img src="logo/aProfiler_Logo_v1.png" alt="aProfiler Logo" width="260"/>
</p>

# **aProfiler**
## MSA Statistics & Visualization Toolkit

aProfiler examines **Multiple Sequence Alignments (MSAs)** and emits **useful statistics, publication-grade plots, codon-aware metrics (RSCU and essential amino acid summaries), embeddings, and CSV tables** from a single command-line call.

[![PyPI version](https://badge.fury.io/py/aprofiler.svg)](https://pypi.org/project/aprofiler/)
<!-- [![Tests](https://github.com/aglucaci/aProfiler/actions/workflows/pytest.yml/badge.svg)](https://github.com/aglucaci/aProfiler/actions/workflows/pytest.yml) -->
[![Python Versions](https://img.shields.io/pypi/pyversions/aprofiler)](https://pypi.org/project/aprofiler/)
[![License](https://img.shields.io/pypi/l/aprofiler.svg)](https://pypi.org/project/aprofiler/)

---

## Installation

```bash
pip install -e .
````

(For a future versioned release: `pip install aprofiler`)

---

## CLI Usage

```bash
aprofiler --input alignment.fasta --mode auto --report
```

### Modes

| Mode    | Description                                                                       |
| ------- | --------------------------------------------------------------------------------- |
| `nt`    | Nucleotide MSA (DNA/RNA, IUPAC tolerated)                                         |
| `aa`    | Amino acid MSA (protein residues)                                                 |
| `codon` | Coding sequence MSA analyzed at codon level (standard genetic code by default)    |
| `auto`  | Auto-detect between NT and AA (codon is never auto-selected and must be explicit) |

### Common Flags

| Flag              | Purpose                                                    |
| ----------------- | ---------------------------------------------------------- |
| `--input`         | Input MSA file (FASTA MSA, A3M, or fixed-column alignment) |
| `--mode`          | `nt`, `aa`, `codon`, or `auto`                             |
| `--report`        | Generate a summary report (`.md` or `.html`)               |
| `--report-format` | `md` or `html` (default: `md`)                             |
| `--no-plots`      | Skip plots, output CSV tables only                         |
| `--max-seqs`      | Downsample sequences for visualizations                    |
| `--max-sites`     | Crop columns for logos/heatmaps/logos                      |
| `--seed`          | Fixed seed for reproducible embeddings/visuals             |
| `--no-embeddings` | Skip dimensionality-reduction steps (PCA/UMAP)             |
| `--summary-card`  | Output condensed tables for social media or talks          |

---

## Outputs (Saved Automatically)

All results are saved under:

```
./results/{alignment_name}/
```

### CSV Tables

| Output                                              | File                         |
| --------------------------------------------------- | ---------------------------- |
| Global NT or AA frequencies                         | `*_global_freqs.csv`         |
| Per-site NT stats + entropy + GC%                   | `*_nt_per_site.csv`          |
| PCA embeddings                                      | `*_pca_embedding.csv`        |
| UMAP embeddings                                     | `*_umap_embedding.csv`       |
| Codon usage table                                   | `*_codon_global.csv`         |
| Relative Synonymous Codon Usage (RSCU)              | `*_codon_rscu.csv`           |
| Amino acid usage derived from codons                | `*_aa_from_codons.csv`       |
| Essential vs non-essential AA summary (from codons) | `*_aa_essential_summary.csv` |

### Plots (On by default unless disabled)

| Plot                        | Purpose                                                |
| --------------------------- | ------------------------------------------------------ |
| Nucleotide logo plot        | Position-wise base enrichment                          |
| GC% per-site                | GC landscape across MSA                                |
| Entropy per-site            | Conservation skyline                                   |
| AA/NT per-site heatmaps     | Residue/base prevalence                                |
| PCA scatter                 | Sequence-space clustering                              |
| UMAP scatter                | Similarity-space embedding                             |
| Pairwise identity histogram | Sequence similarity distribution                       |
| Gap fraction histogram      | Alignment completeness QC                              |
| Codon usage barplot         | Most frequent codons                                   |
| RSCU heatmap                | Synonymous codon bias by AA                            |
| Essential AA barplot        | Essential vs non-essential AA trends (codon mode only) |

---

## Alignment Format Compatibility and Constraints

| Input alignment type  | Supported? | Notes                                       |
| --------------------- | ---------- | ------------------------------------------- |
| FASTA MSA             | Yes        | Sequences must be equal length              |
| A3M                   | Yes        | Lowercase letters denote inserted columns   |
| ALN/Clustal           | Yes        | Must be in fixed columns or converted first |
| Codon FASTA           | Yes        | Requires explicit `--mode codon`            |
| Mixed NT+AA alphabets | No         | Alphabet must be uniform per file           |

> **All alignments are treated as fixed, rectangular matrices; sequences must have equal alignment length.**

---

## Output Guarantees

> **All profiling artifacts are written into `results/` without overwriting unrelated files. Ambiguous input characters are tolerated but tracked, not silently discarded. Codon mode metrics are only computed when explicitly requested.**

---

## Example Output Directory Tree

```
results/
  TP53_alignment/
    TP53_global_freqs.csv
    TP53_nt_per_site.csv
    TP53_pca_embedding.csv
    TP53_umap_embedding.csv
    TP53_report.md
    plots/
      TP53_nt_logo.png
      TP53_entropy.png
      TP53_gc.png
      TP53_pca.png
      TP53_umap.png
      TP53_pairwise_identity.png
      TP53_gap_fraction.png
      TP53_rscu_heatmap.png
      TP53_aa_essential_bar.png
```

---

## Scalability Notes

> **Optimized for MSAs up to ~20k sequences × 10k columns on standard hardware. Larger inputs may require downsampling for logos/heatmaps (`--max-seqs`, `--max-sites`).**

---

## Reproducibility

> **All embeddings (PCA/UMAP) and stochastic plots use fixed `random_state` when `--seed` is passed for reproducibility.**

---

## Testing

```bash
pip install pytest
pytest -q
```

Tests validate:

* equal-length enforcement
* stable fallback embedding behavior
* non-empty CSV outputs
* plot and artifact creation without crashes

---

## Python API Example

```python
from aprofiler.profiler import AlignmentProfiler

prof = AlignmentProfiler("alignment.fasta", mode="auto", out_dir="results")
prof.load_alignment()
outputs = prof.run_full_profile()
report_path = prof.generate_report(outputs, fmt="md")

print("Outputs generated:", outputs)
print("Report saved to:", report_path)
```

---

## Citation

If you use aProfiler in a publication, please cite:

```
Lucaci, Alexander G., *aProfiler: MSA Statistics & Visualization Toolkit*, 2025.
```

For formal reproducible citation, you can later replace this with a Zenodo or PyPI DOI once released.

---

## Contributing

Issues, discussions, and pull requests are welcome. Ensure contributions are:

* statistically useful
* plot-rich by design
* free of silent failures
* non-destructive to unrelated files
* aligned with package philosophy and constraints

---

## Brand Identity

aProfiler is **alignment-useful, metrics-first, codon-aware when requested, embedding-stable, plot-default, science-fun.**
