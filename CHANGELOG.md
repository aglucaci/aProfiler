```md
# Changelog

All notable changes to **aProfiler** will be documented in this file.  
The format follows *Keep a Changelog* and *Semantic Versioning* principles.

---

## [0.1.0] — 2025-11-26
### Added
- Initial release of the `MSA Statistics & Visualization Toolkit`
- Command-line program accessible via:
```

aprofiler --input alignment.fasta --mode auto

```
- Core MSA parsing via **BioPython**
- Auto-detected modes: `nt` and `aa` (with bias toward nucleotide sequences containing IUPAC ambiguity codes)
- Explicit **codon mode** analysis pipeline (never auto-detected)
- CSV table outputs:
- Global nucleotide or amino acid frequencies
- Per-site nucleotide counts, entropy, and GC%
- Codon usage table
- Relative Synonymous Codon Usage (RSCU) table
- Amino-acid usage derived from codons
- Essential vs non-essential amino-acid summary from codons
- PCA embedding matrix
- UMAP embedding matrix
- Plot outputs enabled by default:
- Nucleotide sequence logo
- Per-site entropy skyline
- Per-site GC% landscape
- PCA sequence clustering scatterplot
- UMAP sequence-space scatterplot
- Pairwise identity histogram
- Gap fraction histogram
- RSCU heatmap
- Essential AA barplot
- Sequence-derived feature heatmaps
- Minimal pytest coverage for all 3 pipeline modes (launch coverage)

---

## [0.1.1] — 2025-11-27
### Fixed
- Improved nucleotide vs amino-acid **mode detection** to safely support extended IUPAC alphabet
- Eliminated crashes in UMAP on small or sparse feature matrices by:
- Clamping `n_neighbors`
- Providing **PCA-derived fallback scatter column aliases (`UMAP1/UMAP2`)** for tiny inputs
- Ensured padding of ragged sequences at load instead of failing on indexing

### Added
- Internal accounting of skipped ambiguous codons in codon mode
- Optional report generator scaffolding for `.md` and `.html` outputs
- CLI flags prepared for future:
```

--report
--max-seqs
--max-sites
--seed
--no-embeddings

```

---

## [0.1.2] — 2025-11-28
### Added
- Richer alignment-scoped plots:
- Covariance matrix heatmap
- Per-amino-acid essential vs non-essential comparisons
- Expanded embedding visuals
- Artifact organization guarantees toward reproducibility by naming convention

---

## [0.2.0] — 2025-12-10
### Planned
- Stabilize protein mode to include:
- Per-site amino-acid entropy
- Amino-acid class fingerprints (hydrophobic, polar, charged, aromatic)
- Per-position amino-acid prevalence heatmaps
- Extend codon mode to include:
- Codon-position stratified summary metrics
- 2-fold vs 4-fold degenerate site statistics
- Configurable genetic code support
- Include automatic QC filtering presets for long MSAs
- Introduce `--summary-card` outputs for talk-ready or social-ready tables

---

### Notes
- Earlier versions focused on **robust personal profiling output**, not destructive alignment modification.
- Future versions will balance **scientific completeness**, **visual usability**, and **runtime stability at scale**.
```