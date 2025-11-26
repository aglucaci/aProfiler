# -*- coding: utf-8 -*-
"""

@author: Alexander G. Lucaci
"""

# =============================================================================
# Imports
# =============================================================================

from __future__ import annotations

from collections import Counter
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import csv
import itertools
import random

import numpy as np
import pandas as pd

from Bio import SeqIO
from scipy.stats import entropy
from tqdm.auto import tqdm

import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from umap import UMAP

import logomaker as lm
from plotnine import ggplot, aes, geom_point, theme_minimal, labs

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="umap")
warnings.filterwarnings("ignore", category=UserWarning, module="plotnine")
warnings.filterwarnings("ignore", category=UserWarning, module="plotnine.ggplot")
#warnings.filterwarnings("ignore", category=PlotnineWarning, module="plotnine")
warnings.filterwarnings("ignore", message="Tight layout not applied")

# =============================================================================
# Declares
# =============================================================================

# Standard genetic code (DNA) – single-letter AAs, '*' for stop
CODON_TABLE = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",

    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",

    "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",

    "TGT": "C", "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Essential amino acids (generic mammalian set).
# You can tweak if needed.
ESSENTIAL_AA = set(["H", "I", "L", "K", "M", "F", "T", "W", "V"])

# =============================================================================
# Main Class
# =============================================================================

class AlignmentProfiler:
    """
    Core engine for aProfiler.

    Modes:
      - 'nt'    : nucleotide alignment (A/C/G/T/U + IUPAC)
      - 'aa'    : amino-acid alignment (treated as symbols for now)
      - 'codon' : codon-aware nucleotide alignment (length ~ multiple of 3)
      - 'auto'  : auto-detect nt vs aa (never 'codon' automatically)

    Outputs:
      - CSV tables (results/)
      - PNG plots (results/)
    """

    def __init__(
        self,
        input_path: str,
        mode: str = "auto",
        out_dir: Optional[str | Path] = None,
        max_seqs: Optional[int] = None,
        max_sites: Optional[int] = None,
        seed: Optional[int] = None,
        embeddings: bool = True,
    ) -> None:
        self.path = Path(input_path)
        self.mode = mode
        self.out_dir = Path(out_dir) if out_dir is not None else Path("results")

        self.max_seqs = max_seqs
        self.max_sites = max_sites
        self.random_state = seed
        self.embeddings_enabled = embeddings

        self.headers = []
        self.sequences = []
        self.alignment_length = 0

        self.out_dir.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------------------
    # Mode detection (IUPAC-aware)
    # -------------------------------------------------------------------------
    def detect_sequence_mode(self) -> str:
        """
        Infer nt vs aa based on composition.

        Heuristic:
        - Uses full IUPAC nucleotide alphabet: A C G T U R Y K M S W B D H V N
        - Compares fraction of those vs std 20 amino acids.
        - Biases toward 'nt' in ambiguous cases.
        """
        if not self.sequences:
            raise RuntimeError("Alignment not loaded; cannot auto-detect mode.")

        combined = "".join(self.sequences).upper()
        letters = [ch for ch in combined if ch.isalpha()]
        if not letters:
            print("[aProfiler] No alphabetic characters found; defaulting to 'nt'.")
            return "nt"

        nuc_extended = set("ACGTURYKMSWBDHVN")   # DNA/RNA + ambiguity
        aa_set = set("ACDEFGHIKLMNPQRSTVWY")    # standard 20 AAs

        nuc_count = sum(1 for ch in letters if ch in nuc_extended)
        aa_count = sum(1 for ch in letters if ch in aa_set)
        total = len(letters)

        frac_nuc = nuc_count / total
        frac_aa = aa_count / total

        if frac_nuc >= 0.7 and frac_nuc >= 2 * frac_aa:
            mode = "nt"
        elif frac_aa >= 0.7 and frac_aa >= 2 * frac_nuc:
            mode = "aa"
        else:
            # Ambiguous; bias toward nt but log
            mode = "nt" if frac_nuc >= frac_aa else "aa"
            print(
                "[aProfiler] Warning: ambiguous composition for mode detection.\n"
                f"  frac_nuc={frac_nuc:.3f}, frac_aa={frac_aa:.3f}, total={total}\n"
                f"  → falling back to '{mode}'. Use --mode nt/aa to override."
            )

        print(
            f"[aProfiler] Auto-detected mode: {mode} "
            f"(frac_nuc={frac_nuc:.3f}, frac_aa={frac_aa:.3f}, total={total})"
        )
        return mode

    # -------------------------------------------------------------------------
    # I/O
    # -------------------------------------------------------------------------
    def load_alignment(self) -> None:
        """
        Load an MSA from a FASTA file into memory.

        - Reads all sequences with BioPython.
        - Forces a rectangular alignment by trimming to the minimum length
          if there is minor length variation.
        - Sets `self.alignment_length` to the resulting alignment width.
        """

        if not self.path.exists():
            raise FileNotFoundError(f"Input file not found: {self.path}")

        records = list(SeqIO.parse(str(self.path), "fasta"))

        if not records:
            raise ValueError(f"No sequences found in FASTA file: {self.path}")

        headers: List[str] = []
        seqs_raw: List[str] = []

        for rec in records:
            headers.append(rec.id)
            seqs_raw.append(str(rec.seq).upper())

        lengths = [len(s) for s in seqs_raw]
        min_len = min(lengths)
        max_len = max(lengths)

        if min_len == 0:
            raise ValueError(
                f"At least one sequence in {self.path} has length 0. "
                "Alignment length cannot be zero."
            )

        if min_len != max_len:
            print(
                f"[aProfiler] Warning: alignment has variable sequence lengths "
                f"(min={min_len}, max={max_len}). "
                f"Truncating all sequences to min length {min_len}."
            )
            seqs = [s[:min_len] for s in seqs_raw]
            aln_len = min_len
        else:
            seqs = seqs_raw
            aln_len = max_len

        self.headers = headers
        self.sequences = seqs
        self.alignment_length = aln_len

        # Auto-detect mode if requested
        if self.mode == "auto":
            self.mode = self.detect_sequence_mode()
        else:
            print(f"[aProfiler] Mode set by user: {self.mode}")

        print("[aProfiler] Loaded alignment")
        print(f"  • File      : {self.path}")
        print(f"  • Sequences : {len(self.sequences)}")
        print(f"  • Length    : {self.alignment_length}")
        print(f"  • Mode      : {self.mode}")

    # -------------------------------------------------------------------------
    # Core representations
    # -------------------------------------------------------------------------
    def to_array(self) -> np.ndarray:
        """Return alignment as (n_seq, L) array of single-character strings."""
        if not self.sequences:
            raise RuntimeError("Alignment not loaded.")
        return np.array([list(s) for s in self.sequences], dtype="U1")

    # -------------------------------------------------------------------------
    # GLOBAL FREQUENCIES (nt or generic symbol)
    # -------------------------------------------------------------------------
    def global_frequencies(
        self,
        include_gaps: bool = True,
        nt_mode: bool = True,
    ) -> Tuple[Dict[str, int], Dict[str, float], int]:
        """
        Global symbol frequencies.

        If nt_mode=True, treat A/C/G/T specially and group others as N/other.
        Otherwise, just count everything as generic symbols.
        """
        if not self.sequences:
            raise RuntimeError("Alignment not loaded.")

        all_seq = "".join(self.sequences)
        counts_raw = Counter(all_seq)
        counts: Dict[str, int] = {}
        total = 0

        if nt_mode:
            canonical = ["A", "C", "G", "T"]
            for base in canonical:
                c = counts_raw.pop(base, 0)
                counts[base] = c
                total += c

            gap_count = counts_raw.pop("-", 0)
            if include_gaps and gap_count > 0:
                counts["-"] = gap_count
                total += gap_count

            other_count = sum(counts_raw.values())
            if other_count > 0:
                counts["N/other"] = other_count
                total += other_count
        else:
            # For AA or generic symbols, just count everything including '-'
            for sym, c in counts_raw.items():
                counts[sym] = c
                total += c

        freqs = {k: (v / total if total else 0.0) for k, v in counts.items()}
        return counts, freqs, total

    def global_freqs_dataframe(self) -> pd.DataFrame:
        """Return global frequencies as a DataFrame (mode-aware)."""
        nt_mode = (self.mode == "nt")
        counts, freqs, total = self.global_frequencies(include_gaps=True, nt_mode=nt_mode)
        rows = []
        for sym in counts:
            rows.append(
                {
                    "symbol": sym,
                    "count": counts[sym],
                    "frequency": freqs[sym],
                    "total_alignment_chars": total,
                }
            )
        return pd.DataFrame(rows)

    def save_global_freqs_csv(self) -> Path:
        """Save global frequencies (mode-aware) to CSV."""
        df = self.global_freqs_dataframe()
        suffix = "nucfreq_global" if self.mode == "nt" else "symfreq_global"
        out_file = self.out_dir / f"{self.path.stem}_{suffix}.csv"
        df.to_csv(out_file, index=False)
        return out_file

    def plot_global_freqs(self) -> Path:
        """Barplot of global frequencies (mode-aware)."""
        df = self.global_freqs_dataframe()
        df = df.sort_values("symbol")

        fig, ax = plt.subplots(figsize=(7, 4))
        sns.barplot(data=df, x="symbol", y="frequency", ax=ax)
        ax.set_xlabel("Symbol")
        ax.set_ylabel("Frequency")
        ax.set_title(f"Global frequencies ({self.mode}) – {self.path.stem}")
        fig.tight_layout()

        suffix = "nucfreq_global" if self.mode == "nt" else "symfreq_global"
        out_file = self.out_dir / f"{self.path.stem}_{suffix}.png"
        fig.savefig(out_file, dpi=300)
        plt.close(fig)
        return out_file

    # -------------------------------------------------------------------------
    # PER-SITE NUCLEOTIDE METRICS (nt mode)
    # -------------------------------------------------------------------------
    def per_site_nt_metrics(self, include_gaps: bool = True) -> pd.DataFrame:
        """
        Per-site nucleotide metrics: A/C/G/T counts, gap, other, entropy, GC.

        Only meaningful in nt mode; for aa mode we currently skip.
        """
        if not self.sequences:
            raise RuntimeError("Alignment not loaded.")
        if self.alignment_length is None:
            raise RuntimeError("Alignment length not set.")

        L = self.alignment_length
        bases = ["A", "C", "G", "T"]
        rows = []

        print("[aProfiler] Computing per-site nt metrics (counts, entropy, GC)...")
        for pos in tqdm(range(L), desc="Sites", unit="site"):
            col = [seq[pos] for seq in self.sequences]
            counts = Counter(col)

            row = {"position": pos + 1}
            total_used = 0

            for b in bases:
                c = counts.pop(b, 0)
                row[b] = c
                total_used += c

            gap_c = counts.pop("-", 0)
            row["gap"] = gap_c
            if include_gaps:
                total_used += gap_c

            other_c = sum(counts.values())
            row["other"] = other_c
            total_used += other_c

            # entropy over A/C/G/T only
            probs = []
            for b in bases:
                p = (row[b] / total_used) if total_used else 0.0
                probs.append(p)
            row["entropy_bits"] = entropy(probs, base=2) if any(probs) else 0.0

            # GC fraction (among A/C/G/T only)
            atgc = row["A"] + row["C"] + row["G"] + row["T"]
            gc = row["G"] + row["C"]
            row["gc_fraction"] = (gc / atgc) if atgc else 0.0

            rows.append(row)

        df = pd.DataFrame(rows)
        return df

    def save_per_site_nt_csv(self) -> Path:
        """Save per-site nt metrics to CSV."""
        df = self.per_site_nt_metrics(include_gaps=True)
        out_file = self.out_dir / f"{self.path.stem}_per_site_nt_metrics.csv"
        df.to_csv(out_file, index=False)
        return out_file

    def plot_entropy_track(self) -> Optional[Path]:
        """
        Plot per-site entropy (in bits) across the alignment.

        Expects a per-site NT metrics dataframe, but tolerates:
        - empty dataframes
        - missing or differently-named columns (e.g. `site` or `entropy`)
        """

        df = self.per_site_nt_metrics()

        # If completely empty, skip
        if df is None or len(df) == 0:
            print("[aProfiler] Warning: per-site NT metrics empty; skipping entropy plot.")
            return None

        # Ensure we have a positional axis
        if "position" not in df.columns:
            if "site" in df.columns:
                df = df.rename(columns={"site": "position"})
            else:
                # Fallback: synthesize 1..N
                df["position"] = np.arange(1, len(df) + 1)

        # Ensure we have an entropy column named entropy_bits
        if "entropy_bits" not in df.columns:
            if "entropy" in df.columns:
                df = df.rename(columns={"entropy": "entropy_bits"})
            else:
                print(
                    "[aProfiler] Warning: no entropy column found in per-site NT metrics; "
                    "skipping entropy plot."
                )
                return None

        fig, ax = plt.subplots(figsize=(10, 3))
        sns.lineplot(data=df, x="position", y="entropy_bits", ax=ax)
        ax.set_xlabel("Position")
        ax.set_ylabel("Entropy (bits)")
        ax.set_title(f"Per-site entropy – {self.path.stem}")
        fig.tight_layout()

        out_path = self.out_dir / f"{self.path.stem}_entropy.png"
        fig.savefig(out_path, dpi=300)
        plt.close(fig)
        return out_path


    # -------------------------------------------------------------------------
    # SEQUENCE LOGO (nt mode via logomaker)
    # -------------------------------------------------------------------------
    def plot_sequence_logo(self, max_len: int | None = 200) -> Path:
        """Sequence logo from per-site nt counts (first max_len sites)."""
        df = self.per_site_nt_metrics(include_gaps=False)
        bases = ["A", "C", "G", "T"]
        counts_df = df[bases].copy()

        if max_len is not None:
            counts_df = counts_df.iloc[:max_len, :]

        # probabilities per column
        probs_df = counts_df.div(counts_df.sum(axis=1), axis=0).fillna(0.0)

        fig, ax = plt.subplots(figsize=(10, 3))
        lm.Logo(probs_df, ax=ax)
        ax.set_xlabel("Position")
        ax.set_ylabel("Relative height")
        ax.set_title(f"Sequence logo (first {len(probs_df)} sites) – {self.path.stem}")
        fig.tight_layout()

        out_file = self.out_dir / f"{self.path.stem}_logo.png"
        fig.savefig(out_file, dpi=300)
        plt.close(fig)
        return out_file

    # -------------------------------------------------------------------------
    # CODON MODE
    # -------------------------------------------------------------------------
    def _codon_matrix(self) -> np.ndarray:
        """
        Return a (n_seq, n_codon) array of codon strings.

        If length not multiple of 3, truncate trailing chars.
        """
        if not self.sequences:
            raise RuntimeError("Alignment not loaded.")
        if self.alignment_length is None:
            raise RuntimeError("Alignment length not set.")

        L = self.alignment_length
        if L % 3 != 0:
            trim_to = (L // 3) * 3
            print(
                f"[aProfiler] Codon mode: length {L} not multiple of 3; truncating to {trim_to}."
            )
        else:
            trim_to = L

        n_codons = trim_to // 3
        codon_rows = []
        for seq in self.sequences:
            s = seq[:trim_to]
            codons = [s[i : i + 3] for i in range(0, trim_to, 3)]
            codon_rows.append(codons)

        return np.array(codon_rows, dtype="U3")

    def codon_global_frequencies(self, include_gaps: bool = True) -> pd.DataFrame:
        """
        Global codon usage frequencies.
        """
        codon_mat = self._codon_matrix()
        flat = codon_mat.ravel()
        counts = Counter(flat)

        rows = []
        total = sum(counts.values())
        for codon, c in counts.items():
            if not include_gaps and codon == "---":
                continue
            freq = c / total if total else 0.0
            rows.append({"codon": codon, "count": c, "frequency": freq})

        df = pd.DataFrame(rows).sort_values("count", ascending=False)
        return df

    def save_codon_global_csv(self) -> Path:
        df = self.codon_global_frequencies(include_gaps=True)
        out_file = self.out_dir / f"{self.path.stem}_codon_global.csv"
        df.to_csv(out_file, index=False)
        return out_file

    def codon_per_site_freqs(self) -> pd.DataFrame:
        """
        Per-codon-site codon frequencies + entropy.
        """
        codon_mat = self._codon_matrix()
        n_seq, n_codons = codon_mat.shape

        records = []
        print("[aProfiler] Computing per-codon-site frequencies and entropy...")
        for pos in tqdm(range(n_codons), desc="Codon sites", unit="codon-site"):
            col = codon_mat[:, pos]
            counts = Counter(col)
            total = sum(counts.values())
            probs = [c / total for c in counts.values()] if total else []
            H = entropy(probs, base=2) if probs else 0.0

            for codon, c in counts.items():
                freq = c / total if total else 0.0
                records.append(
                    {
                        "position_codon": pos + 1,
                        "codon": codon,
                        "count": c,
                        "frequency": freq,
                        "entropy_bits": H,
                    }
                )

        return pd.DataFrame(records)

    def save_codon_per_site_csv(self) -> Path:
        df = self.codon_per_site_freqs()
        out_file = self.out_dir / f"{self.path.stem}_codon_per_site.csv"
        df.to_csv(out_file, index=False)
        return out_file

    def plot_codon_usage(self, top_n: int = 30) -> Path:
        """
        Barplot of top-N codons by frequency.
        """
        df = self.codon_global_frequencies(include_gaps=True)
        df_top = df.head(top_n)

        fig, ax = plt.subplots(figsize=(10, 4))
        sns.barplot(data=df_top, x="codon", y="frequency", ax=ax)
        ax.set_xlabel("Codon")
        ax.set_ylabel("Frequency")
        ax.set_title(f"Global codon usage (top {top_n}) – {self.path.stem}")
        plt.xticks(rotation=90)
        fig.tight_layout()

        out_file = self.out_dir / f"{self.path.stem}_codon_usage_top{top_n}.png"
        fig.savefig(out_file, dpi=300)
        plt.close(fig)
        return out_file
    
    def codon_rscu_table(self) -> pd.DataFrame:
        """
        Compute Relative Synonymous Codon Usage (RSCU) for codons in the alignment.
    
        RSCU_i = (observed count of codon i for aa a) /
                 (expected count if all synonymous codons for a were used equally)
    
        Stops (*) and codons with non-ACGT chars are ignored.
        """
        # Get codon matrix & flatten
        codon_mat = self._codon_matrix()
        flat = codon_mat.ravel()
    
        # Count raw codons
        counts = Counter(flat)
    
        # Build synonyms table by AA (excluding stops)
        synonyms_by_aa = defaultdict(list)
        for codon, aa in CODON_TABLE.items():
            if aa == "*":
                continue
            synonyms_by_aa[aa].append(codon)
    
        # Build AA-level totals from observed codon counts
        aa_totals = defaultdict(int)
        for codon, count in counts.items():
            if len(codon) != 3:
                continue
            if any(b not in "ACGT" for b in codon):
                continue
            aa = CODON_TABLE.get(codon)
            if aa is None or aa == "*":
                continue
            aa_totals[aa] += count
    
        rows = []
        for codon, count in counts.items():
            # Only consider "clean" codons for RSCU
            if len(codon) != 3:
                continue
            if any(b not in "ACGT" for b in codon):
                continue
    
            aa = CODON_TABLE.get(codon)
            if aa is None or aa == "*":
                continue
    
            syn_codons = synonyms_by_aa.get(aa, [])
            if not syn_codons:
                continue
    
            n_syn = len(syn_codons)
            aa_total = aa_totals.get(aa, 0)
    
            if aa_total == 0:
                rscu = np.nan
            else:
                expected = aa_total / n_syn
                rscu = count / expected if expected > 0 else np.nan
    
            rows.append(
                {
                    "codon": codon,
                    "aa": aa,
                    "count": count,
                    "aa_total": aa_total,
                    "n_syn": n_syn,
                    "rscu": rscu,
                }
            )
    
        df = pd.DataFrame(rows)
        if not df.empty:
            df = df.sort_values(["aa", "codon"])
        return df

    def save_codon_rscu_csv(self) -> Path:
        """
        Save RSCU table to CSV.
        """
        df = self.codon_rscu_table()
        out_file = self.out_dir / f"{self.path.stem}_codon_rscu.csv"
        df.to_csv(out_file, index=False)
        return out_file

    def plot_rscu_heatmap(self) -> Path:
        """
        Heatmap of RSCU values (AA x codon). Codons on x, AAs on y.
        """
        df = self.codon_rscu_table()
        if df.empty:
            raise RuntimeError("No RSCU values to plot (empty codon table).")

        # pivot: rows=aa, columns=codon, values=rscu
        mat = df.pivot(index="aa", columns="codon", values="rscu")

        fig, ax = plt.subplots(figsize=(12, max(4, 0.4 * mat.shape[0])))
        sns.heatmap(mat, ax=ax, cmap="viridis", annot=False)
        ax.set_xlabel("Codon")
        ax.set_ylabel("Amino acid")
        ax.set_title(f"RSCU heatmap – {self.path.stem}")
        fig.tight_layout()

        out_file = self.out_dir / f"{self.path.stem}_codon_rscu_heatmap.png"
        fig.savefig(out_file, dpi=300)
        plt.close(fig)
        return out_file

    def aa_usage_from_codons(self) -> pd.DataFrame:
        """
        Aggregate codon counts to amino acid usage (from codon alignment).

        Returns columns:
          - aa
          - count
          - frequency
          - essential (True/False)
        """
        codon_mat = self._codon_matrix()
        flat = codon_mat.ravel()

        aa_counts = defaultdict(int)
        total_codons = 0

        for codon in flat:
            if len(codon) != 3:
                continue
            if any(b not in "ACGT" for b in codon):
                continue

            aa = CODON_TABLE.get(codon)
            if aa is None or aa == "*" or aa == "-":
                continue

            aa_counts[aa] += 1
            total_codons += 1

        rows = []
        for aa, c in aa_counts.items():
            freq = c / total_codons if total_codons else 0.0
            rows.append(
                {
                    "aa": aa,
                    "count": c,
                    "frequency": freq,
                    "essential": aa in ESSENTIAL_AA,
                }
            )

        df = pd.DataFrame(rows)
        if not df.empty:
            df = df.sort_values("aa")
        return df

    def save_aa_usage_from_codons_csv(self) -> Path:
        """
        Save amino acid usage (derived from codons) to CSV.
        """
        df = self.aa_usage_from_codons()
        out_file = self.out_dir / f"{self.path.stem}_aa_usage_from_codons.csv"
        df.to_csv(out_file, index=False)
        return out_file

    def save_essential_aa_summary_csv(self) -> Path:
        """
        Save essential vs non-essential amino acid summary (counts + frequencies).
        """
        df = self.aa_usage_from_codons()
        if df.empty:
            raise RuntimeError("No amino acid usage data to summarize.")

        total = df["count"].sum()
        df["class"] = df["essential"].map({True: "essential", False: "non_essential"})

        summary = (
            df.groupby("class")["count"]
            .sum()
            .reset_index()
            .rename(columns={"count": "total_count"})
        )
        summary["fraction"] = summary["total_count"] / total if total else 0.0

        out_file = self.out_dir / f"{self.path.stem}_aa_essential_summary.csv"
        summary.to_csv(out_file, index=False)
        return out_file

    def plot_essential_aa_bar(self) -> Path:
        """
        Barplot of total codon-derived usage for essential vs non-essential AAs.
        """
        df = self.aa_usage_from_codons()
        if df.empty:
            raise RuntimeError("No amino acid usage data to plot.")

        df["class"] = df["essential"].map({True: "essential", False: "non_essential"})
        summary = df.groupby("class")["count"].sum().reset_index()

        fig, ax = plt.subplots(figsize=(4, 4))
        sns.barplot(data=summary, x="class", y="count", ax=ax)
        ax.set_xlabel("AA class")
        ax.set_ylabel("Total codon count")
        ax.set_title(f"Essential vs non-essential AA usage – {self.path.stem}")
        fig.tight_layout()

        out_file = self.out_dir / f"{self.path.stem}_aa_essential_bar.png"
        fig.savefig(out_file, dpi=300)
        plt.close(fig)
        return out_file


    # -------------------------------------------------------------------------
    # SEQUENCE-LEVEL FEATURES + PCA / UMAP EMBEDDINGS
    # -------------------------------------------------------------------------
    def _sequence_level_features(self) -> pd.DataFrame:
        """
        Simple per-sequence features:
        - fraction A/C/G/T
        - gap fraction
        """
        if not self.sequences:
            raise RuntimeError("Alignment not loaded.")

        rows = []
        bases = ["A", "C", "G", "T"]
        for header, seq in zip(self.headers, self.sequences):
            counts = Counter(seq)
            L = len(seq)
            row = {"id": header}
            for b in bases:
                row[f"frac_{b}"] = counts.get(b, 0) / L if L else 0.0
            row["frac_gap"] = counts.get("-", 0) / L if L else 0.0
            rows.append(row)
        return pd.DataFrame(rows)

    def sequence_embedding_pca(self, n_components: int = 2) -> pd.DataFrame:
        df_feat = self._sequence_level_features()
        feature_cols = [c for c in df_feat.columns if c.startswith("frac_")]
        X = df_feat[feature_cols].values
        X_scaled = StandardScaler().fit_transform(X)

        pca = PCA(n_components=n_components)
        emb = pca.fit_transform(X_scaled)
        for i in range(n_components):
            df_feat[f"PC{i+1}"] = emb[:, i]
        return df_feat

    def sequence_embedding_umap(self, n_components: int = 2) -> pd.DataFrame:
        df_feat = self._sequence_level_features()
        feature_cols = [c for c in df_feat.columns if c.startswith("frac_")]
        X = df_feat[feature_cols].values
        X_scaled = StandardScaler().fit_transform(X)

        umap_model = UMAP(n_components=n_components, random_state=42)
        emb = umap_model.fit_transform(X_scaled)
        for i in range(n_components):
            df_feat[f"UMAP{i+1}"] = emb[:, i]
        return df_feat

    def save_embedding_pca_csv(self) -> Path:
        df = self.sequence_embedding_pca(n_components=2)
        out_file = self.out_dir / f"{self.path.stem}_pca_embedding.csv"
        df.to_csv(out_file, index=False)
        return out_file

    def save_embedding_umap_csv(self) -> Path:
        df = self.sequence_embedding_umap(n_components=2)
        out_file = self.out_dir / f"{self.path.stem}_umap_embedding.csv"
        df.to_csv(out_file, index=False)
        return out_file

    def plot_embedding_pca_plotnine(self) -> Path:
        """
        PCA embedding scatter plot using plotnine.
        """
        df_emb = self.sequence_embedding_pca(n_components=2)
        p = (
            ggplot(df_emb, aes(x="PC1", y="PC2"))
            + geom_point(alpha=0.7)
            + theme_minimal()
            + labs(
                title=f"PCA embedding – {self.path.stem}",
                x="PC1",
                y="PC2",
            )
        )

        out_file = self.out_dir / f"{self.path.stem}_pca_embedding.png"
        p.save(out_file, dpi=300)
        return out_file

    def plot_gc_track(self) -> Path:
        """
        Line plot of per-site GC fraction (nt mode).
        """
        if self.mode != "nt":
            raise RuntimeError("GC track is only meaningful in nt mode.")
    
        df = self.per_site_nt_metrics(include_gaps=True)
    
        fig, ax = plt.subplots(figsize=(9, 3))
        sns.lineplot(data=df, x="position", y="gc_fraction", ax=ax)
        ax.set_xlabel("Position")
        ax.set_ylabel("GC fraction")
        ax.set_title(f"Per-site GC fraction – {self.path.stem}")
        fig.tight_layout()
    
        out_file = self.out_dir / f"{self.path.stem}_gc_track.png"
        fig.savefig(out_file, dpi=300)
        plt.close(fig)
        return out_file

    def plot_gap_track(self) -> Path:
        """
        Line plot of per-site gap fraction (nt mode).
        """
        if self.mode != "nt":
            raise RuntimeError("Gap track is currently implemented for nt mode.")

        df = self.per_site_nt_metrics(include_gaps=True)

        # Compute gap fraction across all characters at that column
        total_counts = df[["A", "C", "G", "T", "gap", "other"]].sum(axis=1)
        gap_frac = df["gap"] / total_counts.replace(0, np.nan)
        df_plot = df.copy()
        df_plot["gap_fraction"] = gap_frac.fillna(0.0)

        fig, ax = plt.subplots(figsize=(9, 3))
        sns.lineplot(data=df_plot, x="position", y="gap_fraction", ax=ax)
        ax.set_xlabel("Position")
        ax.set_ylabel("Gap fraction")
        ax.set_title(f"Per-site gap fraction – {self.path.stem}")
        fig.tight_layout()

        out_file = self.out_dir / f"{self.path.stem}_gap_track.png"
        fig.savefig(out_file, dpi=300)
        plt.close(fig)
        return out_file

    def plot_nt_heatmap(self, max_len: int = 200) -> Path:
        """
        Heatmap of per-site A/C/G/T frequency across the alignment
        (first max_len sites to keep it readable).
        """
        if self.mode != "nt":
            raise RuntimeError("Nucleotide heatmap is only meaningful in nt mode.")

        df = self.per_site_nt_metrics(include_gaps=False)
        bases = ["A", "C", "G", "T"]

        if max_len is not None:
            df = df.iloc[:max_len, :]

        # convert counts to frequencies per site
        counts = df[bases]
        row_sums = counts.sum(axis=1)
        freq_df = counts.div(row_sums.replace(0, np.nan), axis=0).fillna(0.0)
        freq_df.index = df["position"]

        # Heatmap: rows = bases, columns = positions
        mat = freq_df.T  # shape (4, positions)

        fig, ax = plt.subplots(figsize=(min(12, 0.06 * mat.shape[1] + 2), 4))
        sns.heatmap(mat, ax=ax, cmap="viridis")  # default colormap
        ax.set_xlabel("Position")
        ax.set_ylabel("Base")
        ax.set_title(f"A/C/G/T frequency heatmap (first {mat.shape[1]} sites) – {self.path.stem}")
        fig.tight_layout()

        out_file = self.out_dir / f"{self.path.stem}_nt_heatmap.png"
        fig.savefig(out_file, dpi=300)
        plt.close(fig)
        return out_file

    def plot_sequence_gap_fraction_bar(self) -> Path:
        """
        Barplot of per-sequence gap fraction.
        """
        if not self.sequences:
            raise RuntimeError("Alignment not loaded.")

        rows = []
        for header, seq in zip(self.headers, self.sequences):
            L = len(seq)
            gap_frac = seq.count("-") / L if L else 0.0
            rows.append({"id": header, "gap_fraction": gap_frac})

        df = pd.DataFrame(rows)

        fig, ax = plt.subplots(figsize=(max(6, 0.15 * len(df)), 4))
        sns.barplot(data=df, x="id", y="gap_fraction", ax=ax)
        ax.set_xlabel("Sequence")
        ax.set_ylabel("Gap fraction")
        ax.set_title(f"Per-sequence gap fraction – {self.path.stem}")
        ax.tick_params(axis="x", rotation=90)
        fig.tight_layout()

        out_file = self.out_dir / f"{self.path.stem}_seq_gap_fraction.png"
        fig.savefig(out_file, dpi=300)
        plt.close(fig)
        return out_file

    def plot_pairwise_identity_histogram(self, max_pairs: int = 5000) -> Path:
        """
        Histogram of pairwise sequence identity (all pairs or a random subset if many).
        Identity is computed over the full alignment length, including gaps.
        """
        if not self.sequences:
            raise RuntimeError("Alignment not loaded.")
        if self.alignment_length is None:
            raise RuntimeError("Alignment length not set.")

        n = len(self.sequences)
        if n < 2:
            raise RuntimeError("Need at least 2 sequences for pairwise identity.")

        seq_arr = self.to_array()  # shape: (n_seq, L)

        # All unique pairs
        all_pairs = list(itertools.combinations(range(n), 2))
        if len(all_pairs) > max_pairs:
            all_pairs = random.sample(all_pairs, max_pairs)

        id_vals = []
        for i, j in all_pairs:
            s1 = seq_arr[i]
            s2 = seq_arr[j]
            matches = (s1 == s2)
            ident = matches.mean()  # fraction identity
            id_vals.append(ident)

        df = pd.DataFrame({"pairwise_identity": id_vals})

        fig, ax = plt.subplots(figsize=(6, 4))
        sns.histplot(df["pairwise_identity"], bins=30, kde=True, ax=ax)
        ax.set_xlabel("Pairwise identity")
        ax.set_ylabel("Count")
        ax.set_title(f"Pairwise identity distribution (n_pairs={len(all_pairs)}) – {self.path.stem}")
        fig.tight_layout()

        out_file = self.out_dir / f"{self.path.stem}_pairwise_identity_hist.png"
        fig.savefig(out_file, dpi=300)
        plt.close(fig)
        return out_file

    def plot_embedding_umap_plotnine(self) -> Path:
        """
        UMAP embedding scatter plot using plotnine.
        Requires that UMAP can be computed from simple sequence-level features.
        """
        df_emb = self.sequence_embedding_umap(n_components=2)

        p = (
            ggplot(df_emb, aes(x="UMAP1", y="UMAP2"))
            + geom_point(alpha=0.7)
            + theme_minimal()
            + labs(
                title=f"UMAP embedding – {self.path.stem}",
                x="UMAP1",
                y="UMAP2",
            )
        )

        out_file = self.out_dir / f"{self.path.stem}_umap_embedding.png"
        p.save(out_file, dpi=300)
        return out_file


    # -------------------------------------------------------------------------
    # Master run
    # -------------------------------------------------------------------------
    def generate_report(
        self,
        outputs: Dict[str, Path],
        fmt: str = "md",
        title: str | None = None,
    ) -> Path:
        """
        Generate a simple report (Markdown or HTML) summarizing all outputs.
    
        Parameters
        ----------
        outputs : dict
            Mapping from logical keys to filesystem Paths (from run_full_profile).
        fmt : {"md", "html"}
            Report format.
        title : str, optional
            Title for the report; defaults to file stem.
    
        Returns
        -------
        Path
            Path to the generated report file.
        """
        title = title or f"aProfiler report – {self.path.stem}"
    
        # Separate tables (csv) from plots (images)
        table_items = []
        plot_items = []
        other_items = []
    
        for key, path in outputs.items():
            suffix = path.suffix.lower()
            if suffix == ".csv":
                table_items.append((key, path))
            elif suffix in {".png", ".jpg", ".jpeg", ".svg"}:
                plot_items.append((key, path))
            else:
                other_items.append((key, path))
    
        # Use relative paths in report
        def rel(p: Path) -> str:
            return str(p.relative_to(self.out_dir))
    
        if fmt == "md":
            lines: list[str] = []
            lines.append(f"# {title}")
            lines.append("")
            lines.append(f"**Input file:** `{self.path}`  ")
            lines.append(f"**Mode:** `{self.mode}`  ")
            lines.append("")
            lines.append("## Tables")
            lines.append("")
            if table_items:
                for key, p in table_items:
                    lines.append(f"- **{key}** – `{rel(p)}`")
            else:
                lines.append("_No tables generated._")
            lines.append("")
    
            lines.append("## Plots")
            lines.append("")
            if plot_items:
                for key, p in plot_items:
                    lines.append(f"### {key}")
                    lines.append("")
                    lines.append(f"![{key}]({rel(p)})")
                    lines.append("")
            else:
                lines.append("_No plots generated._")
                lines.append("")
    
            if other_items:
                lines.append("## Other outputs")
                lines.append("")
                for key, p in other_items:
                    lines.append(f"- **{key}** – `{rel(p)}`")
    
            report_path = self.out_dir / f"{self.path.stem}_{self.mode}_report.md"
            report_path.write_text("\n".join(lines), encoding="utf-8")
            return report_path
    
        elif fmt == "html":
            from html import escape
    
            def esc(x: str) -> str:
                return escape(x, quote=True)
    
            html_lines: list[str] = []
            html_lines.append("<!DOCTYPE html>")
            html_lines.append("<html><head><meta charset='utf-8'>")
            html_lines.append(f"<title>{esc(title)}</title>")
            html_lines.append("</head><body>")
            html_lines.append(f"<h1>{esc(title)}</h1>")
            html_lines.append(f"<p><b>Input file:</b> {esc(str(self.path))}<br>")
            html_lines.append(f"<b>Mode:</b> {esc(self.mode)}</p>")
    
            html_lines.append("<h2>Tables</h2><ul>")
            if table_items:
                for key, p in table_items:
                    html_lines.append(
                        f"<li><b>{esc(key)}</b> – <code>{esc(rel(p))}</code></li>"
                    )
            else:
                html_lines.append("<li><i>No tables generated.</i></li>")
            html_lines.append("</ul>")
    
            html_lines.append("<h2>Plots</h2>")
            if plot_items:
                for key, p in plot_items:
                    html_lines.append(f"<h3>{esc(key)}</h3>")
                    html_lines.append(
                        f"<img src='{esc(rel(p))}' alt='{esc(key)}' style='max-width:100%;'>"
                    )
            else:
                html_lines.append("<p><i>No plots generated.</i></p>")
    
            if other_items:
                html_lines.append("<h2>Other outputs</h2><ul>")
                for key, p in other_items:
                    html_lines.append(
                        f"<li><b>{esc(key)}</b> – <code>{esc(rel(p))}</code></li>"
                    )
                html_lines.append("</ul>")
    
            html_lines.append("</body></html>")
    
            report_path = self.out_dir / f"{self.path.stem}_{self.mode}_report.html"
            report_path.write_text("\n".join(html_lines), encoding="utf-8")
            return report_path
    
        else:
            raise ValueError(f"Unsupported report format: {fmt}")

    
    def run_full_profile(self, make_plots: bool = True) -> Dict[str, Path]:
        """
        High-level driver:
          - For 'nt' and 'aa': global frequencies, per-site nt metrics (nt only), logo, embeddings
          - For 'codon': codon usage + per-codon-site metrics
          - Plots ON by default.
        """
        outputs: Dict[str, Path] = {}

        if self.mode == "codon":
            # Codon-aware pipeline
            outputs["codon_global_csv"] = self.save_codon_global_csv()
            outputs["codon_per_site_csv"] = self.save_codon_per_site_csv()

            # NEW: RSCU + AA usage + essential vs non-essential summary
            outputs["codon_rscu_csv"] = self.save_codon_rscu_csv()
            outputs["aa_usage_from_codons_csv"] = self.save_aa_usage_from_codons_csv()
            outputs["aa_essential_summary_csv"] = self.save_essential_aa_summary_csv()

            if make_plots:
                outputs["codon_usage_plot"] = self.plot_codon_usage(top_n=30)
                outputs["codon_rscu_heatmap_plot"] = self.plot_rscu_heatmap()
                outputs["aa_essential_bar_plot"] = self.plot_essential_aa_bar()

        else:
            # nt / aa pipeline
            outputs["global_freqs_csv"] = self.save_global_freqs_csv()

            if self.mode == "nt":
                # Tables
                outputs["per_site_nt_csv"] = self.save_per_site_nt_csv()

                if make_plots:
                    # Existing plots
                    outputs["global_freqs_plot"] = self.plot_global_freqs()
                    entropy_path = self.plot_entropy_track()
                    if entropy_path is not None:
                        outputs["entropy_plot"] = entropy_path

                    outputs["logo_plot"] = self.plot_sequence_logo()

                    # NEW nt-specific plots
                    outputs["gc_track_plot"] = self.plot_gc_track()
                    outputs["gap_track_plot"] = self.plot_gap_track()
                    outputs["nt_heatmap_plot"] = self.plot_nt_heatmap(max_len=200)

            else:
                # aa: still global symbol freqs (and below, embeddings)
                if make_plots:
                    outputs["global_freqs_plot"] = self.plot_global_freqs()

            # Embeddings (mode-agnostic)
            outputs["pca_embedding_csv"] = self.save_embedding_pca_csv()
            outputs["umap_embedding_csv"] = self.save_embedding_umap_csv()
            if make_plots:
                outputs["pca_embedding_plot"] = self.plot_embedding_pca_plotnine()
                outputs["umap_embedding_plot"] = self.plot_embedding_umap_plotnine()

            # Sequence-level structure plots (all modes)
            if make_plots:
                outputs["seq_gap_fraction_plot"] = self.plot_sequence_gap_fraction_bar()
                outputs["pid_hist_plot"] = self.plot_pairwise_identity_histogram()

        return outputs

    


def run_profile(
    input_path: str,
    mode: str = "auto",
    make_plots: bool = True,
    make_report: bool = False,
    report_format: str = "md",
    max_seqs: Optional[int] = None,
    max_sites: Optional[int] = None,
    seed: Optional[int] = None,
    embeddings: bool = True,
    make_summary_card: bool = False,
) -> None:
    """
    CLI-facing helper for running a full profile on an alignment.

    Parameters
    ----------
    input_path : str
        Path to the input MSA file.
    mode : {"nt","aa","codon","auto"}
        Alignment mode.
    make_plots : bool
        Whether to generate plot files (PNG/SVG) as well as CSV tables.
    make_report : bool
        Whether to generate a summary report file.
    report_format : {"md","html"}
        Format for report if make_report is True.
    max_seqs : int, optional
        Max number of sequences to use for heavy visualizations/embeddings.
    max_sites : int, optional
        Reserved for future: max number of sites (columns) for heavy plots.
    seed : int, optional
        Random seed for subsampling and embedding reproducibility.
    embeddings : bool
        If False, skip PCA/UMAP computations.
    make_summary_card : bool
        If True, generate a condensed summary card Markdown file.
    """
    profiler = AlignmentProfiler(
        input_path=input_path,
        mode=mode,
        max_seqs=max_seqs,
        max_sites=max_sites,
        seed=seed,
        embeddings=embeddings,
    )
    profiler.load_alignment()
    outputs: Dict[str, Path] = profiler.run_full_profile(make_plots=make_plots)

    print("[aProfiler] Profiling complete. Outputs:")
    for key, path in outputs.items():
        print(f"  • {key}: {path}")

    if make_report:
        if hasattr(profiler, "generate_report"):
            report_path = profiler.generate_report(outputs, fmt=report_format)
            print(f"[aProfiler] Report generated: {report_path}")
        else:
            print(
                "[aProfiler] Warning: --report requested, "
                "but generate_report(...) is not implemented in this version."
            )

    if make_summary_card:
        if hasattr(profiler, "generate_summary_card"):
            card_path = profiler.generate_summary_card(outputs, fmt="md")
            print(f"[aProfiler] Summary card generated: {card_path}")
        else:
            print(
                "[aProfiler] Warning: --summary-card requested, "
                "but generate_summary_card(...) is not implemented in this version."
            )
            

# =============================================================================
# End of file
# =============================================================================
