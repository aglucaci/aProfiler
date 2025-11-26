# alignmentprofiler/cli.py

import argparse
import json
import csv
from .core import AlignmentProfiler
from . import plots
#from alignmentprofiler.codon_plot import plot_codon_diversity

def main():
    parser = argparse.ArgumentParser(description="AlignmentProfiler: summarize a multiple sequence alignment.")
    parser.add_argument("input", help="Path to alignment file (FASTA, PHYLIP, etc.)")
    parser.add_argument("--format", default="fasta", help="Alignment format (default: fasta)")
    parser.add_argument("--output", default=None, help="Optional output JSON file")
    parser.add_argument("--csv", default=None, help="Optional output CSV file for sequence-level stats")
    parser.add_argument("--plot", choices=["entropy", "gaps", "codon_diversity", "protein_diversity"], help="Optional plot type to display or save")
    parser.add_argument("--plot-out", default=None, help="Path to save the plot instead of displaying it")
    parser.add_argument("--rscu-csv", default=None, help="Optional output CSV file for RSCU_")
    parser.add_argument("--aa-csv", default=None, help="Optional CSV file to export amino acid frequencies")
    args = parser.parse_args()

    profiler = AlignmentProfiler(args.input, args.format)
    summary = profiler.summarize()
    summary["nucleotide"] = profiler.nucleotide_summary()
    summary["codon"] = profiler.codon_summary()

    # JSON (Standard output) --------------------------------------------------
    
    if args.output:
        with open(args.output, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"[+] Summary written to {args.output}")
    else:
        print(json.dumps(summary, indent=2))
    # end if

    # TABLES (CSV-Formatted output) -------------------------------------------
    
    if args.csv:
        with open(args.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=summary["sequences"][0].keys())
            writer.writeheader()
            writer.writerows(summary["sequences"])
        print(f"[+] CSV summary written to {args.csv}")
    # end if
    
    if args.rscu_csv:
        with open(args.rscu_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Codon", "RSCU"])
            for codon, score in summary["codon"]["rscu"].items():
                writer.writerow([codon, score])
        print(f"[+] RSCU table written to {args.rscu_csv}")
    # end if
    
    if args.aa_csv:
        with open(args.aa_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Amino Acid", "Frequency"])
            for aa, freq in summary["protein"]["frequencies"].items():
                writer.writerow([aa, freq])
        print(f"[+] Amino acid frequency table written to {args.aa_csv}")

    
    # ALL PLOTS ---------------------------------------------------------------
    if args.plot:
        if args.plot == "entropy":
            plots.plot_entropy(profiler.alignment, save_path=args.plot_out)
        elif args.plot == "gaps":
            plots.plot_gap_heatmap(profiler.alignment, save_path=args.plot_out)
        elif args.plot == "codon_diversity":
            plots.plot_codon_diversity(summary["codon"]["diversity"], save_path=args.plot_out)
        elif args.plot == "protein_diversity":
            plots.plot_protein_diversity(summary["protein"]["diversity"], save_path=args.plot_out)
        # end if
    # end if
    
# end method

if __name__ == "__main__":
    main()

# End of file
