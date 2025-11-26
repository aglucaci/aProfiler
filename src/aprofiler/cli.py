# -*- coding: utf-8 -*-
"""

@author: Alexander G. Lucaci
"""

import argparse
import sys
from .profiler import run_profile


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="aprofiler",
        description="aProfiler — MSA statistics & visualization.",
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input MSA file (FASTA/A3M/ALN/etc.)",
    )
    parser.add_argument(
        "--mode",
        choices=["nt", "aa", "codon", "auto"],
        default="auto",
        help="Sequence type mode (nt, aa, codon, or auto-detect).",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Disable plot generation (CSV only).",
    )
    parser.add_argument(
        "--report",
        action="store_true",
        help="Generate a summary report (HTML by default).",
    )
    parser.add_argument(
        "--report-format",
        choices=["md", "html"],
        default="html",
        help="Report format if --report is set (default: html).",
    )

    return parser
# end method


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = build_parser()
    if not argv:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args(argv)

    make_plots = not args.no_plots

    run_profile(
        input_path=args.input,
        mode=args.mode,
        make_plots=make_plots,
        make_report=args.report,
        report_format=args.report_format,
    )


if __name__ == "__main__":
    main()

