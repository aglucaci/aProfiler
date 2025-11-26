# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 13:50:25 2025

@author: Alexander G. Lucaci
"""

from pathlib import Path
from aprofiler.profiler import AlignmentProfiler


def test_codon_profile_rscu_and_essential(tmp_path: Path):
    data_dir = Path(__file__).parent / "data"
    fasta = data_dir / "toy_codon.fasta"

    profiler = AlignmentProfiler(str(fasta), mode="codon", out_dir=tmp_path)
    profiler.load_alignment()
    outputs = profiler.run_full_profile(make_plots=False)

    # core codon outputs
    assert "codon_global_csv" in outputs
    assert "codon_per_site_csv" in outputs
    assert "codon_rscu_csv" in outputs
    assert "aa_usage_from_codons_csv" in outputs
    assert "aa_essential_summary_csv" in outputs

    for key in [
        "codon_global_csv",
        "codon_per_site_csv",
        "codon_rscu_csv",
        "aa_usage_from_codons_csv",
        "aa_essential_summary_csv",
    ]:
        assert outputs[key].exists()
