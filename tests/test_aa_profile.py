# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 13:50:08 2025

@author: Alexander G. Lucaci
"""

from pathlib import Path
from aprofiler.profiler import AlignmentProfiler


def test_aa_profile_basics(tmp_path: Path):
    data_dir = Path(__file__).parent / "data"
    fasta = data_dir / "toy_aa.fasta"

    profiler = AlignmentProfiler(str(fasta), mode="aa", out_dir=tmp_path)
    profiler.load_alignment()
    outputs = profiler.run_full_profile(make_plots=False)

    assert "global_freqs_csv" in outputs
    assert outputs["global_freqs_csv"].exists()
