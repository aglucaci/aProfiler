# -*- coding: utf-8 -*-
"""

@author: Alexander G. Lucaci
"""

import os
from pathlib import Path

from aprofiler.profiler import AlignmentProfiler, run_profile


def test_nt_profile_basic(tmp_path: Path):
    data_dir = Path(__file__).parent / "data"
    fasta = data_dir / "toy_nt.fasta"

    # run with plots off to make tests lighter
    profiler = AlignmentProfiler(str(fasta), mode="nt", out_dir=tmp_path)
    profiler.load_alignment()
    outputs = profiler.run_full_profile(make_plots=False)

    # basic checks
    assert "global_freqs_csv" in outputs
    assert outputs["global_freqs_csv"].exists()

    if "per_site_nt_csv" in outputs:
        assert outputs["per_site_nt_csv"].exists()

    # make sure CSV has some rows
    import pandas as pd
    df = pd.read_csv(outputs["global_freqs_csv"])
    assert not df.empty
