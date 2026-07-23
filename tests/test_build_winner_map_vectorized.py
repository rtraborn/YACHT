#!/usr/bin/env python
"""
Correctness tests for the vectorized build_winner_map().

The vectorized implementation must produce a winner_map identical to the original
serial reduction, including tie-breaking (earliest organism in row order wins on
equal ANI) and skipping of NaN-ANI / empty-sketch organisms.

Run with:
    pytest tests/test_build_winner_map_vectorized.py -v
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from yacht.hypothesis_recovery_src import build_winner_map


class _FakeMinHash:
    def __init__(self, hashes):
        # mirror sourmash: .hashes is a dict of {hash: abundance}
        self.hashes = {int(h): 1 for h in hashes}


class _FakeSketch:
    def __init__(self, hashes):
        self.minhash = _FakeMinHash(hashes)


def _serial_reference(final_stats_df):
    """Original serial winner-map reduction, kept here as ground truth."""
    winner_map = {}
    for idx in range(len(final_stats_df)):
        row = final_stats_df.iloc[idx]
        organism_name = row["organism_name"]
        ani = row["final_est_ani"]
        if pd.isna(ani):
            continue
        for kmer in row["genome_sketch"].minhash.hashes.keys():
            if kmer not in winner_map or ani > winner_map[kmer][0]:
                winner_map[kmer] = (ani, organism_name)
    return winner_map


def _make_df(rows):
    """rows: list of (organism_name, ani, [hashes])"""
    return pd.DataFrame(
        {
            "organism_name": [r[0] for r in rows],
            "final_est_ani": [r[1] for r in rows],
            "genome_sketch": [_FakeSketch(r[2]) for r in rows],
        }
    )


def _assert_maps_equal(df):
    expected = _serial_reference(df)
    got = build_winner_map(df, path_to_genome_temp_dir="", ksize=31)
    assert got == expected, f"\nexpected={expected}\ngot={got}"


def test_basic_disjoint():
    df = _make_df([
        ("org_A", 0.99, [10, 20, 30]),
        ("org_B", 0.97, [40, 50]),
    ])
    _assert_maps_equal(df)


def test_shared_kmer_higher_ani_wins():
    df = _make_df([
        ("org_A", 0.95, [10, 20]),
        ("org_B", 0.99, [20, 30]),  # k-mer 20 shared; org_B has higher ANI
    ])
    _assert_maps_equal(df)
    got = build_winner_map(df, "", 31)
    assert got[20][1] == "org_B"


def test_tie_earliest_row_wins():
    # Equal ANI on shared k-mer 20: the earlier row (org_A) must keep it,
    # matching the serial strict-'>' semantics.
    df = _make_df([
        ("org_A", 0.98, [10, 20]),
        ("org_B", 0.98, [20, 30]),
    ])
    _assert_maps_equal(df)
    got = build_winner_map(df, "", 31)
    assert got[20][1] == "org_A"


def test_nan_ani_skipped():
    df = _make_df([
        ("org_A", np.nan, [10, 20]),
        ("org_B", 0.96, [20, 30]),
    ])
    _assert_maps_equal(df)
    got = build_winner_map(df, "", 31)
    assert got[20][1] == "org_B"  # org_A skipped entirely


def test_empty_sketch_skipped():
    df = _make_df([
        ("org_A", 0.99, []),
        ("org_B", 0.96, [40, 50]),
    ])
    _assert_maps_equal(df)


def test_all_nan_returns_empty():
    df = _make_df([
        ("org_A", np.nan, [10, 20]),
        ("org_B", np.nan, [30]),
    ])
    assert build_winner_map(df, "", 31) == {}


def test_randomized_fuzz():
    rng = np.random.default_rng(1234)
    for _ in range(200):
        n_org = int(rng.integers(1, 12))
        rows = []
        for i in range(n_org):
            # small hash universe forces frequent sharing and ties
            n_h = int(rng.integers(0, 8))
            hashes = rng.integers(0, 25, size=n_h).tolist()
            # occasionally inject NaN ANI
            ani = np.nan if rng.random() < 0.15 else round(float(rng.integers(90, 100)) / 100, 2)
            rows.append((f"org_{i}", ani, hashes))
        _assert_maps_equal(_make_df(rows))
