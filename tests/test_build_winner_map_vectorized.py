#!/usr/bin/env python
"""
Correctness tests for the vectorized winner-takes-all back-half (Phase 2).

The array-based implementation must reproduce the original serial reduction and the
original consumers bit-for-bit:
  * build_winner_map                -> won-k-mers-per-organism (was a {kmer:(ani,org)} dict)
  * recalculate_ani_from_winner_map -> final_est_ani / reassignment_status
  * estimate_relative_abundance     -> rel_abund / kmers_lost

Each test compares the new functions against inline serial references (the pre-Phase-2
logic), including ties, NaN ANIs, empty sketches, elimination, and both one-pass and
two-pass orderings.

Run with:
    pytest tests/test_build_winner_map_vectorized.py -v
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from yacht.hypothesis_recovery_src import (
    build_winner_map,
    recalculate_ani_from_winner_map,
    estimate_relative_abundance,
)
from yacht.utils import (
    ratio_lambda,
    ani_from_lambda,
    SAMPLE_SIZE_CUTOFF,
    MIN_COUNT_THRESH,
)

KSIZE = 31
SCALE = 1000
MIN_ANI = 0.90


# --------------------------------------------------------------------------- #
# Fake sourmash sketch (only the attributes the code touches)
# --------------------------------------------------------------------------- #
class _FakeMinHash:
    def __init__(self, hashes):
        self.hashes = {int(h): 1 for h in hashes}
        self.scaled = SCALE


class _FakeSketch:
    def __init__(self, hashes):
        self.minhash = _FakeMinHash(hashes)


def _make_df(rows):
    """rows: list of (organism_name, ani, [hashes])"""
    return pd.DataFrame(
        {
            "organism_name": [r[0] for r in rows],
            "final_est_ani": [r[1] for r in rows],
            "genome_sketch": [_FakeSketch(r[2]) for r in rows],
        }
    )


def _sample_arrays(sample_hashes):
    sk = np.fromiter(sample_hashes.keys(), dtype=np.uint64)
    sa = np.fromiter(sample_hashes.values(), dtype=np.int64)
    order = np.argsort(sk)
    return sk[order], sa[order]


# --------------------------------------------------------------------------- #
# Serial references = the original (pre-Phase-2) implementations
# --------------------------------------------------------------------------- #
def _serial_build(df):
    winner_map = {}
    for idx in range(len(df)):
        row = df.iloc[idx]
        ani = row["final_est_ani"]
        if pd.isna(ani):
            continue
        org = row["organism_name"]
        for kmer in row["genome_sketch"].minhash.hashes.keys():
            if kmer not in winner_map or ani > winner_map[kmer][0]:
                winner_map[kmer] = (ani, org)
    return winner_map


def _serial_recalc(df, winner_map, sample_hashes, ksize, min_ani):
    df["reassignment_status"] = "active"
    df["original_ani"] = df["final_est_ani"].copy()
    for idx in range(len(df)):
        row = df.iloc[idx]
        org = row["organism_name"]
        won_in_sample = []
        total_won = 0
        for kmer in row["genome_sketch"].minhash.hashes.keys():
            if kmer in winner_map and winner_map[kmer][1] == org:
                total_won += 1
                if kmer in sample_hashes and sample_hashes[kmer] > 0:
                    won_in_sample.append(sample_hashes[kmer])
        if total_won == 0:
            df.at[idx, "reassignment_status"] = "eliminated"
            df.at[idx, "final_est_ani"] = float("nan")
            continue
        if len(won_in_sample) < SAMPLE_SIZE_CUTOFF:
            naive = (len(won_in_sample) / total_won) ** (1 / ksize)
            if naive >= min_ani:
                df.at[idx, "final_est_ani"] = naive
            df.at[idx, "reassignment_status"] = "lambda_failed"
            continue
        full_cov = [0] * (total_won - len(won_in_sample)) + won_in_sample
        new_lambda = ratio_lambda(full_cov, MIN_COUNT_THRESH)
        if new_lambda is None:
            df.at[idx, "reassignment_status"] = "lambda_failed"
            continue
        mean_cov = sum(full_cov) / len(full_cov)
        new_ani = ani_from_lambda(new_lambda, mean_cov, ksize, full_cov)
        if new_ani is not None:
            df.at[idx, "final_est_ani"] = new_ani
        else:
            df.at[idx, "reassignment_status"] = "ani_failed"
    return df


def _serial_estimate(df, winner_map, sample_hashes):
    df["kmers_lost"] = 0
    df["rel_abund"] = 0.0
    has_status = "reassignment_status" in df.columns
    for idx in range(len(df)):
        row = df.iloc[idx]
        org = row["organism_name"]
        if has_status and row["reassignment_status"] == "eliminated":
            continue
        sig = row["genome_sketch"]
        lost = 0
        cov = 0.0
        for kmer in sig.minhash.hashes.keys():
            if kmer in winner_map:
                if winner_map[kmer][1] != org:
                    lost += 1
                elif kmer in sample_hashes:
                    cov += sample_hashes[kmer]
        df.at[idx, "kmers_lost"] = lost
        gsize = len(sig.minhash.hashes) * sig.minhash.scaled
        df.at[idx, "rel_abund"] = cov / gsize if gsize > 0 else 0.0
    total = df["rel_abund"].sum()
    if total > 0:
        df["rel_abund"] = df["rel_abund"] / total
    return df


# --------------------------------------------------------------------------- #
# Comparison helpers
# --------------------------------------------------------------------------- #
def _won_to_map(df, won_by_row):
    """Reconstruct the old {kmer:(ani,org)} dict from the new array output."""
    names = df["organism_name"].to_numpy()
    anis = df["final_est_ani"].to_numpy()
    m = {}
    for r, arr in won_by_row.items():
        for k in arr:
            m[int(k)] = (float(anis[r]), names[r])
    return m


def _assert_frames_equivalent(ds, dn):
    assert list(dn["reassignment_status"]) == list(ds["reassignment_status"])
    assert np.array_equal(dn["kmers_lost"].to_numpy(), ds["kmers_lost"].to_numpy())
    assert np.allclose(
        dn["final_est_ani"].to_numpy(dtype=float),
        ds["final_est_ani"].to_numpy(dtype=float),
        equal_nan=True,
    )
    assert np.allclose(
        dn["rel_abund"].to_numpy(dtype=float),
        ds["rel_abund"].to_numpy(dtype=float),
        equal_nan=True,
    )


# --------------------------------------------------------------------------- #
# build_winner_map: assignment equivalence
# --------------------------------------------------------------------------- #
def _assert_build_equal(df):
    expected = _serial_build(df.copy())
    won = build_winner_map(df.copy(), KSIZE)
    got = _won_to_map(df, won)
    assert got == expected, f"\nexpected={expected}\ngot={got}"


def test_build_disjoint():
    _assert_build_equal(_make_df([("A", 0.99, [10, 20, 30]), ("B", 0.97, [40, 50])]))


def test_build_higher_ani_wins():
    df = _make_df([("A", 0.95, [10, 20]), ("B", 0.99, [20, 30])])
    won = build_winner_map(df, KSIZE)
    assert 20 in won[1] and 20 not in won.get(0, np.array([]))


def test_build_tie_earliest_row_wins():
    df = _make_df([("A", 0.98, [10, 20]), ("B", 0.98, [20, 30])])
    _assert_build_equal(df)
    won = build_winner_map(df, KSIZE)
    assert 20 in won[0] and 20 not in won.get(1, np.array([]))


def test_build_nan_and_empty_skipped():
    _assert_build_equal(_make_df([("A", np.nan, [10, 20]), ("B", 0.96, [20, 30])]))
    _assert_build_equal(_make_df([("A", 0.99, []), ("B", 0.96, [40, 50])]))


def test_build_all_nan_returns_empty():
    df = _make_df([("A", np.nan, [10]), ("B", np.nan, [20])])
    assert build_winner_map(df, KSIZE) == {}


# --------------------------------------------------------------------------- #
# End-to-end: consumers vs serial references
# --------------------------------------------------------------------------- #
def _run_serial_two_pass(df, sample_hashes):
    wm = _serial_build(df)
    _serial_recalc(df, wm, sample_hashes, KSIZE, MIN_ANI)
    wm = _serial_build(df)
    _serial_estimate(df, wm, sample_hashes)
    return df


def _run_new_two_pass(df, sample_hashes):
    sk, sa = _sample_arrays(sample_hashes)
    won = build_winner_map(df, KSIZE)
    recalculate_ani_from_winner_map(df, won, sk, sa, KSIZE, min_ani=MIN_ANI, num_threads=2)
    won = build_winner_map(df, KSIZE)
    estimate_relative_abundance(df, won, sk, sa, SCALE)
    return df


def _run_serial_one_pass(df, sample_hashes):
    wm = _serial_build(df)
    df["reassignment_status"] = "one_pass"
    df["original_ani"] = df["final_est_ani"].copy()
    _serial_estimate(df, wm, sample_hashes)
    return df


def _run_new_one_pass(df, sample_hashes):
    sk, sa = _sample_arrays(sample_hashes)
    won = build_winner_map(df, KSIZE)
    df["reassignment_status"] = "one_pass"
    df["original_ani"] = df["final_est_ani"].copy()
    estimate_relative_abundance(df, won, sk, sa, SCALE)
    return df


def test_two_pass_small():
    rows = [
        ("A", 0.99, list(range(0, 40))),
        ("B", 0.97, list(range(30, 70))),
        ("C", np.nan, list(range(60, 65))),
        ("D", 0.93, [200, 201]),
    ]
    sample = {h: (h % 4) for h in range(0, 70)}  # some zero abundances
    df = _make_df(rows)
    ds = _run_serial_two_pass(df.copy(), dict(sample))
    dn = _run_new_two_pass(df.copy(), dict(sample))
    _assert_frames_equivalent(ds, dn)


def test_one_pass_small():
    rows = [
        ("A", 0.99, list(range(0, 40))),
        ("B", 0.97, list(range(30, 70))),
        ("C", np.nan, list(range(60, 65))),
    ]
    sample = {h: (h % 3) for h in range(0, 70)}
    df = _make_df(rows)
    ds = _run_serial_one_pass(df.copy(), dict(sample))
    dn = _run_new_one_pass(df.copy(), dict(sample))
    _assert_frames_equivalent(ds, dn)


def test_two_pass_randomized_fuzz():
    rng = np.random.default_rng(20240724)
    for _ in range(25):
        n_org = int(rng.integers(1, 8))
        rows = []
        for i in range(n_org):
            n_h = int(rng.integers(0, 80))
            hashes = rng.integers(0, 200, size=n_h).tolist()
            ani = np.nan if rng.random() < 0.15 else round(float(rng.integers(90, 100)) / 100, 2)
            rows.append((f"org_{i}", ani, hashes))
        universe = int(rng.integers(50, 220))
        sample = {int(h): int(rng.integers(0, 5)) for h in range(universe)}
        df = _make_df(rows)
        ds = _run_serial_two_pass(df.copy(), dict(sample))
        dn = _run_new_two_pass(df.copy(), dict(sample))
        _assert_frames_equivalent(ds, dn)
