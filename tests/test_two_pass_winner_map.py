#!/usr/bin/env python
"""
Test script to compare one-pass vs two-pass winner-takes-all approaches.

This script compares:
1. One-pass: Original approach using initial ANI estimates
2. Two-pass: Sylph-aligned approach that recalculates ANI using won k-mers

USAGE - Compare one-pass vs two-pass on the same sample:
=========================================================

Run YACHT twice with the same data, once with two-pass (default) and once without:

    # Two-pass (default, sylph-aligned)
    python -m yacht run_YACHT --json <config.json> --sample_file <sample.sig.zip> \\
        --winner_takes_all --out results_two_pass.xlsx

    # One-pass (original method)
    python -m yacht run_YACHT --json <config.json> --sample_file <sample.sig.zip> \\
        --winner_takes_all --no_two_pass --out results_one_pass.xlsx

Then compare the results:
    python tests/test_two_pass_winner_map.py --one_pass results_one_pass.xlsx \\
        --two_pass results_two_pass.xlsx

Or run with pytest for the unit tests:
    pytest tests/test_two_pass_winner_map.py -v
"""

import argparse
import sys
import os
import pandas as pd
import numpy as np
from typing import Dict, Tuple
from copy import deepcopy

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sourmash
from yacht.utils import load_signature_with_ksize
from yacht.hypothesis_recovery_src import (
    build_winner_map,
    recalculate_ani_from_winner_map,
    estimate_relative_abundance,
)


def run_one_pass(
    final_stats_df: pd.DataFrame,
    path_to_genome_temp_dir: str,
    ksize: int,
    sample_sig: sourmash.SourmashSignature,
    batch_size: int = 1000
) -> pd.DataFrame:
    """
    Run the original one-pass winner-takes-all approach.

    This is the approach YACHT used before implementing sylph's two-pass method.
    """
    df = final_stats_df.copy()

    # Single pass: build winner map and calculate abundance
    print("ONE-PASS: Building winner map with initial ANI estimates...")
    winner_map = build_winner_map(df, path_to_genome_temp_dir, ksize, batch_size)

    print("ONE-PASS: Calculating relative abundance...")
    df = estimate_relative_abundance(df, winner_map, sample_sig, batch_size)

    # Add columns for consistency
    df['reassignment_status'] = 'one_pass'
    df['original_ani'] = df['final_est_ani'].copy()

    return df


def run_two_pass(
    final_stats_df: pd.DataFrame,
    path_to_genome_temp_dir: str,
    ksize: int,
    sample_sig: sourmash.SourmashSignature,
    batch_size: int = 1000
) -> pd.DataFrame:
    """
    Run the new two-pass winner-takes-all approach (sylph-aligned).

    Pass 1: Build winner map using initial ANI estimates
    Recalculate ANI using only won k-mers
    Pass 2: Rebuild winner map with refined ANI estimates
    """
    df = final_stats_df.copy()

    # Pass 1
    print("TWO-PASS: Pass 1 - Building initial winner map...")
    winner_map = build_winner_map(df, path_to_genome_temp_dir, ksize, batch_size)

    # Recalculate ANI
    print("TWO-PASS: Recalculating ANI using won k-mers...")
    df = recalculate_ani_from_winner_map(df, winner_map, sample_sig, ksize, batch_size)

    # Pass 2
    print("TWO-PASS: Pass 2 - Rebuilding winner map with refined ANI...")
    winner_map = build_winner_map(df, path_to_genome_temp_dir, ksize, batch_size)

    # Calculate relative abundance
    print("TWO-PASS: Calculating relative abundance...")
    df = estimate_relative_abundance(df, winner_map, sample_sig, batch_size)

    return df


def compare_results(one_pass_df: pd.DataFrame, two_pass_df: pd.DataFrame) -> Dict:
    """
    Compare results between one-pass and two-pass approaches.

    Returns a dictionary with comparison metrics.
    """
    results = {}

    # Merge on organism_name for comparison
    merged = one_pass_df.merge(
        two_pass_df,
        on='organism_name',
        suffixes=('_1pass', '_2pass')
    )

    # Count organisms
    results['total_organisms'] = len(merged)

    # Check for eliminated organisms in two-pass
    if 'reassignment_status_2pass' in merged.columns:
        eliminated = merged[merged['reassignment_status_2pass'] == 'eliminated']
        results['eliminated_in_two_pass'] = len(eliminated)
        results['eliminated_organisms'] = eliminated['organism_name'].tolist()
    else:
        results['eliminated_in_two_pass'] = 0
        results['eliminated_organisms'] = []

    # Compare ANI values (only for non-eliminated organisms)
    active_mask = merged['reassignment_status_2pass'] != 'eliminated'
    active = merged[active_mask]

    if len(active) > 0 and 'final_est_ani_1pass' in active.columns:
        ani_diff = active['final_est_ani_2pass'] - active['final_est_ani_1pass']
        results['ani_mean_diff'] = ani_diff.mean()
        results['ani_std_diff'] = ani_diff.std()
        results['ani_max_diff'] = ani_diff.abs().max()

        # Organisms with significant ANI change (>0.01)
        significant_change = active[ani_diff.abs() > 0.01]
        results['organisms_with_significant_ani_change'] = len(significant_change)

        if len(significant_change) > 0:
            results['significant_ani_changes'] = significant_change[[
                'organism_name', 'final_est_ani_1pass', 'final_est_ani_2pass'
            ]].to_dict('records')

    # Compare relative abundance
    if 'rel_abund_1pass' in merged.columns and 'rel_abund_2pass' in merged.columns:
        # Only compare non-NaN values
        valid_mask = ~(merged['rel_abund_1pass'].isna() | merged['rel_abund_2pass'].isna())
        valid = merged[valid_mask]

        if len(valid) > 0:
            abund_diff = valid['rel_abund_2pass'] - valid['rel_abund_1pass']
            results['rel_abund_mean_diff'] = abund_diff.mean()
            results['rel_abund_std_diff'] = abund_diff.std()
            results['rel_abund_max_diff'] = abund_diff.abs().max()

            # Top differences
            valid_sorted = valid.copy()
            valid_sorted['abund_diff'] = abund_diff.values
            valid_sorted = valid_sorted.sort_values('abund_diff', key=abs, ascending=False)

            results['top_abundance_differences'] = valid_sorted[[
                'organism_name', 'rel_abund_1pass', 'rel_abund_2pass', 'abund_diff'
            ]].head(10).to_dict('records')

    return results


def print_comparison_report(results: Dict):
    """Print a formatted comparison report."""
    print("\n" + "=" * 80)
    print("COMPARISON REPORT: One-Pass vs Two-Pass Winner-Takes-All")
    print("=" * 80)

    print(f"\nTotal organisms analyzed: {results.get('total_organisms', 'N/A')}")
    print(f"Organisms eliminated in two-pass: {results.get('eliminated_in_two_pass', 0)}")

    if results.get('eliminated_organisms'):
        print("\nEliminated organisms:")
        for org in results['eliminated_organisms'][:10]:  # Show first 10
            print(f"  - {org}")
        if len(results['eliminated_organisms']) > 10:
            print(f"  ... and {len(results['eliminated_organisms']) - 10} more")

    print("\n--- ANI Comparison ---")
    if 'ani_mean_diff' in results:
        print(f"Mean ANI difference (2pass - 1pass): {results['ani_mean_diff']:.6f}")
        print(f"Std ANI difference: {results['ani_std_diff']:.6f}")
        print(f"Max absolute ANI difference: {results['ani_max_diff']:.6f}")
        print(f"Organisms with significant ANI change (>0.01): {results.get('organisms_with_significant_ani_change', 0)}")

        if results.get('significant_ani_changes'):
            print("\nSignificant ANI changes:")
            for change in results['significant_ani_changes'][:5]:
                print(f"  {change['organism_name'][:50]}...")
                print(f"    1-pass: {change['final_est_ani_1pass']:.4f} -> 2-pass: {change['final_est_ani_2pass']:.4f}")

    print("\n--- Relative Abundance Comparison ---")
    if 'rel_abund_mean_diff' in results:
        print(f"Mean rel_abund difference: {results['rel_abund_mean_diff']:.6f}")
        print(f"Std rel_abund difference: {results['rel_abund_std_diff']:.6f}")
        print(f"Max absolute rel_abund difference: {results['rel_abund_max_diff']:.6f}")

        if results.get('top_abundance_differences'):
            print("\nTop abundance differences:")
            for diff in results['top_abundance_differences'][:5]:
                print(f"  {diff['organism_name'][:50]}...")
                print(f"    1-pass: {diff['rel_abund_1pass']:.6f} -> 2-pass: {diff['rel_abund_2pass']:.6f} (diff: {diff['abund_diff']:.6f})")

    print("\n" + "=" * 80)


# ============================================================================
# Unit Tests (run with pytest)
# ============================================================================

def test_recalculate_ani_handles_empty_winner_map():
    """Test that recalculate_ani handles organisms with no won k-mers."""
    # This would need mock data - placeholder for now
    pass


def test_two_pass_produces_different_results():
    """Test that two-pass produces different results than one-pass for related organisms."""
    # This would need mock data with closely related organisms - placeholder for now
    pass


def test_eliminated_organisms_have_zero_abundance():
    """Test that eliminated organisms have rel_abund = 0."""
    # This would need mock data - placeholder for now
    pass


# ============================================================================
# Main CLI
# ============================================================================

def load_yacht_results(filepath: str) -> pd.DataFrame:
    """Load YACHT results from Excel file."""
    # YACHT outputs multiple sheets; we want the one with coverage data
    # Try to find the sheet with rel_abund column
    xlsx = pd.ExcelFile(filepath)
    for sheet_name in xlsx.sheet_names:
        df = pd.read_excel(filepath, sheet_name=sheet_name)
        if 'rel_abund' in df.columns or 'organism_name' in df.columns:
            return df
    # Default to first sheet
    return pd.read_excel(filepath, sheet_name=0)


def compare_excel_results(one_pass_file: str, two_pass_file: str) -> Dict:
    """Compare results from two YACHT Excel output files."""
    print(f"Loading one-pass results from {one_pass_file}...")
    one_pass_df = load_yacht_results(one_pass_file)

    print(f"Loading two-pass results from {two_pass_file}...")
    two_pass_df = load_yacht_results(two_pass_file)

    return compare_results(one_pass_df, two_pass_df)


def main():
    parser = argparse.ArgumentParser(
        description="Compare one-pass vs two-pass winner-takes-all approaches"
    )
    parser.add_argument(
        "--one_pass",
        help="Path to YACHT results Excel file from --no_two_pass run"
    )
    parser.add_argument(
        "--two_pass",
        help="Path to YACHT results Excel file from default (two-pass) run"
    )
    parser.add_argument(
        "--output",
        help="Output CSV file for detailed comparison (optional)"
    )

    args = parser.parse_args()

    if args.one_pass and args.two_pass:
        # Compare two result files
        results = compare_excel_results(args.one_pass, args.two_pass)
        print_comparison_report(results)

        if args.output:
            # Save detailed comparison
            print(f"\nSaving detailed comparison to {args.output}")
            comparison_df = pd.DataFrame([results])
            comparison_df.to_csv(args.output, index=False)
    else:
        print(__doc__)
        print("\nTo compare results, run YACHT twice:")
        print("  1. With two-pass (default):  --winner_takes_all")
        print("  2. With one-pass:            --winner_takes_all --no_two_pass")
        print("\nThen compare:")
        print("  python tests/test_two_pass_winner_map.py --one_pass <file1.xlsx> --two_pass <file2.xlsx>")


if __name__ == "__main__":
    main()
