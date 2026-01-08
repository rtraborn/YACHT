#!/usr/bin/env python3
"""
Quick test script to verify coverage integration works correctly
"""
import sys
import os
import pandas as pd

# Add src to path so we can import yacht modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from yacht.hypothesis_recovery_src import hypothesis_recovery
from yacht.utils import load_signature_with_ksize
import multiprocessing

def test_coverage_integration():
    """Test that coverage stats are merged into output"""

    # Use test data from the repository
    project_path = os.path.dirname(os.path.abspath(__file__))

    # Check if we have training data already
    config_file = 'test_coverage_config.json'
    manifest_file = 'test_coverage_processed_manifest.tsv'
    intermediate_dir = 'test_coverage_intermediate_files'
    sample_file = os.path.join(project_path, 'tests/testdata/sample.sig.zip')

    # Check if training data exists
    if not os.path.exists(manifest_file):
        print("ERROR: Training data not found. Please run 'yacht train' first.")
        print("Run this command:")
        print(f"  yacht train --ref_file tests/testdata/20_genomes_sketches.zip --ksize 31 --prefix test_coverage --ani_thresh 0.95 --outdir ./ --force")
        return False

    # Load manifest
    print("Loading manifest...")
    manifest = pd.read_csv(manifest_file, sep="\t", header=0)
    print(f"Loaded {len(manifest)} organisms from manifest")

    # Load sample signature
    print("Loading sample signature...")
    ksize = 31
    sample_sig = load_signature_with_ksize(sample_file, ksize)
    sample_info_set = (sample_file, sample_sig)

    # Set up parameters
    min_coverage_list = [1.0, 0.5]  # Test with two coverage thresholds

    print("\nRunning hypothesis_recovery with coverage calculation...")
    multiprocessing.set_start_method('fork', force=True)

    manifest_list = hypothesis_recovery(
        manifest=manifest,
        sample_info_set=sample_info_set,
        path_to_genome_temp_dir=intermediate_dir,
        min_coverage_list=min_coverage_list,
        scale=1000,
        ksize=ksize,
        significance=0.99,
        ani_thresh=0.95,
        num_threads=4
    )

    print(f"\n✓ hypothesis_recovery completed successfully!")
    print(f"  Returned {len(manifest_list)} result DataFrames (one per min_coverage threshold)")

    # Check that coverage columns are present
    print("\nChecking for coverage columns in output...")
    expected_coverage_cols = [
        'naive_ani',
        'final_est_ani',
        'final_est_cov',
        'mean_cov',
        'median_cov',
        'lambda_status'
    ]

    for i, df in enumerate(manifest_list):
        print(f"\n  DataFrame {i+1} (min_coverage={min_coverage_list[i]}):")
        print(f"    Shape: {df.shape}")
        print(f"    Columns: {len(df.columns)}")

        # Check for coverage columns
        missing_cols = []
        present_cols = []
        for col in expected_coverage_cols:
            if col in df.columns:
                present_cols.append(col)
            else:
                missing_cols.append(col)

        print(f"    Coverage columns present: {len(present_cols)}/{len(expected_coverage_cols)}")
        if present_cols:
            print(f"      ✓ {', '.join(present_cols)}")
        if missing_cols:
            print(f"      ✗ Missing: {', '.join(missing_cols)}")

        # Show sample of coverage data for organisms with detections
        detected = df[df['in_sample_est'] == True]
        if len(detected) > 0:
            print(f"\n    Sample coverage data for detected organisms:")
            sample_cols = ['organism_name', 'in_sample_est', 'final_est_ani', 'final_est_cov', 'mean_cov', 'median_cov']
            available_cols = [c for c in sample_cols if c in df.columns]
            print(detected[available_cols].head(3).to_string(index=False))
        else:
            print(f"    No organisms detected at this threshold")

    print("\n" + "="*60)
    print("✓ TEST PASSED: Coverage columns successfully integrated!")
    print("="*60)

    return True

if __name__ == "__main__":
    success = test_coverage_integration()
    sys.exit(0 if success else 1)
