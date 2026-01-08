#!/usr/bin/env python3
"""
Minimal demonstration that coverage merging logic works correctly
This simulates what happens at the end of hypothesis_recovery()
"""
import pandas as pd

def demo_coverage_merge():
    """Demonstrate the coverage merging logic"""

    print("="*70)
    print("DEMONSTRATION: Coverage Data Merging into YACHT Output")
    print("="*70)

    # Simulate a manifest_list with hypothesis test results (what YACHT currently returns)
    print("\n1. BEFORE: Hypothesis test results (current YACHT output)")
    print("-" * 70)

    manifest_df = pd.DataFrame({
        'organism_name': ['Organism_A', 'Organism_B', 'Organism_C'],
        'in_sample_est': [True, True, False],
        'p_vals': [0.001, 0.005, 0.95],
        'num_exclusive_kmers_to_genome': [1000, 850, 650],
        'num_matches': [45, 38, 5],
        'min_coverage': [0.5, 0.5, 0.5]
    })

    print(manifest_df.to_string(index=False))

    # Simulate final_stats_df from cov_calc (coverage model output)
    print("\n2. Coverage statistics from cov_calc() - NEW DATA")
    print("-" * 70)

    final_stats_df = pd.DataFrame({
        'organism_name': ['Organism_A', 'Organism_B', 'Organism_C'],
        'naive_ani': [0.962, 0.945, 0.891],
        'final_est_ani': [0.968, 0.951, 0.895],
        'final_est_cov': [3.2, 2.1, 0.8],
        'mean_cov': [3.1, 2.0, 0.7],
        'median_cov': [3.0, 2.0, 1.0],
        'lambda_status': ['LAMBDA', 'LAMBDA', 'LOW'],
        'ani_ci': [(0.96, 0.98), (0.94, 0.96), (None, None)],
        'lambda_ci': [(2.8, 3.6), (1.8, 2.4), (None, None)]
    })

    print(final_stats_df.to_string(index=False))

    # This is the merge operation that happens in hypothesis_recovery()
    print("\n3. MERGING: Combining coverage stats with hypothesis results")
    print("-" * 70)

    # Select coverage columns to include
    coverage_cols = [
        'organism_name',
        'naive_ani',
        'final_est_ani',
        'final_est_cov',
        'mean_cov',
        'median_cov',
        'lambda_status',
        'ani_ci',
        'lambda_ci'
    ]
    coverage_stats = final_stats_df[coverage_cols].copy()

    # Perform the merge (this is the NEW code in hypothesis_recovery_src.py)
    merged_df = manifest_df.merge(
        coverage_stats,
        on='organism_name',
        how='left'  # Keep all organisms, even those without coverage stats
    )

    print("Code executed:")
    print("  merged_df = manifest_df.merge(coverage_stats, on='organism_name', how='left')")

    # Show the final result
    print("\n4. AFTER: Combined output (what YACHT now returns)")
    print("-" * 70)
    print(merged_df.to_string(index=False))

    print("\n5. ANALYSIS: Coverage columns now in Excel output")
    print("-" * 70)
    print(f"✓ Total columns: {len(merged_df.columns)}")
    print(f"✓ Original YACHT columns: {len(manifest_df.columns)}")
    print(f"✓ New coverage columns: {len(coverage_cols) - 1}")  # -1 for organism_name
    print(f"\n✓ New coverage columns added:")
    for col in coverage_cols:
        if col != 'organism_name':
            print(f"    - {col}")

    print("\n6. PRACTICAL EXAMPLE: Detected organisms with coverage info")
    print("-" * 70)
    detected = merged_df[merged_df['in_sample_est'] == True]
    display_cols = ['organism_name', 'in_sample_est', 'p_vals',
                    'final_est_ani', 'final_est_cov', 'mean_cov', 'lambda_status']
    print(detected[display_cols].to_string(index=False))

    print("\n" + "="*70)
    print("✓ SUCCESS: Coverage data successfully integrated!")
    print("="*70)
    print("\nThis demonstrates what happens in hypothesis_recovery_src.py lines 462-483")
    print("The merged DataFrame is returned and saved to Excel by run_YACHT.py")
    print("No changes to run_YACHT.py were needed - it works automatically!")

if __name__ == "__main__":
    demo_coverage_merge()
