# Coverage Model Integration- Summary

Integrate the coverage model (cov_calc) output into YACHT's final Excel output without modifying collaborator's code.

Merging the coverage data inside `hypothesis_recovery_src.py`
- All changes were confined to `hypothesis_recovery_src.py` (already modified on superyacht branch)
- Zero changes were to `run_YACHT.py` (sought to avoid editing existing yacht code)
- Coverage statistics automatically intergrated into existing output pipeline

## Details

### File Modified
- `src/yacht/hypothesis_recovery_src.py` (31 insertions, 31 deletions)

### Changes

1. Added organism_name to coverage DataFrame (Line 214)
```python
# Add organism_name to final_stats_df for merging (stats are in same order as in the sub_manifest)
final_stats_df['organism_name'] = sub_manifest['organism_name'].values
```

2. Merged coverage stats into output (Lines 462-483)
```python
# Merge coverage statistics into each manifest DataFrame
# Select key coverage columns to include in output
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

# Merge coverage stats into each manifest in the list
for i in range(len(manifest_list)):
    manifest_list[i] = manifest_list[i].merge(
        coverage_stats,
        on='organism_name',
        how='left'  # Keeps all organisms, even those without coverage data
    )
```

#### 3. Reverted return signature
```python
return manifest_list # no longer a tuple as it was before
```

4. Code Cleanup
- Removed unused debugging variables (e.g. sample_hashes_keys, samp_dict, etc.)
- Removed commented-out print statements
- Removed unused summary statistics calculation
- Added additional comments 

## New Coverage Columns in Output

The Excel output now includes these additional columns for each organism:

| Column | Description |
|--------|-------------|
| `naive_ani` | Simple containment-based ANI estimate |
| `final_est_ani` | Coverage-adjusted ANI estimate (Shaw & Yu, 2024) |
| `final_est_cov` | Lambda - expected coverage value |
| `mean_cov` | Mean coverage across matched k-mers |
| `median_cov` | Median coverage across matched k-mers |
| `lambda_status` | Adjustment status (LAMBDA/HIGH/LOW enum) |
| `ani_ci` | ANI confidence interval (low, high) |
| `lambda_ci` | Lambda confidence interval (low, high) |

## Demonstration of the output

```
BEFORE (Original YACHT output):
organism_name  in_sample_est  p_vals  num_exclusive_kmers_to_genome  num_matches
   Organism_A           True   0.001                           1000           45
   Organism_B           True   0.005                            850           38

AFTER (i.e. with coverage integration):
organism_name  in_sample_est  p_vals  ...  final_est_ani  final_est_cov  mean_cov  lambda_status
   Organism_A           True   0.001  ...          0.968            3.2       3.1         LAMBDA
   Organism_B           True   0.005  ...          0.951            2.1       2.0         LAMBDA
```

## Testing

### Demonstration Script
- `demo_coverage_merge.py` - Demonstrates the merge logic with sample data
- Shows before/after DataFrames
- Confirms coverage columns are properly integrated

### Syntax Validation
- ✓ Python syntax check passed (`python3 -m py_compile`)
- ✓ No import errors in modified file

## Future Considerations

- Coverage stats appear as `NaN` for organisms with zero overlap (correct behavior): should this become "0"?
- If bootstrap confidence intervals are disabled (`ci_int=False`), the `ani_ci` and `lambda_ci` columns will show (None, None)
- The `lambda_status` column uses your new Pythonic enum (`AdjustStatusType`)

## Files Created

1. `demo_coverage_merge.py` - Demonstration of the merge logic
2. `test_coverage_integration.py` - Integration test script (requires YACHT installation)
3. `coverage_integration_summary.md` - This file

## Next Steps

1. Test with real data when YACHT environment is available
2. Verify that the coverage columns are included in the Excel
3. Consider if additional coverage columns should be included
4. Update documentation to describe new output columns?

---

