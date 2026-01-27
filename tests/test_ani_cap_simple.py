#!/usr/bin/env python
"""
Simple inline test for ANI capping logic
"""
import math

LAMBDA_EPSILON = 1e-10

def ani_from_lambda_with_cap(lambda_val, mean_cov, k_value, full_cov):
    """Copy of the fixed ani_from_lambda function"""
    if lambda_val == None:
        return None

    # Check if lambda is too close to zero so that we avoid dividing by zero
    if abs(lambda_val) < LAMBDA_EPSILON:
        return None

    contain_count = 0
    zero_count = 0
    for x in full_cov:
        if x != 0:
            contain_count += 1
        else:
            zero_count += 1

    if not full_cov:
        return None

    adj_index = contain_count / (1.0 - math.exp(-lambda_val)) / len(full_cov)

    # Calculating ani using math.pow for clarity (and to maintain type correctness)
    ani = math.pow(adj_index, 1.0 / k_value)

    # Cap ANI at 1.0 (100% identity) - biologically impossible to exceed this
    # This can happen with very low lambda values where the adjustment over-corrects
    if ani > 1.0:
        print(f"  → Capping: ANI {ani:.6f} exceeds 1.0 (lambda={lambda_val:.4f}), capping at 1.0")
        ani = 1.0

    if ani < 0.0 or math.isnan(ani):
        ret_ani = None
    else:
        ret_ani = ani

    return ret_ani

# Test scenario
contain_count = 30
full_cov_length = 100
full_cov = [1] * 30 + [0] * 70
mean_cov = sum(full_cov) // len(full_cov)

print("=" * 80)
print("Testing ANI Capping Fix")
print("=" * 80)
print(f"\nScenario: {contain_count} k-mers with coverage out of {full_cov_length} total")
print(f"Naive containment: {contain_count/full_cov_length:.4f}\n")
print(f"{'Lambda':<10} {'ANI Result':<12} {'Status'}")
print("-" * 80)

test_cases = [0.01, 0.05, 0.10, 0.20, 0.50, 1.00, 2.00, 5.00]

all_passed = True
for lambda_val in test_cases:
    ani_result = ani_from_lambda_with_cap(lambda_val, mean_cov, 31, full_cov)

    if ani_result is None:
        status = "✗ FAILED (None)"
        all_passed = False
    elif ani_result > 1.0:
        status = f"✗ FAILED ({ani_result:.6f} > 1.0)"
        all_passed = False
    elif ani_result == 1.0:
        status = "✓ CAPPED at 1.0"
    else:
        status = "✓ OK"

    print(f"{lambda_val:<10.2f} {ani_result:<12.6f} {status}")

print("\n" + "=" * 80)
if all_passed:
    print("✓ ALL TESTS PASSED - No ANI values exceed 1.0")
else:
    print("✗ SOME TESTS FAILED")
print("=" * 80)
