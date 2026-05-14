load("random_polynomials.sage")
load("test_R_matrix.sage")

# Re-analyze: Can det(R^*R)^(1/(2d)) be bounded by poly(d)?

print("="*70)
print("RE-ANALYSIS: Investigating poly(d) bound for det(R^*R)^(1/(2d))")
print("="*70)

# For f = x^d - 2:
# All |alpha_s| = 2^(1/d)
# u_s = (|alpha|^(2d) - 1) / (|alpha|^2 - 1)
#     = (2^2 - 1) / (2^(2/d) - 1)
#     = 3 / (2^(2/d) - 1)

# As d -> infinity, 2^(2/d) -> 1, so denominator -> 0
# Taylor expansion: 2^(2/d) ≈ 1 + (2ln2)/d
# So u_s ≈ 3 / ((2ln2)/d) = (3d)/(2ln2) ≈ 2.16 * d

print("\nAnalysis of u_s for |alpha| = 2^(1/d):")
print("u_s = 3 / (2^(2/d) - 1)")
print("\nd     |alpha|      2^(2/d)    u_s         u_s/d")

for d in [2, 3, 5, 10, 20, 50]:
    alpha = float(2^(1/d))
    alpha_2d = float(2^(2/d))
    u_s = float(3 / (alpha_2d - 1))
    print(f"{d:2d}    {alpha:.6f}   {alpha_2d:.6f}   {u_s:9.3f}   {u_s/d:.3f}")

print("\n" + "="*70)
print("Key observation: u_s ≈ const * d (linear in d)")
print("="*70)

# For this polynomial:
# - There are d roots, all with |alpha| = 2^(1/d)
# - Vandermonde term: depends on root separation
# - Geometric term: prod_s u_s ≈ (c*d)^d for some constant c

# Therefore:
# det(R^*R) ≈ V * (c*d)^d
# det(R^*R)^(1/(2d)) ≈ V^(1/(2d)) * (c*d)^(1/2) ≈ const * sqrt(d)

print("\nExpected growth of det(R^*R)^(1/(2d)):")
print("If det ≈ V * (c*d)^d, then")
print("det^(1/(2d)) ≈ V^(1/(2d)) * (c*d)^(1/2)")
print("            ≈ O(sqrt(d)) or O(d) depending on V")

# Let's verify with actual computation
print("\n" + "="*70)
print("Verification with f = x^d - 2:")
print("="*70)

R.<x> = QQ[]
results = []

for d in range(2, 15):
    f = x^d - 2
    det_val = compute_det_R_star_R_formula(f, d, precision=100)
    det_root = det_val^(1/(2*d))
    
    # Check growth rate
    ratio_d = det_root / d
    ratio_d2 = det_root / (d^2)
    ratio_sqrtd = det_root / sqrt(d)
    
    results.append((d, det_root, ratio_d, ratio_d2, ratio_sqrtd))
    
    print(f"d={d:2d}: det^(1/(2d))={det_root:8.3f}, "
          f"ratio to d={ratio_d:.3f}, "
          f"ratio to sqrt(d)={ratio_sqrtd:.3f}")

print("\n" + "="*70)
print("Analysis of growth patterns:")
print("="*70)

# Check which polynomial fits best
avg_ratio_d = numpy.mean([r[2] for r in results])
avg_ratio_sqrtd = numpy.mean([r[4] for r in results])

print(f"Average det^(1/(2d)) / d       = {avg_ratio_d:.3f}")
print(f"Average det^(1/(2d)) / sqrt(d) = {avg_ratio_sqrtd:.3f}")

if avg_ratio_d < 3:  # If ratio to d is roughly constant
    print("\nCONCLUSION: det^(1/(2d)) ≈ O(d) = poly(d) ✓")
    print("The bound IS polynomial!")
else:
    print("\nNeed to check for super-polynomial growth...")

# Try to find a truly bad case
print("\n" + "="*70)
print("Looking for worse counterexamples...")
print("="*70)

# Try polynomials with one very large root
for B in [10, 100, 1000]:
    print(f"\nTrying f = x^5 - {B}:")
    d = 5
    f = x^d - B
    alphas = f.roots(ring=ComplexField(100), multiplicities=False)
    max_alpha = max(abs(alpha) for alpha in alphas)
    
    det_val = compute_det_R_star_R_formula(f, d, precision=100)
    det_root = det_val^(1/(2*d))
    
    print(f"  max|alpha| = {max_alpha:.3f}")
    print(f"  det^(1/(2d)) = {det_root:.3f}")
    print(f"  det^(1/(2d)) / d = {det_root/d:.3f}")
