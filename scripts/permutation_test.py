"""
Permutation test on T-type → shell-bearing-fraction trend.

For §3.3 of Bibb (v6.2): tests whether the observed Spearman rho = -0.81
across T = 2..9 could have arisen from a random labeling of T-types.

Usage:
    python permutation_test.py [path_to_canonical_fits.csv]

The default path is ../data/sparc_T2-T9_canonical_fits.csv relative to this
script. Output is printed to stdout: observed rho, permutation p-values
(one and two-sided, N=10000 shuffles), and a non-parametric bootstrap 95%
confidence interval on rho.
"""
import sys
import os
import numpy as np
import pandas as pd
from scipy import stats

CSV = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(__file__), '..', 'data', 'sparc_T2-T9_canonical_fits.csv'
)

def per_bin_fractions(T_arr, s_arr):
    """Return Nx3 array of (T, n_galaxies, shell-bearing fraction) for each unique T."""
    out = []
    for t in sorted(np.unique(T_arr)):
        mask = T_arr == t
        out.append((t, mask.sum(), s_arr[mask].sum() / mask.sum()))
    return np.array(out)

def main():
    df = pd.read_csv(CSV)
    T = df['T'].values
    shell = (df['fw_best_n_shells'] > 0).astype(int).values
    print(f"N galaxies = {len(T)}, shell-bearing = {shell.sum()}")

    # Observed per-bin Spearman
    bins = per_bin_fractions(T, shell)
    obs_rho, obs_p_analytic = stats.spearmanr(bins[:,0], bins[:,2])
    print(f"\nObserved per-bin Spearman: rho = {obs_rho:.4f}, "
          f"analytic p = {obs_p_analytic:.4f}")

    # Permutation test
    rng = np.random.default_rng(42)
    N = 10000
    null_rhos = np.empty(N)
    for i in range(N):
        T_shuf = rng.permutation(T)
        b = per_bin_fractions(T_shuf, shell)
        r, _ = stats.spearmanr(b[:,0], b[:,2])
        null_rhos[i] = r if not np.isnan(r) else 0.0
    p_two = (np.abs(null_rhos) >= np.abs(obs_rho)).sum() / N
    p_one = (null_rhos <= obs_rho).sum() / N
    print(f"\nPermutation test (N = {N}, T labels shuffled):")
    print(f"  one-sided P(rho_shuf <= {obs_rho:.3f}) = {p_one:.4f}")
    print(f"  two-sided P(|rho_shuf| >= {abs(obs_rho):.3f}) = {p_two:.4f}")
    print(f"  null mean rho = {null_rhos.mean():+.4f}, std = {null_rhos.std():.4f}")

    # Bootstrap 95% CI
    rng = np.random.default_rng(123)
    N_boot = 10000
    boot_rhos = []
    n = len(T)
    for _ in range(N_boot):
        idx = rng.integers(0, n, n)
        b = per_bin_fractions(T[idx], shell[idx])
        if len(b) < 4:
            continue
        r, _ = stats.spearmanr(b[:,0], b[:,2])
        if not np.isnan(r):
            boot_rhos.append(r)
    boot_rhos = np.array(boot_rhos)
    ci_low, ci_high = np.percentile(boot_rhos, [2.5, 97.5])
    print(f"\nBootstrap 95% CI on rho: [{ci_low:.3f}, {ci_high:.3f}] "
          f"(from {len(boot_rhos)} valid resamples)")

if __name__ == "__main__":
    main()
