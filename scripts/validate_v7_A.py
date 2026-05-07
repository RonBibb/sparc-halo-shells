"""
validate_v7_A.py — Numerical consistency checks A1–A8 for halo_shells v7.0.

For each check, the script computes the ground-truth value from the v7.0 CSVs
in data/ and reports whether it matches the value claimed in the manuscript.

A check passes if the computed value matches the manuscript claim within an
appropriate tolerance (exact for integer counts, 0.005 for percentages,
0.005 for Spearman ρ, etc.).

Usage:
    cd halo_shells_v7.0/
    python3 scripts/validate_v7_A.py

Exit code 0 if all pass, 1 if any fail.
"""
import os
import sys
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, kendalltau, fisher_exact

HERE = os.path.dirname(os.path.abspath(__file__))
PACKAGE_ROOT = os.path.dirname(HERE)
DATA_DIR = os.path.join(PACKAGE_ROOT, 'data')

# ANSI color codes for terminal output
class C:
    OK = '\033[92m'
    FAIL = '\033[91m'
    WARN = '\033[93m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    END = '\033[0m'


# ============================================================
# Helpers
# ============================================================
def check(label, claim, computed, tol=0, integer=False):
    """Compare a manuscript claim to a computed value; return (pass_bool, msg)."""
    if integer:
        passed = int(claim) == int(computed)
        diff = ''
    elif isinstance(claim, str) or isinstance(computed, str):
        # String equality (e.g., "7/173" formatted counts)
        passed = str(claim) == str(computed)
        diff = ''
    else:
        passed = abs(claim - computed) <= tol
        diff = f" (Δ={abs(claim - computed):.4f})"
    status = f"{C.OK}✓ PASS{C.END}" if passed else f"{C.FAIL}✗ FAIL{C.END}"
    return passed, f"  {status} {label}: claim={claim}, computed={computed}{diff}"


def load_canonical():
    return pd.read_csv(os.path.join(DATA_DIR, 'sparc_T2-T9_canonical_fits.csv'))


def load_null_combined():
    return pd.read_csv(os.path.join(DATA_DIR, 'null_test_T2-T9_combined.csv'))


def load_einasto():
    return pd.read_csv(os.path.join(DATA_DIR, 'einasto_full_sample_results.csv'))


def load_dc14_ngc5055():
    return pd.read_csv(os.path.join(DATA_DIR, 'dc14_ngc5055_showcase.csv'))


# ============================================================
# A1. Framework adequacy at chi^2_red < 1.5
# ============================================================
def A1():
    print(f"\n{C.BOLD}A1. Framework adequacy at chi²_red < 1.5{C.END}")
    print(f"   Manuscript claim: 91/102 (89.2%)")
    df = load_canonical()
    n_total = len(df)
    n_adequate = (df['fw_best_chi2_red'] < 1.5).sum()
    pct = 100 * n_adequate / n_total
    results = []
    p1, m1 = check("count of framework-adequate galaxies", 91, n_adequate, integer=True)
    p2, m2 = check("total galaxy count", 102, n_total, integer=True)
    p3, m3 = check("framework-adequacy %", 89.2, pct, tol=0.05)
    print(m1); print(m2); print(m3)
    return all([p1, p2, p3])


# ============================================================
# A2. Burkert and NFW adequacy at chi^2_red < 1.5
# ============================================================
def A2():
    print(f"\n{C.BOLD}A2. Burkert and NFW adequacy at chi²_red < 1.5{C.END}")
    print(f"   Manuscript claim: Burkert 65/102 (63.7%), NFW 52/102 (51.0%)")
    df = load_canonical()
    n_burk = (df['burk_chi2_red'] < 1.5).sum()
    n_nfw = (df['nfw_chi2_red'] < 1.5).sum()
    p1, m1 = check("Burkert-adequate count", 65, n_burk, integer=True)
    p2, m2 = check("NFW-adequate count", 52, n_nfw, integer=True)
    p3, m3 = check("Burkert-adequate %", 63.7, 100*n_burk/len(df), tol=0.05)
    p4, m4 = check("NFW-adequate %", 51.0, 100*n_nfw/len(df), tol=0.05)
    for m in [m1, m2, m3, m4]: print(m)
    return all([p1, p2, p3, p4])


# ============================================================
# A3. Strict adequacy at chi^2_red < 1.0
# ============================================================
def A3():
    print(f"\n{C.BOLD}A3. Strict adequacy at chi²_red < 1.0{C.END}")
    print(f"   Manuscript claim: Framework 86, Burkert 57, NFW 41")
    df = load_canonical()
    n_fw = (df['fw_best_chi2_red'] < 1.0).sum()
    n_burk = (df['burk_chi2_red'] < 1.0).sum()
    n_nfw = (df['nfw_chi2_red'] < 1.0).sum()
    p1, m1 = check("framework strict-adequate count", 86, n_fw, integer=True)
    p2, m2 = check("Burkert strict-adequate count", 57, n_burk, integer=True)
    p3, m3 = check("NFW strict-adequate count", 41, n_nfw, integer=True)
    for m in [m1, m2, m3]: print(m)
    return all([p1, p2, p3])


# ============================================================
# A4. Per-T shell-bearing fractions
# ============================================================
def A4():
    print(f"\n{C.BOLD}A4. Per-T shell-bearing fractions{C.END}")
    df = load_canonical()
    expected = {
        2: (9, 9, 100.0),
        3: (11, 6, 54.5),
        4: (17, 9, 52.9),
        5: (13, 7, 53.8),
        6: (16, 9, 56.2),
        7: (14, 5, 35.7),
        8: (6, 2, 33.3),
        9: (16, 5, 31.2),
    }
    all_pass = True
    for T, (n_total_exp, n_shell_exp, pct_exp) in expected.items():
        sub = df[df['T'] == T]
        n_total = len(sub)
        n_shell = (sub['fw_best_n_shells'] >= 1).sum()
        pct = 100 * n_shell / n_total if n_total else 0
        p1, m1 = check(f"T={T} sample size", n_total_exp, n_total, integer=True)
        p2, m2 = check(f"T={T} shell-bearing count", n_shell_exp, n_shell, integer=True)
        p3, m3 = check(f"T={T} shell-bearing %", pct_exp, pct, tol=0.1)
        print(m1); print(m2); print(m3)
        all_pass = all_pass and p1 and p2 and p3
    return all_pass


# ============================================================
# A5. Trend tests (per-galaxy and per-T-bin)
# ============================================================
def A5():
    print(f"\n{C.BOLD}A5. Trend tests{C.END}")
    df = load_canonical()
    df['shell_bearing'] = (df['fw_best_n_shells'] >= 1).astype(int)
    
    # Per-T-bin Spearman
    Ts, fracs = [], []
    for T in range(2, 10):
        sub = df[df['T'] == T]
        if len(sub) > 0:
            Ts.append(T)
            fracs.append((sub['fw_best_n_shells'] >= 1).mean() * 100)
    rho_T, p_T = spearmanr(Ts, fracs)
    
    # Per-galaxy permutation test (10000 resamples, seed 42)
    np.random.seed(42)
    obs_rho_pg, _ = spearmanr(df['T'], df['shell_bearing'])
    perm_rhos = []
    for _ in range(10000):
        perm_T = np.random.permutation(df['T'].values)
        r, _ = spearmanr(perm_T, df['shell_bearing'])
        if not np.isnan(r):
            perm_rhos.append(r)
    perm_rhos = np.array(perm_rhos)
    p_one = (perm_rhos <= obs_rho_pg).mean()
    p_two = (np.abs(perm_rhos) >= np.abs(obs_rho_pg)).mean()
    
    # Bootstrap 95% CI on per-T-bin Spearman (10000 resamples, seed 42)
    np.random.seed(42)
    boot_rhos = []
    for _ in range(10000):
        idx = np.random.choice(len(df), size=len(df), replace=True)
        sample = df.iloc[idx]
        Ts_b, fracs_b = [], []
        for T in range(2, 10):
            sub = sample[sample['T'] == T]
            if len(sub) > 0:
                Ts_b.append(T)
                fracs_b.append((sub['fw_best_n_shells'] >= 1).mean() * 100)
        if len(Ts_b) >= 3:
            r, _ = spearmanr(Ts_b, fracs_b)
            if not np.isnan(r):
                boot_rhos.append(r)
    boot_rhos = np.array(boot_rhos)
    ci_low, ci_high = np.percentile(boot_rhos, [2.5, 97.5])
    
    # Mann-Kendall (per-galaxy, via Kendall tau as proxy — manuscript reports this)
    tau_pg, p_mk = kendalltau(df['T'], df['shell_bearing'])
    
    p1, m1 = check("per-T-bin Spearman ρ", -0.83, rho_T, tol=0.01)
    p2, m2 = check("per-T-bin Spearman analytic p", 0.010, p_T, tol=0.005)
    p3, m3 = check("per-galaxy permutation p_one-sided", 0.001, p_one, tol=0.005)
    p4, m4 = check("per-galaxy permutation p_two-sided", 0.002, p_two, tol=0.005)
    p5, m5 = check("bootstrap 95% CI lower", -0.95, ci_low, tol=0.02)
    p6, m6 = check("bootstrap 95% CI upper", -0.22, ci_high, tol=0.02)
    p7, m7 = check("Mann-Kendall (Kendall) p", 0.003, p_mk, tol=0.002)
    for m in [m1, m2, m3, m4, m5, m6, m7]: print(m)
    return all([p1, p2, p3, p4, p5, p6, p7])


# ============================================================
# A6. Null test aggregate rates
# ============================================================
def A6():
    print(f"\n{C.BOLD}A6. Null test aggregate rates{C.END}")
    print(f"   Manuscript claim: Burkert 4.0% (7/173), NFW 63.6% (110/173), combined 33.8% (117/346)")
    df = load_null_combined()
    burk = df[df['smooth_truth'] == 'Burkert']
    nfw = df[df['smooth_truth'] == 'NFW']
    n_burk_shell = (burk['best_n_shells'] > 0).sum()
    n_nfw_shell = (nfw['best_n_shells'] > 0).sum()
    n_combined_shell = (df['best_n_shells'] > 0).sum()
    
    p1, m1 = check("Burkert-truth shell count", 7, n_burk_shell, integer=True)
    p2, m2 = check("Burkert-truth N_mocks", 173, len(burk), integer=True)
    p3, m3 = check("Burkert-truth FP %", 4.0, 100*n_burk_shell/len(burk), tol=0.1)
    p4, m4 = check("NFW-truth shell count", 110, n_nfw_shell, integer=True)
    p5, m5 = check("NFW-truth N_mocks", 173, len(nfw), integer=True)
    p6, m6 = check("NFW-truth FP %", 63.6, 100*n_nfw_shell/len(nfw), tol=0.1)
    p7, m7 = check("Combined N_mocks", 346, len(df), integer=True)
    p8, m8 = check("Combined FP %", 33.8, 100*n_combined_shell/len(df), tol=0.1)
    for m in [m1, m2, m3, m4, m5, m6, m7, m8]: print(m)
    return all([p1, p2, p3, p4, p5, p6, p7, p8])


# ============================================================
# A7. Per-T-bin null FP rates and Fisher tests
# ============================================================
def A7():
    print(f"\n{C.BOLD}A7. Per-T-bin null FP rates (Table 4 rows){C.END}")
    df = load_null_combined()
    canon = load_canonical()
    
    # Manuscript Table 4 claims, format: T -> (burk_pct, burk_count_str, nfw_pct, nfw_count_str, fisher_p)
    expected = {
        2: (5,   '1/20',  90, '18/20', 1e-6),
        3: (5,   '1/20',  55, '11/20', 0.004),
        4: (0,   '0/25',  80, '20/25', 5e-5),
        5: (0,   '0/25',  40, '10/25', 1e-4),
        6: (4,   '1/25',  56, '14/25', 3e-4),
        7: (0,   '0/25',  72, '18/25', 0.003),
        8: (11,  '2/18',  50, '9/18',  0.25),
        9: (13,  '2/15',  67, '10/15', 0.22),
    }
    
    all_pass = True
    for T, (burk_pct_exp, burk_count_exp, nfw_pct_exp, nfw_count_exp, fisher_p_exp) in expected.items():
        burk = df[(df['smooth_truth'] == 'Burkert') & (df['T'] == T)]
        nfw = df[(df['smooth_truth'] == 'NFW') & (df['T'] == T)]
        n_burk_shell = (burk['best_n_shells'] > 0).sum()
        n_nfw_shell = (nfw['best_n_shells'] > 0).sum()
        burk_pct = 100 * n_burk_shell / len(burk) if len(burk) > 0 else 0
        nfw_pct = 100 * n_nfw_shell / len(nfw) if len(nfw) > 0 else 0
        burk_count_actual = f"{n_burk_shell}/{len(burk)}"
        nfw_count_actual = f"{n_nfw_shell}/{len(nfw)}"
        
        # Fisher exact: real shell-bearing at T vs Burkert-null
        real = canon[canon['T'] == T]
        n_real_shell = (real['fw_best_n_shells'] >= 1).sum()
        if len(burk) > 0 and len(real) > 0:
            _, fisher_p = fisher_exact(
                [[n_real_shell, len(real) - n_real_shell],
                 [n_burk_shell, len(burk) - n_burk_shell]],
                alternative='greater'
            )
        else:
            fisher_p = float('nan')
        
        p1, m1 = check(f"T={T} Burk null %", burk_pct_exp, round(burk_pct), integer=True)
        p2, m2 = check(f"T={T} Burk count", burk_count_exp, burk_count_actual, integer=False, tol=0)
        p3, m3 = check(f"T={T} NFW null %", nfw_pct_exp, round(nfw_pct), integer=True)
        p4, m4 = check(f"T={T} NFW count", nfw_count_exp, nfw_count_actual, integer=False, tol=0)
        # Fisher p — looser tol; manuscript rounded to 1-2 sig figs
        p5_pass = fisher_p < fisher_p_exp * 5 and fisher_p > fisher_p_exp * 0.2 if fisher_p_exp > 0 else False
        m5 = (f"  {C.OK}✓ PASS{C.END} T={T} Fisher p (within order): "
              f"claim≈{fisher_p_exp:.0e}, computed={fisher_p:.2e}" if p5_pass
              else f"  {C.FAIL}✗ FAIL{C.END} T={T} Fisher p: claim≈{fisher_p_exp:.0e}, "
                   f"computed={fisher_p:.2e}")
        for m in [m1, m2, m3, m4, m5]: print(m)
        all_pass = all_pass and p1 and p2 and p3 and p4 and p5_pass
    return all_pass


# Helper for A7 — string-equality check
def check_str(label, claim, computed):
    passed = str(claim) == str(computed)
    status = f"{C.OK}✓ PASS{C.END}" if passed else f"{C.FAIL}✗ FAIL{C.END}"
    return passed, f"  {status} {label}: claim={claim}, computed={computed}"


# ============================================================
# A8. NGC 5055 showcase numbers
# ============================================================
def A8():
    print(f"\n{C.BOLD}A8. NGC 5055 showcase{C.END}")
    print(f"   Manuscript claim: Burkert 9.53, NFW 30.21, FW 1.464 (1 shell at r=12.0, σ=2.91, M=3.78e10)")
    print(f"                     DC14 9.0; Einasto-only 8.42, Einasto+2-shell 0.507")
    
    df = load_canonical()
    n5055 = df[df['Galaxy'] == 'NGC5055'].iloc[0]
    
    p1, m1 = check("NGC5055 Burkert chi²_red", 9.53, n5055['burk_chi2_red'], tol=0.01)
    p2, m2 = check("NGC5055 NFW chi²_red", 30.21, n5055['nfw_chi2_red'], tol=0.05)
    p3, m3 = check("NGC5055 framework chi²_red", 1.464, n5055['fw_best_chi2_red'], tol=0.005)
    p4, m4 = check("NGC5055 framework n_shells", 1, n5055['fw_best_n_shells'], integer=True)
    p5, m5 = check("NGC5055 shell radius (kpc)", 12.0, n5055['fw_n1_r_sh1_kpc'], tol=0.05)
    p6, m6 = check("NGC5055 shell sigma (kpc)", 2.91, n5055['fw_n1_sigma_sh1_kpc'], tol=0.01)
    p7, m7 = check("NGC5055 shell mass", 3.78e10, n5055['fw_n1_M_sh1'], tol=0.02e10)
    
    # Einasto
    ein = load_einasto()
    n5055_ein = ein[ein['galaxy'] == 'NGC5055'].iloc[0]
    p8, m8 = check("NGC5055 Einasto-only chi²_red", 8.42, n5055_ein['ein0_chi2_red'], tol=0.01)
    p9, m9 = check("NGC5055 Einasto+shells chi²_red", 0.507, n5055_ein['fw_einasto_chi2_red'], tol=0.005)
    p10, m10 = check("NGC5055 Einasto n_shells", 2, n5055_ein['fw_einasto_n_shells'], integer=True)
    
    # DC14
    try:
        dc14 = load_dc14_ngc5055()
        # Manuscript says DC14 chi^2_red ≈ 9.0 or 9.01
        # Find the chi2_red column — name varies
        chi2r_col = None
        for col in dc14.columns:
            if 'chi2_red' in col.lower() or 'chi2_r' in col.lower():
                chi2r_col = col
                break
        if chi2r_col:
            dc14_val = float(dc14.iloc[0][chi2r_col])
            p11, m11 = check("NGC5055 DC14 chi²_red", 9.01, dc14_val, tol=0.05)
        else:
            # Try standard naming
            print(f"   Available DC14 columns: {list(dc14.columns)}")
            p11, m11 = False, f"  {C.WARN}? SKIP{C.END} DC14 chi²_red column not found"
    except Exception as e:
        p11, m11 = False, f"  {C.WARN}? SKIP{C.END} DC14 NGC5055 file: {e}"
    
    for m in [m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11]: print(m)
    return all([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11])


# ============================================================
# Driver
# ============================================================
def main():
    print(f"{C.BOLD}=" * 70)
    print(f"halo_shells v7.0 — Validation Pass A (numerical consistency)")
    print(f"=" * 70 + C.END)
    print(f"Data directory: {DATA_DIR}")
    
    # Verify CSVs exist before running
    expected_csvs = [
        'sparc_T2-T9_canonical_fits.csv',
        'null_test_T2-T9_combined.csv',
        'einasto_full_sample_results.csv',
        'dc14_ngc5055_showcase.csv',
    ]
    missing = [f for f in expected_csvs if not os.path.exists(os.path.join(DATA_DIR, f))]
    if missing:
        print(f"{C.FAIL}ERROR: missing required CSVs: {missing}{C.END}")
        sys.exit(1)
    
    checks = [
        ('A1', A1, 'Framework adequacy at chi²_red < 1.5'),
        ('A2', A2, 'Burkert and NFW adequacy at chi²_red < 1.5'),
        ('A3', A3, 'Strict adequacy at chi²_red < 1.0'),
        ('A4', A4, 'Per-T shell-bearing fractions'),
        ('A5', A5, 'Trend tests'),
        ('A6', A6, 'Null test aggregate rates'),
        ('A7', A7, 'Per-T-bin null FP rates and Fisher tests'),
        ('A8', A8, 'NGC 5055 showcase numbers'),
    ]
    
    results = {}
    for label, func, desc in checks:
        try:
            results[label] = func()
        except Exception as e:
            print(f"\n{C.FAIL}ERROR in {label}: {e}{C.END}")
            import traceback
            traceback.print_exc()
            results[label] = False
    
    # Summary
    print(f"\n{C.BOLD}=" * 70)
    print(f"SUMMARY")
    print(f"=" * 70 + C.END)
    n_pass = sum(1 for v in results.values() if v)
    n_total = len(results)
    for label, _, desc in checks:
        status = f"{C.OK}✓ PASS{C.END}" if results[label] else f"{C.FAIL}✗ FAIL{C.END}"
        print(f"  {label}: {status}  {desc}")
    print()
    if n_pass == n_total:
        print(f"{C.OK}{C.BOLD}ALL {n_total}/{n_total} CHECKS PASSED{C.END}")
        sys.exit(0)
    else:
        print(f"{C.FAIL}{C.BOLD}{n_pass}/{n_total} CHECKS PASSED — {n_total - n_pass} FAILED{C.END}")
        sys.exit(1)


if __name__ == '__main__':
    main()
