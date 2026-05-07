"""
null_test.py — Synthetic smooth-halo null test (v7.0)

Tests whether the framework's BIC-based shell selection produces false-positive
shells when fed Gaussian-noise mocks generated from a smooth (Burkert or NFW)
halo. If shell selection happens rarely on smooth-halo mocks, the morphology
gradient observed on real data cannot be dismissed as overfitting.

This is the v7.0 revision — uses the strict sigma/r_shell <= 0.4 constraint
matching run_canonical_fits.py, via the sigma_frac reparameterization. The
v6.5 version of this script used the looser sigma <= r_max * 0.4 bound.

Procedure:
  1. For each real SPARC galaxy in T=2..7, take its actual radial sampling
     and errorbars
  2. Use the galaxy's BEST-FIT Burkert (or NFW) parameters from the canonical
     fits CSV
  3. Generate mock V_obs by sampling Gaussian noise around V_smooth
  4. Fit the framework (Burkert + 0/1/2 BIC shells, sigma/r<=0.4 strict) to
     the mock
  5. Record whether BIC selects shells (and how many)
  6. Repeat for many noise realizations
  7. Output: false-positive shell selection rate, by smooth-truth and T-type

Usage:
  Place this script in a directory with:
    - Rotmod_LTG/                                  (SPARC rotmod files)
    - sparc_T2-T9_canonical_fits.csv               (output of run_canonical_fits.py)
  
  Then:
    python3 null_test.py
  
  Output:
    - null_test_results.csv  (one row per (galaxy, smooth-truth, realization))
    - Console summary

  Expected runtime: 10-30 minutes (depends on N_REALIZATIONS).

Requires: numpy, scipy, pandas. No other dependencies.
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.special import erf
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
# Path resolution: prefer package layout (PACKAGE_ROOT/Rotmod_LTG/),
# fall back to current directory (./Rotmod_LTG/) for legacy v6.5 layouts.
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_PACKAGE_DATA = _os.path.join(_os.path.dirname(_HERE), 'Rotmod_LTG')
_LOCAL_DATA = _os.path.join(_os.getcwd(), 'Rotmod_LTG')
DATA_DIR = _PACKAGE_DATA if _os.path.isdir(_PACKAGE_DATA) else _LOCAL_DATA
_PACKAGE_CSV = _os.path.join(_os.path.dirname(_HERE), 'data', 'sparc_T2-T9_canonical_fits.csv')
_LOCAL_CSV = _os.path.join(_os.getcwd(), 'sparc_T2-T9_canonical_fits.csv')
CANONICAL_CSV = _PACKAGE_CSV if _os.path.isfile(_PACKAGE_CSV) else _LOCAL_CSV
OUTPUT_CSV = './null_test_results.csv'

G = 4.302e-6   # kpc · (km/s)^2 / M_sun
SIGMA_V_FLOOR = 1.0
UPSILON_DISK = 0.5
UPSILON_BULGE = 0.7

SHELL_R_MAX_GRID = [3.0, 6.0, 12.0]   # multi-restart upper bounds on r_shell
SHELL_WIDTH_MAX_FRAC = 0.4             # sigma_frac upper bound (strict)

# T-bins to test in the main null run; extension script handles T=8/9
SAMPLE_PER_T = {2: 4, 3: 4, 4: 5, 5: 5, 6: 5, 7: 5}
N_REALIZATIONS = 5
RNG_SEED = 12345


# ============================================================
# Physics
# ============================================================
def baryon_squared(vg, vd, vb):
    return (vg * np.abs(vg)
            + UPSILON_DISK * vd * np.abs(vd)
            + UPSILON_BULGE * vb * np.abs(vb))


def burkert_v(r, rho_0, a):
    a = max(a, 1e-6)
    x = r / a
    M = np.pi * rho_0 * a**3 * (np.log(1+x**2) + 2*np.log(1+x) - 2*np.arctan(x))
    return np.sqrt(np.maximum(G * M / np.maximum(r, 1e-6), 0))


def nfw_v(r, rho_s, r_s):
    r_s = max(r_s, 1e-6)
    x = r / r_s
    M = 4*np.pi * rho_s * r_s**3 * (np.log(1+x) - x/(1+x))
    return np.sqrt(np.maximum(G * M / np.maximum(r, 1e-6), 0))


def shell_v2(r, M, r_sh, sigma):
    sigma = max(sigma, 1e-6)
    M_enc = 0.5 * M * (1 + erf((r - r_sh) / (np.sqrt(2) * sigma)))
    return np.maximum(G * M_enc / np.maximum(r, 1e-6), 0)


# ============================================================
# Models (curve_fit compatible)
# ============================================================
def model_fw0(V2_bar):
    def vt(r, rho_0, a):
        v2 = V2_bar + burkert_v(r, rho_0, a)**2
        return np.sqrt(np.maximum(v2, 0))
    return vt


def model_fw1(V2_bar):
    """Framework + 1 shell. sigma_frac = sigma/r_shell, bounded [0.01, 0.4]."""
    def vt(r, rho_0, a, M_sh, r_sh, sigma_frac):
        sigma = sigma_frac * r_sh
        v2 = V2_bar + burkert_v(r, rho_0, a)**2 + shell_v2(r, M_sh, r_sh, sigma)
        return np.sqrt(np.maximum(v2, 0))
    return vt


def model_fw2(V2_bar):
    """Framework + 2 shells. Each shell uses sigma_frac in [0.01, 0.4]."""
    def vt(r, rho_0, a, M1, r1, sf1, M2, r2, sf2):
        s1 = sf1 * r1
        s2 = sf2 * r2
        v2 = (V2_bar
              + burkert_v(r, rho_0, a)**2
              + shell_v2(r, M1, r1, s1)
              + shell_v2(r, M2, r2, s2))
        return np.sqrt(np.maximum(v2, 0))
    return vt


# ============================================================
# Fitters (matches run_canonical_fits.py conventions)
# ============================================================
def fit_fw_n_shells(rd, vd, ed, V2_bar, n_shells, time_budget=30):
    sig = np.maximum(ed, SIGMA_V_FLOOR)
    
    if n_shells == 0:
        vt = model_fw0(V2_bar)
        starts = [[3e7, 8.0], [1e7, 15.0], [1e8, 4.0], [3e6, 30.0], [3e8, 2.0]]
        bounds = ([1e3, 0.1], [1e10, 200.0])
        best_chi2, best_p = np.inf, None
        t0 = time.time()
        for p0 in starts:
            try:
                p, _ = curve_fit(vt, rd, vd, p0=p0, bounds=bounds, sigma=sig, maxfev=10000)
                c = float(np.sum(((vd - vt(rd, *p)) / sig)**2))
                if c < best_chi2:
                    best_chi2, best_p = c, p
            except Exception:
                continue
            if time.time() - t0 > time_budget:
                break
        return best_p, best_chi2
    
    # n_shells = 1 or 2
    if n_shells == 1:
        vt = model_fw1(V2_bar)
        starts = []
        for r_max in SHELL_R_MAX_GRID:
            for rho0 in [1e7, 3e7, 1e8]:
                for a_init in [5.0, 10.0, 20.0]:
                    for r_frac in [0.3, 0.5, 0.7]:
                        starts.append(([rho0, a_init, 3e9, r_frac*r_max, 0.25], r_max))
    else:  # n_shells == 2
        vt = model_fw2(V2_bar)
        starts = []
        for r_max in SHELL_R_MAX_GRID:
            for rho0 in [1e7, 3e7, 1e8]:
                for a_init in [5.0, 10.0, 20.0]:
                    starts.append(([rho0, a_init, 3e9, r_max/3, 0.25, 3e9, 2*r_max/3, 0.25], r_max))
    
    best_chi2, best_p_physical = np.inf, None
    t0 = time.time()
    for p0, r_max in starts:
        if n_shells == 1:
            lb = [1e3, 0.1, 1e6, 0.2, 0.01]
            ub = [1e10, 200.0, 5e10, r_max, SHELL_WIDTH_MAX_FRAC]
        else:
            lb = [1e3, 0.1, 1e6, 0.2, 0.01, 1e6, 0.2, 0.01]
            ub = [1e10, 200.0, 5e10, r_max, SHELL_WIDTH_MAX_FRAC,
                  5e10, r_max, SHELL_WIDTH_MAX_FRAC]
        try:
            p, _ = curve_fit(vt, rd, vd, p0=p0, bounds=(lb, ub), sigma=sig, maxfev=15000)
            c = float(np.sum(((vd - vt(rd, *p)) / sig)**2))
            if c < best_chi2:
                best_chi2 = c
                if n_shells == 1:
                    rho0, a, M, rsh, sf = p
                    best_p_physical = np.array([rho0, a, M, rsh, sf*rsh])
                else:
                    rho0, a, M1, r1, sf1, M2, r2, sf2 = p
                    best_p_physical = np.array([rho0, a, M1, r1, sf1*r1, M2, r2, sf2*r2])
        except Exception:
            continue
        if time.time() - t0 > time_budget:
            break
    return best_p_physical, best_chi2


def best_n_via_bic(rd, vd, ed, V2_bar):
    """Return (best_n, best_chi2, best_bic) using BIC selection over n=0,1,2."""
    n_pts = len(rd)
    best_n, best_bic, best_chi2 = 0, np.inf, np.inf
    for ns in [0, 1, 2]:
        p, c = fit_fw_n_shells(rd, vd, ed, V2_bar, ns)
        if p is None:
            continue
        k = 2 + 3*ns   # 2 / 5 / 8
        bic = c + k * np.log(n_pts)
        if bic < best_bic:
            best_n, best_bic, best_chi2 = ns, bic, c
    return best_n, best_chi2, best_bic


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 72)
    print("NULL TEST (v7.0, strict sigma/r ≤ 0.4)")
    print("=" * 72)
    
    if not os.path.isdir(DATA_DIR):
        print(f"ERROR: {DATA_DIR} not found")
        sys.exit(1)
    if not os.path.exists(CANONICAL_CSV):
        print(f"ERROR: {CANONICAL_CSV} not found")
        sys.exit(1)
    
    canon = pd.read_csv(CANONICAL_CSV)
    
    # Stratified sample by T
    galaxies_to_test = []
    rng_sample = np.random.default_rng(seed=RNG_SEED)
    for T, n_sample in SAMPLE_PER_T.items():
        sub_T = canon[canon['T'] == T]
        if len(sub_T) == 0:
            continue
        n = min(n_sample, len(sub_T))
        # Use deterministic sampling for reproducibility
        sub = sub_T.sample(n=n, random_state=42)
        galaxies_to_test.extend(sub.to_dict('records'))
    
    print(f"\nSampling {len(galaxies_to_test)} galaxies, "
          f"{N_REALIZATIONS} realizations × 2 smooth models = "
          f"{len(galaxies_to_test) * 2 * N_REALIZATIONS} total mock fits")
    
    rng = np.random.default_rng(seed=RNG_SEED)
    results = []
    t_start = time.time()
    
    for i, gal in enumerate(galaxies_to_test):
        g = gal['Galaxy']
        path = os.path.join(DATA_DIR, f'{g}_rotmod.dat')
        if not os.path.exists(path):
            print(f"  [SKIP] {g}: rotmod file not found")
            continue
        
        d = np.loadtxt(path, comments='#')
        r, vobs, evobs = d[:,0], d[:,1], d[:,2]
        vgas, vdisk, vbul = d[:,3], d[:,4], d[:,5]
        Vbar2 = baryon_squared(vgas, vdisk, vbul)
        mask = vobs**2 > Vbar2
        n_pts = int(mask.sum())
        if n_pts < 5:
            continue
        
        rd = r[mask]
        ed = evobs[mask]
        V2_bd = Vbar2[mask]
        
        # Smooth predictions
        V_burk_smooth = burkert_v(rd, gal['burk_rho0'], gal['burk_a_kpc'])
        V_burk_total = np.sqrt(np.maximum(V2_bd + V_burk_smooth**2, 0))
        
        V_nfw_smooth = nfw_v(rd, gal['nfw_rho_s'], gal['nfw_r_s_kpc'])
        V_nfw_total = np.sqrt(np.maximum(V2_bd + V_nfw_smooth**2, 0))
        
        for smooth_label, V_truth in [('Burkert', V_burk_total), ('NFW', V_nfw_total)]:
            for k in range(N_REALIZATIONS):
                noise = rng.normal(0, np.maximum(ed, SIGMA_V_FLOOR))
                v_mock = np.maximum(V_truth + noise, 0.1)
                best_n, best_chi2, best_bic = best_n_via_bic(rd, v_mock, ed, V2_bd)
                results.append({
                    'Galaxy': g,
                    'T': int(gal['T']),
                    'V_flat': float(gal['V_flat']),
                    'smooth_truth': smooth_label,
                    'realization': k,
                    'n_pts': n_pts,
                    'best_n_shells': best_n,
                    'best_chi2': best_chi2,
                    'best_bic': best_bic,
                })
        
        elapsed = time.time() - t_start
        eta = elapsed / (i+1) * (len(galaxies_to_test) - i - 1)
        df_so_far = pd.DataFrame(results)
        rate = (df_so_far['best_n_shells'] > 0).mean() * 100 if len(df_so_far) > 0 else 0
        print(f"  [{i+1}/{len(galaxies_to_test)}] {g:<12} T={int(gal['T'])} n_pts={n_pts}  "
              f"shell-rate-so-far={rate:.1f}%  ({elapsed:.0f}s, ETA {eta:.0f}s)")
    
    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSaved: {OUTPUT_CSV}")
    
    # Summary
    print("\n" + "="*70)
    print("NULL TEST SUMMARY (v7.0)")
    print("="*70)
    print(f"\nTotal mock fits: {len(df)}")
    print(f"\nFalse shell selection rate (overall):")
    print(f"  Selected ≥1 shell:  {(df['best_n_shells']>0).sum()}/{len(df)} ({100*(df['best_n_shells']>0).mean():.1f}%)")
    print(f"  Selected 1 shell:   {(df['best_n_shells']==1).sum()}/{len(df)} ({100*(df['best_n_shells']==1).mean():.1f}%)")
    print(f"  Selected 2 shells:  {(df['best_n_shells']==2).sum()}/{len(df)} ({100*(df['best_n_shells']==2).mean():.1f}%)")
    
    print(f"\nBy smooth-halo truth:")
    for label in ['Burkert', 'NFW']:
        sub = df[df['smooth_truth'] == label]
        rate = (sub['best_n_shells'] > 0).mean() * 100
        print(f"  {label:<8}: shell-bearing fraction in mocks = {rate:.1f}% (N={len(sub)})")
    
    print(f"\nBy T-type (Burkert-truth mocks only):")
    for T in [2,3,4,5,6,7]:
        sub = df[(df['smooth_truth'] == 'Burkert') & (df['T'] == T)]
        if len(sub) > 0:
            rate = (sub['best_n_shells'] > 0).mean() * 100
            print(f"  T={T}: false-positive shell rate = {rate:.1f}% (N={len(sub)})")


if __name__ == '__main__':
    main()
