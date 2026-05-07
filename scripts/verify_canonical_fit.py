"""
verify_canonical_fit.py — Single-galaxy verification for the canonical SPARC
                          framework fits.

Reproduces the columns of sparc_T2-T9_canonical_fits.csv for one galaxy by
re-fitting from raw SPARC rotmod data, then prints a side-by-side comparison
against the shipped CSV. Intended as a Category-C spot-check verifier.

Usage:
  Place this script in a directory containing:
    - Rotmod_LTG/                       (SPARC rotmod .dat files)
    - sparc_T2-T9_canonical_fits.csv    (the canonical CSV to verify against)
    - sparc_sample123.csv               (SPARC catalog with T, V_flat etc.)

  Then:
    python3 verify_canonical_fit.py NGC5055
    python3 verify_canonical_fit.py NGC2403
    python3 verify_canonical_fit.py UGC06787

  Output: comparison table and pass/fail flags per quantity.

Conventions (reverse-engineered from the canonical CSV; any differences from
the original producer script may surface in the spot-check as discrepancies):

  - Mass-to-light: Upsilon_disk = 0.5, Upsilon_bulge = 0.7 (SPARC default).
  - Exclusion rule: drop points with V_obs^2 <= V_bar^2.
  - Sigma_V floor: 1 km/s.
  - BIC formula: chi^2 + k * ln(n_pts_used), with k = 2 / 2 / 5 / 8 for
    Burkert / NFW / FW 1-shell / FW 2-shell.
  - chi^2_red = chi^2 / max(n_pts_used - k, 1).
  - is_clean = (n_excluded == 0).
  - BIC selection: argmin over n0/n1/n2 BIC values.

Parameter bounds (reverse-engineered):

  - Burkert: rho_0 in [1e3, 1e11] M_sun/kpc^3; a in [1.0, 200] kpc.
  - NFW (free-c): rho_s in [1e3, 1e11] M_sun/kpc^3; r_s in [0.1, 500] kpc.
  - Shell M in [1e6, 5e10] M_sun; r in [0.2, 12.0] kpc; sigma in [0.05, ?] kpc
    with the structural constraint sigma/r <= 0.4 enforced.

Multi-restart starts: matches dc14_universal_failures.py (same producer set):
  rho ∈ {1e6, 1e7, 5e7, 1e8, 1e9}, length-scale ∈ {3, 8, 15, 30, 50}.

Requires: numpy, scipy, pandas. No other dependencies.
"""

import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
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
SAMPLE_CSV = './sparc_sample123.csv'

G = 4.302e-6   # kpc (km/s)^2 / M_sun
SIGMA_V_FLOOR = 1.0   # km/s
UPSILON_DISK = 0.5
UPSILON_BULGE = 0.7

# Multi-restart starting grids
RHO_STARTS = [1e6, 1e7, 5e7, 1e8, 1e9]
LENGTH_STARTS = [3.0, 8.0, 15.0, 30.0, 50.0]
SHELL_R_STARTS = [1.0, 3.0, 6.0, 9.0]      # within [0.2, 12.0] bound
SHELL_M_STARTS = [1e8, 1e9, 1e10]           # within [1e6, 5e10] bound
SHELL_SIGMA_FRAC_STARTS = [0.1, 0.25, 0.4]  # as fraction of r (within sigma/r ≤ 0.4)

TIME_BUDGET_PER_FIT = 30   # seconds; FW 2-shell can be slow; cap to keep verifier fast

# Parameter bounds
BURK_BOUNDS = ([1e3, 1.0], [1e11, 200.0])              # (rho_0, a)
NFW_BOUNDS  = ([1e3, 0.1], [1e11, 500.0])              # (rho_s, r_s)


# ============================================================
# Physics
# ============================================================
def baryon_squared(vg, vd, vb):
    """SPARC signed-square baryonic V^2 (km/s)^2."""
    return (vg * np.abs(vg)
            + UPSILON_DISK  * vd * np.abs(vd)
            + UPSILON_BULGE * vb * np.abs(vb))


def burkert_v2(r, rho_0, a):
    """Burkert profile circular velocity squared. Returns V^2 in (km/s)^2."""
    x = r / max(a, 1e-6)
    M_enc = np.pi * rho_0 * a**3 * (
        np.log(1 + x**2) + 2 * np.log(1 + x) - 2 * np.arctan(x)
    )
    v2 = G * M_enc / np.maximum(r, 1e-6)
    return np.maximum(v2, 0)


def nfw_v2(r, rho_s, r_s):
    """NFW profile circular velocity squared (free-c)."""
    x = r / max(r_s, 1e-6)
    M_enc = 4 * np.pi * rho_s * r_s**3 * (
        np.log(1 + x) - x / (1 + x)
    )
    v2 = G * M_enc / np.maximum(r, 1e-6)
    return np.maximum(v2, 0)


def shell_v2(r, M_sh, r_sh, sigma_sh):
    """
    Gaussian-shell circular velocity squared. Implementation: shell mass M
    distributed in a Gaussian shell of mean radius r_sh and width sigma_sh.
    Enclosed mass at radius r is M times the cumulative Gaussian to r.

    Approximation: treat the Gaussian shell as adding mass with cumulative
    M_enc(r) = M * 0.5 * (1 + erf((r - r_sh) / (sqrt(2) * sigma_sh))).
    This is the standard analytic form used for Gaussian-shell components.
    """
    from scipy.special import erf
    sigma = max(sigma_sh, 1e-6)
    M_enc = 0.5 * M_sh * (1 + erf((r - r_sh) / (np.sqrt(2) * sigma)))
    v2 = G * M_enc / np.maximum(r, 1e-6)
    return np.maximum(v2, 0)


# ============================================================
# Models
# ============================================================
def model_burkert(r, V2_bar, rho_0, a):
    """Total V (km/s) = sqrt(V_bar^2 + V_burkert^2)."""
    return np.sqrt(V2_bar + burkert_v2(r, rho_0, a))


def model_nfw(r, V2_bar, rho_s, r_s):
    return np.sqrt(V2_bar + nfw_v2(r, rho_s, r_s))


def model_fw1(r, V2_bar, rho_0, a, M_sh, r_sh, sigma_sh):
    return np.sqrt(V2_bar + burkert_v2(r, rho_0, a) + shell_v2(r, M_sh, r_sh, sigma_sh))


def model_fw2(r, V2_bar, rho_0, a, M1, r1, s1, M2, r2, s2):
    return np.sqrt(V2_bar
                   + burkert_v2(r, rho_0, a)
                   + shell_v2(r, M1, r1, s1)
                   + shell_v2(r, M2, r2, s2))


# ============================================================
# Fitters
# ============================================================
def chi2_of(model_fn, params, r, vobs, sig, V2_bar):
    """Compute chi^2 for a model evaluated at (r, V2_bar) with given params."""
    pred = model_fn(r, V2_bar, *params)
    return float(np.sum(((vobs - pred) / sig)**2))


def fit_burkert(r, vobs, sig, V2_bar, time_budget=TIME_BUDGET_PER_FIT):
    """Multi-restart Burkert fit. Returns (best_params, best_chi2)."""
    def model(r, rho_0, a):
        return model_burkert(r, V2_bar, rho_0, a)

    starts = [[rho, a] for rho in RHO_STARTS for a in LENGTH_STARTS]
    best_chi2 = np.inf
    best_p = None
    t0 = time.time()
    for p0 in starts:
        try:
            p, _ = curve_fit(model, r, vobs, p0=p0, bounds=BURK_BOUNDS,
                             sigma=sig, maxfev=15000)
            chi2 = float(np.sum(((vobs - model(r, *p))/sig)**2))
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_p = p
        except Exception:
            continue
        if time.time() - t0 > time_budget:
            break
    return best_p, best_chi2


def fit_nfw(r, vobs, sig, V2_bar, time_budget=TIME_BUDGET_PER_FIT):
    def model(r, rho_s, r_s):
        return model_nfw(r, V2_bar, rho_s, r_s)
    starts = [[rho, rs] for rho in RHO_STARTS for rs in LENGTH_STARTS]
    best_chi2 = np.inf
    best_p = None
    t0 = time.time()
    for p0 in starts:
        try:
            p, _ = curve_fit(model, r, vobs, p0=p0, bounds=NFW_BOUNDS,
                             sigma=sig, maxfev=15000)
            chi2 = float(np.sum(((vobs - model(r, *p))/sig)**2))
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_p = p
        except Exception:
            continue
        if time.time() - t0 > time_budget:
            break
    return best_p, best_chi2


def fit_fw1(r, vobs, sig, V2_bar, time_budget=TIME_BUDGET_PER_FIT):
    """Fit Burkert + 1 shell. Bounds: shell M ∈ [1e6, 5e10], r ∈ [0.2, 12], σ ∈ [0.05, 0.4*r]."""
    # We enforce sigma <= 0.4 * r by mapping sigma -> sigma_frac (in [0, 0.4]) then
    # multiplying by r. But curve_fit needs box bounds. Easiest: bound sigma to a
    # global max of, say, 5 kpc, then post-filter / penalize sigma > 0.4*r. The
    # canonical CSV shows max sigma = 4.69, so 5 is a safe cap. The constraint
    # is enforced via residual penalty inside the model.
    SIGMA_MAX = 5.0   # global cap; effective constraint applied via penalty
    SH_BOUNDS = ([1e3, 1.0,   1e6, 0.2, 0.05],
                 [1e11, 200., 5e10, 12.0, SIGMA_MAX])

    def model(r, rho_0, a, M_sh, r_sh, sigma_sh):
        v = model_fw1(r, V2_bar, rho_0, a, M_sh, r_sh, sigma_sh)
        # Soft penalty for sigma/r > 0.4 (force to upper bound)
        if sigma_sh > 0.4 * r_sh:
            penalty = 1e6 * (sigma_sh / r_sh - 0.4)**2
            v = v + penalty
        return v

    starts = []
    for rho in [1e7, 1e8]:
        for a_l in [3.0, 15.0]:
            for r_sh in SHELL_R_STARTS:
                for M_sh in SHELL_M_STARTS:
                    for sf in SHELL_SIGMA_FRAC_STARTS:
                        starts.append([rho, a_l, M_sh, r_sh, sf * r_sh])

    best_chi2 = np.inf
    best_p = None
    t0 = time.time()
    for p0 in starts:
        try:
            p, _ = curve_fit(model, r, vobs, p0=p0, bounds=SH_BOUNDS,
                             sigma=sig, maxfev=15000)
            # Enforce sigma/r <= 0.4 by clipping if needed and recomputing
            if p[4] > 0.4 * p[3]:
                p[4] = 0.4 * p[3]
            chi2 = chi2_of(model_fw1, p, r, vobs, sig, V2_bar)
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_p = p
        except Exception:
            continue
        if time.time() - t0 > time_budget:
            break
    return best_p, best_chi2


def fit_fw2(r, vobs, sig, V2_bar, time_budget=TIME_BUDGET_PER_FIT * 2):
    """Fit Burkert + 2 shells. 8 free parameters."""
    SIGMA_MAX = 5.0
    SH2_BOUNDS = ([1e3, 1.0,   1e6, 0.2, 0.05,  1e6, 0.2, 0.05],
                  [1e11, 200., 5e10, 12.0, SIGMA_MAX, 5e10, 12.0, SIGMA_MAX])

    def model(r, rho_0, a, M1, r1, s1, M2, r2, s2):
        v = model_fw2(r, V2_bar, rho_0, a, M1, r1, s1, M2, r2, s2)
        # Soft penalties
        pen = 0
        if s1 > 0.4 * r1: pen += 1e6 * (s1/r1 - 0.4)**2
        if s2 > 0.4 * r2: pen += 1e6 * (s2/r2 - 0.4)**2
        return v + pen

    # Coarser starting grid for 2-shell to keep runtime manageable, but
    # with enough diversity to catch most local optima.
    starts = []
    for rho in [1e7, 5e7, 1e8]:
        for a_l in [3.0, 8.0, 15.0]:
            for r1 in [1.0, 3.0, 6.0]:
                for r2 in [4.0, 8.0, 11.0]:
                    if r2 <= r1:
                        continue
                    for M_pair in [(1e8, 1e9), (1e9, 5e9), (5e9, 1e10), (1e10, 5e10)]:
                        starts.append([rho, a_l, M_pair[0], r1, 0.25*r1, M_pair[1], r2, 0.25*r2])

    best_chi2 = np.inf
    best_p = None
    t0 = time.time()
    for p0 in starts:
        try:
            p, _ = curve_fit(model, r, vobs, p0=p0, bounds=SH2_BOUNDS,
                             sigma=sig, maxfev=15000)
            if p[4] > 0.4 * p[3]: p[4] = 0.4 * p[3]
            if p[7] > 0.4 * p[6]: p[7] = 0.4 * p[6]
            chi2 = chi2_of(model_fw2, p, r, vobs, sig, V2_bar)
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_p = p
        except Exception:
            continue
        if time.time() - t0 > time_budget:
            break
    return best_p, best_chi2


# ============================================================
# Verification helpers
# ============================================================
def chi2_red_of(chi2, n_pts_used, k):
    """Reproduce the canonical CSV's chi^2_red formula."""
    return chi2 / max(n_pts_used - k, 1)


def bic_of(chi2, n_pts_used, k):
    return chi2 + k * np.log(n_pts_used)


def cmp(actual, expected, tol_rel=0.05, tol_abs=None):
    """Compare two numbers, return formatted pass/fail string."""
    if expected == 0:
        diff = abs(actual)
        ok = diff < (tol_abs if tol_abs is not None else 1e-6)
    else:
        diff = abs(actual - expected) / abs(expected)
        ok = diff < tol_rel
    if tol_abs is not None and abs(actual - expected) < tol_abs:
        ok = True
    flag = "PASS" if ok else "DIFF"
    return f"[{flag}] act={actual:.6g}  exp={expected:.6g}  rel_diff={diff:.4f}"


# ============================================================
# Main
# ============================================================
def verify_galaxy(galaxy_name):
    print("=" * 78)
    print(f"VERIFYING: {galaxy_name}")
    print("=" * 78)

    # ---- Load data ----
    rotmod_path = os.path.join(DATA_DIR, f'{galaxy_name}_rotmod.dat')
    if not os.path.exists(rotmod_path):
        print(f"ERROR: {rotmod_path} not found")
        sys.exit(1)

    if not os.path.exists(CANONICAL_CSV):
        print(f"ERROR: {CANONICAL_CSV} not found")
        sys.exit(1)

    df_can = pd.read_csv(CANONICAL_CSV)
    expected = df_can[df_can['Galaxy'] == galaxy_name]
    if len(expected) == 0:
        print(f"ERROR: galaxy {galaxy_name} not in canonical CSV")
        sys.exit(1)
    expected = expected.iloc[0]

    d = np.loadtxt(rotmod_path, comments='#')
    r_all, vobs_all, evobs_all = d[:, 0], d[:, 1], d[:, 2]
    vgas, vdisk, vbul = d[:, 3], d[:, 4], d[:, 5]
    Vbar2_all = baryon_squared(vgas, vdisk, vbul)

    n_total = len(r_all)
    mask = vobs_all**2 > Vbar2_all
    n_used = int(mask.sum())
    n_excluded = n_total - n_used

    r = r_all[mask]
    vobs = vobs_all[mask]
    evobs = evobs_all[mask]
    Vbar2 = Vbar2_all[mask]
    sig = np.maximum(evobs, SIGMA_V_FLOOR)

    print(f"\n--- Data point counts ---")
    print(f"  n_pts_total: {n_total}     " + cmp(n_total, expected['n_pts_total']))
    print(f"  n_pts_used:  {n_used}      " + cmp(n_used, expected['n_pts_used']))
    print(f"  n_excluded:  {n_excluded}  " + cmp(n_excluded, expected['n_excluded']))

    # ---- Burkert ----
    print(f"\n--- Burkert (k=2) ---")
    t0 = time.time()
    p_burk, chi2_burk = fit_burkert(r, vobs, sig, Vbar2)
    elapsed = time.time() - t0
    if p_burk is None:
        print("  FAILED")
        return
    chi2r_burk = chi2_red_of(chi2_burk, n_used, 2)
    bic_burk = bic_of(chi2_burk, n_used, 2)
    print(f"  Fit time: {elapsed:.1f}s")
    print(f"  rho_0:     " + cmp(p_burk[0], expected['burk_rho0'], tol_rel=0.10))
    print(f"  a:         " + cmp(p_burk[1], expected['burk_a_kpc'], tol_rel=0.10))
    print(f"  chi^2:     " + cmp(chi2_burk, expected['burk_chi2'], tol_rel=0.05))
    print(f"  chi^2_red: " + cmp(chi2r_burk, expected['burk_chi2_red'], tol_rel=0.05))
    print(f"  BIC:       " + cmp(bic_burk, expected['burk_bic'], tol_rel=0.02))

    # ---- NFW ----
    print(f"\n--- NFW (k=2) ---")
    t0 = time.time()
    p_nfw, chi2_nfw = fit_nfw(r, vobs, sig, Vbar2)
    elapsed = time.time() - t0
    chi2r_nfw = chi2_red_of(chi2_nfw, n_used, 2)
    bic_nfw = bic_of(chi2_nfw, n_used, 2)
    print(f"  Fit time: {elapsed:.1f}s")
    print(f"  rho_s:     " + cmp(p_nfw[0], expected['nfw_rho_s'], tol_rel=0.10))
    print(f"  r_s:       " + cmp(p_nfw[1], expected['nfw_r_s_kpc'], tol_rel=0.10))
    print(f"  chi^2:     " + cmp(chi2_nfw, expected['nfw_chi2'], tol_rel=0.05))
    print(f"  chi^2_red: " + cmp(chi2r_nfw, expected['nfw_chi2_red'], tol_rel=0.05))
    print(f"  BIC:       " + cmp(bic_nfw, expected['nfw_bic'], tol_rel=0.02))

    # ---- Framework 1-shell ----
    print(f"\n--- Framework 1-shell (k=5) ---")
    t0 = time.time()
    p_fw1, chi2_fw1 = fit_fw1(r, vobs, sig, Vbar2)
    elapsed = time.time() - t0
    chi2r_fw1 = chi2_red_of(chi2_fw1, n_used, 5)
    bic_fw1 = bic_of(chi2_fw1, n_used, 5)
    print(f"  Fit time: {elapsed:.1f}s")
    if p_fw1 is not None:
        print(f"  shell M:    " + cmp(p_fw1[2], expected['fw_n1_M_sh1'], tol_rel=0.20))
        print(f"  shell r:    " + cmp(p_fw1[3], expected['fw_n1_r_sh1_kpc'], tol_rel=0.10))
        print(f"  shell σ:    " + cmp(p_fw1[4], expected['fw_n1_sigma_sh1_kpc'], tol_rel=0.20))
    print(f"  chi^2:     " + cmp(chi2_fw1, expected['fw_n1_chi2'], tol_rel=0.05))
    print(f"  chi^2_red: " + cmp(chi2r_fw1, expected['fw_n1_chi2_red'], tol_rel=0.05))
    print(f"  BIC:       " + cmp(bic_fw1, expected['fw_n1_bic'], tol_rel=0.02))

    # ---- Framework 2-shell ----
    print(f"\n--- Framework 2-shell (k=8) ---")
    t0 = time.time()
    p_fw2, chi2_fw2 = fit_fw2(r, vobs, sig, Vbar2)
    elapsed = time.time() - t0
    chi2r_fw2 = chi2_red_of(chi2_fw2, n_used, 8)
    bic_fw2 = bic_of(chi2_fw2, n_used, 8)
    print(f"  Fit time: {elapsed:.1f}s")
    print(f"  chi^2:     " + cmp(chi2_fw2, expected['fw_n2_chi2'], tol_rel=0.10))
    print(f"  chi^2_red: " + cmp(chi2r_fw2, expected['fw_n2_chi2_red'], tol_rel=0.10))
    print(f"  BIC:       " + cmp(bic_fw2, expected['fw_n2_bic'], tol_rel=0.05))

    # ---- BIC selection ----
    bics = {0: bic_burk, 1: bic_fw1, 2: bic_fw2}
    best_n = int(min(bics, key=bics.get))
    print(f"\n--- BIC selection ---")
    print(f"  best_n_shells: actual={best_n}  expected={int(expected['fw_best_n_shells'])}  "
          f"{'PASS' if best_n == int(expected['fw_best_n_shells']) else 'DIFF'}")

    # Summary
    print(f"\n{'=' * 78}")
    print(f"SUMMARY: see PASS/DIFF flags above. Brief acceptance criterion:")
    print(f"  - chi^2_red and BIC agree within ~5% → fit pipeline reproduces")
    print(f"  - Larger discrepancies → check multi-restart starting grid or")
    print(f"    parameter-bound details before declaring a real difference.")
    print(f"{'=' * 78}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('galaxy', help='Galaxy name, e.g. NGC5055')
    args = parser.parse_args()
    verify_galaxy(args.galaxy)
