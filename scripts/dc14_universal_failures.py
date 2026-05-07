"""
DC14 fits for the 10 universal-failure galaxies.

Tests whether DC14 (Di Cintio et al. 2014, feedback-modified halo profile) can
adequately fit (chi^2_red < 1.5) any of the galaxies where Burkert, NFW, and
the Burkert+shell framework all failed.

Usage:
  Place this script in a directory with a Rotmod_LTG/ subfolder containing
  the SPARC rotmod .dat files. Then:
    python3 dc14_universal_failures.py
  
  Output:
    - Console summary table
    - dc14_universal_failures_results.csv (machine-readable)
    - dc14_universal_failures_summary.txt (human-readable)
  
  Expected runtime: 5-10 minutes depending on machine.

Requires: numpy, scipy. No other dependencies.
"""

import os
import time
import csv
import numpy as np
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
# Path resolution: prefer package layout (PACKAGE_ROOT/Rotmod_LTG/),
# fall back to current directory (./Rotmod_LTG/) for legacy v6.5 layouts.
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_PACKAGE_DATA = _os.path.join(_os.path.dirname(_HERE), 'Rotmod_LTG')
_LOCAL_DATA = _os.path.join(_os.getcwd(), 'Rotmod_LTG')
DATA_DIR = _PACKAGE_DATA if _os.path.isdir(_PACKAGE_DATA) else _LOCAL_DATA

# Universal-failure galaxies from the paper (T=2-9, 10 galaxies)
# 9 SMBH-class galaxies in T=2-7 + UGC 00128 in T=8
UNIVERSAL_FAILURES = [
    ('NGC6674',  3, 241, 'SMBH-class'),
    ('UGC06787', 2, 248, 'SMBH-class'),
    ('UGC09133', 2, 227, 'SMBH-class'),
    ('NGC5907',  5, 215, 'SMBH-class'),
    ('NGC5033',  5, 194, 'SMBH-class'),
    ('NGC2841',  3, 285, 'SMBH-class'),
    ('UGC02953', 2, 265, 'SMBH-class'),
    ('NGC5985',  3, 294, 'SMBH-class'),
    ('NGC0289',  4, 163, 'SMBH-class'),
    ('UGC00128', 8,  73, 'late-type'),
]

# Multi-restart grid for DC14
RHO_S_STARTS = [1e6, 1e7, 5e7, 1e8, 1e9]
R_S_STARTS = [3.0, 8.0, 15.0, 30.0, 50.0]
X_STARTS = [-3.5, -3.0, -2.5, -2.0]
TIME_BUDGET_PER_GALAXY = 60   # seconds; stops restart loop if exceeded

G = 4.302e-6   # kpc (km/s)^2 / M_sun

# ------------------------------------------------------------------
# Physics
# ------------------------------------------------------------------
def baryon_squared(vg, vd, vb):
    """SPARC signed-square baryonic velocity. Standard Upsilon_disk=0.5, bulge=0.7."""
    return vg*np.abs(vg) + 0.5*vd*np.abs(vd) + 0.7*vb*np.abs(vb)


def dc14_shape_params(X):
    """
    DC14 alpha-beta-gamma profile parameters as functions of X = log10(M*/M_halo).
    From Di Cintio, Brook, Maccio et al. 2014.
    Valid range: X in [-4.1, -1.3]. Outside, behavior reduces to NFW.
    """
    Xc = np.clip(X, -4.1, -1.3)
    a = 2.94 - np.log10(10**((Xc + 2.33)*(-1.08)) + 10**((Xc + 2.33)*2.29))
    b = 4.23 + 1.34*Xc + 0.26*Xc**2
    g = -0.06 + np.log10(10**((Xc + 2.56)*(-0.68)) + 10**(Xc + 2.56))
    # Physical bounds
    a = np.clip(a, 0.5, 5.0)
    b = np.clip(b, 2.5, 6.0)
    g = np.clip(g, 0.0, 1.5)
    return a, b, g


def dc14_v_circular(r_grid_for_integ, r_eval, rho_s, r_s, X):
    """
    Compute DC14 circular velocity at r_eval points by integrating density on
    a fine logarithmically-spaced grid (cumulative trapezoid).
    """
    a, b, g = dc14_shape_params(X)
    x = r_grid_for_integ / max(r_s, 1e-3)
    with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
        rho = rho_s / (x**g * (1 + x**a)**((b-g)/a))
    rho = np.nan_to_num(rho, nan=0.0, posinf=0.0, neginf=0.0)
    integrand = 4 * np.pi * r_grid_for_integ**2 * rho
    M = np.zeros_like(r_grid_for_integ)
    M[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r_grid_for_integ))
    M_at_eval = np.interp(r_eval, r_grid_for_integ, M)
    v2 = G * M_at_eval / np.maximum(r_eval, 1e-3)
    return np.sqrt(np.maximum(v2, 0))


def fit_dc14_one_galaxy(rd, vd, ed, V2, time_budget=60):
    """
    Multi-restart DC14 fit. Returns (best_params, best_chi2, n_restarts_completed).
    """
    sig = np.maximum(ed, 1.0)
    r_grid = np.geomspace(0.01, max(rd.max()*2, 200), 500)

    def model(r, rho_s, r_s, X):
        v_dm = dc14_v_circular(r_grid, r, rho_s, r_s, X)
        return np.sqrt(np.maximum(V2 + v_dm**2, 0))

    starts = []
    for rho_s in RHO_S_STARTS:
        for r_s in R_S_STARTS:
            for X in X_STARTS:
                starts.append([rho_s, r_s, X])

    best_chi2 = np.inf
    best_p = None
    t0 = time.time()
    n_completed = 0

    for p0 in starts:
        try:
            p, _ = curve_fit(
                model, rd, vd, p0=p0,
                bounds=([1e3, 0.1, -4.0], [1e11, 200, -1.5]),
                sigma=sig, maxfev=10000
            )
            chi2 = float(np.sum(((vd - model(rd, *p))/sig)**2))
            n_completed += 1
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_p = p
        except Exception:
            continue
        if time.time() - t0 > time_budget:
            break

    return best_p, best_chi2, n_completed


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def main():
    print("="*72)
    print("DC14 fits for 10 universal-failure galaxies")
    print(f"Data directory: {DATA_DIR}")
    print(f"Time budget per galaxy: {TIME_BUDGET_PER_GALAXY}s")
    print("="*72)

    if not os.path.isdir(DATA_DIR):
        print(f"ERROR: Cannot find data directory '{DATA_DIR}'")
        print("Please adjust DATA_DIR at the top of this script.")
        return

    results = []
    t_total_start = time.time()

    for galaxy, T, V_flat, gclass in UNIVERSAL_FAILURES:
        path = os.path.join(DATA_DIR, f'{galaxy}_rotmod.dat')
        if not os.path.exists(path):
            print(f"\n  [SKIP] {galaxy}: file not found at {path}")
            continue

        d = np.loadtxt(path, comments='#')
        r, vobs, evobs = d[:,0], d[:,1], d[:,2]
        vgas, vdisk, vbul = d[:,3], d[:,4], d[:,5]
        Vbar2 = baryon_squared(vgas, vdisk, vbul)
        mask = vobs**2 > Vbar2
        n_pts = mask.sum()
        if n_pts < 5:
            print(f"\n  [SKIP] {galaxy}: only {n_pts} usable data points")
            continue
        rd = r[mask]; vd = vobs[mask]; ed = evobs[mask]; Vbd = Vbar2[mask]

        print(f"\n  Fitting {galaxy} (T={T}, V_flat={V_flat}, n_pts={n_pts})...")
        t0 = time.time()
        p, chi2, n_restarts = fit_dc14_one_galaxy(
            rd, vd, ed, Vbd, time_budget=TIME_BUDGET_PER_GALAXY
        )
        elapsed = time.time() - t0

        if p is None:
            print(f"    FAILED to converge in {elapsed:.0f}s, {n_restarts} restarts")
            continue

        n_free = 3
        chi2_red = chi2 / max(n_pts - n_free, 1)
        bic = chi2 + n_free * np.log(n_pts)
        a, b, gamma = dc14_shape_params(p[2])
        verdict = 'PASS' if chi2_red < 1.5 else 'FAIL'

        print(f"    Done in {elapsed:.0f}s, {n_restarts} restarts")
        print(f"    rho_s={p[0]:.2e}  r_s={p[1]:.1f}kpc  X={p[2]:.2f}  "
              f"(alpha={a:.2f}, beta={b:.2f}, gamma={gamma:.2f})")
        print(f"    chi^2={chi2:.1f}  chi^2_red={chi2_red:.2f}  BIC={bic:.1f}  [{verdict}]")

        results.append({
            'galaxy': galaxy, 'T': T, 'V_flat': V_flat, 'class': gclass,
            'n_pts': n_pts, 'n_restarts': n_restarts,
            'rho_s': float(p[0]), 'r_s': float(p[1]), 'X': float(p[2]),
            'alpha': float(a), 'beta': float(b), 'gamma': float(gamma),
            'chi2': float(chi2), 'chi2_red': float(chi2_red), 'bic': float(bic),
            'verdict': verdict,
            'elapsed_s': float(elapsed),
        })

    # -------- Summary --------
    total_time = time.time() - t_total_start
    n_pass = sum(1 for r in results if r['verdict'] == 'PASS')
    n_fail = sum(1 for r in results if r['verdict'] == 'FAIL')

    print("\n" + "="*72)
    print(f"SUMMARY: {n_pass} pass / {n_fail} fail (out of {len(results)} fits)")
    print(f"Total runtime: {total_time:.0f}s")
    print("="*72)
    print(f"{'Galaxy':<10} {'T':>3} {'V_flat':>7} {'n_pts':>6} "
          f"{'gamma':>6} {'chi^2_red':>10} {'verdict':>8}")
    print('-'*72)
    for r in results:
        print(f"{r['galaxy']:<10} {r['T']:>3} {r['V_flat']:>7} {r['n_pts']:>6} "
              f"{r['gamma']:>6.2f} {r['chi2_red']:>10.2f} {r['verdict']:>8}")

    print("\nReference values (from paper, T=2-9 sample):")
    print("  Burkert chi^2_red on these galaxies: all >= 1.5 (universal fail)")
    print("  NFW chi^2_red on these galaxies: all >= 1.5 (universal fail)")
    print("  Framework chi^2_red on these galaxies: all >= 1.5 (universal fail)")
    print()
    if n_pass == 0:
        print(">> RESULT: DC14 also fails on all universal-failure galaxies.")
        print("   This strengthens the paper's claim that no smooth profile,")
        print("   even feedback-modified, captures the structure these galaxies")
        print("   contain. Adds to defensibility of section 4.1 and 5.3.")
    elif n_pass <= 3:
        print(f">> RESULT: DC14 rescues {n_pass} of {len(results)} universal failures.")
        print("   Worth investigating which ones. Paper section 4.1 should be")
        print("   revised to acknowledge DC14 captures part (but not most) of")
        print("   the unmodeled structure.")
    else:
        print(f">> RESULT: DC14 rescues {n_pass} of {len(results)} universal failures.")
        print("   This is a meaningful finding. Paper sections 4.1, 5.3, and 5.4")
        print("   need revision to account for DC14's reach. Recommend re-running")
        print("   DC14 across the full 102-galaxy sample before submission.")

    # -------- Save results --------
    if results:
        csv_path = 'dc14_universal_failures_results.csv'
        with open(csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
            writer.writeheader()
            writer.writerows(results)
        print(f"\nSaved: {csv_path}")

        with open('dc14_universal_failures_summary.txt', 'w') as f:
            f.write(f"DC14 universal-failure test\n")
            f.write(f"Run on {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total runtime: {total_time:.0f}s\n\n")
            f.write(f"Pass: {n_pass}/{len(results)}\n")
            f.write(f"Fail: {n_fail}/{len(results)}\n\n")
            for r in results:
                f.write(f"{r['galaxy']:<10} T={r['T']} V={r['V_flat']:>4} "
                        f"chi^2_red={r['chi2_red']:.2f}  [{r['verdict']}]\n")
        print(f"Saved: dc14_universal_failures_summary.txt")


if __name__ == '__main__':
    main()
