"""
DC14 fit for NGC 5055 — single-galaxy showcase (§3.4).

This is a focused single-galaxy version of the universal-failure DC14 analysis,
provided as a separate script because §3.4's single-galaxy showcase is a
distinct scope from §3.6's universal-failure subset (Table 7). The numerical
results from this fit are reported in §3.4 (line 174) and Table 5.

Pipeline matches dc14_universal_failures.py exactly: same DC14 implementation,
same multi-restart optimization, same SPARC mass-to-light convention
(Upsilon_disk=0.5, Upsilon_bulge=0.7), same parameter bounds.

Usage:
  Place this script in a directory with a Rotmod_LTG/ subfolder containing
  NGC5055_rotmod.dat. Then:
    python3 dc14_ngc5055.py

  Output:
    - Console summary
    - dc14_ngc5055_showcase.csv (machine-readable, single-row table)

  Expected runtime: ~30-60 seconds.

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
GALAXY = 'NGC5055'
T_TYPE = 4
V_FLAT = 179.0
GALAXY_CLASS = 'showcase'   # tagged for clarity vs the universal-failure CSV

# Multi-restart grid (matches dc14_universal_failures.py)
RHO_S_STARTS = [1e6, 1e7, 5e7, 1e8, 1e9]
R_S_STARTS = [3.0, 8.0, 15.0, 30.0, 50.0]
X_STARTS = [-3.5, -3.0, -2.5, -2.0]
TIME_BUDGET = 120   # seconds; one galaxy gets a bit more budget than the bulk script

G = 4.302e-6   # kpc (km/s)^2 / M_sun

# ------------------------------------------------------------------
# Physics (identical to dc14_universal_failures.py)
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
    a = np.clip(a, 0.5, 5.0)
    b = np.clip(b, 2.5, 6.0)
    g = np.clip(g, 0.0, 1.5)
    return a, b, g


def dc14_v_circular(r_grid_for_integ, r_eval, rho_s, r_s, X):
    """DC14 circular velocity at r_eval (km/s)."""
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


def fit_dc14(rd, vd, ed, V2, time_budget=120):
    """Multi-restart DC14 fit. Returns (best_params, best_chi2, n_restarts_completed)."""
    sig = np.maximum(ed, 1.0)
    r_grid = np.geomspace(0.01, max(rd.max()*2, 200), 500)

    def model(r, rho_s, r_s, X):
        v_dm = dc14_v_circular(r_grid, r, rho_s, r_s, X)
        return np.sqrt(np.maximum(V2 + v_dm**2, 0))

    starts = [[rho_s, r_s, X] for rho_s in RHO_S_STARTS
                                for r_s in R_S_STARTS
                                for X in X_STARTS]

    best_chi2 = np.inf
    best_p = None
    t0 = time.time()
    n_completed = 0

    for p0 in starts:
        try:
            p, _ = curve_fit(model, rd, vd, p0=p0,
                             bounds=([1e3, 0.1, -4.0], [1e11, 200, -1.5]),
                             sigma=sig, maxfev=10000)
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
    print(f"DC14 fit for {GALAXY} (single-galaxy showcase, §3.4)")
    print(f"Data directory: {DATA_DIR}")
    print(f"Time budget: {TIME_BUDGET}s")
    print("="*72)

    path = os.path.join(DATA_DIR, f'{GALAXY}_rotmod.dat')
    if not os.path.exists(path):
        print(f"ERROR: Cannot find {path}")
        print("Adjust DATA_DIR at the top of this script.")
        return

    d = np.loadtxt(path, comments='#')
    r, vobs, evobs = d[:,0], d[:,1], d[:,2]
    vgas, vdisk, vbul = d[:,3], d[:,4], d[:,5]
    Vbar2 = baryon_squared(vgas, vdisk, vbul)
    mask = vobs**2 > Vbar2
    n_pts = mask.sum()
    rd = r[mask]; vd = vobs[mask]; ed = evobs[mask]; Vbd = Vbar2[mask]

    print(f"\nLoaded {GALAXY}: {len(r)} total points, {n_pts} usable (V_obs^2 > V_bar^2)")
    print(f"\nFitting DC14 with multi-restart...")

    t0 = time.time()
    best_p, best_chi2, n_completed = fit_dc14(rd, vd, ed, Vbd, TIME_BUDGET)
    elapsed = time.time() - t0

    if best_p is None:
        print(f"FAILED: no successful fit after {n_completed} restarts")
        return

    rho_s, r_s, X = best_p
    n_params = 3
    dof = n_pts - n_params
    chi2_red = best_chi2 / dof
    bic = best_chi2 + n_params * np.log(n_pts)

    a, b, g = dc14_shape_params(X)

    boundary = '*BOUNDARY*' if abs(X - (-1.5)) < 0.01 else 'interior'
    verdict = 'PASS' if chi2_red < 1.5 else 'FAIL'

    print(f"\nResults ({n_completed} restarts, {elapsed:.1f}s):")
    print(f"  rho_s        = {rho_s:.3e} M_sun/kpc^3")
    print(f"  r_s          = {r_s:.2f} kpc")
    print(f"  X            = {X:+.2f} ({boundary})")
    print(f"  alpha,beta,gamma = {a:.2f}, {b:.2f}, {g:.2f}")
    print(f"  chi^2        = {best_chi2:.2f}")
    print(f"  chi^2_red    = {chi2_red:.2f} ({verdict} adequacy <1.5)")
    print(f"  BIC          = {bic:.2f}")
    print(f"  n_pts        = {n_pts}")

    # Write CSV (single row, same column structure as dc14_universal_failures_results.csv)
    out_path = 'dc14_ngc5055_showcase.csv'
    with open(out_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['galaxy', 'T', 'V_flat', 'class', 'n_pts', 'n_restarts',
                         'rho_s', 'r_s', 'X', 'alpha', 'beta', 'gamma',
                         'chi2', 'chi2_red', 'bic', 'verdict', 'elapsed_s'])
        writer.writerow([GALAXY, T_TYPE, V_FLAT, GALAXY_CLASS, int(n_pts), int(n_completed),
                         f'{rho_s:.6e}', f'{r_s:.4f}', f'{X:.4f}',
                         f'{a:.4f}', f'{b:.4f}', f'{g:.4f}',
                         f'{best_chi2:.4f}', f'{chi2_red:.4f}', f'{bic:.4f}',
                         verdict, f'{elapsed:.2f}'])

    print(f"\nSaved: {out_path}")


if __name__ == '__main__':
    main()
