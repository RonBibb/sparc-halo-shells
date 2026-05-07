"""
einasto_backbone.py — Einasto-backbone framework fit (v7.0)

This is the v7.0 revision — uses strict sigma/r_shell <= 0.4 constraint
enforcement matching run_canonical_fits.py, via the sigma_frac
reparameterization. The v6.5 version of this library used the looser
sigma <= rmax * 0.4 bound, which permitted shells whose width could exceed
the shell radius itself.

External API is unchanged:
  - fit_galaxy(name, rotmod_path, verbose=True) returns the same dict shape
  - 'fw_einasto_popt' contains physical sigma values (sigma_frac resolved)
  - All other keys identical to v6.5

Pipeline conventions (matching v6.2 / v6.5):
  - Fit V_obs(r) directly (not V_DM)
  - Model = sqrt(V2_DM_model + V2_baryonic_data)
  - Use only points with V2_DM > 0 (filter baryonic overshoot)
  - Multi-restart with curve_fit, BIC selection over N_shells in {0, 1, 2}

Profiles tested:
  Burkert (sanity reproduction of v6.2 numbers)
  Einasto-only (k=3)
  Einasto + 1 shell (k=6)
  Einasto + 2 shells (k=9)

For comparison, v6.2 NGC 5055:
  Burkert chi^2 = 190.58 (chi^2_red = 9.53)
  Framework (Burkert + 1 shell) chi^2 = 24.88 (chi^2_red = 1.46)
"""
import numpy as np
import pandas as pd
import os
from scipy.optimize import curve_fit
from scipy.special import gammainc, gamma, erf

G = 4.30091e-6  # (km/s)^2 kpc / Msun

# Strict constraint: sigma/r_shell <= SHELL_WIDTH_MAX_FRAC
SHELL_WIDTH_MAX_FRAC = 0.4

# ============================================================ profile primitives

def M_einasto(r, rho_s, r_s, alpha):
    M_inf = 2*np.pi * rho_s * r_s**3 * (alpha/2)**(3/alpha - 1) * np.exp(2/alpha) * gamma(3/alpha)
    s = (2.0/alpha) * (r/r_s)**alpha
    return M_inf * gammainc(3/alpha, s)

def V2_einasto(r, log_rho_s, log_r_s, alpha):
    return G * M_einasto(r, 10**log_rho_s, 10**log_r_s, alpha) / np.maximum(r, 1e-6)

def M_burkert(r, rho0, a):
    x = r/a
    return np.pi * rho0 * a**3 * (np.log(1+x**2) + 2*np.log(1+x) - 2*np.arctan(x))

def V2_burkert(r, log_rho0, log_a):
    return G * M_burkert(r, 10**log_rho0, 10**log_a) / np.maximum(r, 1e-6)

def V2_shell(r, M, rc, sigma):
    M_at_r = 0.5*M*(erf((r-rc)/(sigma*np.sqrt(2))) - erf((-rc)/(sigma*np.sqrt(2))))
    return G * M_at_r / np.maximum(r, 1e-6)

def V2_baryonic(V_gas, V_disk, V_bulge, Y_disk=0.5, Y_bulge=0.7):
    return V_gas*np.abs(V_gas) + Y_disk*V_disk*np.abs(V_disk) + Y_bulge*V_bulge*np.abs(V_bulge)

# ============================================================ closure factories
# v7.0: shell width parameter is sigma_frac (= sigma/rc), bounded [0.01, 0.4]

def make_models(backbone, V2_bar):
    """backbone in {'burkert', 'einasto'}. Returns (m0, m1, m2) total-V models.
    
    For models with shells, the shell width is parameterized as sigma_frac
    (= sigma/rc), with bounds [0.01, 0.4] hard-enforced via curve_fit's box
    bounds.
    """
    if backbone == 'einasto':
        def m0(r, log_rho_s, log_r_s, alpha):
            return np.sqrt(np.maximum(V2_einasto(r, log_rho_s, log_r_s, alpha) + V2_bar, 1e-6))
        def m1(r, log_rho_s, log_r_s, alpha, log_M, rc, sigma_frac):
            sigma = sigma_frac * rc
            v2 = V2_einasto(r, log_rho_s, log_r_s, alpha) + V2_bar
            v2 += V2_shell(r, 10**log_M, rc, sigma)
            return np.sqrt(np.maximum(v2, 1e-6))
        def m2(r, log_rho_s, log_r_s, alpha,
               log_M1, rc1, sf1, log_M2, rc2, sf2):
            s1 = sf1 * rc1
            s2 = sf2 * rc2
            v2 = V2_einasto(r, log_rho_s, log_r_s, alpha) + V2_bar
            v2 += V2_shell(r, 10**log_M1, rc1, s1)
            v2 += V2_shell(r, 10**log_M2, rc2, s2)
            return np.sqrt(np.maximum(v2, 1e-6))
        return m0, m1, m2
    elif backbone == 'burkert':
        def m0(r, log_rho0, log_a):
            return np.sqrt(np.maximum(V2_burkert(r, log_rho0, log_a) + V2_bar, 1e-6))
        def m1(r, log_rho0, log_a, log_M, rc, sigma_frac):
            sigma = sigma_frac * rc
            v2 = V2_burkert(r, log_rho0, log_a) + V2_bar
            v2 += V2_shell(r, 10**log_M, rc, sigma)
            return np.sqrt(np.maximum(v2, 1e-6))
        def m2(r, log_rho0, log_a, log_M1, rc1, sf1, log_M2, rc2, sf2):
            s1 = sf1 * rc1
            s2 = sf2 * rc2
            v2 = V2_burkert(r, log_rho0, log_a) + V2_bar
            v2 += V2_shell(r, 10**log_M1, rc1, s1)
            v2 += V2_shell(r, 10**log_M2, rc2, s2)
            return np.sqrt(np.maximum(v2, 1e-6))
        return m0, m1, m2

# ============================================================ multi-restart fit

def fit_with_restarts(model, r, V_obs, sigma_v, p0_list, bounds, maxfev=15000):
    best = None
    for p0 in p0_list:
        try:
            popt, _ = curve_fit(model, r, V_obs, p0=p0, sigma=sigma_v,
                                bounds=bounds, maxfev=maxfev, absolute_sigma=False)
            V_pred = model(r, *popt)
            chi2 = np.sum(((V_obs - V_pred)/sigma_v)**2)
            if best is None or chi2 < best['chi2']:
                best = {'popt': popt, 'chi2': chi2}
        except Exception:
            continue
    return best

def bic(chi2, n, k):
    return chi2 + k*np.log(n)

# ============================================================ helpers to convert popt back to physical

def _resolve_popt_e1(popt):
    """Convert (..., rc, sigma_frac) -> (..., rc, sigma)."""
    log_rho_s, log_r_s, alpha, log_M, rc, sf = popt
    return np.array([log_rho_s, log_r_s, alpha, log_M, rc, sf * rc])

def _resolve_popt_e2(popt):
    """Convert (..., rc1, sf1, ..., rc2, sf2) -> (..., rc1, s1, ..., rc2, s2)."""
    log_rho_s, log_r_s, alpha, log_M1, rc1, sf1, log_M2, rc2, sf2 = popt
    return np.array([log_rho_s, log_r_s, alpha,
                     log_M1, rc1, sf1 * rc1,
                     log_M2, rc2, sf2 * rc2])

# ============================================================ per-galaxy driver

def fit_galaxy(name, rotmod_path, verbose=True):
    data = np.loadtxt(rotmod_path)
    r = data[:,0]; V_obs = data[:,1]; sig_V = np.maximum(data[:,2], 1.0)
    V_gas = data[:,3]; V_disk = data[:,4]; V_bulge = data[:,5]
    V2_bar = V2_baryonic(V_gas, V_disk, V_bulge)

    # Filter overshoot
    keep = (V_obs**2 - V2_bar) > 0
    r = r[keep]; V_obs = V_obs[keep]; sig_V = sig_V[keep]; V2_bar = V2_bar[keep]
    n = len(r)
    if n < 5:
        return None
    rmax = r.max()

    # ============= Burkert sanity (k=2) =============
    mb = make_models('burkert', V2_bar)
    p0_burk = [[lr, la] for lr in [6.5, 7.5, 8.5] for la in [0.0, 0.7, 1.3, 2.0]]
    fit_burk = fit_with_restarts(mb[0], r, V_obs, sig_V, p0_burk, ([5,-1],[11,2.5]))
    bic_burk = bic(fit_burk['chi2'], n, 2)

    # ============= Einasto-only (k=3) ===============
    me = make_models('einasto', V2_bar)
    p0_ein = []
    for log_rho in [6.5, 7.5, 8.5, 9.5]:
        for log_rs in [0.0, 0.7, 1.3, 2.0]:
            for a in [0.12, 0.20, 0.35, 0.6]:
                p0_ein.append([log_rho, log_rs, a])
    bounds_ein = ([4, -1, 0.05], [11.5, 3, 1.2])
    fit_e0 = fit_with_restarts(me[0], r, V_obs, sig_V, p0_ein, bounds_ein)
    bic_e0 = bic(fit_e0['chi2'], n, 3)

    # ============= Einasto + 1 shell (k=6) ==========
    # v7.0: 6th param is sigma_frac (= sigma/rc), bounded [0.01, 0.4]
    p0_e1 = []
    for s in [p0_ein[i] for i in [21, 25, 41, 45, 53]]:  # mid-range Einasto seeds
        for rc_init in [rmax*0.2, rmax*0.45, rmax*0.7]:
            for logM in [9.0, 9.7, 10.3]:
                p0_e1.append(s + [logM, rc_init, 0.25])  # sigma_frac init = 0.25
    bounds_e1 = (list(bounds_ein[0]) + [6.0, 0.5, 0.01],
                 list(bounds_ein[1]) + [10.7, rmax, SHELL_WIDTH_MAX_FRAC])
    fit_e1 = fit_with_restarts(me[1], r, V_obs, sig_V, p0_e1, bounds_e1)
    bic_e1 = bic(fit_e1['chi2'], n, 6)

    # ============= Einasto + 2 shells (k=9) =========
    p0_e2 = []
    for s in [p0_ein[25], p0_ein[41]]:
        for rc1 in [rmax*0.2, rmax*0.4]:
            for rc2 in [rmax*0.55, rmax*0.8]:
                if rc2 > rc1*1.4:
                    p0_e2.append(s + [9.5, rc1, 0.25,
                                      9.5, rc2, 0.25])
    bounds_e2 = (list(bounds_ein[0]) + [6.0, 0.5, 0.01, 6.0, 0.5, 0.01],
                 list(bounds_ein[1]) + [10.7, rmax, SHELL_WIDTH_MAX_FRAC,
                                         10.7, rmax, SHELL_WIDTH_MAX_FRAC])
    fit_e2 = fit_with_restarts(me[2], r, V_obs, sig_V, p0_e2, bounds_e2) if p0_e2 else None
    bic_e2 = bic(fit_e2['chi2'], n, 9) if fit_e2 else None

    # BIC selection (over Einasto-family options only, matching v6.5 behavior)
    bics = [bic_e0, bic_e1] + ([bic_e2] if bic_e2 is not None else [])
    chi2s = [fit_e0['chi2'], fit_e1['chi2']] + ([fit_e2['chi2']] if fit_e2 else [])
    n_shells_best = int(np.argmin(bics))
    chi2_best = chi2s[n_shells_best]
    k_best = 3 + 3*n_shells_best
    chi2_red_best = chi2_best / max(n - k_best, 1)

    # Resolve sigma_frac back to physical sigma in popt for output
    if n_shells_best == 1:
        popt_physical = _resolve_popt_e1(fit_e1['popt']).tolist()
    elif n_shells_best == 2 and fit_e2:
        popt_physical = _resolve_popt_e2(fit_e2['popt']).tolist()
    else:
        popt_physical = fit_e0['popt'].tolist()

    if verbose:
        print(f"\n=== {name} (n={n} usable, rmax={rmax:.1f} kpc) ===")
        print(f"  Burkert (k=2):         chi^2 = {fit_burk['chi2']:7.2f}  chi^2_red = {fit_burk['chi2']/(n-2):6.2f}  BIC = {bic_burk:7.2f}")
        print(f"  Einasto only (k=3):    chi^2 = {fit_e0['chi2']:7.2f}  chi^2_red = {fit_e0['chi2']/(n-3):6.2f}  BIC = {bic_e0:7.2f}")
        print(f"  Einasto + 1sh (k=6):   chi^2 = {fit_e1['chi2']:7.2f}  chi^2_red = {fit_e1['chi2']/(n-6):6.2f}  BIC = {bic_e1:7.2f}")
        if fit_e2:
            print(f"  Einasto + 2sh (k=9):   chi^2 = {fit_e2['chi2']:7.2f}  chi^2_red = {fit_e2['chi2']/(n-9):6.2f}  BIC = {bic_e2:7.2f}")
        print(f"  -> BIC-selected: N_shells = {n_shells_best}, chi^2_red = {chi2_red_best:.2f}")
        if n_shells_best == 1:
            log_rho_s, log_r_s, alpha, log_M, rc, sigma = popt_physical
            print(f"     Shell: M = {10**log_M:.2e} Msun, rc = {rc:.2f} kpc, sigma = {sigma:.2f} kpc, sigma/rc = {sigma/rc:.2f}")
            print(f"     Einasto: rho_s = {10**log_rho_s:.2e}, r_s = {10**log_r_s:.2f} kpc, alpha = {alpha:.3f}")

    return {
        'galaxy': name, 'n': n, 'rmax': rmax,
        'burk_chi2': fit_burk['chi2'], 'burk_chi2_red': fit_burk['chi2']/(n-2), 'burk_bic': bic_burk,
        'ein0_chi2': fit_e0['chi2'], 'ein0_chi2_red': fit_e0['chi2']/(n-3), 'ein0_bic': bic_e0,
        'ein1_chi2': fit_e1['chi2'], 'ein1_chi2_red': fit_e1['chi2']/(n-6), 'ein1_bic': bic_e1,
        'ein2_chi2': fit_e2['chi2'] if fit_e2 else None,
        'ein2_chi2_red': (fit_e2['chi2']/max(n-9, 1)) if fit_e2 else None,
        'ein2_bic': bic_e2,
        'fw_einasto_n_shells': n_shells_best,
        'fw_einasto_chi2': chi2_best,
        'fw_einasto_chi2_red': chi2_red_best,
        'fw_einasto_popt': popt_physical,
    }

if __name__ == "__main__":
    import sys
    # Quick smoke test on NGC 5055; usage:
    #   python3 einasto_backbone.py [path/to/NGC5055_rotmod.dat]
    # Defaults to ./Rotmod_LTG/NGC5055_rotmod.dat
    rotmod_path = sys.argv[1] if len(sys.argv) > 1 else './Rotmod_LTG/NGC5055_rotmod.dat'
    res = fit_galaxy('NGC 5055', rotmod_path)
