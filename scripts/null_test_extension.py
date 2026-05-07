"""
Extend the synthetic null test to T=8 and T=9.

Reads the canonical fits CSV (for the input fit parameters used to seed each
mock) and writes a 66-row T=8/9 extension CSV plus a combined T=2-9 CSV.

Run from the package root:
    cd halo_shells_v7.0/
    python3 scripts/null_test_extension.py

NOTE: in v7.0 this CSV is preserved unchanged from v6.5 because the v6.5
T=8/9 producer already enforced sigma/r <= 0.4 via clamping. Re-running this
script is supported but not required.
"""
import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.special import erf

# Package-relative paths
HERE = os.path.dirname(os.path.abspath(__file__))
PACKAGE_ROOT = os.path.dirname(HERE)
DATA_DIR = os.path.join(PACKAGE_ROOT, 'Rotmod_LTG')        # SPARC rotmod files
CANONICAL_CSV = os.path.join(PACKAGE_ROOT, 'data', 'sparc_T2-T9_canonical_fits.csv')
EXISTING_NULL_CSV = os.path.join(PACKAGE_ROOT, 'data', 'null_test_results.csv')
OUTPUT_DIR = os.path.join(PACKAGE_ROOT, 'data')
EXTENSION_OUT = os.path.join(OUTPUT_DIR, 'null_test_T8_T9_extension.csv')
COMBINED_OUT = os.path.join(OUTPUT_DIR, 'null_test_T2-T9_combined.csv')

G_KPC = 4.30091e-6
SIGMA_FLOOR = 1.0


def v2_burkert(r, rho0, a):
    r = np.maximum(r, 1e-6)
    a = max(a, 1e-6)
    x = r / a
    M_enc = np.pi * rho0 * a**3 * (
        np.log(1 + x**2) + 2*np.log(1 + x) - 2*np.arctan(x)
    )
    return G_KPC * M_enc / r


def v2_nfw(r, rho_s, r_s):
    r = np.maximum(r, 1e-6)
    r_s = max(r_s, 1e-6)
    x = r / r_s
    M_enc = 4*np.pi*rho_s*r_s**3 * (np.log(1+x) - x/(1+x))
    return G_KPC * M_enc / r


def v2_shell(r, M, r_c, sigma):
    r = np.maximum(r, 1e-6)
    sigma = max(sigma, 1e-6)
    M_enc = 0.5 * M * (
        erf((r - r_c)/(sigma*np.sqrt(2))) - erf(-r_c/(sigma*np.sqrt(2)))
    )
    return G_KPC * M_enc / r


def read_rotmod(galaxy):
    path = os.path.join(DATA_DIR, f"{galaxy}_rotmod.dat")
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            try:
                rad = float(parts[0])
                v_obs = float(parts[1])
                v_err = float(parts[2])
                v_gas = float(parts[3])
                v_disk = float(parts[4])
                v_bulge = float(parts[5])
            except ValueError:
                continue
            rows.append([rad, v_obs, v_err, v_gas, v_disk, v_bulge])
    arr = np.array(rows)
    if len(arr) == 0:
        return None
    return {
        'r': arr[:, 0],
        'V_obs': arr[:, 1],
        'sV': np.maximum(arr[:, 2], SIGMA_FLOOR),
        'V_gas': arr[:, 3],
        'V_disk': arr[:, 4],
        'V_bulge': arr[:, 5]
    }


def baryonic_v2(d):
    return (d['V_gas']*np.abs(d['V_gas'])
            + 0.5*d['V_disk']*np.abs(d['V_disk'])
            + 0.7*d['V_bulge']*np.abs(d['V_bulge']))


def fit_n0(r, V2_DM, sV2):
    seeds = [(1e7, 2.0), (5e7, 5.0), (2e7, 3.0), (1e8, 8.0), (3e7, 1.5)]
    best_chi2 = np.inf
    for s in seeds:
        try:
            popt, _ = curve_fit(
                lambda r, rho0, a: v2_burkert(r, abs(rho0), abs(a)),
                r, V2_DM, p0=s, sigma=sV2, maxfev=15000, absolute_sigma=False
            )
            chi2 = np.sum(((V2_DM - v2_burkert(r, abs(popt[0]), abs(popt[1])))/sV2)**2)
            if chi2 < best_chi2:
                best_chi2 = chi2
        except Exception:
            continue
    return best_chi2


def fit_n1(r, V2_DM, sV2):
    best_chi2 = np.inf
    for r_max in (3.0, 6.0, 12.0):
        for rho_seed in [1e7, 5e7]:
            for r1_init in [r_max*0.4, r_max*0.7]:
                try:
                    def model(r, rho0, a, M1, r1, s1):
                        s1c = min(abs(s1), 0.4*max(abs(r1), 1e-3))
                        return v2_burkert(r, abs(rho0), abs(a)) + v2_shell(r, abs(M1), abs(r1), s1c)
                    p0 = [rho_seed, 3.0, 1e9, r1_init, r1_init*0.3]
                    popt, _ = curve_fit(
                        model, r, V2_DM, p0=p0, sigma=sV2, maxfev=5000,
                        absolute_sigma=False,
                        bounds=([1e4, 0.1, 1e6, 0.1, 0.05],
                                [1e10, 50, 5e10, r_max, r_max*0.5])
                    )
                    chi2 = np.sum(((V2_DM - model(r, *popt))/sV2)**2)
                    if chi2 < best_chi2:
                        best_chi2 = chi2
                except Exception:
                    continue
    return best_chi2


def fit_n2(r, V2_DM, sV2):
    best_chi2 = np.inf
    for r_max in (6.0, 12.0):
        for rho_seed in [1e7, 5e7]:
            for r1_init, r2_init in [(r_max*0.25, r_max*0.6), (r_max*0.4, r_max*0.8)]:
                try:
                    def model(r, rho0, a, M1, r1, s1, M2, r2, s2):
                        s1c = min(abs(s1), 0.4*max(abs(r1), 1e-3))
                        s2c = min(abs(s2), 0.4*max(abs(r2), 1e-3))
                        return (v2_burkert(r, abs(rho0), abs(a))
                                + v2_shell(r, abs(M1), abs(r1), s1c)
                                + v2_shell(r, abs(M2), abs(r2), s2c))
                    p0 = [rho_seed, 3.0, 1e9, r1_init, r1_init*0.3,
                          1e9, r2_init, r2_init*0.3]
                    popt, _ = curve_fit(
                        model, r, V2_DM, p0=p0, sigma=sV2, maxfev=8000,
                        absolute_sigma=False,
                        bounds=([1e4, 0.1, 1e6, 0.1, 0.05, 1e6, 0.1, 0.05],
                                [1e10, 50, 5e10, r_max, r_max*0.5, 5e10, r_max, r_max*0.5])
                    )
                    chi2 = np.sum(((V2_DM - model(r, *popt))/sV2)**2)
                    if chi2 < best_chi2:
                        best_chi2 = chi2
                except Exception:
                    continue
    return best_chi2


def select_best_n(r, V2_DM, sV2, n_pts):
    chi2_0 = fit_n0(r, V2_DM, sV2)
    chi2_1 = fit_n1(r, V2_DM, sV2)
    chi2_2 = fit_n2(r, V2_DM, sV2)

    def bic(chi2, p):
        return chi2 + p*np.log(n_pts)

    bics = [bic(chi2_0, 2), bic(chi2_1, 5), bic(chi2_2, 8)]
    chi2s = [chi2_0, chi2_1, chi2_2]
    best_n = int(np.argmin(bics))
    return best_n, chi2s[best_n], bics[best_n]


def generate_mock(d, V2_bary, truth, params, rng):
    r = d['r']
    if truth == 'Burkert':
        rho0, a = params
        V2_DM_true = v2_burkert(r, rho0, a)
    else:
        rho_s, r_s = params
        V2_DM_true = v2_nfw(r, rho_s, r_s)
    V2_total = V2_bary + V2_DM_true
    V_total = np.sqrt(np.maximum(V2_total, 1e-3))
    V_obs_mock = V_total + rng.normal(0, d['sV'])
    return V_obs_mock


# Load canonical fits
CANON = pd.read_csv(CANONICAL_CSV).set_index('Galaxy')


def truth_params(galaxy, truth):
    row = CANON.loc[galaxy]
    if truth == 'Burkert':
        return (float(row['burk_rho0']), float(row['burk_a_kpc']))
    else:
        return (float(row['nfw_rho_s']), float(row['nfw_r_s_kpc']))


# T=8: all 6 galaxies (small bin); T=9: 5 stratified across V_flat
T8_GALAXIES = ['UGC01281', 'UGC04499', 'UGC07261', 'UGC02259', 'UGC07399', 'UGC00128']
T9_GALAXIES = ['UGCA442', 'UGC07125', 'UGC10310', 'NGC0055', 'UGC05986']


def main():
    results = []
    master_seed = 20260501

    for galaxy_list, T in [(T8_GALAXIES, 8), (T9_GALAXIES, 9)]:
        for galaxy in galaxy_list:
            d = read_rotmod(galaxy)
            if d is None:
                print(f"  SKIP: {galaxy} (no rotmod)")
                continue
            V2_bary = baryonic_v2(d)
            V2_obs = d['V_obs']**2
            good = V2_obs > V2_bary
            if good.sum() < 4:
                print(f"  SKIP: {galaxy} ({good.sum()} points)")
                continue

            r_use = d['r'][good]
            sV_use = d['sV'][good]
            V2_bary_use = V2_bary[good]
            V_flat = float(CANON.loc[galaxy, 'V_flat'])

            for truth in ['Burkert', 'NFW']:
                params = truth_params(galaxy, truth)
                for realization in range(3):
                    seed = (master_seed + hash((galaxy, truth, realization))) % (2**32)
                    rng = np.random.default_rng(seed)
                    V_obs_mock = generate_mock(d, V2_bary, truth, params, rng)
                    V2_obs_mock = V_obs_mock[good]**2
                    V2_DM_mock = V2_obs_mock - V2_bary_use
                    ok = V2_DM_mock > 0
                    if ok.sum() < 4:
                        continue
                    sV2_use = 2 * V_obs_mock[good][ok] * sV_use[ok]
                    sV2_use = np.maximum(sV2_use, 1.0)

                    best_n, best_chi2, best_bic = select_best_n(
                        r_use[ok], V2_DM_mock[ok], sV2_use, ok.sum()
                    )

                    results.append({
                        'Galaxy': galaxy,
                        'T': T,
                        'V_flat': V_flat,
                        'smooth_truth': truth,
                        'realization': realization,
                        'n_pts': int(ok.sum()),
                        'best_n_shells': int(best_n),
                        'best_chi2': float(best_chi2),
                        'best_bic': float(best_bic),
                    })
                    print(f"  {galaxy} T={T} truth={truth} real={realization}: "
                          f"best_n={best_n}, chi2={best_chi2:.1f}")

    new = pd.DataFrame(results)
    new.to_csv(EXTENSION_OUT, index=False)
    print(f"\nNew rows: {len(new)}  →  {EXTENSION_OUT}")

    existing = pd.read_csv(EXISTING_NULL_CSV)
    combined = pd.concat([existing, new], ignore_index=True)
    combined.to_csv(COMBINED_OUT, index=False)
    print(f"Combined rows: {len(combined)}  →  {COMBINED_OUT}")

    print("\n--- T=8/T=9 Summary ---")
    for T in [8, 9]:
        for truth in ['Burkert', 'NFW']:
            sub = new[(new['T'] == T) & (new['smooth_truth'] == truth)]
            if len(sub) == 0:
                continue
            n_shell = (sub['best_n_shells'] > 0).sum()
            rate = n_shell / len(sub)
            print(f"  T={T} {truth}-truth: {n_shell}/{len(sub)} = {rate*100:.1f}%")


if __name__ == "__main__":
    main()
