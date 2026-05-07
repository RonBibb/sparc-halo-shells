"""
figure_panels.py — Reproduce the rotation-curve figure panels (Figures 1, 2, 3,
                   and the supplementary T2-T9 grid) from the canonical fits CSV.

Modes:
  python3 figure_panels.py t-grid [T]    # T-type grid panel (T=2..9)
                                          # If T omitted, all 8 panels generated.
  python3 figure_panels.py ngc5055        # NGC 5055 showcase (Figure 2)
  python3 figure_panels.py universal      # Universal-failure grid (Figure 3)
  python3 figure_panels.py all            # All of the above

Inputs (read from PACKAGE_ROOT, with auto-fallback to ./Rotmod_LTG):
  data/sparc_T2-T9_canonical_fits.csv     # All fitted parameters per galaxy
  Rotmod_LTG/{Galaxy}_rotmod.dat          # SPARC observed rotation curves

Outputs (written to PACKAGE_ROOT/figures/):
  T{N}_burkert_vs_framework.png           # 8 T-grid panels
  figure1_T4_example_grid.png             # Same as T4 panel, copy for manuscript
  figure2_ngc5055_showcase.png            # NGC 5055 showcase + residuals
  figure3_universal_failures.png          # 9-galaxy universal-failure grid

Conventions match run_canonical_fits.py:
  Upsilon_disk = 0.5, Upsilon_bulge = 0.7
  V_bar^2 = Vgas|Vgas| + 0.5*Vdisk|Vdisk| + 0.7*Vbulge|Vbulge|
  Excluded points: where V_bar^2 >= V_obs^2

Run from package root:
  cd halo_shells_v7.0/
  python3 scripts/figure_panels.py all
"""

import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.special import erf

# ============================================================================
# Paths (package-relative with fallback)
# ============================================================================
HERE = os.path.dirname(os.path.abspath(__file__))
PACKAGE_ROOT = os.path.dirname(HERE)

_PRIMARY_DATA = os.path.join(PACKAGE_ROOT, 'Rotmod_LTG')
_LOCAL_DATA = os.path.join(os.getcwd(), 'Rotmod_LTG')
DATA_DIR = _PRIMARY_DATA if os.path.isdir(_PRIMARY_DATA) else _LOCAL_DATA

CANONICAL_CSV = os.path.join(PACKAGE_ROOT, 'data', 'sparc_T2-T9_canonical_fits.csv')
OUTPUT_DIR = os.path.join(PACKAGE_ROOT, 'figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)

G_KPC = 4.302e-6  # kpc (km/s)^2 / Msun  (matches run_canonical_fits.py)
SIGMA_FLOOR = 1.0  # km/s

# Morphological-type labels
T_LABEL = {2: 'Sab', 3: 'Sb', 4: 'Sb', 5: 'Sbc', 6: 'Sc', 7: 'Scd', 8: 'Sd', 9: 'Sdm'}


def refit_framework_backbone(galaxy_row, rotmod):
    """The canonical CSV stores per-galaxy shell parameters (M, r, sigma) but
    NOT the corresponding (rho_0, a) of the framework's Burkert backbone — the
    framework re-fits the backbone simultaneously with the shells, so its
    backbone differs from the stored 'burk_only' (rho_0, a). To accurately
    reproduce framework curves, we re-fit (rho_0, a) here with the stored
    shell parameters held fixed. This is a fast 2-parameter fit.

    Returns (rho_0, a) for the framework's Burkert backbone, or None if
    the BIC-selected configuration is zero shells (in which case the
    framework reduces exactly to burk-only and burk_rho0/burk_a_kpc apply)."""
    from scipy.optimize import curve_fit

    n_shells = int(galaxy_row['fw_best_n_shells'])
    if n_shells == 0:
        return None  # Use burk-only params

    r = rotmod['r']
    Vobs = rotmod['Vobs']
    errV = rotmod['errV']
    v2_bar = v2_baryonic(rotmod)
    mask = v2_bar < Vobs**2

    if n_shells == 1:
        M = galaxy_row['fw_n1_M_sh1']
        r_c = galaxy_row['fw_n1_r_sh1_kpc']
        sigma = galaxy_row['fw_n1_sigma_sh1_kpc']
        shells_v2 = v2_shell(r[mask], M, r_c, sigma)
    else:  # n_shells == 2
        M1 = galaxy_row['fw_n2_M_sh1']
        r1 = galaxy_row['fw_n2_r_sh1_kpc']
        s1 = galaxy_row['fw_n2_sigma_sh1_kpc']
        M2 = galaxy_row['fw_n2_M_sh2']
        r2 = galaxy_row['fw_n2_r_sh2_kpc']
        s2 = galaxy_row['fw_n2_sigma_sh2_kpc']
        shells_v2 = v2_shell(r[mask], M1, r1, s1) + v2_shell(r[mask], M2, r2, s2)

    v2_bar_m = v2_bar[mask]

    def model(r_arr, rho0, a):
        v2 = v2_bar_m + v2_burkert(r_arr, rho0, a) + shells_v2
        return np.sqrt(np.maximum(v2, 0))

    try:
        p, _ = curve_fit(model, r[mask], Vobs[mask],
                         p0=[float(galaxy_row['burk_rho0']),
                             float(galaxy_row['burk_a_kpc'])],
                         bounds=([1e3, 0.1], [1e10, 200]),
                         sigma=errV[mask], maxfev=20000)
        return float(p[0]), float(p[1])
    except Exception as e:
        # Fallback: use burk-only params
        return float(galaxy_row['burk_rho0']), float(galaxy_row['burk_a_kpc'])


# ============================================================================
# Halo profile evaluators (matching run_canonical_fits.py)
# ============================================================================

def v2_burkert(r, rho0, a):
    """Burkert profile V^2 contribution."""
    r = np.maximum(r, 1e-6)
    a = max(a, 1e-6)
    x = r / a
    M_enc = np.pi * rho0 * a**3 * (
        np.log(1 + x**2) + 2*np.log(1 + x) - 2*np.arctan(x)
    )
    return G_KPC * M_enc / r


def v2_nfw(r, rho_s, r_s):
    """NFW profile V^2 contribution."""
    r = np.maximum(r, 1e-6)
    r_s = max(r_s, 1e-6)
    x = r / r_s
    M_enc = 4 * np.pi * rho_s * r_s**3 * (np.log(1 + x) - x / (1 + x))
    return G_KPC * M_enc / r


def v2_shell(r, M, r_c, sigma):
    """Gaussian shell V^2 contribution. Matches run_canonical_fits.py:
    M_enc(r) = 0.5 * M * (1 + erf((r - r_c) / (sqrt(2) * sigma)))
    This integrates the Gaussian from -infinity to r."""
    r = np.maximum(r, 1e-6)
    sigma = max(sigma, 1e-6)
    M_enc = 0.5 * M * (1 + erf((r - r_c) / (np.sqrt(2) * sigma)))
    return G_KPC * M_enc / r


def v2_baryonic(rotmod):
    """V_bar^2 from rotmod data (Upsilon_disk=0.5, Upsilon_bulge=0.7)."""
    Vgas = rotmod['Vgas']
    Vdisk = rotmod['Vdisk']
    Vbulge = rotmod['Vbulge']
    return Vgas * np.abs(Vgas) + 0.5 * Vdisk * np.abs(Vdisk) + 0.7 * Vbulge * np.abs(Vbulge)


# ============================================================================
# Data loading
# ============================================================================

def read_rotmod(galaxy):
    """Read a SPARC rotmod file. Returns dict with arrays of r, Vobs, errV,
    Vgas, Vdisk, Vbulge — None if file missing."""
    # Try a few naming variations (SPARC files use both with and without dashes/underscores)
    candidates = [
        os.path.join(DATA_DIR, f"{galaxy}_rotmod.dat"),
        os.path.join(DATA_DIR, f"{galaxy.replace(' ', '')}_rotmod.dat"),
        os.path.join(DATA_DIR, f"{galaxy.replace(' ', '_')}_rotmod.dat"),
    ]
    for path in candidates:
        if os.path.exists(path):
            break
    else:
        return None

    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                rows.append([float(p) for p in parts[:6]])
            except ValueError:
                continue

    if not rows:
        return None

    arr = np.array(rows)
    return {
        'r': arr[:, 0],
        'Vobs': arr[:, 1],
        'errV': np.maximum(arr[:, 2], SIGMA_FLOOR),
        'Vgas': arr[:, 3],
        'Vdisk': arr[:, 4],
        'Vbulge': arr[:, 5],
    }


def evaluate_models(galaxy_row, r_eval, v2_bar_eval=None, fw_backbone=None):
    """Given a row from the canonical fits CSV and an array of evaluation radii,
    compute V_burkert(r), V_nfw(r), V_framework(r), and shell parameters list.
    If v2_bar_eval is provided, the returned velocities are the TOTAL observed-velocity
    predictions sqrt(V_bar^2 + V_DM^2). If None, returns just the V_DM contributions.

    fw_backbone: optional (rho_0, a) tuple for the framework's backbone (refit).
    If None, uses the burk-only stored (rho_0, a). For shell-bearing galaxies this
    matters because the framework re-fits the backbone alongside the shells.

    Returns dict with keys: v_burkert, v_nfw, v_framework, shells (list of (M,r_c,sigma))."""
    # Burkert (always present) — backbone for "burkert-only" comparison
    rho0 = galaxy_row['burk_rho0']
    a = galaxy_row['burk_a_kpc']
    v2_b = v2_burkert(r_eval, rho0, a)

    # NFW
    rho_s = galaxy_row['nfw_rho_s']
    r_s = galaxy_row['nfw_r_s_kpc']
    v2_n = v2_nfw(r_eval, rho_s, r_s)

    # Framework: Burkert backbone + BIC-selected shells
    n_shells = int(galaxy_row['fw_best_n_shells'])
    shells = []
    if n_shells == 1:
        shells.append((
            galaxy_row['fw_n1_M_sh1'],
            galaxy_row['fw_n1_r_sh1_kpc'],
            galaxy_row['fw_n1_sigma_sh1_kpc'],
        ))
    elif n_shells == 2:
        shells.append((
            galaxy_row['fw_n2_M_sh1'],
            galaxy_row['fw_n2_r_sh1_kpc'],
            galaxy_row['fw_n2_sigma_sh1_kpc'],
        ))
        shells.append((
            galaxy_row['fw_n2_M_sh2'],
            galaxy_row['fw_n2_r_sh2_kpc'],
            galaxy_row['fw_n2_sigma_sh2_kpc'],
        ))

    # Framework backbone — re-fit if provided, otherwise reuse burk-only
    if fw_backbone is not None and n_shells > 0:
        rho0_fw, a_fw = fw_backbone
        v2_fw_back = v2_burkert(r_eval, rho0_fw, a_fw)
    else:
        v2_fw_back = v2_b  # zero-shell case: framework == burk-only

    v2_fw = v2_fw_back.copy()
    for M, r_c, sigma in shells:
        v2_fw = v2_fw + v2_shell(r_eval, M, r_c, sigma)

    # If baryonic V^2 provided, return TOTAL V predictions (matches V_obs scale)
    if v2_bar_eval is not None:
        v2_b_total = v2_bar_eval + v2_b
        v2_n_total = v2_bar_eval + v2_n
        v2_fw_total = v2_bar_eval + v2_fw
    else:
        v2_b_total = v2_b
        v2_n_total = v2_n
        v2_fw_total = v2_fw

    return {
        'v_burkert': np.sqrt(np.maximum(v2_b_total, 0)),
        'v_nfw': np.sqrt(np.maximum(v2_n_total, 0)),
        'v_framework': np.sqrt(np.maximum(v2_fw_total, 0)),
        'shells': shells,
        'n_shells': n_shells,
    }


# ============================================================================
# Panel plotter — single galaxy, one subplot
# ============================================================================

def plot_galaxy_panel(ax, galaxy_row, rotmod, *, show_nfw=False, dense_only=False):
    """Plot a single galaxy's rotation curve and fits onto an Axes object."""
    galaxy = galaxy_row['Galaxy']
    r_obs = rotmod['r']
    Vobs = rotmod['Vobs']
    errV = rotmod['errV']

    # Compute V_baryonic
    v2_bar = v2_baryonic(rotmod)
    V_bar = np.sign(v2_bar) * np.sqrt(np.abs(v2_bar))

    # Identify excluded points (V_bar^2 >= V_obs^2)
    excluded_mask = v2_bar >= Vobs**2

    # Re-fit framework backbone given stored shell params (for shell-bearing galaxies)
    fw_backbone = refit_framework_backbone(galaxy_row, rotmod)

    # Dense radial grid for smooth model curves
    r_dense = np.linspace(r_obs.min(), r_obs.max(), 200)
    # Interpolate v2_bar onto dense grid for total-V model curves
    v2_bar_dense = np.interp(r_dense, r_obs, v2_bar)
    models = evaluate_models(galaxy_row, r_dense, v2_bar_eval=v2_bar_dense,
                             fw_backbone=fw_backbone)

    # --- Plot V_obs ---
    if not dense_only:
        included_mask = ~excluded_mask
        if included_mask.any():
            ax.errorbar(r_obs[included_mask], Vobs[included_mask],
                        yerr=errV[included_mask],
                        fmt='o', color='black', markersize=3.5,
                        elinewidth=0.8, capsize=2, zorder=10)
        if excluded_mask.any():
            ax.scatter(r_obs[excluded_mask], Vobs[excluded_mask],
                       marker='x', color='red', s=40, linewidth=1.5,
                       zorder=11, label='_nolegend_')

    # --- Plot V_baryonic ---
    ax.plot(r_dense, np.interp(r_dense, r_obs, V_bar), color='gray', lw=1.0, alpha=0.7)

    # --- Plot model curves ---
    ax.plot(r_dense, models['v_burkert'], '--', color='tab:orange', lw=1.4, alpha=0.85)
    if show_nfw:
        ax.plot(r_dense, models['v_nfw'], '-.', color='tab:green', lw=1.4, alpha=0.85)
    ax.plot(r_dense, models['v_framework'], '-', color='tab:blue', lw=1.8)

    # --- Shade shell regions ---
    for M, r_c, sigma in models['shells']:
        ax.axvspan(r_c - sigma, r_c + sigma, alpha=0.15, color='steelblue', zorder=1)

    # --- Title bar ---
    n_exc = int(excluded_mask.sum())
    Vflat = galaxy_row['V_flat']
    burk_chi2r = galaxy_row['burk_chi2_red']
    fw_chi2r = galaxy_row['fw_best_chi2_red']
    n_shells = models['n_shells']

    is_shell_bearing = n_shells >= 1
    title_color = '#d4e7d0' if is_shell_bearing else '#f5f5f5'

    if show_nfw:
        nfw_chi2r = galaxy_row['nfw_chi2_red']
        title = (f"{galaxy}  T={int(galaxy_row['T'])}  V={int(Vflat)}  "
                 f"({'CLEAN' if n_exc == 0 else 'dirty'}, n_excl={n_exc})\n"
                 f"Burk χ²r={burk_chi2r:.1f}  |  "
                 f"NFW χ²r={nfw_chi2r:.1f}  |  "
                 f"FW(n={n_shells}) χ²r={fw_chi2r:.1f}")
    else:
        title = (f"{galaxy}  V={int(Vflat)}  excl={n_exc}\n"
                 f"Burk χ²r={burk_chi2r:.1f}  |  "
                 f"FW(n={n_shells}) χ²r={fw_chi2r:.1f}")

    ax.set_title(title, fontsize=8.5, pad=3,
                 bbox=dict(facecolor=title_color, edgecolor='gray',
                           linewidth=0.5, pad=2))
    ax.set_xlabel('r (kpc)', fontsize=8)
    ax.set_ylabel('V (km/s)', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.2)


# ============================================================================
# T-type grid generator
# ============================================================================

def make_t_grid(T, df, output_path=None):
    """Make a multi-panel grid for all galaxies in a single T-bin."""
    sub = df[df['T'] == T].sort_values('Galaxy').reset_index(drop=True)
    n = len(sub)
    if n == 0:
        print(f"  No galaxies for T={T}")
        return

    # Layout: 3 columns; rows = ceil((n+1)/3) so there's room for legend cell
    n_cols = 3
    n_rows = math.ceil((n + 1) / n_cols)

    fig = plt.figure(figsize=(n_cols * 5.0, n_rows * 3.0))
    gs = GridSpec(n_rows, n_cols, figure=fig, hspace=0.55, wspace=0.30)

    plotted = 0
    for i, row in sub.iterrows():
        rotmod = read_rotmod(row['Galaxy'])
        if rotmod is None:
            print(f"  Warning: rotmod file missing for {row['Galaxy']}")
            continue
        ax = fig.add_subplot(gs[i // n_cols, i % n_cols])
        plot_galaxy_panel(ax, row, rotmod, show_nfw=False)
        plotted += 1

    # Legend cell (bottom-right corner)
    legend_ax = fig.add_subplot(gs[n_rows - 1, n_cols - 1])
    legend_ax.axis('off')

    # Synthetic legend handles
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    legend_handles = [
        Line2D([0], [0], marker='o', color='black', linestyle='',
               markersize=5, label='V_obs'),
        Line2D([0], [0], marker='x', color='red', linestyle='',
               markersize=8, label='V_DM clamped (excl.)'),
        Line2D([0], [0], color='gray', lw=1.0, label='V_baryonic'),
        Line2D([0], [0], color='tab:orange', lw=1.5, ls='--', label='Burkert-only'),
        Line2D([0], [0], color='tab:blue', lw=1.8, label='Framework (Burk + N shells)'),
        Patch(facecolor='steelblue', alpha=0.15, label='shell ±σ'),
    ]
    legend_ax.legend(handles=legend_handles, loc='upper center', fontsize=9,
                     framealpha=0.95)

    # Summary stats box
    n_shell_bearing = (sub['fw_best_n_shells'] >= 1).sum()
    sum_chi2_burk = sub['burk_chi2'].sum()
    sum_chi2_fw = sub['fw_best_chi2'].sum()
    reduction = 100 * (1 - sum_chi2_fw / sum_chi2_burk) if sum_chi2_burk > 0 else 0
    summary = (f"Reduction: {reduction:.1f}%\n"
               f"Shell-bearing: {n_shell_bearing}/{n}")
    legend_ax.text(0.5, 0.20, summary, transform=legend_ax.transAxes,
                   ha='center', va='top', fontsize=10,
                   bbox=dict(facecolor='#fff8e1', edgecolor='gray',
                             linewidth=0.5, pad=4))

    fig.suptitle(f"T = {T} galaxies — Burkert vs Framework (Burkert + BIC shells)",
                 fontsize=13, fontweight='bold', y=0.995)

    if output_path is None:
        output_path = os.path.join(OUTPUT_DIR, f"T{T}_burkert_vs_framework.png")
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}  ({plotted}/{n} galaxies plotted)")


# ============================================================================
# NGC 5055 showcase (Figure 2)
# ============================================================================

def make_ngc5055_showcase(df):
    """Two-panel: main rotation curve + residuals. NGC 5055 only."""
    sub = df[df['Galaxy'] == 'NGC5055']
    if len(sub) == 0:
        print("  NGC5055 not found in canonical CSV")
        return
    row = sub.iloc[0]

    rotmod = read_rotmod('NGC5055')
    if rotmod is None:
        print("  NGC5055 rotmod file missing")
        return

    r_obs = rotmod['r']
    Vobs = rotmod['Vobs']
    errV = rotmod['errV']
    v2_bar = v2_baryonic(rotmod)
    V_bar = np.sign(v2_bar) * np.sqrt(np.abs(v2_bar))
    excluded_mask = v2_bar >= Vobs**2

    # Re-fit framework backbone given stored shell params
    fw_backbone = refit_framework_backbone(row, rotmod)

    # Models at observed radii (for residuals) and dense grid (for curves)
    r_dense = np.linspace(r_obs.min(), r_obs.max(), 400)
    v2_bar_dense = np.interp(r_dense, r_obs, v2_bar)
    models_dense = evaluate_models(row, r_dense, v2_bar_eval=v2_bar_dense,
                                   fw_backbone=fw_backbone)
    models_obs = evaluate_models(row, r_obs, v2_bar_eval=v2_bar,
                                 fw_backbone=fw_backbone)

    fig = plt.figure(figsize=(13, 7.5))
    gs = GridSpec(2, 1, figure=fig, height_ratios=[3, 1], hspace=0.05)

    # --- Top panel: curves ---
    ax1 = fig.add_subplot(gs[0])

    included = ~excluded_mask
    ax1.errorbar(r_obs[included], Vobs[included], yerr=errV[included],
                 fmt='o', color='black', markersize=4, elinewidth=0.9,
                 capsize=2, label='V_obs', zorder=10)
    if excluded_mask.any():
        ax1.scatter(r_obs[excluded_mask], Vobs[excluded_mask],
                    marker='x', color='red', s=60, linewidth=1.8,
                    label=f'V_DM clamped to 0 ({excluded_mask.sum()} pts)', zorder=11)

    ax1.plot(r_dense, np.interp(r_dense, r_obs, V_bar), color='gray',
             lw=1.0, alpha=0.8, label='V_baryonic')

    burk_chi2r = row['burk_chi2_red']
    nfw_chi2r = row['nfw_chi2_red']
    fw_chi2r = row['fw_best_chi2_red']
    n_shells = int(row['fw_best_n_shells'])

    ax1.plot(r_dense, models_dense['v_burkert'], '--', color='tab:orange',
             lw=2, label=f"Burkert family best:  χ²_red = {burk_chi2r:.1f}")
    ax1.plot(r_dense, models_dense['v_nfw'], '-.', color='tab:green',
             lw=2, label=f"NFW family best:  χ²_red = {nfw_chi2r:.1f}")
    ax1.plot(r_dense, models_dense['v_framework'], '-', color='tab:blue',
             lw=2.5, label=f"Burkert + {n_shells} shell:  χ²_red = {fw_chi2r:.1f}")

    # Shade shells
    for M, r_c, sigma in models_dense['shells']:
        ax1.axvspan(r_c - sigma, r_c + sigma, alpha=0.18, color='steelblue', zorder=1)

    Vflat = row['V_flat']
    ax1.set_title(f"NGC 5055 — single-galaxy showcase\n"
                  f"V_flat = {Vflat:.0f} km/s, T = {int(row['T'])}",
                  fontsize=12, fontweight='bold')
    ax1.set_ylabel('V (km/s)', fontsize=11)
    ax1.tick_params(labelbottom=False)
    ax1.grid(alpha=0.3)
    ax1.legend(loc='lower right', fontsize=9, framealpha=0.95)

    # --- Bottom panel: residuals ---
    ax2 = fig.add_subplot(gs[1], sharex=ax1)

    # Residuals at observed radii (only for included points)
    resid_burk = Vobs[included] - models_obs['v_burkert'][included]
    resid_nfw = Vobs[included] - models_obs['v_nfw'][included]
    resid_fw = Vobs[included] - models_obs['v_framework'][included]

    ax2.errorbar(r_obs[included] - 0.1, resid_burk, yerr=errV[included],
                 fmt='s', color='tab:orange', markersize=4, elinewidth=0.7,
                 capsize=1.5, label='Burkert resid', alpha=0.8)
    ax2.errorbar(r_obs[included], resid_nfw, yerr=errV[included],
                 fmt='^', color='tab:green', markersize=4, elinewidth=0.7,
                 capsize=1.5, label='NFW resid', alpha=0.8)
    ax2.errorbar(r_obs[included] + 0.1, resid_fw, yerr=errV[included],
                 fmt='o', color='tab:blue', markersize=4, elinewidth=0.7,
                 capsize=1.5, label='Framework resid', alpha=0.85)

    ax2.axhline(0, color='black', lw=0.5)
    ax2.set_xlabel('r (kpc)', fontsize=11)
    ax2.set_ylabel('V_obs - V_model\n(km/s)', fontsize=10)
    ax2.legend(loc='lower right', fontsize=8, ncol=3, framealpha=0.95)
    ax2.grid(alpha=0.3)

    output_path = os.path.join(OUTPUT_DIR, 'figure2_ngc5055_showcase.png')
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


# ============================================================================
# Universal-failure grid (Figure 3)
# ============================================================================

# The 9 universal-failure galaxies in the T=2-7 range, sorted by descending
# framework chi2_red (worst first). Names match the canonical CSV.
UNIVERSAL_FAILURES = [
    'NGC6674', 'UGC06787', 'UGC09133',
    'NGC5907', 'NGC5033', 'NGC2841',
    'UGC02953', 'NGC5985', 'NGC0289',
]


def make_universal_failures(df):
    """3x3 grid of universal-failure galaxies, sorted by worst framework chi2_red."""
    # Filter and order
    found = df[df['Galaxy'].isin(UNIVERSAL_FAILURES)].copy()
    if len(found) == 0:
        print("  No universal-failure galaxies found in CSV")
        return

    # Sort by framework chi2_red descending (worst first)
    found = found.sort_values('fw_best_chi2_red', ascending=False).reset_index(drop=True)

    n = len(found)
    n_cols = 3
    n_rows = math.ceil(n / n_cols)

    fig = plt.figure(figsize=(n_cols * 5.0, n_rows * 4.0))
    gs = GridSpec(n_rows, n_cols, figure=fig, hspace=0.55, wspace=0.30)

    for i, row in found.iterrows():
        rotmod = read_rotmod(row['Galaxy'])
        if rotmod is None:
            print(f"  Warning: rotmod missing for {row['Galaxy']}")
            continue
        ax = fig.add_subplot(gs[i // n_cols, i % n_cols])
        plot_galaxy_panel(ax, row, rotmod, show_nfw=True)

    # Legend at bottom
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [
        Line2D([0], [0], marker='o', color='black', linestyle='',
               markersize=5, label='V_obs'),
        Line2D([0], [0], marker='x', color='red', linestyle='',
               markersize=8, label='V_DM clamped (excl.)'),
        Line2D([0], [0], color='gray', lw=1.0, label='V_baryonic'),
        Line2D([0], [0], color='tab:orange', lw=1.5, ls='--', label='Burkert-only'),
        Line2D([0], [0], color='tab:green', lw=1.5, ls='-.', label='NFW-only'),
        Line2D([0], [0], color='tab:blue', lw=1.8, label='Framework (Burk + N shells)'),
        Patch(facecolor='steelblue', alpha=0.15, label='shell ±σ region'),
    ]
    fig.legend(handles=handles, loc='lower center', ncol=4,
               fontsize=10, framealpha=0.95, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle("Universal-failure galaxies: where all three models give χ²_red ≥ 1.5\n"
                 "Sorted worst→best by framework χ²_red",
                 fontsize=12, fontweight='bold', y=0.995)

    output_path = os.path.join(OUTPUT_DIR, 'figure3_universal_failures.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}  ({n} galaxies)")


# ============================================================================
# Main entry point
# ============================================================================

def main():
    if not os.path.exists(CANONICAL_CSV):
        sys.exit(f"ERROR: canonical fits CSV not found at {CANONICAL_CSV}")
    if not os.path.isdir(DATA_DIR):
        sys.exit(f"ERROR: SPARC rotmod directory not found at {DATA_DIR}\n"
                 "  Set DATA_DIR or place rotmod files in PACKAGE_ROOT/Rotmod_LTG/")

    df = pd.read_csv(CANONICAL_CSV)
    print(f"Loaded {len(df)} galaxies from {CANONICAL_CSV}")
    print(f"Reading rotmod files from {DATA_DIR}")

    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'

    if mode == 't-grid':
        if len(sys.argv) > 2:
            T = int(sys.argv[2])
            print(f"\nGenerating T={T} grid panel...")
            make_t_grid(T, df)
        else:
            print("\nGenerating all 8 T-grid panels (T=2..9)...")
            for T in range(2, 10):
                make_t_grid(T, df)
            # Copy T=4 as figure1
            t4_path = os.path.join(OUTPUT_DIR, 'T4_burkert_vs_framework.png')
            fig1_path = os.path.join(OUTPUT_DIR, 'figure1_T4_example_grid.png')
            if os.path.exists(t4_path):
                import shutil
                shutil.copy(t4_path, fig1_path)
                print(f"  Copied: {fig1_path} (= T4 panel, used as Figure 1 in manuscript)")

    elif mode == 'ngc5055':
        print("\nGenerating NGC 5055 showcase (Figure 2)...")
        make_ngc5055_showcase(df)

    elif mode == 'universal':
        print("\nGenerating universal-failure grid (Figure 3)...")
        make_universal_failures(df)

    elif mode == 'all':
        print("\n=== All T-grid panels ===")
        for T in range(2, 10):
            make_t_grid(T, df)
        # Copy T=4 as figure1
        t4_path = os.path.join(OUTPUT_DIR, 'T4_burkert_vs_framework.png')
        fig1_path = os.path.join(OUTPUT_DIR, 'figure1_T4_example_grid.png')
        if os.path.exists(t4_path):
            import shutil
            shutil.copy(t4_path, fig1_path)
            print(f"  Copied: {fig1_path}")

        print("\n=== NGC 5055 showcase ===")
        make_ngc5055_showcase(df)

        print("\n=== Universal-failure grid ===")
        make_universal_failures(df)

    else:
        sys.exit(f"Unknown mode: {mode}\n"
                 "Usage: figure_panels.py [t-grid [T] | ngc5055 | universal | all]")

    print("\nDone.")


if __name__ == '__main__':
    main()
