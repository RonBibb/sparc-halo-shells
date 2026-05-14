#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ngc5055_marginalized_upsilon.py

Empirical test of whether the BIC shell preference at NGC 5055 (the §3.4
single-galaxy showcase in sparc_shells.pdf, v7.1.0) survives joint
marginalization of the stellar mass-to-light ratios Υ_disk and Υ_bulge
under Gaussian log-priors at the Li et al. 2020 (ApJS 247, 31) convention:

    log10(Υ_disk)  ~ N(log10(0.5), 0.1)   # canonical Υ_disk  = 0.5, 1σ = 0.1 dex
    log10(Υ_bulge) ~ N(log10(0.7), 0.1)   # canonical Υ_bulge = 0.7, 1σ = 0.1 dex

This script directly addresses the §4.4 disclosure gap flagged in the PASP
resubmission review: the canonical T-dependent steelman pushes Υ_disk upward
at early types (T=2 anchored at 0.65), which makes shells more necessary, not
less. It does not adversarially stress-test the showcase galaxy in the
direction that would dissolve its shell preference (Υ_disk pulled down toward
the late-type-blue end). A statistically literate referee will identify this
asymmetry on a first read.

Three possible outcomes:

    (1) Shell preference holds with marginalized Υ
        → §3.4 strengthens; §4.4 disclosure shortens to a single sentence
        → Update §3.4 last paragraph with parenthetical confirming the
          marginalized result and ΔBIC vs Burkert-only-marg.
    (2) Shell preference reverses
        → Replace NGC 5055 as the showcase (NGC 2841 is the obvious
          candidate from the universal-failure set), or reframe §3.4
          around the four-smooth-profile failure pattern alone.
    (3) Marginal result (|ΔBIC| < 6 between marginalized 0-shell and 1-shell)
        → Report honestly; reposition §3.4 as "shells are competitive with
          smooth profiles under marginalized Υ" rather than "decisively win."

REFERENCE NUMBERS (canonical Table 5, v7.1.0):

    Model                    χ²       χ²_red   BIC      ΔBIC vs Framework
    Burkert-only            190.6     9.53    196.8    +156.5
    NFW (free c)            604.1    30.21    610.3    +570.0
    DC14 (3 free)           171.1     9.01    180.4    +140.1
    Framework (1 shell)      24.9     1.46     40.3    ---  (best)
    Shell parameters: M = 3.78e10 M_sun, r = 12.0 kpc, σ = 2.91 kpc

For NGC 5055 specifically, V_bulge = 0 throughout the SPARC rotmod file
(bulgeless in the SPARC photometric decomposition), so Υ_bulge has zero
leverage on V_bar² and will return to its prior mean. The script keeps
Υ_bulge in the parameter list so the same code applies cleanly to galaxies
with non-zero bulge contributions (e.g., NGC 2841 if a replacement showcase
is needed).

OUTPUT
    data/ngc5055_marginalized_upsilon_results.csv

USAGE
    python scripts/ngc5055_marginalized_upsilon.py
    python scripts/ngc5055_marginalized_upsilon.py --data-dir ./Rotmod_LTG
    python scripts/ngc5055_marginalized_upsilon.py --n-restarts-per-rmax 16
    python scripts/ngc5055_marginalized_upsilon.py --galaxy NGC2841

Author: Ronald Bibb
Created: 2026-05-13 for PASP resubmission review
"""
from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
from scipy import optimize
from scipy.special import erf


# =============================================================================
# CONSTANTS AND DEFAULTS
# =============================================================================

# Gravitational constant in units of (km/s)² · kpc / M_sun
G = 4.30091e-6

# Default galaxy and reference values (NGC 5055, from §3.4 / Table 5)
DEFAULT_GALAXY = "NGC5055"
DEFAULT_DATA_DIR = "./Rotmod_LTG"

# Canonical fixed-Υ BIC values (Table 5, v7.1.0) for NGC 5055.
# These are used only for printing the comparison report; the script does not
# refit the canonical case.
CANONICAL_REFERENCE = {
    "NGC5055": {
        "n_data": 22,
        "burkert_only_chi2": 190.6,
        "burkert_only_chi2_red": 9.53,
        "burkert_only_bic": 196.8,
        "framework_1shell_chi2": 24.9,
        "framework_1shell_chi2_red": 1.46,
        "framework_1shell_bic": 40.3,
        "framework_shell_M": 3.78e10,
        "framework_shell_r": 12.0,
        "framework_shell_sigma": 2.91,
    },
}

# Υ priors (Li 2020 convention; widths in dex)
PRIOR_UPS_DISK_MEAN = 0.5
PRIOR_UPS_DISK_SIGMA_DEX = 0.1
PRIOR_UPS_BULGE_MEAN = 0.7
PRIOR_UPS_BULGE_SIGMA_DEX = 0.1

# Shell parameter bounds (match canonical pipeline §2.2)
SHELL_M_MIN = 1.0e6        # M_sun
SHELL_M_MAX = 5.0e10       # M_sun
SHELL_F_MIN = 0.01         # σ/r minimum
SHELL_F_MAX = 0.4          # σ/r maximum (strict architectural constraint)
SHELL_R_MIN = 0.5          # kpc
SHELL_R_MAX = 50.0         # kpc

# Burkert backbone bounds
RHO0_MIN = 1.0e3           # M_sun / kpc³
RHO0_MAX = 1.0e10          # M_sun / kpc³
A_MIN = 0.1                # kpc
A_MAX = 50.0               # kpc

# Υ search box (well outside the prior range to let the prior do the work)
UPS_DISK_BOX_MIN = 0.1
UPS_DISK_BOX_MAX = 2.0
UPS_BULGE_BOX_MIN = 0.1
UPS_BULGE_BOX_MAX = 2.0

# Multi-restart configuration. Extended rmax options give r_init coverage
# spanning small to large radii, so restarts cover both narrow inner shells
# and broad outer shells. NGC 5055's canonical shell sits at r = 12 kpc,
# which is well-covered by the r_max = 15 and r_max = 25 options.
DEFAULT_N_RESTARTS_PER_RMAX = 8       # → 32 total per fit across 4 r_max values
RMAX_OPTIONS = (3.0, 8.0, 15.0, 25.0) # kpc
MAXFEV = 20000
FTOL = 1e-10
GTOL = 1e-8

# Reproducibility
DEFAULT_RANDOM_SEED = 20260513

# σ_V floor (per §2.3)
SIGMA_V_FLOOR_KMS = 1.0

# Adequacy threshold
CHI2_RED_ADEQUATE = 1.5

# Kass & Raftery 1995 ΔBIC interpretation thresholds
DBIC_THRESHOLDS = {
    "very_strong": 10,
    "strong": 6,
    "positive": 2,
}


# =============================================================================
# DATA LOADING
# =============================================================================

def load_rotmod(filepath: str) -> dict:
    """Load a SPARC rotmod file. Returns dict of arrays.

    Columns (per SPARC convention, Lelli et al. 2016):
        Rad(kpc) Vobs(km/s) errV(km/s) Vgas Vdisk Vbul SBdisk SBbul

    The disk and bulge velocity columns are computed at Υ_disk = Υ_bulge = 1
    (i.e., they represent V² contributions per unit Υ in the V_bar² sum of
    equation 1 in §2.1).
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(
            f"Cannot find rotmod file: {filepath}\n"
            f"  Set --data-dir to the directory containing *_rotmod.dat files."
        )

    # Read the distance line from the comment header
    distance_mpc = None
    with open(filepath, "r") as f:
        for line in f:
            if line.lower().startswith("# distance"):
                # Format: "# Distance = X.X Mpc"
                parts = line.split("=")
                if len(parts) >= 2:
                    distance_mpc = float(parts[1].strip().split()[0])
                break

    data = np.loadtxt(filepath, comments="#")
    if data.ndim == 1:
        data = data[np.newaxis, :]

    sigma_V = np.maximum(data[:, 2], SIGMA_V_FLOOR_KMS)

    return {
        "r_kpc": data[:, 0],
        "V_obs": data[:, 1],
        "sigma_V": sigma_V,
        "V_gas": data[:, 3],
        "V_disk": data[:, 4],
        "V_bulge": data[:, 5],
        "distance_mpc": distance_mpc,
    }


# =============================================================================
# MODEL COMPONENTS
# =============================================================================

def v_baryon_squared(Y_disk: float, Y_bulge: float,
                     V_gas: np.ndarray, V_disk: np.ndarray, V_bulge: np.ndarray
                     ) -> np.ndarray:
    """Baryonic V² with free Υ (Equation 1 of §2.1).

    The sign-preserving formulation V_i |V_i| handles cases where the SPARC
    catalog reports negative rotation contributions (inward pressure-gradient
    regions in some galaxies).
    """
    return (V_gas * np.abs(V_gas)
            + Y_disk * V_disk * np.abs(V_disk)
            + Y_bulge * V_bulge * np.abs(V_bulge))


def v_burkert_squared(r_kpc: np.ndarray, rho_0: float, a_kpc: float) -> np.ndarray:
    """Burkert profile V²(r) via enclosed mass (Equation 4 of §2.2).

    ρ(r) = ρ_0 a³ / [(r + a)(r² + a²)]
    M(r) = π ρ_0 a³ [2 ln(1 + r/a) + ln(1 + (r/a)²) - 2 arctan(r/a)]
    V² = G M / r
    """
    x = r_kpc / a_kpc
    M_enc = np.pi * rho_0 * a_kpc**3 * (
        2.0 * np.log(1.0 + x)
        + np.log(1.0 + x**2)
        - 2.0 * np.arctan(x)
    )
    return G * M_enc / r_kpc


def v_shell_squared(r_kpc: np.ndarray, M_shell: float, r_shell: float,
                    sigma_shell: float) -> np.ndarray:
    """Gaussian shell V²(r) via enclosed mass (Equations 5-6 of §2.2)."""
    inv_sigma_sqrt2 = 1.0 / (sigma_shell * np.sqrt(2.0))
    M_enc = 0.5 * M_shell * (
        erf((r_kpc - r_shell) * inv_sigma_sqrt2)
        - erf(-r_shell * inv_sigma_sqrt2)
    )
    return G * M_enc / r_kpc


# =============================================================================
# FIT MODEL
# =============================================================================

# Parameter layout:
#   params[0] = Y_disk
#   params[1] = Y_bulge
#   params[2] = rho_0      (M_sun / kpc³)
#   params[3] = a_kpc      (kpc)
#   For each shell i (0-indexed) in n_shells:
#     params[4 + 3i + 0] = M_shell_i        (M_sun)
#     params[4 + 3i + 1] = r_shell_i        (kpc)
#     params[4 + 3i + 2] = f_shell_i = σ_i / r_i   (dimensionless, in [0.01, 0.4])
#
# Reparameterizing σ = f * r enforces the σ/r ≤ 0.4 constraint as a box bound
# on f, mirroring the canonical pipeline (§2.2).


def n_params(n_shells: int, fix_upsilon: bool = False) -> int:
    # Marginalized: Υ_disk, Υ_bulge, ρ_0, a, plus 3 per shell
    # Fix-Υ: ρ_0, a, plus 3 per shell (Υ values are constants, not free params)
    n_upsilon = 0 if fix_upsilon else 2
    return n_upsilon + 2 + 3 * n_shells


def unpack_params(params: np.ndarray, n_shells: int) -> tuple:
    """Return (Y_disk, Y_bulge, rho_0, a, list_of_(M_i, r_i, sigma_i))."""
    Y_disk, Y_bulge, rho_0, a_kpc = params[:4]
    shells = []
    for i in range(n_shells):
        offset = 4 + 3 * i
        M_i = params[offset]
        r_i = params[offset + 1]
        f_i = params[offset + 2]
        sigma_i = f_i * r_i
        shells.append((M_i, r_i, sigma_i))
    return Y_disk, Y_bulge, rho_0, a_kpc, shells


def v_model_total(params: np.ndarray, n_shells: int, data: dict) -> np.ndarray:
    """V_obs model = sqrt(V_bar² + V_DM²)."""
    Y_disk, Y_bulge, rho_0, a_kpc, shells = unpack_params(params, n_shells)
    V_bar_sq = v_baryon_squared(Y_disk, Y_bulge,
                                 data["V_gas"], data["V_disk"], data["V_bulge"])
    V_dm_sq = v_burkert_squared(data["r_kpc"], rho_0, a_kpc)
    for M_i, r_i, sigma_i in shells:
        V_dm_sq += v_shell_squared(data["r_kpc"], M_i, r_i, sigma_i)
    V_total_sq = V_bar_sq + V_dm_sq
    # Clamp tiny negative values that can arise from numerical noise
    V_total_sq = np.maximum(V_total_sq, 1e-6)
    return np.sqrt(V_total_sq)


def chi2_data(params: np.ndarray, n_shells: int, data: dict,
              mask: np.ndarray) -> float:
    """χ² on velocity residuals only (no prior). Used for BIC."""
    V_pred = v_model_total(params, n_shells, data)
    residuals = (data["V_obs"][mask] - V_pred[mask]) / data["sigma_V"][mask]
    return float(np.sum(residuals ** 2))


def prior_penalty(params: np.ndarray, fix_upsilon: bool = False) -> float:
    """Gaussian log-prior penalty on Υ_disk and Υ_bulge.

    Returns the −2 log P contribution from the prior, i.e., the χ²-equivalent
    additive term: ((log10(Y) - log10(μ)) / σ_dex)²

    Returns 0 when fix_upsilon=True (canonical-recheck mode: Υ values are
    pinned at canonical and the prior is irrelevant).
    """
    if fix_upsilon:
        return 0.0
    Y_disk = max(params[0], 1e-6)
    Y_bulge = max(params[1], 1e-6)
    log_d = np.log10(Y_disk)
    log_b = np.log10(Y_bulge)
    return (
        ((log_d - np.log10(PRIOR_UPS_DISK_MEAN)) / PRIOR_UPS_DISK_SIGMA_DEX) ** 2
        + ((log_b - np.log10(PRIOR_UPS_BULGE_MEAN)) / PRIOR_UPS_BULGE_SIGMA_DEX) ** 2
    )


def loss_for_fitting(params: np.ndarray, n_shells: int, data: dict,
                     mask: np.ndarray, fix_upsilon: bool = False) -> float:
    """Total loss for MAP estimation: χ²_data + prior_penalty."""
    return chi2_data(params, n_shells, data, mask) + prior_penalty(params, fix_upsilon)


# =============================================================================
# MULTI-RESTART FIT
# =============================================================================

def make_bounds(n_shells: int, fix_upsilon: bool = False) -> list:
    if fix_upsilon:
        # Pin Υ to canonical via tight bounds. L-BFGS-B needs nonzero width
        # for finite-difference gradient computation, so use ±1e-6 instead
        # of a true equality constraint.
        eps = 1.0e-6
        lower = [PRIOR_UPS_DISK_MEAN - eps, PRIOR_UPS_BULGE_MEAN - eps,
                 RHO0_MIN, A_MIN]
        upper = [PRIOR_UPS_DISK_MEAN + eps, PRIOR_UPS_BULGE_MEAN + eps,
                 RHO0_MAX, A_MAX]
    else:
        lower = [UPS_DISK_BOX_MIN, UPS_BULGE_BOX_MIN, RHO0_MIN, A_MIN]
        upper = [UPS_DISK_BOX_MAX, UPS_BULGE_BOX_MAX, RHO0_MAX, A_MAX]
    for _ in range(n_shells):
        lower.extend([SHELL_M_MIN, SHELL_R_MIN, SHELL_F_MIN])
        upper.extend([SHELL_M_MAX, SHELL_R_MAX, SHELL_F_MAX])
    return list(zip(lower, upper))


def draw_initial(n_shells: int, rmax_kpc: float, rng: np.random.Generator,
                 fix_upsilon: bool = False,
                 residual_seed: Optional[tuple] = None) -> np.ndarray:
    """Random initial-condition vector. Υ start near the prior mean to keep
    optimizer in a sensible region; the prior penalty handles the rest.

    When fix_upsilon=True, Υ initials are pinned at canonical values exactly
    so they sit inside the tight ±1e-6 bounds from make_bounds.

    When residual_seed is provided as (r_seed_kpc, M_seed_Msun), the FIRST
    shell's initial conditions are set deterministically from that seed and
    only its σ/r and the Burkert/Υ parameters are randomized. Used to bias
    one restart per rmax toward the radius of maximum 0-shell residual.
    """
    if fix_upsilon:
        Y_disk_init = PRIOR_UPS_DISK_MEAN
        Y_bulge_init = PRIOR_UPS_BULGE_MEAN
    else:
        Y_disk_init = 10.0 ** (np.log10(PRIOR_UPS_DISK_MEAN) + rng.uniform(-0.15, 0.15))
        Y_bulge_init = 10.0 ** (np.log10(PRIOR_UPS_BULGE_MEAN) + rng.uniform(-0.15, 0.15))
    rho_0_init = 10.0 ** rng.uniform(5.0, 9.0)
    a_init = rng.uniform(1.0, 20.0)
    p0 = [Y_disk_init, Y_bulge_init, rho_0_init, a_init]
    # Linear-uniform M and r initial conditions over wider ranges so the
    # initial-condition prior covers the SHELL_M_MAX scale and r up to the
    # full rmax_kpc. (Log-uniform biases toward small r, which is the wrong
    # direction when the basin of interest sits at large r.)
    M_init_log_lo, M_init_log_hi = 7.0, np.log10(SHELL_M_MAX)
    for i in range(n_shells):
        if i == 0 and residual_seed is not None:
            r_seed, M_seed = residual_seed
            # Wide log perturbation on the seed: r * 10^[-0.3,+0.3] covers
            # factor 0.5x-2x in radius, and M * 10^[-0.7,+0.7] covers factor
            # 0.2x-5x in mass. This lets one seeded restart escape to a
            # nearby basin even if the residual peak isn't exactly at the
            # mass center.
            r_init = r_seed * 10.0 ** rng.uniform(-0.3, 0.3)
            r_init = float(np.clip(r_init, SHELL_R_MIN, SHELL_R_MAX))
            M_init = M_seed * 10.0 ** rng.uniform(-0.7, 0.7)
            M_init = float(np.clip(M_init, SHELL_M_MIN, SHELL_M_MAX))
            f_init = rng.uniform(0.10, 0.35)
        else:
            M_init = 10.0 ** rng.uniform(M_init_log_lo, M_init_log_hi)
            r_init = rng.uniform(SHELL_R_MIN, max(rmax_kpc, SHELL_R_MIN * 2))
            f_init = rng.uniform(0.05, 0.35)
        p0.extend([M_init, r_init, f_init])
    return np.array(p0)


def _compute_residual_seed(zero_shell_params: np.ndarray, data: dict,
                           mask: np.ndarray) -> Optional[tuple]:
    """From a 0-shell best fit, find the radius of maximum |V_obs - V_pred|/σ_V
    residual and estimate the shell mass needed to absorb it.

    Returns (r_seed_kpc, M_seed_Msun) or None if no clean residual is found.

    Mass estimate uses the V² residual at the peak radius:
        ΔV² ≈ G · M_shell / r_seed   →   M_seed ≈ ΔV² · r_seed / G
    where ΔV² = V_obs² − V_pred² at the peak-residual radius.
    """
    V_pred = v_model_total(zero_shell_params, 0, data)
    r = data["r_kpc"]
    V_obs = data["V_obs"]
    sigma_V = data["sigma_V"]
    # Work in velocity-residual space normalized by σ_V
    normalized_residual = np.abs(V_obs - V_pred) / sigma_V
    # Only consider unmasked points
    normalized_residual = np.where(mask, normalized_residual, -np.inf)
    if not np.any(np.isfinite(normalized_residual)):
        return None
    idx = int(np.argmax(normalized_residual))
    r_seed = float(r[idx])
    # V² residual at the peak. Use signed residual to estimate mass direction:
    # if V_obs > V_pred, need extra mass (positive shell); if V_obs < V_pred,
    # the smooth profile is over-predicting and a "negative" shell isn't in
    # our parameterization, so skip the seed in that case.
    dV_sq = V_obs[idx] ** 2 - V_pred[idx] ** 2
    if dV_sq <= 0:
        return None
    M_seed = dV_sq * r_seed / G
    M_seed = float(np.clip(M_seed, SHELL_M_MIN, SHELL_M_MAX))
    return r_seed, M_seed


def fit_n_shells_marginalized(n_shells: int, data: dict, mask: np.ndarray,
                              n_restarts_per_rmax: int, rng: np.random.Generator,
                              fix_upsilon: bool = False,
                              residual_seed: Optional[tuple] = None,
                              use_differential_evolution: bool = True,
                              verbose: bool = True) -> Optional[dict]:
    """Multi-restart fit for one n_shells configuration with marginalized Υ.

    Returns a dict with the best-fit parameters, χ², BIC, and metadata, or
    None if all restarts failed.

    When fix_upsilon=True, runs the canonical-recheck mode: Υ pinned at
    (0.5, 0.7) via tight bounds, no prior penalty, and parameter count in
    BIC adjusted to match canonical Table 5 conventions.

    When residual_seed=(r_seed, M_seed) is provided (for n_shells ≥ 1), a
    fraction of restarts use the residual-derived seed for the first shell's
    (r, M) initial values; the rest are random.

    When use_differential_evolution=True and n_shells ≥ 1, scipy's
    differential_evolution global optimizer is run before the multi-restart
    L-BFGS-B passes. DE is robust to local minima but slower; for shell-
    bearing fits its result is then polished with L-BFGS-B. The best of all
    approaches is retained.
    """
    bounds = make_bounds(n_shells, fix_upsilon)
    lower_arr = np.array([b[0] for b in bounds])
    upper_arr = np.array([b[1] for b in bounds])

    best_loss = np.inf
    best_params = None
    n_success = 0
    n_attempts = 0
    n_seeded = 0
    used_de = False

    # --- Stage 1: differential_evolution global pass for ALL fits ---
    # Enabled for n_shells=0 too so 0-shell and shell-bearing fits get the
    # same optimizer treatment. Without this, 0-shell underexplores while
    # shell-bearing fits use DE — an asymmetry that inflates the apparent
    # shell preference because the additional shell parameters give DE more
    # dimensions to search even when the resulting shell is inert.
    if use_differential_evolution:
        try:
            # DE explores the full bounded space globally with a stochastic
            # population search; insensitive to initial-condition basin choice.
            # Use a deterministic seed derived from the parent rng for
            # reproducibility.
            de_seed = int(rng.integers(0, 2**31 - 1))
            de_result = optimize.differential_evolution(
                loss_for_fitting,
                bounds=bounds,
                args=(n_shells, data, mask, fix_upsilon),
                seed=de_seed,
                popsize=20,
                maxiter=500,
                tol=1e-8,
                mutation=(0.5, 1.5),
                recombination=0.8,
                polish=True,   # internal L-BFGS-B polish
                workers=1,
                updating="deferred",
            )
            if de_result.success or de_result.fun < np.inf:
                n_success += 1
                used_de = True
                if de_result.fun < best_loss:
                    best_loss = de_result.fun
                    best_params = de_result.x.copy()
        except Exception:
            pass

    # --- Stage 2: residual-seeded + random multi-restart L-BFGS-B ---
    # Use residual-seed for ~half of restarts when available
    seeded_restart_indices = (set(range(n_restarts_per_rmax // 2))
                              if residual_seed is not None else set())

    for rmax_kpc in RMAX_OPTIONS:
        for restart_idx in range(n_restarts_per_rmax):
            n_attempts += 1
            use_seed = (restart_idx in seeded_restart_indices
                        and residual_seed is not None)
            seed_arg = residual_seed if use_seed else None
            if use_seed:
                n_seeded += 1
            p0 = draw_initial(n_shells, rmax_kpc, rng, fix_upsilon,
                              residual_seed=seed_arg)
            p0 = np.clip(p0, lower_arr, upper_arr)

            try:
                result = optimize.minimize(
                    loss_for_fitting,
                    p0,
                    args=(n_shells, data, mask, fix_upsilon),
                    method="L-BFGS-B",
                    bounds=bounds,
                    options={
                        "maxiter": MAXFEV,
                        "ftol": FTOL,
                        "gtol": GTOL,
                    },
                )
                if result.success:
                    n_success += 1
                    if result.fun < best_loss:
                        best_loss = result.fun
                        best_params = result.x.copy()
            except Exception:
                continue

    if best_params is None:
        if verbose:
            print(f"    WARNING: no successful fit for n_shells={n_shells}")
        return None

    # Recompute χ²_data separately (without prior) for BIC reporting
    chi2_d = chi2_data(best_params, n_shells, data, mask)
    prior_val = prior_penalty(best_params, fix_upsilon)
    n_data = int(mask.sum())
    p = n_params(n_shells, fix_upsilon)
    dof = max(n_data - p, 1)
    chi2_red = chi2_d / dof
    bic = chi2_d + p * np.log(n_data)

    Y_disk, Y_bulge, rho_0, a_kpc, shells = unpack_params(best_params, n_shells)

    out = {
        "n_shells": n_shells,
        "n_params": p,
        "n_data": n_data,
        "dof": dof,
        "chi2_data": chi2_d,
        "chi2_red": chi2_red,
        "prior_penalty": prior_val,
        "bic": bic,
        "Y_disk_fit": float(Y_disk),
        "Y_bulge_fit": float(Y_bulge),
        "rho_0_Msun_kpc3": float(rho_0),
        "a_kpc": float(a_kpc),
        "n_attempts": n_attempts,
        "n_success": n_success,
        "n_seeded": n_seeded,
        "used_differential_evolution": used_de,
        "fix_upsilon": fix_upsilon,
        "_best_params": best_params,  # internal: used to seed higher-shell fits
    }
    for i, (M_i, r_i, sigma_i) in enumerate(shells, start=1):
        out[f"M_shell_{i}_Msun"] = float(M_i)
        out[f"r_shell_{i}_kpc"] = float(r_i)
        out[f"sigma_shell_{i}_kpc"] = float(sigma_i)
        out[f"sigma_over_r_{i}"] = float(sigma_i / r_i)

    return out


# =============================================================================
# REPORTING
# =============================================================================

def interpret_dbic(delta_bic: float) -> str:
    """Kass & Raftery 1995 verdict for ΔBIC = BIC_alt − BIC_best."""
    abs_d = abs(delta_bic)
    if abs_d > DBIC_THRESHOLDS["very_strong"]:
        return "very strong"
    if abs_d > DBIC_THRESHOLDS["strong"]:
        return "strong"
    if abs_d > DBIC_THRESHOLDS["positive"]:
        return "positive"
    return "inconclusive"


def print_summary(galaxy: str, results: list[dict],
                  canonical: Optional[dict], fix_upsilon: bool = False) -> None:
    mode_label = "fixed-Υ canonical recheck" if fix_upsilon else "marginalized-Υ"
    bic_col_label = "BIC" if fix_upsilon else "BIC_marg"
    print()
    print("=" * 78)
    print(f"  {galaxy}: {mode_label} framework fits")
    print("=" * 78)
    if fix_upsilon:
        print(f"  Υ_disk pinned at canonical 0.5; Υ_bulge pinned at canonical 0.7")
        print(f"  Prior penalty disabled. Parameter count matches Table 5.")
    else:
        print(f"  Priors: log10(Υ_disk)  ~ N(log10(0.5), 0.1)")
        print(f"          log10(Υ_bulge) ~ N(log10(0.7), 0.1)")
    if results:
        print(f"  n_data = {results[0]['n_data']}  (mask: V_bar²(canonical Υ) < V_obs²)")
    print()
    print(f"  {'n_shells':>8}  {'χ²_data':>10}  {'χ²_red':>8}  {'p':>3}  "
          f"{bic_col_label:>10}  {'Υ_disk':>7}  {'Υ_bulge':>7}")
    print("  " + "-" * 74)
    for r in results:
        print(f"  {r['n_shells']:>8}  {r['chi2_data']:>10.2f}  "
              f"{r['chi2_red']:>8.3f}  {r['n_params']:>3}  "
              f"{r['bic']:>10.2f}  "
              f"{r['Y_disk_fit']:>7.3f}  {r['Y_bulge_fit']:>7.3f}")

    # ΔBIC table relative to best fit in this mode
    if len(results) >= 2:
        bics = [r["bic"] for r in results]
        bic_best = min(bics)
        best_n = results[bics.index(bic_best)]["n_shells"]
        print()
        print(f"  Best {mode_label} fit: n_shells = {best_n}, "
              f"BIC = {bic_best:.2f}")
        for r in results:
            d = r["bic"] - bic_best
            verdict = interpret_dbic(d) if d > 0 else "best"
            print(f"    n_shells={r['n_shells']}: ΔBIC = {d:+8.2f}  "
                  f"({verdict})")

    # Comparison to canonical Table 5
    if canonical is not None and results:
        print()
        if fix_upsilon:
            print("  Canonical Table 5 (reference) vs this run's fixed-Υ recheck:")
            bic_run_label = "BIC_run"
        else:
            print("  Canonical (fixed-Υ Table 5) vs marginalized BIC:")
            bic_run_label = "BIC_marg"
        print(f"    {'Model':<22}  {'BIC_canon':>10}  {bic_run_label:>10}  "
              f"{'ΔBIC':>8}")
        print("    " + "-" * 56)
        for r in results:
            if r["n_shells"] == 0:
                bic_canon = canonical["burkert_only_bic"]
                label = "Burkert-only"
            elif r["n_shells"] == 1:
                bic_canon = canonical["framework_1shell_bic"]
                label = "Framework (1 shell)"
            else:
                continue
            d = r["bic"] - bic_canon
            print(f"    {label:<22}  {bic_canon:>10.2f}  "
                  f"{r['bic']:>10.2f}  {d:>+8.2f}")
        if fix_upsilon:
            print()
            print("    SANITY-CHECK INTERPRETATION:")
            print("    With Υ pinned at canonical, parameter count matches Table 5")
            print("    exactly, so BIC_run should approximately equal BIC_canon.")
            print("    A large positive ΔBIC means our optimizer is finding a")
            print("    worse minimum than canonical (multi-restart problem).")
            print("    A large negative ΔBIC means the canonical Table 5 number")
            print("    is overstated (canonical pipeline didn't find this minimum).")
            print("    |ΔBIC| < ~5 is the expected sanity-check pass band.")
        else:
            print()
            print("    Expected BIC cost of marginalization (per extra parameter):")
            print(f"      ln(n)  = ln({canonical['n_data']}) = "
                  f"{np.log(canonical['n_data']):.3f}")
            print(f"      2·ln(n) for two extra Υ parameters = "
                  f"{2*np.log(canonical['n_data']):.3f}")
            print("    A marginalized BIC that exceeds canonical by ≈ 2·ln(n) and")
            print("    no more indicates the fit prefers Υ near the prior mean")
            print("    (no information gained from freeing Υ).")

    # Verdict for §3.4 — only show in marginalized mode; sanity check shows
    # a different conclusion bucket.
    if len(results) >= 2:
        bic_0 = next((r["bic"] for r in results
                      if r["n_shells"] == 0), None)
        bic_1 = next((r["bic"] for r in results
                      if r["n_shells"] == 1), None)
        if bic_0 is not None and bic_1 is not None:
            d = bic_0 - bic_1
            print()
            print("  -" * 38)
            if fix_upsilon:
                print(f"  Sanity-check verdict for {galaxy} (fixed-Υ recheck):")
                print(f"  ΔBIC(Burkert-only − Framework-1shell) = {d:+.2f}")
                print(f"  Canonical Table 5 value: {canonical['burkert_only_bic'] - canonical['framework_1shell_bic']:+.2f}")
                print(f"  If these two values agree to within a few units, the")
                print(f"  script reproduces the canonical pipeline. If not,")
                print(f"  investigate model / mask / multi-restart conventions.")
            else:
                print(f"  §3.4 verdict for {galaxy}:")
                print(f"  ΔBIC(Burkert-only_marg − Framework-1shell_marg) = "
                      f"{d:+.2f}")
                if d > DBIC_THRESHOLDS["very_strong"]:
                    print(f"  -> Framework still VERY STRONGLY preferred under")
                    print(f"     marginalized Υ. §3.4 showcase holds.")
                elif d > DBIC_THRESHOLDS["strong"]:
                    print(f"  -> Framework STRONGLY preferred under marginalized")
                    print(f"     Υ. §3.4 showcase holds.")
                elif d > DBIC_THRESHOLDS["positive"]:
                    print(f"  -> Framework POSITIVELY preferred under marginalized")
                    print(f"     Υ. §3.4 showcase weakened but standing.")
                elif d > -DBIC_THRESHOLDS["positive"]:
                    print(f"  -> INCONCLUSIVE between Burkert-only and Framework")
                    print(f"     under marginalized Υ. §3.4 needs reframing.")
                else:
                    print(f"  -> Burkert-only PREFERRED under marginalized Υ.")
                    print(f"     §3.4 showcase needs replacement or reframe.")
    print("=" * 78)
    print()


# =============================================================================
# MAIN
# =============================================================================

def run(galaxy: str, data_dir: str, out_dir: str,
        n_restarts_per_rmax: int, random_seed: int,
        fix_upsilon: bool = False, verbose: bool = True
        ) -> pd.DataFrame:
    rng = np.random.default_rng(random_seed)
    mode_label = "fixed-Υ canonical recheck" if fix_upsilon else "marginalized Υ"

    filepath = os.path.join(data_dir, f"{galaxy}_rotmod.dat")
    if verbose:
        print(f"Loading {galaxy} from {filepath}")
    data = load_rotmod(filepath)
    if verbose:
        print(f"  Distance = {data['distance_mpc']} Mpc, "
              f"{len(data['r_kpc'])} data points before masking")
        print(f"  Mode: {mode_label}")

    # Mask: exclude points where V_bar²(canonical Υ) > V_obs²
    # Fixed mask preserves direct comparability to canonical Table 5 BIC values.
    V_bar_canon_sq = v_baryon_squared(
        PRIOR_UPS_DISK_MEAN, PRIOR_UPS_BULGE_MEAN,
        data["V_gas"], data["V_disk"], data["V_bulge"]
    )
    mask = V_bar_canon_sq < data["V_obs"] ** 2
    if verbose:
        n_kept = int(mask.sum())
        n_drop = len(data["r_kpc"]) - n_kept
        print(f"  Mask at canonical Υ: {n_kept} points kept, "
              f"{n_drop} dropped (V_bar² > V_obs²)")

    results = []
    residual_seed: Optional[tuple] = None
    for n_shells in (0, 1, 2):
        if verbose:
            print()
            print(f"Fitting {galaxy} with {n_shells} shell(s), {mode_label} ...")
            seed_msg = ""
            if residual_seed is not None and n_shells >= 1:
                seed_msg = (f"  Residual-seeded restarts active: "
                            f"r_seed = {residual_seed[0]:.2f} kpc, "
                            f"M_seed = {residual_seed[1]:.2e} M☉")
            print(f"  Multi-restart: {n_restarts_per_rmax} per r_max × "
                  f"{len(RMAX_OPTIONS)} r_max values = "
                  f"{n_restarts_per_rmax * len(RMAX_OPTIONS)} total restarts")
            if seed_msg:
                print(seed_msg)
        res = fit_n_shells_marginalized(
            n_shells, data, mask, n_restarts_per_rmax, rng,
            fix_upsilon=fix_upsilon,
            residual_seed=residual_seed,
            verbose=verbose,
        )
        if res is None:
            continue
        if verbose:
            bic_label = "BIC" if fix_upsilon else "BIC_marg"
            print(f"  χ²_data = {res['chi2_data']:.2f}, "
                  f"χ²_red = {res['chi2_red']:.3f}, "
                  f"{bic_label} = {res['bic']:.2f}")
            print(f"  Υ_disk = {res['Y_disk_fit']:.3f}, "
                  f"Υ_bulge = {res['Y_bulge_fit']:.3f}")
            for i in range(1, n_shells + 1):
                print(f"  Shell {i}: M = {res[f'M_shell_{i}_Msun']:.2e} M☉, "
                      f"r = {res[f'r_shell_{i}_kpc']:.2f} kpc, "
                      f"σ = {res[f'sigma_shell_{i}_kpc']:.2f} kpc, "
                      f"σ/r = {res[f'sigma_over_r_{i}']:.3f}")
        # After 0-shell fit lands, compute residual-based seed for n_shells ≥ 1
        if n_shells == 0 and "_best_params" in res:
            try:
                residual_seed = _compute_residual_seed(
                    res["_best_params"], data, mask
                )
                if verbose and residual_seed is not None:
                    print(f"  Residual seed for higher-shell fits: "
                          f"r = {residual_seed[0]:.2f} kpc, "
                          f"M ≈ {residual_seed[1]:.2e} M☉")
            except Exception:
                residual_seed = None
        results.append(res)

    if not results:
        print("ERROR: no successful fits.", file=sys.stderr)
        sys.exit(2)

    # Strip the internal _best_params key (not for CSV output)
    for r in results:
        r.pop("_best_params", None)

    # Append ΔBIC columns
    bic_best = min(r["bic"] for r in results)
    delta_best_col = ("delta_bic_vs_best_fixed" if fix_upsilon
                      else "delta_bic_vs_best_marg")
    for r in results:
        r[delta_best_col] = r["bic"] - bic_best
        canonical = CANONICAL_REFERENCE.get(galaxy)
        if canonical:
            if r["n_shells"] == 0:
                r["canonical_bic"] = canonical["burkert_only_bic"]
                r["canonical_chi2_red"] = canonical["burkert_only_chi2_red"]
            elif r["n_shells"] == 1:
                r["canonical_bic"] = canonical["framework_1shell_bic"]
                r["canonical_chi2_red"] = canonical["framework_1shell_chi2_red"]
            else:
                r["canonical_bic"] = np.nan
                r["canonical_chi2_red"] = np.nan
            r["delta_bic_vs_canonical"] = r["bic"] - r["canonical_bic"]
        r["galaxy"] = galaxy

    df = pd.DataFrame(results)
    # Move 'galaxy' to the front
    cols = ["galaxy"] + [c for c in df.columns if c != "galaxy"]
    df = df[cols]

    os.makedirs(out_dir, exist_ok=True)
    out_basename = ("{}_fixed_upsilon_recheck_results.csv".format(galaxy.lower())
                    if fix_upsilon
                    else "{}_marginalized_upsilon_results.csv".format(galaxy.lower()))
    out_path = os.path.join(out_dir, out_basename)
    df.to_csv(out_path, index=False)
    if verbose:
        print(f"\nResults written to {out_path}")

    if verbose:
        canonical = CANONICAL_REFERENCE.get(galaxy)
        print_summary(galaxy, results, canonical, fix_upsilon=fix_upsilon)

    return df


def main():
    parser = argparse.ArgumentParser(
        description=("NGC 5055 framework fit with marginalized Υ_disk and "
                     "Υ_bulge under Gaussian log-priors (Li 2020 convention). "
                     "Addresses §4.4 marginalization gap raised in PASP "
                     "resubmission review."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--galaxy", default=DEFAULT_GALAXY,
                        help=(f"Galaxy name to fit (default: {DEFAULT_GALAXY}). "
                              "Must have a matching <galaxy>_rotmod.dat file "
                              "in --data-dir."))
    parser.add_argument("--data-dir", default=DEFAULT_DATA_DIR,
                        help=(f"Directory containing rotmod files "
                              f"(default: {DEFAULT_DATA_DIR})"))
    parser.add_argument("--out-dir", default="./data",
                        help="Output directory for CSV (default: ./data)")
    parser.add_argument("--n-restarts-per-rmax", type=int,
                        default=DEFAULT_N_RESTARTS_PER_RMAX,
                        help=(f"Restarts per r_max value (default: "
                              f"{DEFAULT_N_RESTARTS_PER_RMAX}; total restarts "
                              f"= 3× this)"))
    parser.add_argument("--random-seed", type=int, default=DEFAULT_RANDOM_SEED,
                        help=(f"NumPy RNG seed for reproducibility "
                              f"(default: {DEFAULT_RANDOM_SEED})"))
    parser.add_argument("--quiet", action="store_true",
                        help="Suppress verbose progress output")
    parser.add_argument("--fix-upsilon", action="store_true",
                        help=("SANITY-CHECK MODE: pin Υ_disk=0.5 and "
                              "Υ_bulge=0.7 (canonical values), disable the "
                              "Gaussian log-prior, and use the canonical "
                              "Table 5 parameter count for BIC. Use this to "
                              "verify the script reproduces canonical BIC "
                              "values before trusting the marginalized run. "
                              "Output filename becomes <galaxy>_fixed_"
                              "upsilon_recheck_results.csv to avoid "
                              "overwriting the marginalized output."))
    args = parser.parse_args()

    run(
        galaxy=args.galaxy,
        data_dir=args.data_dir,
        out_dir=args.out_dir,
        n_restarts_per_rmax=args.n_restarts_per_rmax,
        random_seed=args.random_seed,
        fix_upsilon=args.fix_upsilon,
        verbose=not args.quiet,
    )


if __name__ == "__main__":
    main()
