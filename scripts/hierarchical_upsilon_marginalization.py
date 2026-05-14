#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hierarchical_upsilon_marginalization.py

Empirical Bayes hierarchical marginalization of Υ_disk and Υ_bulge across the
full SPARC T=2--9 sample (102 galaxies, optionally 101 with NGC 6674 excluded
per Paper II §2.3 convention).

THE QUESTION THIS ANSWERS

Paper I §3.4 and §4.4 disclosed that the per-galaxy BIC preference for shells
is sensitive to joint Υ marginalization (NGC 5055: ΔBIC drops from 157 to
+0.5 under Li 2020 priors). The natural follow-up is whether the population-
level headline numbers change when Υ is marginalized hierarchically across
the full sample:

    (a) Does the framework's adequacy advantage (91/102 vs 65/102) compress
        or hold under marginalization?
    (b) Does the morphology gradient survive marginalization (Paper II §3.1)?
    (c) Where does the population-level Υ_disk distribution converge?

WHAT THIS SCRIPT DOES

For each iteration:
  1. Fit each galaxy with framework at N_shells ∈ {0, 1, 2}, with Υ_disk and
     Υ_bulge jointly free under the current hyperprior:
         log10(Υ_disk_g)  ~ N(μ_disk,  τ_disk)
         log10(Υ_bulge_g) ~ N(μ_bulge, τ_bulge)
  2. Select best N per galaxy by BIC; record per-galaxy fitted Υ values.
  3. Aggregate per-galaxy fitted Υ into population-level estimates:
         μ_disk_new = mean(log10(Υ_disk_fit)) across all galaxies
         τ_disk_new = std(log10(Υ_disk_fit))  across all galaxies
         μ_bulge_new, τ_bulge_new = same, but only over bulged galaxies
         (bulgeless galaxies contribute no information about Υ_bulge)
  4. Check convergence: |Δμ| < 0.01 dex AND |Δτ| < 0.01 dex.
  5. Update hyperprior; repeat.

Initial hyperprior is the Li 2020 convention: (μ_disk, τ_disk) = (log10(0.5), 0.1),
(μ_bulge, τ_bulge) = (log10(0.7), 0.1). This matches Paper I's NGC 5055 test
exactly, so the hierarchical result reduces to Paper I's single-galaxy
marginalization on iteration 1 with no shrinkage.

EMPIRICAL BAYES vs FULL HIERARCHICAL BAYES

This is empirical Bayes (iterated MAP with point-estimate hyperparameters),
not full hierarchical Bayes. The point estimates of μ and τ converge to the
mode of the marginal posterior in the limit of weak hyperpriors, which is what
we have here. For full uncertainty quantification on the hyperparameters, a
PyMC or numpyro implementation would be required; that's a separate piece of
work and is not necessary for the headline questions (a)-(c) above.

PARTIAL OUTPUT / RESUMABILITY

Per-galaxy results are appended to `data/hierarchical_marginalization_per_galaxy.csv`
immediately after each galaxy completes. If the script is interrupted, re-run
with `--resume` to continue from where it left off. Hyperparameter history is
written to `data/hierarchical_marginalization_hyperprior_history.csv` after
each iteration.

The per-galaxy CSV is inspectable at any time during the run.

USAGE

    python scripts/hierarchical_upsilon_marginalization.py
    python scripts/hierarchical_upsilon_marginalization.py --data-dir ./Rotmod_LTG
    python scripts/hierarchical_upsilon_marginalization.py --sample-csv ./data/sparc_sample123.csv
    python scripts/hierarchical_upsilon_marginalization.py --resume
    python scripts/hierarchical_upsilon_marginalization.py --max-galaxies 20  # for quick smoke-test
    python scripts/hierarchical_upsilon_marginalization.py --fast              # no DE, faster but less reliable

ESTIMATED RUNTIME (M-series Mac)

    Per galaxy: ~30-60s (DE + L-BFGS-B multi-restart for N=0,1,2)
    Per iteration: ~60-90 minutes for 102 galaxies
    Total (3 iterations to convergence): ~3-5 hours
    With --fast (no DE): ~1-2 hours total, but less reliable on shell-bearing fits

OUTPUTS

    data/hierarchical_marginalization_per_galaxy.csv
        One row per (galaxy, iteration). Schema:
        - galaxy, iteration, T, n_data, status
        - For each N in {0, 1, 2}: chi2_N, chi2_red_N, bic_N,
          Y_disk_N, Y_bulge_N, [shell params if N>=1]
        - best_N, best_bic
        - delta_bic_framework_vs_burkert (best_bic_with_shells - bic_N0)
        - prior_penalty_at_best
        - is_adequate (chi2_red < 1.5 at best_N)
        - is_shell_bearing (best_N >= 1)

    data/hierarchical_marginalization_hyperprior_history.csv
        One row per iteration. Schema:
        - iteration, mu_disk, tau_disk, mu_bulge, tau_bulge,
          n_galaxies_fitted, n_bulged_galaxies, delta_mu_disk, delta_tau_disk,
          delta_mu_bulge, delta_tau_bulge, converged

REPRODUCIBILITY

    Random seed 20260513 by default. With --resume, the rng is re-seeded
    from the same value but the per-galaxy results already in the CSV are
    not re-fit, so the strict reproducibility holds only if you complete a
    full run without interruption.

Author: Ronald Bibb
Created: 2026-05-13 for hierarchical Υ marginalization follow-up to Paper I.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import dataclass
from typing import Optional, Tuple, Dict, List

import numpy as np
import pandas as pd
from scipy import optimize
from scipy.special import erf


# =============================================================================
# CONSTANTS
# =============================================================================

G = 4.30091e-6  # gravitational constant in (km/s)² · kpc / M_sun

# Shell parameter bounds (match Paper I §2.2)
SHELL_M_MIN = 1.0e6
SHELL_M_MAX = 5.0e10
SHELL_F_MIN = 0.01
SHELL_F_MAX = 0.4
SHELL_R_MIN = 0.5
SHELL_R_MAX = 50.0

# Burkert backbone bounds
RHO0_MIN, RHO0_MAX = 1.0e3, 1.0e10
A_MIN, A_MAX = 0.1, 50.0

# Υ search box (well outside any realistic prior)
UPS_DISK_BOX_MIN, UPS_DISK_BOX_MAX = 0.1, 2.0
UPS_BULGE_BOX_MIN, UPS_BULGE_BOX_MAX = 0.1, 2.0

# Hierarchical hyperprior initialization (Li 2020 convention)
INITIAL_MU_DISK = float(np.log10(0.5))   # ≈ -0.301
INITIAL_TAU_DISK_DEX = 0.1
INITIAL_MU_BULGE = float(np.log10(0.7))  # ≈ -0.155
INITIAL_TAU_BULGE_DEX = 0.1

# Convergence
HYPER_CONVERGENCE_TOL_DEX = 0.01     # |Δμ|, |Δτ| convergence threshold
MAX_HYPER_ITERATIONS = 5

# Hyperprior safety bounds (prevent runaway)
TAU_MIN = 0.02
TAU_MAX = 0.3
MU_DISK_MIN, MU_DISK_MAX = np.log10(0.2), np.log10(1.5)
MU_BULGE_MIN, MU_BULGE_MAX = np.log10(0.3), np.log10(1.5)

# Multi-restart fit configuration
DEFAULT_N_RESTARTS_PER_RMAX = 4   # reduced from Paper I's 8 to save time
RMAX_OPTIONS = (3.0, 8.0, 15.0, 25.0)
MAXFEV = 12000
FTOL = 1e-9
GTOL = 1e-7

# DE configuration (used for n_shells >= 1 by default)
DE_POPSIZE = 18
DE_MAXITER = 300
DE_TOL = 1e-7

# Data conventions
SIGMA_V_FLOOR_KMS = 1.0  # per Paper I §2.3
PRIOR_UPS_DISK_MEAN = 0.5
PRIOR_UPS_BULGE_MEAN = 0.7

# Quality cuts (Paper I §2.1, also matches Paper II §2.1)
T_MIN, T_MAX = 2, 9
Q_MAX = 2
INC_MIN_DEG = 30.0
V_FLAT_MIN_KMS = 25.0
MIN_USABLE_POINTS = 5

# Adequacy threshold
CHI2_RED_ADEQUACY = 1.5

# Reproducibility
DEFAULT_RANDOM_SEED = 20260513

# NGC 6674 convention (Paper II §2.3 excludes due to degenerate two-shell fit)
NGC6674_EXCLUDE_DEFAULT = True


# =============================================================================
# DATA LOADING
# =============================================================================

def load_rotmod(filepath: str) -> dict:
    """Load a SPARC rotmod file. Returns dict of arrays.

    Columns (SPARC convention, Lelli et al. 2016):
        Rad(kpc) Vobs(km/s) errV(km/s) Vgas Vdisk Vbul SBdisk SBbul
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Cannot find rotmod: {filepath}")

    distance_mpc = None
    with open(filepath, "r") as f:
        for line in f:
            if line.lower().startswith("# distance"):
                parts = line.split("=")
                if len(parts) >= 2:
                    try:
                        distance_mpc = float(parts[1].strip().split()[0])
                    except (ValueError, IndexError):
                        pass
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


def load_sample(
    data_dir: str,
    sample_csv_path: Optional[str] = None,
    exclude_ngc6674: bool = True,
    max_galaxies: Optional[int] = None,
) -> List[Dict]:
    """Load the SPARC sample, filtered by Paper I quality cuts.

    If `sample_csv_path` is provided and exists, read galaxy list and metadata
    (including T-type and bulge presence) from it. Otherwise fall back to
    scanning data_dir for *_rotmod.dat files (T-type will be unknown in this
    fallback case).
    """
    galaxies = []

    if sample_csv_path and os.path.exists(sample_csv_path):
        df = pd.read_csv(sample_csv_path)
        # Apply Paper I quality cuts
        mask = (
            (df["T"] >= T_MIN) & (df["T"] <= T_MAX)
            & (df["Q"] <= Q_MAX)
            & (df["Inc"] >= INC_MIN_DEG)
            & (df["Vflat"] >= V_FLAT_MIN_KMS)
        )
        df = df[mask].copy()
        if exclude_ngc6674:
            df = df[df["Galaxy"] != "NGC6674"].copy()

        for _, row in df.iterrows():
            gname = str(row["Galaxy"]).strip()
            rotmod_path = os.path.join(data_dir, f"{gname}_rotmod.dat")
            if not os.path.exists(rotmod_path):
                continue
            galaxies.append({
                "name": gname,
                "T": int(row["T"]),
                "rotmod_path": rotmod_path,
                "Vflat": float(row["Vflat"]),
            })
    else:
        # Fallback: scan data_dir
        if sample_csv_path:
            print(f"WARNING: sample_csv_path not found: {sample_csv_path}",
                  file=sys.stderr)
            print("Falling back to filesystem scan; T-type will be unknown.",
                  file=sys.stderr)
        for fname in sorted(os.listdir(data_dir)):
            if not fname.endswith("_rotmod.dat"):
                continue
            gname = fname[:-len("_rotmod.dat")]
            if exclude_ngc6674 and gname == "NGC6674":
                continue
            galaxies.append({
                "name": gname,
                "T": -1,
                "rotmod_path": os.path.join(data_dir, fname),
                "Vflat": -1.0,
            })

    if max_galaxies is not None and max_galaxies > 0:
        galaxies = galaxies[:max_galaxies]

    return galaxies


# =============================================================================
# MODEL COMPONENTS (mirror Paper I)
# =============================================================================

def v_baryon_squared(Y_disk: float, Y_bulge: float,
                     V_gas: np.ndarray, V_disk: np.ndarray, V_bulge: np.ndarray
                     ) -> np.ndarray:
    return (V_gas * np.abs(V_gas)
            + Y_disk * V_disk * np.abs(V_disk)
            + Y_bulge * V_bulge * np.abs(V_bulge))


def v_burkert_squared(r_kpc: np.ndarray, rho_0: float, a_kpc: float) -> np.ndarray:
    x = r_kpc / a_kpc
    M_enc = np.pi * rho_0 * a_kpc**3 * (
        2.0 * np.log(1.0 + x) + np.log(1.0 + x**2) - 2.0 * np.arctan(x)
    )
    return G * M_enc / r_kpc


def v_shell_squared(r_kpc: np.ndarray, M_shell: float, r_shell: float,
                    sigma_shell: float) -> np.ndarray:
    inv = 1.0 / (sigma_shell * np.sqrt(2.0))
    M_enc = 0.5 * M_shell * (
        erf((r_kpc - r_shell) * inv) - erf(-r_shell * inv)
    )
    return G * M_enc / r_kpc


def n_params(n_shells: int) -> int:
    """Total free parameter count for hierarchical-marginalized fit:
    2 Υ + 2 Burkert + 3 per shell."""
    return 4 + 3 * n_shells


def unpack_params(params: np.ndarray, n_shells: int) -> tuple:
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
    Y_disk, Y_bulge, rho_0, a_kpc, shells = unpack_params(params, n_shells)
    V_bar_sq = v_baryon_squared(Y_disk, Y_bulge,
                                 data["V_gas"], data["V_disk"], data["V_bulge"])
    V_dm_sq = v_burkert_squared(data["r_kpc"], rho_0, a_kpc)
    for M_i, r_i, sigma_i in shells:
        V_dm_sq += v_shell_squared(data["r_kpc"], M_i, r_i, sigma_i)
    V_total_sq = V_bar_sq + V_dm_sq
    V_total_sq = np.maximum(V_total_sq, 1e-6)
    return np.sqrt(V_total_sq)


def chi2_data(params: np.ndarray, n_shells: int, data: dict,
              mask: np.ndarray) -> float:
    V_pred = v_model_total(params, n_shells, data)
    res = (data["V_obs"][mask] - V_pred[mask]) / data["sigma_V"][mask]
    return float(np.sum(res ** 2))


# =============================================================================
# HIERARCHICAL PRIOR
# =============================================================================

@dataclass
class Hyperprior:
    mu_disk: float
    tau_disk: float
    mu_bulge: float
    tau_bulge: float

    def as_tuple(self) -> tuple:
        return (self.mu_disk, self.tau_disk, self.mu_bulge, self.tau_bulge)


def prior_penalty_hier(params: np.ndarray, hp: Hyperprior) -> float:
    Y_disk = max(params[0], 1e-6)
    Y_bulge = max(params[1], 1e-6)
    log_d = np.log10(Y_disk)
    log_b = np.log10(Y_bulge)
    return (((log_d - hp.mu_disk) / hp.tau_disk) ** 2
            + ((log_b - hp.mu_bulge) / hp.tau_bulge) ** 2)


def loss_hier(params: np.ndarray, n_shells: int, data: dict,
              mask: np.ndarray, hp: Hyperprior) -> float:
    return chi2_data(params, n_shells, data, mask) + prior_penalty_hier(params, hp)


# =============================================================================
# PER-GALAXY FIT
# =============================================================================

def make_bounds(n_shells: int) -> list:
    lower = [UPS_DISK_BOX_MIN, UPS_BULGE_BOX_MIN, RHO0_MIN, A_MIN]
    upper = [UPS_DISK_BOX_MAX, UPS_BULGE_BOX_MAX, RHO0_MAX, A_MAX]
    for _ in range(n_shells):
        lower.extend([SHELL_M_MIN, SHELL_R_MIN, SHELL_F_MIN])
        upper.extend([SHELL_M_MAX, SHELL_R_MAX, SHELL_F_MAX])
    return list(zip(lower, upper))


def draw_initial(n_shells: int, rmax_kpc: float, rng: np.random.Generator,
                 hp: Hyperprior,
                 residual_seed: Optional[tuple] = None) -> np.ndarray:
    """Random initial-condition vector for L-BFGS-B restart.

    Υ initials are drawn near the current hyperprior mean so the optimizer
    starts in a sensible region; the prior penalty handles the rest.
    """
    Y_disk_init = 10.0 ** (hp.mu_disk + rng.uniform(-hp.tau_disk, hp.tau_disk))
    Y_bulge_init = 10.0 ** (hp.mu_bulge + rng.uniform(-hp.tau_bulge, hp.tau_bulge))
    rho_0_init = 10.0 ** rng.uniform(5.0, 9.0)
    a_init = rng.uniform(1.0, 20.0)
    p0 = [Y_disk_init, Y_bulge_init, rho_0_init, a_init]
    M_init_log_lo, M_init_log_hi = 7.0, np.log10(SHELL_M_MAX)
    for i in range(n_shells):
        if i == 0 and residual_seed is not None:
            r_seed, M_seed = residual_seed
            r_init = float(np.clip(
                r_seed * 10.0 ** rng.uniform(-0.3, 0.3),
                SHELL_R_MIN, SHELL_R_MAX,
            ))
            M_init = float(np.clip(
                M_seed * 10.0 ** rng.uniform(-0.7, 0.7),
                SHELL_M_MIN, SHELL_M_MAX,
            ))
            f_init = rng.uniform(0.10, 0.35)
        else:
            M_init = 10.0 ** rng.uniform(M_init_log_lo, M_init_log_hi)
            r_init = rng.uniform(SHELL_R_MIN, max(rmax_kpc, SHELL_R_MIN * 2))
            f_init = rng.uniform(0.05, 0.35)
        p0.extend([M_init, r_init, f_init])
    return np.array(p0)


def compute_residual_seed(zero_shell_params: np.ndarray, data: dict,
                          mask: np.ndarray) -> Optional[tuple]:
    """From 0-shell best fit, find radius of max |V_obs - V_pred|/σ_V residual
    and estimate seed shell mass from the V² residual at that radius."""
    V_pred = v_model_total(zero_shell_params, 0, data)
    r = data["r_kpc"]
    V_obs = data["V_obs"]
    sigma_V = data["sigma_V"]
    normalized_residual = np.abs(V_obs - V_pred) / sigma_V
    normalized_residual = np.where(mask, normalized_residual, -np.inf)
    if not np.any(np.isfinite(normalized_residual)):
        return None
    idx = int(np.argmax(normalized_residual))
    r_seed = float(r[idx])
    dV_sq = V_obs[idx] ** 2 - V_pred[idx] ** 2
    if dV_sq <= 0:
        return None
    M_seed = float(np.clip(dV_sq * r_seed / G, SHELL_M_MIN, SHELL_M_MAX))
    return r_seed, M_seed


def fit_n_shells(
    n_shells: int,
    data: dict,
    mask: np.ndarray,
    hp: Hyperprior,
    rng: np.random.Generator,
    n_restarts_per_rmax: int,
    use_differential_evolution: bool,
    residual_seed: Optional[tuple] = None,
) -> Optional[dict]:
    """Fit one (n_shells) configuration on one galaxy."""
    bounds = make_bounds(n_shells)
    lower_arr = np.array([b[0] for b in bounds])
    upper_arr = np.array([b[1] for b in bounds])

    best_loss = np.inf
    best_params = None
    n_success = 0
    used_de = False

    # Stage 1: differential_evolution (global) when enabled
    if use_differential_evolution:
        try:
            de_seed = int(rng.integers(0, 2**31 - 1))
            de_result = optimize.differential_evolution(
                loss_hier,
                bounds=bounds,
                args=(n_shells, data, mask, hp),
                seed=de_seed,
                popsize=DE_POPSIZE,
                maxiter=DE_MAXITER,
                tol=DE_TOL,
                mutation=(0.5, 1.5),
                recombination=0.8,
                polish=True,
                workers=1,
                updating="deferred",
            )
            if de_result.fun < np.inf:
                n_success += 1
                used_de = True
                if de_result.fun < best_loss:
                    best_loss = de_result.fun
                    best_params = de_result.x.copy()
        except Exception:
            pass

    # Stage 2: residual-seeded + random multi-restart L-BFGS-B
    seeded_indices = (set(range(max(1, n_restarts_per_rmax // 2)))
                      if residual_seed is not None else set())

    for rmax_kpc in RMAX_OPTIONS:
        for restart_idx in range(n_restarts_per_rmax):
            use_seed = (restart_idx in seeded_indices and residual_seed is not None)
            seed_arg = residual_seed if use_seed else None
            p0 = draw_initial(n_shells, rmax_kpc, rng, hp, residual_seed=seed_arg)
            p0 = np.clip(p0, lower_arr, upper_arr)
            try:
                result = optimize.minimize(
                    loss_hier, p0,
                    args=(n_shells, data, mask, hp),
                    method="L-BFGS-B", bounds=bounds,
                    options={"maxiter": MAXFEV, "ftol": FTOL, "gtol": GTOL},
                )
                if result.success:
                    n_success += 1
                    if result.fun < best_loss:
                        best_loss = result.fun
                        best_params = result.x.copy()
            except Exception:
                continue

    if best_params is None:
        return None

    chi2_d = chi2_data(best_params, n_shells, data, mask)
    prior_val = prior_penalty_hier(best_params, hp)
    n_data = int(mask.sum())
    p = n_params(n_shells)
    dof = max(n_data - p, 1)
    chi2_red = chi2_d / dof
    bic = chi2_d + p * np.log(n_data)
    Y_disk, Y_bulge, rho_0, a_kpc, shells = unpack_params(best_params, n_shells)

    out = {
        "n_shells": n_shells,
        "n_params": p,
        "chi2": chi2_d,
        "chi2_red": chi2_red,
        "bic": bic,
        "prior_penalty": prior_val,
        "Y_disk": float(Y_disk),
        "Y_bulge": float(Y_bulge),
        "rho_0": float(rho_0),
        "a_kpc": float(a_kpc),
        "n_success": n_success,
        "used_de": used_de,
        "_best_params": best_params,
    }
    for i, (M_i, r_i, sigma_i) in enumerate(shells, start=1):
        out[f"M_shell_{i}"] = float(M_i)
        out[f"r_shell_{i}"] = float(r_i)
        out[f"sigma_shell_{i}"] = float(sigma_i)
        out[f"sigma_over_r_{i}"] = float(sigma_i / r_i) if r_i > 0 else np.nan
    return out


def fit_galaxy(
    galaxy: dict,
    hp: Hyperprior,
    rng: np.random.Generator,
    n_restarts_per_rmax: int,
    use_differential_evolution: bool,
) -> dict:
    """Fit one galaxy across N_shells ∈ {0, 1, 2}, BIC-select best."""
    data = load_rotmod(galaxy["rotmod_path"])

    # Mask: V_bar²(canonical Υ) < V_obs² — fixed mask matches Paper I
    V_bar_sq_canon = v_baryon_squared(
        PRIOR_UPS_DISK_MEAN, PRIOR_UPS_BULGE_MEAN,
        data["V_gas"], data["V_disk"], data["V_bulge"]
    )
    mask = V_bar_sq_canon < data["V_obs"] ** 2
    n_data = int(mask.sum())
    if n_data < MIN_USABLE_POINTS:
        return {
            "galaxy": galaxy["name"],
            "T": galaxy.get("T", -1),
            "n_data": n_data,
            "status": "insufficient_data",
        }

    bulgeless = bool(np.all(data["V_bulge"] == 0))

    per_n = {}
    residual_seed = None
    for n_shells in (0, 1, 2):
        res = fit_n_shells(
            n_shells, data, mask, hp, rng,
            n_restarts_per_rmax=n_restarts_per_rmax,
            use_differential_evolution=use_differential_evolution,
            residual_seed=residual_seed,
        )
        if res is None:
            per_n[n_shells] = None
            continue
        per_n[n_shells] = res
        if n_shells == 0 and "_best_params" in res:
            try:
                residual_seed = compute_residual_seed(res["_best_params"], data, mask)
            except Exception:
                residual_seed = None

    valid = {k: v for k, v in per_n.items() if v is not None}
    if not valid:
        return {
            "galaxy": galaxy["name"],
            "T": galaxy.get("T", -1),
            "n_data": n_data,
            "status": "no_successful_fit",
            "bulgeless": bulgeless,
        }

    best_n = min(valid.keys(), key=lambda k: valid[k]["bic"])
    best = valid[best_n]
    bic_n0 = valid.get(0, {}).get("bic", np.nan)
    bic_best_shellbearing = min(
        (valid[k]["bic"] for k in (1, 2) if k in valid), default=np.nan
    )

    out = {
        "galaxy": galaxy["name"],
        "T": galaxy.get("T", -1),
        "n_data": n_data,
        "bulgeless": bulgeless,
        "status": "ok",
        "best_N": best_n,
        "best_bic": best["bic"],
        "best_chi2": best["chi2"],
        "best_chi2_red": best["chi2_red"],
        "best_Y_disk": best["Y_disk"],
        "best_Y_bulge": best["Y_bulge"],
        "prior_penalty_at_best": best["prior_penalty"],
        "is_adequate": bool(best["chi2_red"] < CHI2_RED_ADEQUACY),
        "is_shell_bearing": bool(best_n >= 1),
        "delta_bic_framework_vs_burkert": (
            float(bic_best_shellbearing - bic_n0)
            if (np.isfinite(bic_n0) and np.isfinite(bic_best_shellbearing))
            else np.nan
        ),
    }
    # Flatten per-N for the CSV (keeps inspection easy)
    for n_shells in (0, 1, 2):
        r = valid.get(n_shells)
        if r is None:
            continue
        out[f"chi2_N{n_shells}"] = r["chi2"]
        out[f"chi2_red_N{n_shells}"] = r["chi2_red"]
        out[f"bic_N{n_shells}"] = r["bic"]
        out[f"Y_disk_N{n_shells}"] = r["Y_disk"]
        out[f"Y_bulge_N{n_shells}"] = r["Y_bulge"]
        for i in range(1, n_shells + 1):
            out[f"M_shell_{i}_N{n_shells}"] = r.get(f"M_shell_{i}", np.nan)
            out[f"r_shell_{i}_N{n_shells}"] = r.get(f"r_shell_{i}", np.nan)
            out[f"sigma_shell_{i}_N{n_shells}"] = r.get(f"sigma_shell_{i}", np.nan)
            out[f"sigma_over_r_{i}_N{n_shells}"] = r.get(f"sigma_over_r_{i}", np.nan)
    return out


# =============================================================================
# POPULATION-LEVEL AGGREGATION (Empirical Bayes update)
# =============================================================================

def aggregate_hyperprior(
    per_galaxy_results: List[dict],
    current_hp: Hyperprior,
) -> Tuple[Hyperprior, dict]:
    """Compute empirical Bayes update of hyperparameters from per-galaxy Υ.

    Disk: all galaxies with status='ok' contribute.
    Bulge: only galaxies with bulgeless==False contribute (bulgeless galaxies
    have no information about Υ_bulge and would dominate the average toward
    the prior mean if included).

    Returns the new hyperprior plus a diagnostics dict.
    """
    log_disks, log_bulges = [], []
    n_bulged = 0
    for r in per_galaxy_results:
        if r.get("status") != "ok":
            continue
        if not np.isfinite(r.get("best_Y_disk", np.nan)):
            continue
        log_disks.append(np.log10(r["best_Y_disk"]))
        if not r.get("bulgeless", True):
            log_bulges.append(np.log10(r["best_Y_bulge"]))
            n_bulged += 1

    if len(log_disks) < 3:
        return current_hp, {
            "n_disk_used": len(log_disks),
            "n_bulge_used": len(log_bulges),
            "note": "insufficient data, hyperprior unchanged",
        }

    mu_disk_new = float(np.clip(np.mean(log_disks), MU_DISK_MIN, MU_DISK_MAX))
    tau_disk_new = float(np.clip(np.std(log_disks, ddof=1) if len(log_disks) > 1 else current_hp.tau_disk,
                                  TAU_MIN, TAU_MAX))

    if len(log_bulges) >= 3:
        mu_bulge_new = float(np.clip(np.mean(log_bulges), MU_BULGE_MIN, MU_BULGE_MAX))
        tau_bulge_new = float(np.clip(np.std(log_bulges, ddof=1),
                                       TAU_MIN, TAU_MAX))
    else:
        # Not enough bulged galaxies → keep current
        mu_bulge_new = current_hp.mu_bulge
        tau_bulge_new = current_hp.tau_bulge

    new_hp = Hyperprior(mu_disk_new, tau_disk_new, mu_bulge_new, tau_bulge_new)
    diag = {
        "n_disk_used": len(log_disks),
        "n_bulge_used": len(log_bulges),
        "delta_mu_disk": mu_disk_new - current_hp.mu_disk,
        "delta_tau_disk": tau_disk_new - current_hp.tau_disk,
        "delta_mu_bulge": mu_bulge_new - current_hp.mu_bulge,
        "delta_tau_bulge": tau_bulge_new - current_hp.tau_bulge,
    }
    return new_hp, diag


def is_converged(diag: dict) -> bool:
    deltas = [
        abs(diag.get("delta_mu_disk", 0.0)),
        abs(diag.get("delta_tau_disk", 0.0)),
        abs(diag.get("delta_mu_bulge", 0.0)),
        abs(diag.get("delta_tau_bulge", 0.0)),
    ]
    return all(d < HYPER_CONVERGENCE_TOL_DEX for d in deltas)


# =============================================================================
# CHECKPOINT / RESUME
# =============================================================================

PER_GALAXY_CSV = "hierarchical_marginalization_per_galaxy.csv"
HYPER_HISTORY_CSV = "hierarchical_marginalization_hyperprior_history.csv"


def load_completed(csv_path: str) -> set:
    """Return set of (galaxy, iteration) tuples already present in the CSV."""
    if not os.path.exists(csv_path):
        return set()
    try:
        df = pd.read_csv(csv_path)
        return set(zip(df["galaxy"].astype(str), df["iteration"].astype(int)))
    except Exception:
        return set()


def append_result(csv_path: str, row: dict) -> None:
    """Append one row to the per-galaxy CSV (creates file if absent)."""
    df = pd.DataFrame([row])
    if os.path.exists(csv_path):
        # Use mode='a' append with no header for incremental writes; this means
        # the CSV's column set is fixed by the first write. We pad missing
        # columns with NaN to keep schema stable.
        try:
            existing = pd.read_csv(csv_path, nrows=1)
            all_cols = list(existing.columns)
            for c in all_cols:
                if c not in df.columns:
                    df[c] = np.nan
            for c in df.columns:
                if c not in all_cols:
                    # If new column appears, we need to rewrite the file with
                    # the union of columns. For simplicity, fall through to
                    # full rewrite.
                    full = pd.read_csv(csv_path)
                    full = pd.concat([full, df], ignore_index=True, sort=False)
                    full.to_csv(csv_path, index=False)
                    return
            df = df[all_cols]
            df.to_csv(csv_path, mode="a", header=False, index=False)
        except Exception:
            # Fallback: rewrite
            full = pd.read_csv(csv_path)
            full = pd.concat([full, df], ignore_index=True, sort=False)
            full.to_csv(csv_path, index=False)
    else:
        df.to_csv(csv_path, index=False)


def append_hyper(csv_path: str, row: dict) -> None:
    if os.path.exists(csv_path):
        existing = pd.read_csv(csv_path)
        full = pd.concat([existing, pd.DataFrame([row])], ignore_index=True, sort=False)
        full.to_csv(csv_path, index=False)
    else:
        pd.DataFrame([row]).to_csv(csv_path, index=False)


# =============================================================================
# MAIN
# =============================================================================

def run(
    data_dir: str,
    sample_csv_path: Optional[str],
    out_dir: str,
    exclude_ngc6674: bool,
    max_galaxies: Optional[int],
    max_iterations: int,
    n_restarts_per_rmax: int,
    use_differential_evolution: bool,
    resume: bool,
    random_seed: int,
) -> None:
    os.makedirs(out_dir, exist_ok=True)
    per_galaxy_csv = os.path.join(out_dir, PER_GALAXY_CSV)
    hyper_csv = os.path.join(out_dir, HYPER_HISTORY_CSV)

    rng = np.random.default_rng(random_seed)

    galaxies = load_sample(
        data_dir=data_dir,
        sample_csv_path=sample_csv_path,
        exclude_ngc6674=exclude_ngc6674,
        max_galaxies=max_galaxies,
    )
    print(f"Loaded sample: {len(galaxies)} galaxies after quality cuts "
          f"(NGC 6674 {'excluded' if exclude_ngc6674 else 'included'}).")

    # Initialize hyperprior
    hp = Hyperprior(INITIAL_MU_DISK, INITIAL_TAU_DISK_DEX,
                    INITIAL_MU_BULGE, INITIAL_TAU_BULGE_DEX)

    # If resuming, find the most recent iteration in the per-galaxy CSV and
    # the corresponding hyperprior. Otherwise start from iteration 0.
    completed = load_completed(per_galaxy_csv) if resume else set()
    if resume and os.path.exists(hyper_csv):
        try:
            hyper_df = pd.read_csv(hyper_csv)
            if len(hyper_df) > 0:
                last = hyper_df.iloc[-1]
                hp = Hyperprior(
                    float(last["mu_disk"]), float(last["tau_disk"]),
                    float(last["mu_bulge"]), float(last["tau_bulge"]),
                )
                print(f"Resuming from hyperprior at iteration "
                      f"{int(last['iteration'])}: μ_d={hp.mu_disk:.4f}, "
                      f"τ_d={hp.tau_disk:.4f}, μ_b={hp.mu_bulge:.4f}, "
                      f"τ_b={hp.tau_bulge:.4f}")
        except Exception as e:
            print(f"Could not parse hyper history: {e}; starting from initial prior.")

    print(f"Initial hyperprior:")
    print(f"  μ_disk  = {hp.mu_disk:.4f}  (Υ_disk  = {10**hp.mu_disk:.3f})")
    print(f"  τ_disk  = {hp.tau_disk:.4f} dex")
    print(f"  μ_bulge = {hp.mu_bulge:.4f}  (Υ_bulge = {10**hp.mu_bulge:.3f})")
    print(f"  τ_bulge = {hp.tau_bulge:.4f} dex")
    print()

    converged = False
    for iteration in range(max_iterations):
        print(f"{'='*70}")
        print(f"Iteration {iteration + 1}/{max_iterations}")
        print(f"  hyperprior: μ_disk = {hp.mu_disk:+.4f}, τ_disk = {hp.tau_disk:.4f}, "
              f"μ_bulge = {hp.mu_bulge:+.4f}, τ_bulge = {hp.tau_bulge:.4f}")
        print(f"{'='*70}")

        iteration_results: List[dict] = []
        t_iter_start = time.time()

        for idx, galaxy in enumerate(galaxies, start=1):
            if (galaxy["name"], iteration) in completed:
                # Read back the result so we can use it for aggregation
                try:
                    df_full = pd.read_csv(per_galaxy_csv)
                    row = df_full[
                        (df_full["galaxy"] == galaxy["name"])
                        & (df_full["iteration"] == iteration)
                    ].iloc[0].to_dict()
                    iteration_results.append(row)
                    print(f"  [{idx:>3}/{len(galaxies)}] {galaxy['name']:<14} "
                          f"(skipped, already in CSV)")
                except Exception:
                    pass
                continue

            t0 = time.time()
            try:
                result = fit_galaxy(
                    galaxy, hp, rng,
                    n_restarts_per_rmax=n_restarts_per_rmax,
                    use_differential_evolution=use_differential_evolution,
                )
            except Exception as e:
                result = {
                    "galaxy": galaxy["name"],
                    "T": galaxy.get("T", -1),
                    "n_data": 0,
                    "status": f"exception: {type(e).__name__}: {e}",
                }
            elapsed = time.time() - t0

            # Add bookkeeping
            result["iteration"] = iteration
            result["mu_disk_at_fit"] = hp.mu_disk
            result["tau_disk_at_fit"] = hp.tau_disk
            result["mu_bulge_at_fit"] = hp.mu_bulge
            result["tau_bulge_at_fit"] = hp.tau_bulge
            result["fit_time_sec"] = elapsed
            append_result(per_galaxy_csv, result)
            iteration_results.append(result)

            # Progress line
            if result.get("status") == "ok":
                print(f"  [{idx:>3}/{len(galaxies)}] {galaxy['name']:<14} "
                      f"T={result['T']:>2}  N*={result['best_N']}  "
                      f"χ²_r={result['best_chi2_red']:.2f}  "
                      f"BIC={result['best_bic']:.2f}  "
                      f"Υ_d={result['best_Y_disk']:.3f}  "
                      f"Υ_b={result['best_Y_bulge']:.3f}  "
                      f"ΔBIC(fw-burk)={result['delta_bic_framework_vs_burkert']:+.2f}  "
                      f"[{elapsed:.1f}s]")
            else:
                print(f"  [{idx:>3}/{len(galaxies)}] {galaxy['name']:<14} "
                      f"status={result.get('status')}  [{elapsed:.1f}s]")
            sys.stdout.flush()

        elapsed_iter = time.time() - t_iter_start
        print(f"\nIteration {iteration + 1} complete in {elapsed_iter/60:.1f} min.")

        # Aggregate
        new_hp, diag = aggregate_hyperprior(iteration_results, hp)
        diag_row = {
            "iteration": iteration,
            "mu_disk": new_hp.mu_disk,
            "tau_disk": new_hp.tau_disk,
            "mu_bulge": new_hp.mu_bulge,
            "tau_bulge": new_hp.tau_bulge,
            "n_disk_used": diag.get("n_disk_used", 0),
            "n_bulge_used": diag.get("n_bulge_used", 0),
            "delta_mu_disk": diag.get("delta_mu_disk", 0.0),
            "delta_tau_disk": diag.get("delta_tau_disk", 0.0),
            "delta_mu_bulge": diag.get("delta_mu_bulge", 0.0),
            "delta_tau_bulge": diag.get("delta_tau_bulge", 0.0),
            "converged": is_converged(diag),
        }
        append_hyper(hyper_csv, diag_row)
        print(f"\nPopulation update:")
        print(f"  μ_disk:  {hp.mu_disk:+.4f} → {new_hp.mu_disk:+.4f}   (Δ = {diag_row['delta_mu_disk']:+.4f}, n = {diag_row['n_disk_used']})")
        print(f"  τ_disk:  {hp.tau_disk:.4f}  → {new_hp.tau_disk:.4f}    (Δ = {diag_row['delta_tau_disk']:+.4f})")
        print(f"  μ_bulge: {hp.mu_bulge:+.4f} → {new_hp.mu_bulge:+.4f}   (Δ = {diag_row['delta_mu_bulge']:+.4f}, n = {diag_row['n_bulge_used']})")
        print(f"  τ_bulge: {hp.tau_bulge:.4f}  → {new_hp.tau_bulge:.4f}    (Δ = {diag_row['delta_tau_bulge']:+.4f})")

        if is_converged(diag):
            print(f"\nCONVERGED: all |Δ| < {HYPER_CONVERGENCE_TOL_DEX} dex.")
            converged = True
            hp = new_hp
            break
        hp = new_hp

    print()
    print("=" * 70)
    print("Final hyperprior:")
    print(f"  μ_disk  = {hp.mu_disk:+.4f}   (Υ_disk  population mean = {10**hp.mu_disk:.3f})")
    print(f"  τ_disk  = {hp.tau_disk:.4f}  dex")
    print(f"  μ_bulge = {hp.mu_bulge:+.4f}  (Υ_bulge population mean = {10**hp.mu_bulge:.3f})")
    print(f"  τ_bulge = {hp.tau_bulge:.4f}  dex")
    print(f"Converged: {converged}")
    print(f"\nPer-galaxy results: {per_galaxy_csv}")
    print(f"Hyperprior history: {hyper_csv}")
    print("\nNext step: run `summarize_hierarchical_results.py` to produce the")
    print("adequacy / morphology-gradient / per-galaxy ΔBIC summary tables.")


def main():
    p = argparse.ArgumentParser(
        description=("Empirical Bayes hierarchical Υ marginalization across "
                     "the full SPARC T=2-9 sample. Resumable; writes per-galaxy "
                     "results immediately as a checkpoint-safe CSV."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--data-dir", default="./Rotmod_LTG",
                   help="Directory containing SPARC *_rotmod.dat files")
    p.add_argument("--sample-csv", default="./data/sparc_sample123.csv",
                   help="Path to sparc_sample123.csv for T-type and quality flags")
    p.add_argument("--out-dir", default="./data",
                   help="Output directory for CSVs")
    p.add_argument("--exclude-ngc6674", action="store_true", default=True,
                   help="Exclude NGC 6674 (Paper II §2.3 default)")
    p.add_argument("--include-ngc6674", dest="exclude_ngc6674",
                   action="store_false",
                   help="Include NGC 6674 (matches Paper I 102-galaxy sample)")
    p.add_argument("--max-galaxies", type=int, default=None,
                   help="Limit to first N galaxies (for smoke testing)")
    p.add_argument("--max-iterations", type=int, default=MAX_HYPER_ITERATIONS,
                   help=f"Max empirical-Bayes iterations (default {MAX_HYPER_ITERATIONS})")
    p.add_argument("--n-restarts-per-rmax", type=int,
                   default=DEFAULT_N_RESTARTS_PER_RMAX,
                   help=f"L-BFGS-B restarts per r_max (default {DEFAULT_N_RESTARTS_PER_RMAX}; "
                        f"total restarts = 4× this)")
    p.add_argument("--fast", action="store_true",
                   help="Disable differential_evolution; faster but less reliable")
    p.add_argument("--resume", action="store_true",
                   help="Resume from existing CSV; skip galaxy-iteration pairs already present")
    p.add_argument("--random-seed", type=int, default=DEFAULT_RANDOM_SEED,
                   help=f"Random seed (default {DEFAULT_RANDOM_SEED})")
    args = p.parse_args()

    run(
        data_dir=args.data_dir,
        sample_csv_path=args.sample_csv,
        out_dir=args.out_dir,
        exclude_ngc6674=args.exclude_ngc6674,
        max_galaxies=args.max_galaxies,
        max_iterations=args.max_iterations,
        n_restarts_per_rmax=args.n_restarts_per_rmax,
        use_differential_evolution=(not args.fast),
        resume=args.resume,
        random_seed=args.random_seed,
    )


if __name__ == "__main__":
    main()
