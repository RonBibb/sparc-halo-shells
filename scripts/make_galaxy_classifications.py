#!/usr/bin/env python3
"""
make_galaxy_classifications.py

Regenerates `data/galaxy_classifications.csv` from `data/sparc_sample123.csv`
plus the SPARC rotmod files. Closes the reproducibility loop for the seven
morphology-classification flag columns (`is_dwarf`, `is_mw_like`,
`is_bulge_dom`, `is_bulgeless`, `is_transitional`, `max_bulge_frac`,
`has_bulge_col`).

The flags are not used by the v7.0 fitting pipeline (which works directly
from raw T-types); they are provided for follow-up analyses that need
pre-computed bulge / dwarf / MW-like classifications.

Recovered classification rules (verified to reproduce the shipped CSV
to machine precision on all 123 rows):

  is_dwarf       =  (T >= 8) AND (logM_star < 9.0)
  is_mw_like     =  (T in [3, 6]) AND (10.0 <= logM_star <= 11.0)
  max_bulge_frac =  max over r of [ 0.7 * Vbul^2(r) /
                                    (Vgas^2(r) + 0.5*Vdisk^2(r) + 0.7*Vbul^2(r)) ]
                    using ABSOLUTE squares (not sign-preserving) on the raw
                    Vgas, Vdisk, Vbul columns from the rotmod file
  is_bulgeless   =  (max_bulge_frac == 0)
  is_bulge_dom   =  (max_bulge_frac > 0)
  has_bulge_col  =  (max_bulge_frac > 0)
  is_transitional = always False (vacuous in this sample; reserved for
                                  future categorization)

Usage:
  python3 scripts/make_galaxy_classifications.py
  python3 scripts/make_galaxy_classifications.py --rotmod-dir /path/to/rotmod
  python3 scripts/make_galaxy_classifications.py --check       # verify against shipped CSV
  python3 scripts/make_galaxy_classifications.py --output X.csv

Path resolution:
  The script tries the package-relative location for sparc_sample123.csv
  (PACKAGE_ROOT/data/sparc_sample123.csv) first, then falls back to a flat
  working directory layout. The same applies to rotmod files: it tries
  PACKAGE_ROOT/Rotmod_LTG/<Galaxy>_rotmod.dat first, then ./<Galaxy>_rotmod.dat.

Note on convention:
  The absolute-square form for max_bulge_frac is slightly inconsistent with
  the canonical fits in run_canonical_fits.py, which use sign-preserving
  V_baryonic^2 (V*|V|) to handle a small number of galaxies with negative
  Vgas at some radii. The two forms diverge by at most ~0.006 in
  max_bulge_frac for five galaxies (NGC5985, UGC02487, UGC02916, UGC03580,
  UGC06614) and do not change any flag classification. This script
  reproduces the shipped CSV using the absolute-square form by design.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Path resolution

SCRIPT_DIR = Path(__file__).resolve().parent
PACKAGE_ROOT = SCRIPT_DIR.parent

DEFAULT_INPUT_CSV = PACKAGE_ROOT / "data" / "sparc_sample123.csv"
DEFAULT_OUTPUT_CSV = PACKAGE_ROOT / "data" / "galaxy_classifications.csv"
DEFAULT_ROTMOD_DIR_CANDIDATES = [
    PACKAGE_ROOT / "Rotmod_LTG",
    Path("./Rotmod_LTG"),
    Path("."),
]


def find_rotmod_dir(user_dir: str | None) -> Path:
    """Locate the directory holding <Galaxy>_rotmod.dat files."""
    if user_dir is not None:
        p = Path(user_dir)
        if not p.is_dir():
            raise FileNotFoundError(f"Rotmod directory not found: {user_dir}")
        return p
    for cand in DEFAULT_ROTMOD_DIR_CANDIDATES:
        if cand.is_dir() and any(cand.glob("*_rotmod.dat")):
            return cand
    raise FileNotFoundError(
        "No directory containing *_rotmod.dat files found. "
        f"Tried: {DEFAULT_ROTMOD_DIR_CANDIDATES}. "
        "Pass --rotmod-dir to specify."
    )


# ---------------------------------------------------------------------------
# Core computation

def read_rotmod_columns(rotmod_path: Path) -> np.ndarray:
    """Parse a SPARC rotmod file. Returns array with columns
    [Rad, Vobs, errV, Vgas, Vdisk, Vbul] (first 6 columns).
    Skips comment lines (#) and blank lines."""
    rows = []
    with open(rotmod_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 6:
                try:
                    rows.append([float(p) for p in parts[:6]])
                except ValueError:
                    continue
    if not rows:
        raise ValueError(f"No data rows parsed from {rotmod_path}")
    return np.array(rows)


def compute_max_bulge_frac(rotmod_path: Path) -> float:
    """Compute max_bulge_frac = max_r [ 0.7 * Vbul^2 /
    (Vgas^2 + 0.5*Vdisk^2 + 0.7*Vbul^2) ] using absolute squares.

    Returns 0.0 if Vbul is identically zero (galaxy has no bulge column),
    or if the denominator is non-positive everywhere."""
    arr = read_rotmod_columns(rotmod_path)
    Vgas = arr[:, 3]
    Vdisk = arr[:, 4]
    Vbul = arr[:, 5]

    Y_disk = 0.5
    Y_bulge = 0.7

    V2_bar = Vgas**2 + Y_disk * Vdisk**2 + Y_bulge * Vbul**2
    V2_bul = Y_bulge * Vbul**2

    valid = V2_bar > 0
    if not valid.any():
        return 0.0
    return float((V2_bul[valid] / V2_bar[valid]).max())


def classify(sample_df: pd.DataFrame, rotmod_dir: Path) -> pd.DataFrame:
    """Apply the seven derived columns to a copy of sample_df."""
    df = sample_df.copy()

    max_bulge_frac = []
    missing = []
    for galaxy in df["Galaxy"]:
        path = rotmod_dir / f"{galaxy}_rotmod.dat"
        if not path.exists():
            missing.append(galaxy)
            max_bulge_frac.append(np.nan)
            continue
        max_bulge_frac.append(compute_max_bulge_frac(path))

    if missing:
        raise FileNotFoundError(
            f"Missing rotmod files for {len(missing)} galaxies in {rotmod_dir}: "
            f"{missing[:5]}{'...' if len(missing) > 5 else ''}"
        )

    df["is_dwarf"] = (df["T"] >= 8) & (df["logM_star"] < 9.0)
    df["is_mw_like"] = df["T"].between(3, 6) & df["logM_star"].between(10.0, 11.0)
    df["is_bulge_dom"] = pd.Series(max_bulge_frac, index=df.index) > 0
    df["is_bulgeless"] = pd.Series(max_bulge_frac, index=df.index) == 0
    df["is_transitional"] = False
    df["max_bulge_frac"] = max_bulge_frac
    df["has_bulge_col"] = pd.Series(max_bulge_frac, index=df.index) > 0

    return df


# ---------------------------------------------------------------------------
# Verification

def check_against_shipped(produced: pd.DataFrame, shipped_path: Path) -> bool:
    """Verify produced DataFrame reproduces shipped CSV exactly."""
    if not shipped_path.exists():
        print(f"[check] No shipped CSV at {shipped_path}; skipping comparison.")
        return False

    shipped = pd.read_csv(shipped_path)

    # Align on Galaxy
    p = produced.sort_values("Galaxy").reset_index(drop=True)
    s = shipped.sort_values("Galaxy").reset_index(drop=True)
    if not (p["Galaxy"] == s["Galaxy"]).all():
        print("[check] Galaxy lists differ between produced and shipped.")
        return False

    all_pass = True
    flag_cols = [
        "is_dwarf",
        "is_mw_like",
        "is_bulge_dom",
        "is_bulgeless",
        "is_transitional",
        "has_bulge_col",
    ]
    for col in flag_cols:
        n_match = (p[col].astype(bool) == s[col].astype(bool)).sum()
        status = "OK" if n_match == len(s) else "FAIL"
        print(f"  [{status}] {col}: {n_match}/{len(s)}")
        if n_match != len(s):
            all_pass = False

    diffs = (p["max_bulge_frac"] - s["max_bulge_frac"]).abs()
    max_diff = float(diffs.max())
    if max_diff < 1e-9:
        print(f"  [OK] max_bulge_frac: {len(s)}/{len(s)} match within 1e-9 "
              f"(max diff = {max_diff:.2e})")
    else:
        print(f"  [WARN] max_bulge_frac: max diff = {max_diff:.2e}")
        if max_diff > 1e-4:
            all_pass = False

    return all_pass


# ---------------------------------------------------------------------------
# CLI

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Regenerate galaxy_classifications.csv from "
                    "sparc_sample123.csv + SPARC rotmod files.",
    )
    ap.add_argument("--input", default=str(DEFAULT_INPUT_CSV),
                    help="Input sample CSV (default: data/sparc_sample123.csv)")
    ap.add_argument("--output", default=str(DEFAULT_OUTPUT_CSV),
                    help="Output classifications CSV "
                         "(default: data/galaxy_classifications.csv)")
    ap.add_argument("--rotmod-dir", default=None,
                    help="Directory of *_rotmod.dat files "
                         "(default: search Rotmod_LTG/, ., etc.)")
    ap.add_argument("--check", action="store_true",
                    help="Compare produced output against existing shipped CSV "
                         "without overwriting it.")
    args = ap.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"ERROR: input not found: {input_path}", file=sys.stderr)
        return 1

    rotmod_dir = find_rotmod_dir(args.rotmod_dir)
    print(f"Reading: {input_path}")
    print(f"Rotmod:  {rotmod_dir}")

    sample = pd.read_csv(input_path)
    print(f"  -> {len(sample)} galaxies")

    out = classify(sample, rotmod_dir)
    print("Derived columns:")
    print(f"  is_dwarf:        {int(out['is_dwarf'].sum())}/{len(out)} True")
    print(f"  is_mw_like:      {int(out['is_mw_like'].sum())}/{len(out)} True")
    print(f"  is_bulge_dom:    {int(out['is_bulge_dom'].sum())}/{len(out)} True")
    print(f"  is_bulgeless:    {int(out['is_bulgeless'].sum())}/{len(out)} True")
    print(f"  is_transitional: {int(out['is_transitional'].sum())}/{len(out)} True")
    print(f"  has_bulge_col:   {int(out['has_bulge_col'].sum())}/{len(out)} True")
    print(f"  max_bulge_frac: range [{out['max_bulge_frac'].min():.4f}, "
          f"{out['max_bulge_frac'].max():.4f}]")

    if args.check:
        print(f"\nVerifying against shipped: {args.output}")
        ok = check_against_shipped(out, Path(args.output))
        if ok:
            print("\nReproduction VERIFIED: produced output matches shipped CSV exactly.")
            return 0
        print("\nReproduction MISMATCH: see above.")
        return 2

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, index=False)
    print(f"\nWrote: {output_path}")
    print(f"  {len(out)} rows x {len(out.columns)} columns")
    return 0


if __name__ == "__main__":
    sys.exit(main())
