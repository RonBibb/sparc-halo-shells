"""
run_einasto_robustness.py — Einasto-backbone framework on a 16-galaxy
                            stratified subsample (T=2-9, 2/bin) (v7.0)

This is the v7.0 revision — uses the v7.0 einasto_backbone.py library
which enforces the strict sigma/r_shell <= 0.4 constraint via sigma_frac
reparameterization. The v6.5 version of this script imported `einasto_v2`
and used hardcoded paths to /home/claude and /mnt/project that won't work
on a fresh checkout.

Comparison galaxies and their v6.5 (Burkert-backbone) BIC-selected n_shells:
  T=2: UGC02916, UGC09133
  T=3: NGC0289,  NGC5985
  T=4: NGC5055,  NGC3198
  T=5: NGC3521,  NGC5033
  T=6: NGC3893,  NGC6946
  T=7: NGC2403,  NGC6503
  T=8: UGC01281, NGC3109
  T=9: UGC04305, DDO154

Usage:
  Place this script in a directory with:
    - einasto_backbone.py                    (v7.0 library)
    - sparc_T2-T9_canonical_fits.csv         (output of run_canonical_fits.py)
    - Rotmod_LTG/                             (SPARC rotmod files)
  
  Then:
    python3 run_einasto_robustness.py
  
  Output: einasto_robustness_results.csv with one row per galaxy.
  Expected runtime: 5-15 minutes.

Requires: numpy, scipy, pandas. einasto_backbone.py must be in the same
directory.
"""
import os
import sys
import numpy as np
import pandas as pd

# Import the v7.0 einasto fitter from this directory
HERE = os.path.dirname(os.path.abspath(__file__))
PACKAGE_ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from einasto_backbone import fit_galaxy

# ============================================================
# Configuration
# ============================================================
# Prefer package layout (PACKAGE_ROOT/Rotmod_LTG/), fall back to scripts/Rotmod_LTG/
_PRIMARY_DATA = os.path.join(PACKAGE_ROOT, 'Rotmod_LTG')
_LEGACY_DATA = os.path.join(HERE, 'Rotmod_LTG')
DATA_DIR = _PRIMARY_DATA if os.path.isdir(_PRIMARY_DATA) else _LEGACY_DATA
CANONICAL_CSV = os.path.join(PACKAGE_ROOT, 'data', 'sparc_T2-T9_canonical_fits.csv')
OUTPUT_CSV = os.path.join(PACKAGE_ROOT, 'data', 'einasto_robustness_results.csv')

# Galaxy list with T-types
GALAXIES = [
    ('UGC02916', 2), ('UGC09133', 2),
    ('NGC0289',  3), ('NGC5985',  3),
    ('NGC5055',  4), ('NGC3198',  4),
    ('NGC3521',  5), ('NGC5033',  5),
    ('NGC3893',  6), ('NGC6946',  6),
    ('NGC2403',  7), ('NGC6503',  7),
    ('UGC01281', 8), ('NGC3109',  8),
    ('UGC04305', 9), ('DDO154',   9),
]


def main():
    print("=" * 72)
    print("EINASTO ROBUSTNESS RUN (v7.0, strict sigma/r ≤ 0.4)")
    print(f"Data directory: {DATA_DIR}")
    print(f"Canonical CSV:  {CANONICAL_CSV}")
    print(f"Output:         {OUTPUT_CSV}")
    print("=" * 72)
    
    if not os.path.isdir(DATA_DIR):
        print(f"\nERROR: SPARC rotmod directory not found at '{DATA_DIR}'")
        sys.exit(1)
    
    if not os.path.exists(CANONICAL_CSV):
        print(f"\nERROR: Canonical CSV not found at '{CANONICAL_CSV}'")
        sys.exit(1)
    
    # Load the v7.0 canonical fits for comparison
    canon = pd.read_csv(CANONICAL_CSV).set_index('Galaxy')
    
    results = []
    for name, T in GALAXIES:
        rotmod = os.path.join(DATA_DIR, f'{name}_rotmod.dat')
        if not os.path.exists(rotmod):
            print(f"  MISSING: {rotmod}")
            continue
        
        res = fit_galaxy(name, rotmod, verbose=True)
        if res is None:
            print(f"  FAILED:  {name}")
            continue
        
        res['T'] = T
        # Bring in v7.0 Burkert-backbone framework results for comparison
        if name in canon.index:
            row = canon.loc[name]
            res['canon_burk_chi2_red'] = float(row['burk_chi2_red'])
            res['canon_fw_burkert_n_shells'] = int(row['fw_best_n_shells'])
            res['canon_fw_burkert_chi2_red'] = float(row['fw_best_chi2_red'])
        results.append(res)
    
    # Build summary DataFrame
    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_CSV, index=False)
    
    print(f"\n{'=' * 72}")
    print(f"DONE: saved {len(df)} rows to {OUTPUT_CSV}")
    print(f"{'=' * 72}")
    
    # Quick agreement summary
    if 'canon_fw_burkert_n_shells' in df.columns:
        df['burk_shell'] = df['canon_fw_burkert_n_shells'] >= 1
        df['eina_shell'] = df['fw_einasto_n_shells'] >= 1
        agree = (df['burk_shell'] == df['eina_shell']).sum()
        print(f"\nBurkert/Einasto classification agreement: {agree}/{len(df)} = {100*agree/len(df):.1f}%")
        print(f"\nPer-T-type breakdown:")
        for T in range(2, 10):
            sub = df[df['T'] == T]
            if len(sub) == 0:
                continue
            n_burk = sub['burk_shell'].sum()
            n_eina = sub['eina_shell'].sum()
            print(f"  T={T}: Burkert {n_burk}/{len(sub)} shell-bearing, Einasto {n_eina}/{len(sub)}")


if __name__ == "__main__":
    main()
