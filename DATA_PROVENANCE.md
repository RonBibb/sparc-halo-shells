# Data Provenance — halo_shells v7.0

This document maps every machine-readable data file in `source/data/` to its
producer script in `source/code/`. Every CSV cited in the manuscript should
trace back to a single producer script that takes raw SPARC rotmod files
(plus, in some cases, an upstream CSV) as input and produces the file
deterministically.

## Inventory

| CSV file | Rows × Cols | Producer | Inputs | Manuscript section |
|---|---|---|---|---|
| `sparc_T2-T9_canonical_fits.csv` | 102 × 43 | `run_canonical_fits.py` | SPARC rotmod files, `sparc_sample123.csv` | abstract, §3.1, §3.2, §3.3, §3.6, Tables 1-4, 6 |
| `dc14_ngc5055_showcase.csv` | 1 × 17 | `dc14_ngc5055.py` | SPARC rotmod (NGC 5055 only) | §3.4, Table 6 |
| `dc14_universal_failures_results.csv` | 10 × 17 | `dc14_universal_failures.py` | SPARC rotmod files | §4.1, Table 7 |
| `einasto_robustness_results.csv` | 16 × 23 | `run_einasto_robustness.py` | SPARC rotmod files, `einasto_backbone.py` (library) | §3.6 |
| `einasto_full_sample_results.csv` | 102 × 24 | `run_einasto_full_sample.py` | SPARC rotmod files, `sparc_T2-T9_canonical_fits.csv`, `einasto_backbone.py` | §3.6, Appendix A |
| `null_test_results.csv` | 280 × 9 | `null_test.py` | `sparc_T2-T9_canonical_fits.csv` (for input fits) + SPARC rotmod radii | §3.5 |
| `null_test_T8_T9_extension.csv` | 66 × 9 | `null_test_extension.py` | same as above, T=8/9 only | §3.5 |
| `null_test_T2-T9_combined.csv` | 346 × 9 | combination of above two | merge of the two null_test CSVs | §3.5 |

## Producer-script-only files (no CSV output, but used by the pipeline)

| Script | Role |
|---|---|
| `einasto_backbone.py` | Library defining the Einasto profile and fitter; imported by `run_einasto_robustness.py` and `run_einasto_full_sample.py`. Not run directly. |
| `permutation_test.py` | Reads `sparc_T2-T9_canonical_fits.csv`; produces stdout statistics only (Spearman ρ, permutation p, bootstrap CI). Outputs to terminal, not file. Manuscript §3.3 quotes its values. |
| `verify_canonical_fit.py` | Single-galaxy spot-check verifier. Reads canonical CSV and SPARC rotmod for one galaxy; reports whether the verifier reproduces canonical numbers. Used for Category-C validation. |
| `figure4_300dpi.py` | Reads canonical CSV; produces `figures/figure4.pdf` and `figures/figure4.png`. |
| `figure5_300dpi.py` | Reads canonical CSV; produces `figures/figure5.pdf` and `figures/figure5.png`. |

## Upstream input that is *not* produced by any script in this package

`sparc_sample123.csv` — the SPARC catalog filtered to 123 galaxies, originally
produced by Paper A's pipeline (`01_parse_sparc_master.py` from the Paper A
v3 source). It is treated as a frozen input here. Documented in Paper A.

## Reproducibility chain

```
SPARC rotmod files (Lelli+2016 raw data)
       │
       ├── sparc_sample123.csv  (Paper A pipeline; treated as frozen input)
       │
       ├── run_canonical_fits.py
       │      → sparc_T2-T9_canonical_fits.csv
       │             ↓
       │             ├── permutation_test.py → stdout (§3.3)
       │             ├── figure4_300dpi.py → figures/figure4.{pdf,png}
       │             ├── figure5_300dpi.py → figures/figure5.{pdf,png}
       │             ├── null_test.py → null_test_results.csv (§3.5)
       │             ├── null_test_extension.py → null_test_T8_T9_extension.csv (§3.5)
       │             ├── run_einasto_full_sample.py → einasto_full_sample_results.csv (§3.6)
       │             └── verify_canonical_fit.py → spot-check verification (Category-C)
       │
       ├── dc14_ngc5055.py
       │      → dc14_ngc5055_showcase.csv (§3.4, Table 5)
       │
       ├── dc14_universal_failures.py
       │      → dc14_universal_failures_results.csv (§3.6, Table 7)
       │
       └── run_einasto_robustness.py
              → einasto_robustness_results.csv (§3.6)
```

## Conventions shared across all SPARC-fitting scripts

These conventions are reverse-engineered from the canonical CSV and verified
in the producer scripts:

- **Mass-to-light ratios**: Υ_disk = 0.5, Υ_bulge = 0.7 (SPARC default).
- **Baryonic V²**: signed-square sum: V² = V_gas|V_gas| + Υ_disk·V_disk|V_disk| + Υ_bulge·V_bulge|V_bulge|.
- **Exclusion rule**: drop points where V_obs² ≤ V_bar².
- **σ_V floor**: max(eV_obs, 1.0) km/s.
- **Sample filter for canonical analysis**: T ∈ [2, 9] (102 galaxies of the 123 in `sparc_sample123.csv`).
- **BIC formula**: chi² + k·ln(n_pts_used), with k = 2 / 2 / 5 / 8 for Burkert / NFW / FW 1-shell / FW 2-shell.
- **chi²_red formula**: chi² / max(n_pts_used − k, 1).
- **G**: 4.302 × 10⁻⁶ kpc · (km/s)² / M_☉.

## v6.5 → v7.0 changes

The v6.5 release shipped without a producer script for the canonical fits
CSV; the script existed only in the v6.3.6 development environment and was
not preserved in the release package. The v7.0 release adds the missing
producer (`run_canonical_fits.py`) and corrects a methodology-vs-code
mismatch: the σ/r ≤ 0.4 shell-width constraint stated in the manuscript was
not strictly enforced in the v6.5 producer (which used σ ≤ r_max × 0.4 with
r_max from a multi-restart configuration, not the fitted shell radius). The
v7.0 producer enforces σ/r ≤ 0.4 strictly via reparameterization.

The new constraint affected 6 of 102 BIC verdicts and tightened the
morphology gradient (Spearman ρ −0.81 → −0.83). The headline T=2 = 100%
shell-bearing claim is unchanged. See Appendix A of the manuscript for full
v6.5 → v7.0 comparison.
