# Data Provenance — halo_shells v7.1.1

This document maps every machine-readable data file in `data/` to its
producer script in `scripts/`. Every CSV cited in the manuscript should
trace back to a single producer script that takes raw SPARC rotmod files
(plus, in some cases, an upstream CSV) as input and produces the file
deterministically.

## Inventory

### Canonical fits and robustness data (v7.0; carried forward unchanged in v7.1.x)

| CSV file | Rows × Cols | Producer | Inputs | Manuscript section |
| --- | --- | --- | --- | --- |
| `sparc_T2-T9_canonical_fits.csv` | 102 × 43 | `run_canonical_fits.py` | SPARC rotmod files, `sparc_sample123.csv` | abstract, §3.1–3.3, §3.6, Tables 1-4, 6 |
| `sparc_T2-T9_y_T_fits.csv` | 102 × 43 | `run_canonical_fits_y_T.py` | SPARC rotmod files, `sparc_sample123.csv` | §4.4 (M/L systematics robustness re-fit) |
| `dc14_ngc5055_showcase.csv` | 1 × 17 | `dc14_ngc5055.py` | SPARC rotmod (NGC 5055 only) | §3.4, Table 6 |
| `dc14_universal_failures_results.csv` | 10 × 17 | `dc14_universal_failures.py` | SPARC rotmod files | §4.1, Table 7 |
| `einasto_robustness_results.csv` | 16 × 23 | `run_einasto_robustness.py` | SPARC rotmod files, `einasto_backbone.py` | §3.6 |
| `einasto_full_sample_results.csv` | 102 × 24 | `run_einasto_full_sample.py` | SPARC rotmod files, `sparc_T2-T9_canonical_fits.csv`, `einasto_backbone.py` | §3.6, Appendix A |
| `null_test_results.csv` | 280 × 9 | `null_test.py` | `sparc_T2-T9_canonical_fits.csv` + SPARC rotmod radii | §3.5 |
| `null_test_T8_T9_extension.csv` | 66 × 9 | `null_test_extension.py` | same, T=8/9 only | §3.5 |
| `null_test_T2-T9_combined.csv` | 346 × 9 | combination of above two | merge of the two null_test CSVs | §3.5 |
| `sparc_sample123.csv` | 123 × 25 | upstream sample table (not produced by any script in this package) | — | sample-construction provenance |
| `galaxy_classifications.csv` | 123 × 32 | `make_galaxy_classifications.py` | `sparc_sample123.csv` + SPARC rotmod files | morphology-flag provenance |

### Marginalization analyses (NEW in v7.1.1)

| CSV file | Rows × Cols | Producer | Inputs | Manuscript section |
| --- | --- | --- | --- | --- |
| `ngc5055_marginalized_upsilon_results.csv` | 3 × ~30 | `ngc5055_marginalized_upsilon.py` | SPARC rotmod (NGC 5055 only) | §3.4 (per-galaxy Υ marginalization disclosure) |
| `ngc5055_fixed_upsilon_recheck_results.csv` | 3 × ~30 | `ngc5055_marginalized_upsilon.py --fix-upsilon` | SPARC rotmod (NGC 5055 only) | §3.4 (sanity check: --fix-upsilon mode pins Υ at canonical and matches Burkert-only canonical BIC to within ΔBIC = +0.5, verifying that the marginalization fitter reproduces the canonical pipeline at fixed Υ) |
| `hierarchical_marginalization_per_galaxy.csv` | ~303 × ~50 | `hierarchical_upsilon_marginalization.py` | SPARC rotmod files, `sparc_sample123.csv` | §4.4 (empirical-Bayes hierarchical extension; preliminary results referenced in Paper I §4.4, full treatment in Bibb 2026b, in prep.) |
| `hierarchical_marginalization_hyperprior_history.csv` | 3–5 × 12 | `hierarchical_upsilon_marginalization.py` | per-galaxy fits at each iteration | §4.4 (hyperprior convergence trace; converges in 3 iterations to μ_disk = −0.315, τ_disk = 0.081 dex, μ_bulge = −0.125, τ_bulge = 0.081 dex) |
| `ngc2841_marginalized_upsilon_results.csv` | 3 × ~30 | `ngc5055_marginalized_upsilon.py` (run with NGC 2841 as target) | SPARC rotmod (NGC 2841 only) | Supplementary: documents the prior-escape failure mode (NGC 2841 pushes Υ_disk to ~1.08, ~7σ above the hierarchical hyperprior). Not referenced in v7.1.1 manuscript but retained for reproducibility of the prior-escape pattern discussion in Paper II / Paper III. |

## Producer-script-only files (no CSV output, but used by the pipeline)

### From v7.0 (preserved unchanged)

| Script | Role |
| --- | --- |
| `einasto_backbone.py` | Library defining the Einasto profile and fitter; imported by `run_einasto_robustness.py` and `run_einasto_full_sample.py`. Not run directly. |
| `permutation_test.py` | Reads `sparc_T2-T9_canonical_fits.csv`; produces stdout statistics only (Spearman ρ, permutation p, bootstrap CI). Manuscript §3.3 quotes its values. |
| `verify_canonical_fit.py` | Single-galaxy spot-check verifier. Reads canonical CSV and SPARC rotmod for one galaxy; reports whether the verifier reproduces canonical numbers. |
| `figure4_300dpi.py` | Reads canonical CSV; produces `figures/figure4_shellfrac_vs_T.{pdf,png}` (v7.1.x manuscript Figure 2). |
| `figure5_300dpi.py` | Reads canonical CSV; produces `figures/figure5_dbic_histogram.{pdf,png}` (v7.1.x manuscript Figure 1). |
| `figure_panels.py` | Reads canonical CSV + SPARC rotmod files; produces the rotation-curve panel figures: `figure1_T4_example_grid.png` (v7.1.x manuscript Figure 4), `figure2_ngc5055_showcase.png` (v7.1.x manuscript Figure 3), `figure3_universal_failures.png` (v7.1.x manuscript Figure 5), and the eight supplementary T-bin grids `T2_burkert_vs_framework.png` through `T9_burkert_vs_framework.png`. Modes: `t-grid [T]`, `ngc5055`, `universal`, `all`. |
| `make_galaxy_classifications.py` | Regenerates `galaxy_classifications.csv` from `sparc_sample123.csv` + SPARC rotmod files. |
| `validate_v7_A.py` | A1–A8 numerical validation (manuscript ↔ data consistency). |
| `validate_v7_B.py` | B1–B7 cross-section validation (manuscript ↔ manuscript internal consistency). |

### NEW in v7.1.1

| Script | Role |
| --- | --- |
| `ngc5055_marginalized_upsilon.py` | Joint Υ_disk/Υ_bulge marginalization fitter for NGC 5055 (§3.4 reproduction). Implements `scipy.optimize.differential_evolution` for the global search followed by residual-seeded multi-restart L-BFGS-B, with Gaussian log-priors at the Li 2020 convention (`Y_disk ∼ N(log10(0.5), 0.1 dex)`, `Y_bulge ∼ N(log10(0.7), 0.1 dex)`). Supports `--fix-upsilon` sanity-check mode that pins Υ at canonical and disables priors (matches the canonical pipeline at fixed Υ). Seed-stable across {1, 42, 99999, 20260513}: ΔBIC(0-shell − 1-shell) = +0.4815 to four decimal places. Produces `ngc5055_marginalized_upsilon_results.csv` (or `ngc5055_fixed_upsilon_recheck_results.csv` in sanity-check mode). |
| `hierarchical_upsilon_marginalization.py` | Empirical-Bayes hierarchical Υ marginalization across the full 101-galaxy sample (NGC 6674 excluded by default per Paper II §2.3 convention). Iterates between per-galaxy MAP fits with the framework (Burkert + N ∈ {0, 1, 2} shells) under a population-level Gaussian log-prior on Υ, and population-level updates of the hyperprior from per-galaxy posteriors. Converges in 3 iterations with all |Δhyperparameter| < 0.01 dex. Resumable via `--resume`; per-galaxy CSV appended after each galaxy completes (safe for interruption). Produces `hierarchical_marginalization_per_galaxy.csv` and `hierarchical_marginalization_hyperprior_history.csv`. |
| `summarize_hierarchical_results.py` | Companion summary script. Reads the per-galaxy CSV (partial or complete) and produces headline summary tables: (1) adequacy comparison (marginalized vs canonical 91/102, 65/102), (2) morphology gradient by T-type under marginalization with Spearman tests, (3) per-galaxy ΔBIC distribution, (4) Υ population statistics by T-type and by bulge presence. Headline marginalized results: framework adequate 86/101, Burkert-only adequate 72/101, adequacy gap 14 (vs canonical 26); morphology gradient per-galaxy ρ = −0.245, p = 0.013 (vs canonical p = 0.003); 33.7% of galaxies decisively framework-preferred under marginalization, 53.5% decisively Burkert-only, 11.9% inconclusive. |

## Upstream input shipped with the package but not produced by any script in this package

`sparc_sample123.csv` — the SPARC catalog filtered to 123 candidate galaxies that passed the author's pre-fitting quality screen. The 102-galaxy T=2–9 analysis sample reported in this paper is a subset (T=2 through T=9 only; 21 galaxies at T=10 or with insufficient quality flags excluded). It is treated as a frozen input to the v7.x pipeline and is shipped here for sample-construction provenance.

`galaxy_classifications.csv` — `sparc_sample123.csv` augmented with seven morphology-classification flag columns (`is_dwarf`, `is_mw_like`, `is_bulge_dom`, `is_bulgeless`, `is_transitional`, `max_bulge_frac`, `has_bulge_col`). These flags are not used by the v7.0 fitting pipeline (which works directly from raw T-types) but are shipped for use by follow-up analyses that may need pre-computed bulge/dwarf classification without re-deriving it from the SPARC photometry. Reproducible from `sparc_sample123.csv` and the SPARC rotmod files via `scripts/make_galaxy_classifications.py`; the rules are:

- `is_dwarf` = (T ≥ 8) AND (logM_star < 9.0)
- `is_mw_like` = (T ∈ [3, 6]) AND (10.0 ≤ logM_star ≤ 11.0)
- `max_bulge_frac` = max over r of [0.7·Vbul²(r) / (Vgas²(r) + 0.5·Vdisk²(r) + 0.7·Vbul²(r))] using the absolute-square form (not sign-preserving) on the raw rotmod Vgas/Vdisk/Vbul columns
- `is_bulgeless` = (max_bulge_frac == 0)
- `is_bulge_dom` = `has_bulge_col` = (max_bulge_frac > 0)
- `is_transitional` = always False (vacuous in the v7.0 sample; reserved for a future categorization)

The script reproduces every row to machine precision (max diff 2.22e-16 on `max_bulge_frac`; 123/123 exact on all six flag columns).

`sparc_T2-T9_y_T_fits.csv` — output of `scripts/run_canonical_fits_y_T.py`, identical in column structure to the canonical fits CSV but produced under a T-binned `Upsilon_disk` schedule with anchors `Y(T=2)=0.65`, `Y(T=5)=0.50`, `Y(T=9)=0.40` per Schombert+2022 color calibrations (the strongest realistic gradient discussed in the literature). Input for §4.4's M/L-systematics robustness test (T=2 stays at 9/9 = 100% shell-bearing, full-sample per-galaxy classification stability is 96.1%, morphology-gradient finding survives under T-binned Υ at per-galaxy permutation `p_two_sided = 0.008` vs 0.002 canonical).

## Reproducibility chain

```
SPARC rotmod files (Lelli+2016 raw data)
       │
       ├── sparc_sample123.csv  (frozen upstream input)
       │      ↓
       │      └── make_galaxy_classifications.py
       │             → galaxy_classifications.csv
       │
       ├── run_canonical_fits.py
       │      → sparc_T2-T9_canonical_fits.csv
       │             ↓
       │             ├── permutation_test.py → stdout (§3.3)
       │             ├── figure4_300dpi.py → figures/figure4_shellfrac_vs_T.{pdf,png}
       │             ├── figure5_300dpi.py → figures/figure5_dbic_histogram.{pdf,png}
       │             ├── figure_panels.py → figures/figure[1-3]_*.png and T[2-9]_*.png
       │             ├── null_test.py → null_test_results.csv (§3.5)
       │             ├── null_test_extension.py → null_test_T8_T9_extension.csv (§3.5)
       │             ├── run_einasto_full_sample.py → einasto_full_sample_results.csv (§3.6)
       │             └── verify_canonical_fit.py → spot-check verification
       │
       ├── dc14_ngc5055.py
       │      → dc14_ngc5055_showcase.csv (§3.4, Table 6)
       │
       ├── dc14_universal_failures.py
       │      → dc14_universal_failures_results.csv (§4.1, Table 7)
       │
       ├── run_canonical_fits_y_T.py  (T-binned Upsilon_disk re-fit)
       │      → sparc_T2-T9_y_T_fits.csv (§4.4 M/L systematics robustness)
       │
       ├── run_einasto_robustness.py
       │      → einasto_robustness_results.csv (§3.6)
       │
       │  ── NEW in v7.1.1 ──
       │
       ├── ngc5055_marginalized_upsilon.py
       │      → ngc5055_marginalized_upsilon_results.csv (§3.4 marginalization disclosure)
       │      → ngc5055_fixed_upsilon_recheck_results.csv (--fix-upsilon sanity check)
       │      → ngc2841_marginalized_upsilon_results.csv (supplementary; prior-escape failure mode)
       │
       └── hierarchical_upsilon_marginalization.py
              → hierarchical_marginalization_per_galaxy.csv (§4.4, Paper III in prep.)
              → hierarchical_marginalization_hyperprior_history.csv (hyperprior convergence trace)
                  ↓
                  └── summarize_hierarchical_results.py → headline summary tables (stdout + CSV)
```

## Conventions shared across all SPARC-fitting scripts

These conventions are reverse-engineered from the canonical CSV and verified in the producer scripts:

- **Mass-to-light ratios** (canonical): Υ_disk = 0.5, Υ_bulge = 0.7 (SPARC default).
- **Marginalization priors** (v7.1.1, §3.4 only): log10(Υ_disk) ∼ N(log10(0.5), 0.1 dex), log10(Υ_bulge) ∼ N(log10(0.7), 0.1 dex) per Li 2020 convention.
- **Baryonic V²**: signed-square sum: V² = V_gas|V_gas| + Υ_disk·V_disk|V_disk| + Υ_bulge·V_bulge|V_bulge|.
- **Exclusion rule**: drop points where V_obs² ≤ V_bar².
- **σ_V floor**: max(eV_obs, 1.0) km/s.
- **Sample filter for canonical analysis**: T ∈ [2, 9] (102 galaxies of the 123 in `sparc_sample123.csv`).
- **Sample filter for hierarchical marginalization**: T ∈ [2, 9], NGC 6674 excluded (101 galaxies; matches Paper II §2.3 convention due to NGC 6674's degenerate two-shell fit).
- **BIC formula** (canonical-Υ): chi² + k·ln(n_pts_used), with k = 2 / 2 / 5 / 8 for Burkert / NFW / FW 1-shell / FW 2-shell.
- **BIC formula** (marginalized Υ): chi² + k·ln(n_pts_used), with k = 4 / 4 / 7 / 10 (two additional free parameters Υ_disk, Υ_bulge added to each model).
- **chi²_red formula**: chi² / max(n_pts_used − k, 1).
- **G**: 4.302 × 10⁻⁶ kpc · (km/s)² / M_☉.

## Manuscript figure ordering: v7.0 vs v7.1.x

The v7.1.x PASP manuscript reorders the figure sequence relative to v7.0. **File names retain their v7.0 numbering** for cross-version reproducibility (anyone citing the v7.1.0 Zenodo deposit by figure file name will find the same files in v7.1.x). Paper figure numbers are assigned by `\begin{figure}` placement in the LaTeX source.

| Image file (v7.0 naming, preserved) | v7.0 paper Fig # | v7.1.x paper Fig # | v7.1.x section |
| --- | --- | --- | --- |
| `figure5_dbic_histogram.{pdf,png}` | 5 | 1 | §3.2 win/loss analysis |
| `figure4_shellfrac_vs_T.{pdf,png}` | 4 | 2 | §3.3 morphology gradient |
| `figure2_ngc5055_showcase.png` | 2 | 3 | §3.4 NGC 5055 showcase |
| `figure1_T4_example_grid.png` | 1 | 4 | §3.5 rotation-curve example grid |
| `figure3_universal_failures.png` | 3 | 5 | §4.1 universal failures |
| `T{2-9}_burkert_vs_framework.png` | (supplementary) | (supplementary) | §3.5 supplementary set |

The reordering reflects the v7.1.x editorial reframe: leading with the population-level adequacy comparison (`figure5_dbic_histogram`) and morphology gradient (`figure4_shellfrac_vs_T`) as load-bearing findings, with the single-galaxy NGC 5055 showcase demoted to an architectural demonstration rather than a headline result.

## v7.0 → v7.1.x changes

The v7.1.x release is an **editorial and scope-bounding revision** of v7.0. **No canonical fits, null tests, Einasto results, DC14 results, or Υ(T) steelman were re-run.** All v7.0 numerical results carry forward unchanged.

What v7.1.x adds:

1. **Reframe for PASP target.** Original v7.0/v7.1.0 manuscript was framed as an empirical-observational paper for AJ; received same-day desk-reject. v7.1.x reframes as a methods paper with disclosed scope limitations, suitable for PASP's editorial profile. The reframe demotes per-galaxy interpretation, elevates the population-level 91/102-vs-65/102 adequacy comparison as the load-bearing claim, and reorders figures accordingly.

2. **Per-galaxy joint Υ marginalization disclosure (§3.4).** New analysis on NGC 5055 with Υ_disk and Υ_bulge jointly free under Gaussian log-priors at the Li 2020 convention. Result: ΔBIC drops from +157 (canonical Υ) to +0.5 (inconclusive) under marginalization. Documented as a scope-bounding disclosure rather than a negative result: per-galaxy preferences are Υ-sensitive; population-level adequacy at canonical Υ is not. Reproduction script: `ngc5055_marginalized_upsilon.py`.

3. **§4.4 scope-bounding paragraphs.** Three blue paragraphs added clarifying that (a) the Υ(T) steelman tests population-level robustness, not per-galaxy adversarial direction; (b) per-galaxy joint Υ marginalization for NGC 5055 has been performed and shows the per-galaxy Υ-degeneracy; (c) the population-level 91/102-vs-65/102 claim is computed at canonical Υ across the full sample and is the load-bearing methodological result of the paper.

4. **Empirical-Bayes hierarchical Υ marginalization across the full sample (§4.4 forward reference; full treatment in Paper III in prep.).** New analysis on 101 galaxies (NGC 6674 excluded) with population-level hyperprior on Υ that is iteratively updated from per-galaxy fits. Converges in 3 iterations with hyperprior near canonical SPARC values (μ_disk = −0.315, τ_disk = 0.081 dex, μ_bulge = −0.125, τ_bulge = 0.081 dex). Adequacy gap compresses from 26 to 14 galaxies but is preserved; morphology gradient survives at per-galaxy significance (ρ = −0.245, p = 0.013). Reproduction scripts: `hierarchical_upsilon_marginalization.py` + `summarize_hierarchical_results.py`.

5. **Einasto adequacy paragraph in §3.6.** New "Fourth" finding paragraph quantifying that the framework's adequacy advantage over its smooth-halo counterpart persists under the Einasto backbone, not only the morphology gradient. Numbers derived from `einasto_full_sample_results.csv` (unchanged from v7.0).

6. **Forward references in §6 to Paper II and Paper III in preparation.** Paper II addresses consolidated artifact-channel robustness for the morphology gradient (spatial-coherence destructive nulls, Υ/D/i perturbation analyses, anti-warp clean subsample, backbone-shift coupling). Paper III addresses hierarchical marginalization.

7. **Editorial compression of §4.3 and §4.4.** Reduced defensive prose without changing methodological content; net reduction ~4.4% in word count from v3.1 of the manuscript.

8. **Manuscript figure reordering** (see table above). Files unchanged; only LaTeX `\includegraphics{}` placement reordered.

The v6.5 → v7.0 strict σ/r ≤ 0.4 constraint enforcement (described in the v7.0 release notes) carries forward unchanged into v7.1.x.
