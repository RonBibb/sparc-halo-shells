# Validation Status — halo_shells v7.0

**Package version:** v7.0
**Release date:** May 7, 2026
| Manuscript file:** `source/sparc_shells_body.tex` (527 lines, 81 kB)
**Derived formats:** `source/derived/sparc_shells.md` (markdown), `source/derived/sparc_shells.docx` (Word)

---

## What v7.0 changed relative to v6.5

The v6.5 manuscript stated `σ_i / r_i ≤ 0.4` as a shell-width constraint, but the v6.5 fit code implemented this as `σ_i ≤ r_max × 0.4`, where `r_max ∈ {3, 6, 12} kpc` is the multi-restart starting bound on the shell radius — not the fitted shell radius itself. This let `σ/r` exceed 0.4 by factors of several whenever the fitted `r_i` came in below `r_max`. In the v6.5 canonical fits, 20 of 102 galaxies had at least one σ/r > 0.4 in the n=1 column; among the 37 BIC-selected 1-shell galaxies, 11 violated the stated constraint.

The v7.0 producers reparameterize each shell's width as `σ_i = f_i · r_i` with `f_i ∈ [0.01, 0.4]` as a box-bounded optimizer parameter, so no fitted shell can exceed σ/r = 0.4 regardless of `r_max` or optimizer trajectory. The radial bound `r_max ∈ {3, 6, 12} kpc` is retained as a multi-restart starting bound only.

All v6.5 producer scripts whose σ/r enforcement was affected by this issue were re-run under v7.0. The DC14 fits (`scripts/dc14_*.py`) are constraint-independent and were not re-run; their v6.5 outputs are preserved unchanged. The v6.5 T=8/T=9 null test extension already enforced σ/r ≤ 0.4 via post-fit clamping and is preserved unchanged in `data/null_test_T8_T9_extension.csv`.

---

## Producer status

| Producer script                  | Producer version | Output CSV                                    | Status   |
| -------------------------------- | ---------------- | --------------------------------------------- | -------- |
| `run_canonical_fits.py`          | v7.0             | `sparc_T2-T9_canonical_fits.csv`              | ✅ DONE   |
| `null_test.py` (T=2-7, main)     | v7.0             | `null_test_results.csv`                       | ✅ DONE   |
| `null_test_extension.py` (T=8-9) | v6.5 (preserved) | `null_test_T8_T9_extension.csv`               | ✅ DONE   |
| (combine step)                   | trivial concat   | `null_test_T2-T9_combined.csv`                | ✅ DONE   |
| `run_einasto_full_sample.py`     | v7.0             | `einasto_full_sample_results.csv`             | ✅ DONE   |
| `run_einasto_robustness.py`      | v7.0             | `einasto_robustness_results.csv`              | ✅ DONE   |
| `dc14_universal_failures.py`     | v6.5 (preserved) | `dc14_universal_failures_results.csv`         | ✅ DONE   |
| `dc14_ngc5055.py`                | v6.5 (preserved) | `dc14_ngc5055_showcase.csv`                   | ✅ DONE   |

All 8 producers complete; all 7 output CSVs present in `data/`.

### CSV provenance (MD5 hashes — see `data/MD5SUMS.txt`)

```
sparc_T2-T9_canonical_fits.csv      b5eddfe46563ec194e952c2eb0cc3706    (v7.0, May 7 2026)
null_test_results.csv               4e3fc0fe99aeaccc5d0219e2e38e3ccb    (v7.0, May 7 2026)
null_test_T8_T9_extension.csv       d931208e8167e1315e32da739804f42e    (v6.5, preserved)
null_test_T2-T9_combined.csv        f3d09c474469e4e8258a1bb6c913aa89    (v7.0 + v6.5 concat)
einasto_full_sample_results.csv     1061231c64e9a3e4120edaf41c4a05c8    (v7.0, May 7 2026)
einasto_robustness_results.csv      d6841d1117503abd61a15fe01f5b4755    (v7.0, May 7 2026)
dc14_universal_failures_results.csv fdba51a7f90d3b31bfd743733115a524    (v6.5, preserved)
dc14_ngc5055_showcase.csv           6228c870cb2a9441634bc16c42146aa8    (v6.5, preserved)
```

The full hash file is in `data/MD5SUMS.txt`. To re-verify after any modification: `cd data && md5sum -c MD5SUMS.txt`.

---

## Manuscript update status

| Section / element                 | v7.0 update status |
| --------------------------------- | ------------------ |
| Abstract                          | ✅ DONE             |
| §1 Introduction                   | ✅ DONE             |
| §2 Methods (sample, framework, fitting) | ✅ DONE — §2.2 shell-width constraint description rewritten to match strict reparameterization |
| §3.1 Aggregate fit quality        | ✅ DONE             |
| §3.2 Win/loss analysis            | ✅ DONE — count of NFW-preferring galaxies revised 9→8 (strict KR-positive); detailed per-galaxy paragraph rewritten with v7.0 set |
| §3.3 T-type morphology gradient   | ✅ DONE — T=8 wobble subsubsection deleted (T=8 at 33% is on-trend in v7.0) |
| §3.4 NGC 5055 showcase            | ✅ Reproduces exactly under v7.0 — no edit needed |
| §3.5 Synthetic null test          | ✅ DONE — substantive rewrite of FP-trend argument: "opposite trend" replaced by "magnitude argument" |
| §3.6 Robustness to backbone       | ✅ DONE             |
| §4 Limitations (universal-failure) | ✅ Constraint-independent — no edit needed |
| §5.3 Base-profile shape coupling  | ✅ DONE — substantive rewrite mirroring §3.5 change |
| §6 Discussion                     | ✅ DONE             |
| §7 Conclusion                     | ✅ DONE             |
| Table 1 (aggregate fit quality)   | ✅ DONE             |
| Table 2 (ΔBIC win/loss)           | ✅ DONE             |
| Table 3 (sample composition)      | ✅ DONE — T=8 row 3/6 (50%) → 2/6 (33%); total 53 → 52 |
| Table 4 (morphology + null)       | ✅ DONE — all 8 rows + caption |
| Table 5 (null test summary)       | ✅ DONE — all cells + footnote rewritten |
| Table 6 (NGC 5055)                | ✅ Reproduces exactly — no edit needed |
| Table 7 (DC14 universal failures) | ✅ Constraint-independent — no edit needed |
| Data Availability section         | ✅ DONE             |
| Bibliography reference            | ✅ DONE — v6.5 → v7.0 |
| Appendix A: v7.0 changelog        | ✅ DONE — new appendix added documenting all changes |

---

## Outstanding work

### Required before submission

1. **Figures 4 and 5 — re-run from v7.0 CSVs** ✅ **DONE**
   - `figures/figure4_shellfrac_vs_T.pdf` and `.png` regenerated from v7.0 CSVs
   - `figures/figure5_dbic_histogram.pdf` and `.png` regenerated from v7.0 CSVs
   - Both scripts now compute everything from `data/` rather than hardcoding values

2. **Validation pass** ✅ **DONE**
   - A1-A8 numerical consistency: 8/8 pass, 97 individual claims verified against CSVs
   - B1-B6 cross-section consistency: 6/6 pass, 24 patterns verified across 30+ locations
   - Run via `python3 scripts/validate_v7_A.py` and `python3 scripts/validate_v7_B.py`

3. **CSV hash file** ✅ **DONE** — `data/MD5SUMS.txt` generated

### Optional

4. **PDF compilation check** — render the manuscript and visually verify pagination, table widths, figure references resolve correctly.

5. **Repository upload** — Zenodo deposit, GitHub mirror, ORCID linkage. Manuscript Data and Code Availability section currently says "Code repository URL pending."

---

## Headline result comparison: v6.5 → v7.0

| Result                                          | v6.5                | v7.0                |
| ----------------------------------------------- | ------------------- | ------------------- |
| Framework adequacy at χ²_red < 1.5              | 92/102              | 91/102              |
| Burkert adequacy                                | 65/102              | 65/102 (unchanged)  |
| NFW adequacy                                    | 52/102              | 52/102 (unchanged)  |
| Framework strongly preferred over NFW (ΔBIC<-10)| 53/102              | 53/102 (unchanged)  |
| Median ΔBIC among strongly-preferred subset    | -34                 | -35                 |
| NFW BIC-preferred (strict KR, ΔBIC > 2)         | 9/102               | 8/102               |
| Per-galaxy permutation p_two-sided              | 0.028               | 0.002               |
| Per-T-bin Spearman ρ                            | -0.81               | -0.83               |
| Bootstrap 95% CI on ρ                           | [-0.93, -0.12]      | [-0.95, -0.22]      |
| Mann-Kendall p                                  | 0.035               | 0.003               |
| Burkert-truth aggregate FP rate                 | 4.3% (5/117)        | 4.0% (7/173)        |
| NFW-truth aggregate FP rate                     | 57.3% (67/117)      | 63.6% (110/173)     |
| Null test mock count                            | 234                 | 346                 |
| Burkert/Einasto classification agreement        | 89/102 (87%)        | 90/102 (88%)        |
| Einasto per-galaxy Spearman ρ                   | -0.327              | -0.347              |
| Einasto per-galaxy p                            | 0.0008              | 0.0004              |
| T=8 shell-bearing fraction                      | 50% (3/6)           | 33% (2/6)           |
| NGC 5055 framework χ²_red                       | 1.464               | 1.464 (unchanged)   |
| NGC 5055 Einasto+2-shell χ²_red                 | 0.507               | 0.507 (unchanged)   |

**Net assessment:** the v7.0 numbers are uniformly stronger than the v6.5 numbers on every trend test, with no headline result reversed. The framework's empirical claims are tighter under strict σ/r enforcement.

### Argument changes

One substantive argument was retired: the v6.5 §3.5/§5.3 "Burkert-null FP trend opposite to real trend" argument relied on ρ_FP = +0.87, p = 0.005. Under v7.0 strict enforcement, ρ_FP = +0.33, p = 0.42 — no significant trend. Replaced by the "magnitude argument": real-data shell-bearing fractions exceed Burkert-truth FP rates by factors of 6-20 in every T-bin, ruling out smooth-Burkert truth as the dominant gradient driver at every point in the morphological range.

This is documented in detail in Appendix A of the manuscript.

---

## Verdict changes (per-galaxy)

### Canonical fits (v7.0 vs v6.5): 6 verdicts changed

| Galaxy   | T | v6.5 n_shells | v7.0 n_shells | Notes                              |
| -------- | - | ------------- | ------------- | ---------------------------------- |
| NGC0891  | 3 | 2             | 1             | Shell count refinement             |
| NGC3198  | 5 | 1             | 2             | Shell count refinement             |
| NGC5371  | 4 | 1             | 2             | Shell count refinement             |
| NGC6946  | 6 | 1             | 2             | Shell count refinement             |
| UGC05253 | 2 | 2             | 1             | Shell count refinement             |
| UGC07399 | 8 | 1             | 0             | **Shell-bearing → no-shell flip**  |

Net shell-bearing population: 53 → 52 (one galaxy lost; this is the source of the T=8 fraction change from 50% to 33%).

### Einasto full-sample (v7.0 vs v6.5): 5 verdicts changed

| Galaxy   | T | v6.5 n_shells | v7.0 n_shells |
| -------- | - | ------------- | ------------- |
| IC4202   | 4 | 2             | 1             |
| NGC0289  | 4 | 1             | 2             |
| NGC6195  | 3 | 0             | 1             |
| UGC02885 | 5 | 1             | 0             |
| UGC09133 | 2 | 1             | 2             |

Net Einasto/Burkert classification agreement: 89/102 → 90/102.

### Einasto 16-galaxy stratified (v7.0 vs v6.5): 2 verdicts changed

| Galaxy   | T | v6.5 n_shells | v7.0 n_shells |
| -------- | - | ------------- | ------------- |
| UGC09133 | 2 | 1             | 2             |
| NGC0289  | 3 | 1             | 2             |

Both within shell-bearing population (no shell-bearing/no-shell flips).

---

## Reproducibility

The full v7.0 result set can be regenerated from the package:

```bash
cd halo_shells_v7.0/scripts/
# Place SPARC rotmod files in halo_shells_v7.0/Rotmod_LTG/
python3 run_canonical_fits.py        # ~30 min on Apple Silicon Ultra
python3 null_test.py                  # ~10 min
python3 null_test_extension.py        # ~5 min (preserved v6.5 code, σ/r already enforced)
python3 run_einasto_full_sample.py    # ~30-60 min
python3 run_einasto_robustness.py     # ~10 min
python3 dc14_universal_failures.py    # constraint-independent; v6.5 output preserved
python3 dc14_ngc5055.py               # constraint-independent; v6.5 output preserved
python3 figure4_300dpi.py             # ⚠️ pending — needs running on v7.0 data
python3 figure5_300dpi.py             # ⚠️ pending — needs running on v7.0 data
```

The producer scripts are designed for a flat working directory: drop `data/sparc_T2-T9_canonical_fits.csv` into the same directory as the script and place SPARC rotmod files in `Rotmod_LTG/` alongside.

---

## Known minor issues — fixed in v7.0 release

1. **`einasto_backbone.py:237` — divide-by-zero in `ein2_chi2_red`** ✅ **FIXED in v7.0**
   - Pre-fix: `fit_e2['chi2']/(n-9)` divides by zero when `n ≤ 9`. Affected 17 galaxies.
   - Fix: `fit_e2['chi2']/max(n-9, 1)` — consistent with the canonical-fits convention.
   - Impact: cosmetic only (all 17 affected galaxies selected n_shells=0 or 1, the divide-by-zero never affected any BIC selection or downstream analysis).

2. **Hardcoded paths in `figure4_300dpi.py`, `figure5_300dpi.py`, `null_test_extension.py`, `einasto_backbone.py`** ✅ **FIXED in v7.0**
   - Pre-fix: scripts contained `/mnt/...` and `/home/claude/...` paths from earlier work.
   - Fix: all scripts now use package-relative paths (`HERE/PACKAGE_ROOT`-based) with auto-fallback to local `./Rotmod_LTG` for legacy v6.5 layouts.

3. **`figure4_300dpi.py` had stale v6.5 numbers and annotations** ✅ **FIXED in v7.0**
   - Pre-fix: hardcoded `n_shell = [9, 6, 9, 7, 9, 5, 3, 5]` (v6.5 values; T=8 is now 2 not 3), stale subtitle "Mann-Kendall p = 0.035", stale T=8 outlier annotation.
   - Fix: rewritten to compute all values from `data/sparc_T2-T9_canonical_fits.csv` and `data/null_test_T2-T9_combined.csv`. No hardcoded numbers; figure stays in sync with data.

## Known minor issues remaining

1. **`null_test.py` random_state stability** — the v7.0 main run sampled a different 28-galaxy stratified subset than the v6.5 main run because the input CSV path changed (the `random_state=42` is deterministic, but the input list is not). The v7.0 sample is independent of the v6.5 sample for T=2–7 but yields aggregate FP rates statistically consistent. The T=8/T=9 extension preserves the v6.5 11-galaxy subset.
