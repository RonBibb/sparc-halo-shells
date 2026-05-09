# Validation Plan — halo_shells v7.0

This document lists the consistency checks that should be run before the v7.0 manuscript is submitted. Each check is binary (pass/fail) and points at a specific number or claim in the manuscript that must match a specific value in a specific data file.

**Status legend:** ⏳ not yet run, ✅ passed, ❌ failed

---

## Part A — Numerical consistency (manuscript ↔ data files)

Each Ax check verifies that a number quoted in the manuscript matches the value computed from a v7.0 CSV.

### A1. Framework adequacy at χ²_red < 1.5
- **Manuscript claim:** 91/102 (89.2%)
- **CSV:** `data/sparc_T2-T9_canonical_fits.csv`
- **Check:** `(fw_best_chi2_red < 1.5).sum()` should equal 91
- **Status:** ⏳

### A2. Burkert and NFW adequacy at χ²_red < 1.5
- **Manuscript claim:** Burkert 65/102, NFW 52/102
- **CSV:** `data/sparc_T2-T9_canonical_fits.csv`
- **Check:** `(burk_chi2_red < 1.5).sum() == 65` and `(nfw_chi2_red < 1.5).sum() == 52`
- **Status:** ⏳

### A3. Strict adequacy at χ²_red < 1.0
- **Manuscript claim:** Framework 86, Burkert 57, NFW 41
- **CSV:** `data/sparc_T2-T9_canonical_fits.csv`
- **Check:** corresponding `(* < 1.0).sum()`
- **Status:** ⏳

### A4. Per-T shell-bearing fractions
- **Manuscript claim (Table 3, abstract):** T=2: 100%, T=3: 54.5%, T=4: 52.9%, T=5: 53.8%, T=6: 56.2%, T=7: 35.7%, T=8: 33.3%, T=9: 31.2%
- **CSV:** `data/sparc_T2-T9_canonical_fits.csv`
- **Check:** group by `T`, compute `(fw_best_n_shells >= 1).mean()` per group
- **Status:** ⏳

### A5. Trend tests (per-galaxy and per-T-bin)
- **Manuscript claim:** per-T-bin Spearman ρ = -0.83, p = 0.010; per-galaxy permutation p_two-sided = 0.002, p_one-sided = 0.001; bootstrap 95% CI [-0.95, -0.22]; Mann-Kendall p = 0.003
- **CSV:** `data/sparc_T2-T9_canonical_fits.csv`
- **Check:** recompute Spearman, permutation (10000 resamples, seed=42), bootstrap (10000 resamples, seed=42), Kendall τ
- **Status:** ⏳

### A6. Null test aggregate rates
- **Manuscript claim:** Burkert-truth 4.0% (7/173), NFW-truth 63.6% (110/173), combined 33.8% (117/346)
- **CSV:** `data/null_test_T2-T9_combined.csv`
- **Check:** group by `smooth_truth`, compute `(best_n_shells > 0).sum() / count`
- **Status:** ⏳

### A7. Per-T-bin null FP rates
- **Manuscript claim (Table 4):** Burkert and NFW per-T values (8 rows × 2 columns)
- **CSV:** `data/null_test_T2-T9_combined.csv`
- **Check:** group by `(smooth_truth, T)`, compute `(best_n_shells > 0).mean()`
- **Status:** ⏳

### A8. NGC 5055 showcase numbers
- **Manuscript claim:** Burkert χ²_red = 9.53, NFW χ²_red = 30.21, Framework χ²_red = 1.464 (1 shell at r=12.0 kpc, σ=2.91 kpc, M=3.78×10¹⁰ M☉); Einasto-only χ²_red = 8.42, Einasto+2-shell χ²_red = 0.507
- **CSVs:** `data/sparc_T2-T9_canonical_fits.csv` (NGC5055 row), `data/einasto_full_sample_results.csv` (NGC5055 row), `data/dc14_ngc5055_showcase.csv`
- **Check:** spot-check each value in the corresponding row
- **Status:** ⏳

---

## Part B — Cross-section consistency (manuscript ↔ manuscript)

Each Bx check verifies that a number quoted in two or more places in the manuscript is consistent.

### B1. Framework adequacy count
- **Locations:** Abstract, §1, §3.1, §6.1, §7
- **Common value:** 91/102
- **Check:** all 5 locations say "91" not "92"
- **Status:** ⏳

### B2. NFW BIC-preferring count
- **Locations:** Abstract, §3.2 (twice — count and detailed analysis), §6.1, §7
- **Common value:** 8 (under strict KR, dBIC > 2)
- **Check:** all locations say "8 BIC-prefer NFW", not "9" or "14"
- **Status:** ⏳

### B3. Trend statistics
- **Locations:** Abstract, §1, §3.3, §6.1, Table 4 caption
- **Common values:** ρ = -0.83, CI [-0.95, -0.22], p_two-sided = 0.002, MK p = 0.003
- **Check:** all locations agree
- **Status:** ⏳

### B4. Null test aggregate rates and mock count
- **Locations:** Abstract, §3.5, §5.3, §6.6, Table 5
- **Common values:** 4.0% Burkert FP, 63.6% NFW FP, 346 total mocks
- **Check:** all locations agree
- **Status:** ⏳

### B5. Einasto agreement count
- **Locations:** Abstract, §3.6, §5.3, §6.1, §6.3
- **Common value:** 90/102 (88%)
- **Check:** all locations agree; no leftover "89/102 (87%)"
- **Status:** ⏳

### B6. Einasto morphology gradient stats
- **Locations:** Abstract, §3.6, §7
- **Common values:** Einasto per-galaxy ρ = -0.347, p = 0.0004; Burkert per-galaxy ρ = -0.296, p = 0.003
- **Check:** all locations agree
- **Status:** ⏳

### B7. T-binned Upsilon_disk refit (M/L systematics)
- **Locations:** Abstract, §6.4
- **Common values:** classification stability 96.1%; per-galaxy permutation `p_two = 0.008`; per-galaxy ρ = -0.264; T=2 vs T=3-6 Fisher p = 0.004; bootstrap CI [-0.905, -0.096]; adequacy 90/66/51; Schombert 2022 anchor Y(T=2) = 0.65
- **Check:** all values consistent with `data/sparc_T2-T9_y_T_fits.csv`; the §6.4 refutation paragraph contains all headline numbers
- **Status:** ⏳

---

## Part C — Producer reproducibility (optional but recommended)

### C1. Canonical fits reproducibility
- **Action:** Re-run `scripts/run_canonical_fits.py` against `Rotmod_LTG/` and verify the output CSV byte-matches `data/sparc_T2-T9_canonical_fits.csv` (MD5: `b5eddfe46563ec194e952c2eb0cc3706`)
- **Status:** ⏳

### C2. T-binned Y refit reproducibility
- **Action:** Re-run `scripts/run_canonical_fits_y_T.py` against `Rotmod_LTG/` and verify the output CSV byte-matches `data/sparc_T2-T9_y_T_fits.csv` (MD5: `1dbfabbc134ae8023f99747f578a2d53`)
- **Caveat:** Multi-restart fitting may produce small floating-point differences across machines; sanity-check via the B7 headline numbers if exact byte-match fails
- **Status:** ⏳

### C3. Null test reproducibility
- **Action:** Re-run `scripts/null_test.py` and verify the output CSV byte-matches `data/null_test_results.csv`
- **Caveat:** depends on numpy random_state stability across machines and numpy versions
- **Status:** ⏳

### C4. Einasto reproducibility
- **Action:** Re-run `scripts/run_einasto_full_sample.py` and verify byte-match
- **Status:** ⏳

---

## Running the validation

The validation is implemented as two scripts in `scripts/`:

- `validate_v7_A.py` — automates checks A1-A8 (manuscript prose claims vs CSV).
  Loads all v7.0 CSVs, parses numeric claims from `source/sparc_shells_body.tex`,
  and verifies each claim against the corresponding CSV computation.
- `validate_v7_B.py` — automates checks B1-B7 (cross-section consistency:
  whether the same numbers appear consistently in abstract, body sections,
  and table captions; whether stale v6.5 values were cleanly removed).

Both should pass clean before any submission. Run from the package root:

```
python3 scripts/validate_v7_A.py    # Should print ALL 8/8 CHECKS PASSED
python3 scripts/validate_v7_B.py    # Should print ALL 7/7 CHECKS PASSED
```

These were the last scripts run before v7.0 release lock-in.

---

## Submission checklist (post-validation)

After all A1-A8 and B1-B7 checks pass:

- [ ] Compile manuscript to PDF; visual review of pagination, figure references, table widths
- [ ] Generate `data/MD5SUMS.txt` for reproducibility
- [ ] Update `Code repository URL` in §Data and Code Availability
- [ ] Zenodo deposit (DOI for citation)
- [ ] AJ submission portal upload
