# SPARC Halo-Shells Project — File Inventory

**Active paper:** *Localized Residual Structure in SPARC Rotation Curves: A Burkert-Backbone Model with BIC-Selected Shells* (halo_shells v7.1.0 frozen; **v7.1.1 candidate in progress for PASP resubmission** after AJ desk-reject)

**Author:** Ronald Bibb (Independent Researcher, Lilburn, GA)
**Correspondence:** ronbibb@gmail.com — ORCID 0009-0004-1153-2464

**Reproducibility chain:**
- GitHub: `github.com/RonBibb/sparc-halo-shells` at tag `v7.1.0` (frozen AJ submission)
- Zenodo concept DOI: `10.5281/zenodo.20072882` (always latest)
- Zenodo v7.1.0 DOI: `10.5281/zenodo.20100475` (frozen)
- Compiled manuscript: `sparc_shells.pdf` in this project (v7.1.0)
- v7.1.1 work-in-progress: PASP-targeted methods-paper reframe; Categories A + B body changes integrated; NGC 5055 marginalized-Υ analysis pending local execution (see new script below). Not yet tagged on GitHub; will become release `v7.1.1` once §3.4 update lands.

**Project knowledge scope:** This project contains only files relevant to the v7 halo_shells paper and its forward extensions. All Paper A/B/C/D framework / k_SMBH coupling / two-component decomposition work has been deliberately excluded. If a file appears that doesn't match the inventory below, it is stale and should be flagged.

---

## File inventory

### Foundational SPARC data

| File | Rows | Cols | Description |
|---|---|---|---|
| `*_rotmod.dat` (175 files) | varies | 8 | SPARC rotation curve data per Lelli, McGaugh & Schombert 2016 (AJ 152, 157). Header line gives distance in Mpc; columns are `Rad(kpc) Vobs(km/s) errV(km/s) Vgas Vdisk Vbul SBdisk(L/pc²) SBbul(L/pc²)`. Comment lines start with `#`. The 175 files cover the full SPARC release; the v7 paper uses 102 of them after Q≤2, Inc≥30°, V_flat>25 km/s, and T=2-9 cuts. |
| `sparc_sample123.csv` | 123 | 25 | Master sample table after quality cuts (Q≤2, Inc≥30°, V_flat>25 km/s). Includes structural parameters from the SPARC master (Galaxy, T, D, e_D, f_D, Inc, e_Inc, L36, e_L36, Reff, SBeff, Rdisk, SBdisk, MHI, RHI, Vflat, e_Vflat, Q, Ref) plus derived quantities (M_star, r_vir, M_halo, logM_star, logM_halo, logM_halo_per_star). The 123-galaxy count is the pre-T-cut total; 102 of these enter the main analysis. |
| `galaxy_classifications.csv` | 123 | 32 | Same 25 columns as `sparc_sample123.csv` plus 7 classification flags: `is_dwarf`, `is_mw_like`, `is_bulge_dom`, `is_bulgeless`, `is_transitional`, `max_bulge_frac`, `has_bulge_col`. Drives the morphology-gradient analysis and the bulge-dominated subset definitions. |

### Comparison fits (cited in paper as baselines)

| File | Rows | Cols | Description |
|---|---|---|---|
| `nfw_fixedc_fits.csv` | 123 | 12 | NFW fits with concentration fixed to the Dutton & Macciò (2014) M-c relation. Columns: Galaxy, c_DM14, rms, chi2, chi2_red, bic, n_points + 5 classification flags. The "fixed-c" variant in the paper's adequacy comparisons. |

### Reference manuscript

| File | Description |
|---|---|
| `sparc_shells.pdf` | Compiled v7.1.0 manuscript (~1.4 MB). Authoritative source for all headline numerical claims pending fit-table loadout. |

### Paper I v7.1.1 (PASP resubmission) — work in progress

Scripts and outputs supporting the v7.1.1 candidate manuscript. The v7.1.1 changes are a methods-paper reframe (Categories A + B integrated in `full_manuscript_review_v2.tex`) plus one substantive empirical addition: a joint Υ_disk / Υ_bulge marginalization test on the §3.4 NGC 5055 showcase, addressing the PASP-review-flagged gap that the §4.4 piecewise-linear Υ_T steelman pushes Υ upward (making shells more necessary, not less) and therefore does not adversarially stress-test the showcase galaxy.

| File | Type | Description |
|---|---|---|
| `scripts/ngc5055_marginalized_upsilon.py` | Python script | Fits NGC 5055 (default) under the framework architecture (Burkert backbone + N_shells ∈ {0, 1, 2}) with Υ_disk and Υ_bulge **jointly free** under Gaussian log-priors at the Li et al. 2020 convention: log10(Υ_disk) ~ N(log10(0.5), 0.1) and log10(Υ_bulge) ~ N(log10(0.7), 0.1). Mirrors the canonical pipeline (multi-restart with 8 restarts × 3 r_max values per n_shells; strict σ/r ≤ 0.4 via reparameterization σ_i = f_i · r_i, f_i ∈ [0.01, 0.4]). Reports BIC computed on the data χ² alone (prior contribution excluded from BIC for direct comparability to canonical Table 5). Configurable via `--galaxy` to rerun on NGC 2841 if the NGC 5055 showcase needs replacement. Random seed 20260513 for reproducibility. Default `--data-dir ./Rotmod_LTG` matches the Mac-side convention. **Run locally; not yet executed.** |
| `data/ngc5055_marginalized_upsilon_results.csv` | CSV (generated) | Output of the script above. One row per fitted configuration (`n_shells = 0, 1, 2`). Schema: `galaxy, n_shells, n_params, n_data, dof, chi2_data, chi2_red, prior_penalty, bic_marginalized, Y_disk_fit, Y_bulge_fit, rho_0_Msun_kpc3, a_kpc, n_attempts, n_success, [M_shell_i_Msun, r_shell_i_kpc, sigma_shell_i_kpc, sigma_over_r_i for each shell], delta_bic_vs_best_marg, canonical_bic, canonical_chi2_red, delta_bic_vs_canonical`. The `n_data` value should be 22 (matches canonical Table 5 mask for NGC 5055). The `bic_marginalized` column carries +2·ln(n) ≈ +6.18 expected BIC cost relative to canonical for the two extra Υ parameters. |

**Outcome interpretation gates for §3.4 update (built into the script's printed summary):**
- ΔBIC(Burkert-only_marg − Framework-1shell_marg) > 10 → framework still very strongly preferred; §3.4 holds, add single-sentence parenthetical in §3.4 Para 4 and shorten §4.4 disclosure.
- 2 < ΔBIC ≤ 10 → framework still preferred but weakened; reframe §3.4 as "competitive under marginalized Υ" rather than "decisively wins."
- |ΔBIC| ≤ 2 → inconclusive; consider NGC 2841 as replacement showcase or reframe §3.4 around the four-smooth-profile failure pattern alone (independent of shell preference).
- ΔBIC < −2 → Burkert-only preferred under marginalized Υ; replace NGC 5055 showcase (NGC 2841 is the canonical alternative from the universal-failure list).

---

## Known gaps — files NOT currently loaded

These files exist in the GitHub repo's `data/` directory but are not in project knowledge. If future work needs to reference v7 fit outputs directly (e.g., "show me the chi²_red distribution for the Burkert+shells fits"), these should be loaded:

| Missing file (canonical name in repo) | Why it might be wanted |
|---|---|
| `nfw_freec_fits.csv` | NFW free-concentration comparison; paper's 52/102 adequacy claim depends on this. Loading it back would close the comparison-fits set. |
| Burkert-only baseline fits | Paper's 65/102 adequacy claim |
| Burkert + BIC-selected shells main fits | Primary paper result (91/102 adequacy + per-galaxy shell parameters) |
| Einasto backbone fits | Paper's robustness check (90/102 classification agreement) |
| DC14 universal-failure subset fits | 10/10 fail demonstration |
| Y_T M/L steelman refit table | Paper's strongest robustness defense |
| 346-mock synthetic null test outputs | Validates 4.0% false-positive shell-selection rate |
| Per-T classification table | Per-bin shell-bearing fractions table |
| Shell catalog (per-galaxy shell parameters) | Detailed shell architecture per galaxy |

If a session needs any of these and they're not loaded, ask Ron to upload them rather than re-deriving — the repo at v7.1.0 is the canonical source.

---

## Notes for future-Claude on session start

1. **Read this README first.** Don't guess at file purposes from filenames; this document is more reliable than the file listing.
2. **The 102 vs 123 distinction matters.** `sparc_sample123.csv` and `galaxy_classifications.csv` have all 123 quality-cut galaxies. The paper's main analysis sample is 102 (after T=2-9 morphology cuts). When asked about "the sample" without qualifier, confirm which one Ron means.
3. **For numerical claims about the paper itself, check `sparc_shells.pdf`.** Headline numbers (91/102, ΔBIC distributions, p-values, etc.) are in the manuscript. Don't re-derive them; cite them.
4. **The `_rotmod.dat` format is fixed.** Whitespace-separated, `#` comments, 8 numeric columns after the header. The first comment line carries the distance (`# Distance = X.X Mpc`). For Mac-side scripts, default DATA_DIR is `./Rotmod_LTG`; in this project the files are at the project root.
5. **Excluded work.** Do not invoke or reference the framework / k_SMBH coupling / two-component decomposition / α-β / Paper A/B/C/D content. That work is dead-end and shouldn't appear in current analyses or recommendations.
6. **v7.1.1 candidate is the active editorial state.** Paper text is now methods-framed (PASP target) with Categories A + B body changes integrated. The NGC 5055 marginalized-Υ analysis (`scripts/ngc5055_marginalized_upsilon.py`) is the last remaining substantive item before submission. Once the CSV lands, §3.4 Para 4 gets a parenthetical with the marginalized BIC, and §4.4 gets either a one-sentence reference to that result (if shell preference holds) or a full disclosure paragraph (if it doesn't).

---

*Last updated: 2026-05-13. If the file list at the top of a session disagrees with this README, the file list wins for what's loaded — but flag the discrepancy so the README can be updated.*
