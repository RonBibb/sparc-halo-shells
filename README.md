# halo_shells v7.1.1

**Localized Residual Structure in SPARC Rotation Curves: A Burkert-Backbone Model with BIC-Selected Shells**

Author: Ron Bibb (<ronbibb@gmail.com>)
ORCID: 0009-0004-1153-2464
Version: v7.1.1
Release date: 2026-05-17 (PASP submission state)

---

## What is this?

This package contains the manuscript, data, producer scripts, and figures for the v7.1.1 revision of the SPARC rotation-curves halo-shells paper. It is a self-contained scientific release intended for journal submission and public deposit.

## Version lineage

- **v6.5** (May 2026): preceding state, AJ-style empirical paper.
- **v7.0** (May 7, 2026): strict `σ/r ≤ 0.4` constraint enforced via reparameterization; producer script for canonical fits committed (had been absent in v6.5). Six BIC verdicts changed (~5% of sample); morphology gradient tightened (Spearman ρ −0.81 → −0.83). Headline T=2 = 100% shell-bearing unchanged.
- **v7.1.0** (May 9, 2026): minor frozen snapshot of v7.0, archived at Zenodo `10.5281/zenodo.20100475`. AJ submission attempt; desk-rejected without referee feedback as scope-mismatched (methods-framed paper not the right fit for AJ's empirical-observational profile).
- **v7.1.1** (this release): editorial reframe for PASP, with per-galaxy joint Υ marginalization analysis added as a scope-bounding disclosure. No changes to canonical fits, null tests, Einasto results, DC14 results, or Υ(T) steelman — all v7.0 results carry over verbatim. New content:
  - Single-galaxy joint Υ marginalization fitter for NGC 5055 (`scripts/ngc5055_marginalized_upsilon.py`)
  - Empirical-Bayes hierarchical Υ marginalization across the full 101-galaxy sample (`scripts/hierarchical_upsilon_marginalization.py`, with companion summary script)
  - Three new data products documenting the marginalization disclosure (`data/ngc5055_marginalized_upsilon_results.csv`, `data/ngc5055_fixed_upsilon_recheck_results.csv`, `data/hierarchical_marginalization_per_galaxy.csv`, `data/hierarchical_marginalization_hyperprior_history.csv`)
  - Editorial compression of §4.3 and §4.4 (reduced defensive prose; methodological content unchanged)
  - Manuscript figure reordering — same files as v7.0, different `\includegraphics{}` positions in the LaTeX. File names retain their v7.0 numbering for cross-version reproducibility.

See `VALIDATION_STATUS.md` and `DATA_PROVENANCE.md` for full version-comparison details.

## Package contents

    sparc-halo-shells/
    ├── README.md                       (this file)
    ├── VALIDATION_STATUS.md            current state of v7.1.1 validation
    ├── VALIDATION_PLAN.md              checklist of consistency tests for submission
    ├── DATA_PROVENANCE.md              CSV-to-producer-script trace
    │
    ├── source/
    │   └── sparc_shells_body.tex       v7.1.1 manuscript (LaTeX, PASP target)
    │
    ├── data/
    │   ├── sparc_T2-T9_canonical_fits.csv          v7.0 canonical Burkert+shells fits, all 102 galaxies
    │   ├── null_test_results.csv                   v7.0 main null test (T=2-7, 280 mocks)
    │   ├── null_test_T8_T9_extension.csv           v6.5 extension preserved (T=8-9, 66 mocks)
    │   ├── null_test_T2-T9_combined.csv            combined null test (346 mocks total)
    │   ├── einasto_full_sample_results.csv         v7.0 Einasto-backbone full sample (102 galaxies)
    │   ├── einasto_robustness_results.csv          v7.0 Einasto 16-galaxy stratified subsample
    │   ├── sparc_T2-T9_y_T_fits.csv                v7.0 T-binned Y_disk re-fit (M/L systematics, §6.4)
    │   ├── dc14_universal_failures_results.csv     v6.5 preserved (constraint-independent)
    │   ├── dc14_ngc5055_showcase.csv               v6.5 preserved (constraint-independent)
    │   ├── ngc5055_marginalized_upsilon_results.csv          NEW v7.1.1: NGC 5055 joint-Y marginalization
    │   ├── ngc5055_fixed_upsilon_recheck_results.csv         NEW v7.1.1: --fix-upsilon sanity check
    │   ├── hierarchical_marginalization_per_galaxy.csv       NEW v7.1.1: empirical-Bayes hierarchical fit, all 101 galaxies × 3 iterations
    │   ├── hierarchical_marginalization_hyperprior_history.csv  NEW v7.1.1: hyperprior convergence trace
    │   ├── ngc2841_marginalized_upsilon_results.csv          NEW v7.1.1: NGC 2841 supplementary marginalization (prior-escape failure-mode documentation)
    │   ├── sparc_sample123.csv                     SPARC catalog after quality cuts (123 galaxies)
    │   ├── galaxy_classifications.csv              sparc_sample123 augmented with 7 morphology flags
    │   └── MD5SUMS.txt                             checksums for all data CSVs
    │
    ├── scripts/
    │   ├── run_canonical_fits.py                   v7.0 canonical fitter
    │   ├── run_canonical_fits_y_T.py               v7.0 T-binned Y_disk re-fit (§6.4)
    │   ├── null_test.py                            v7.0 main null test producer
    │   ├── null_test_extension.py                  v6.5 T=8/9 extension (preserved)
    │   ├── einasto_backbone.py                     v7.0 Einasto framework library
    │   ├── run_einasto_full_sample.py              v7.0 full-sample Einasto runner
    │   ├── run_einasto_robustness.py               v7.0 16-galaxy stratified Einasto runner
    │   ├── dc14_universal_failures.py              v6.5 preserved
    │   ├── dc14_ngc5055.py                         v6.5 preserved
    │   ├── permutation_test.py                     stats helper for §3.3
    │   ├── verify_canonical_fit.py                 verification utility
    │   ├── figure4_300dpi.py                       figure generator (regenerates figure4_shellfrac_vs_T.{pdf,png})
    │   ├── figure5_300dpi.py                       figure generator (regenerates figure5_dbic_histogram.{pdf,png})
    │   ├── figure_panels.py                        rotation-curve panel figures (T=4 grid, NGC 5055 showcase, universal failures, T2-T9 supplementary grids)
    │   ├── make_galaxy_classifications.py          regenerator for galaxy_classifications.csv
    │   ├── ngc5055_marginalized_upsilon.py         NEW v7.1.1: joint Y_disk/Y_bulge marginalization fitter for NGC 5055 (§3.4 reproduction). Produces ngc5055_marginalized_upsilon_results.csv. Supports --fix-upsilon sanity-check mode.
    │   ├── hierarchical_upsilon_marginalization.py NEW v7.1.1: empirical-Bayes hierarchical Y marginalization across the full 101-galaxy sample. Produces hierarchical_marginalization_per_galaxy.csv and hierarchical_marginalization_hyperprior_history.csv. Resumable via --resume.
    │   ├── summarize_hierarchical_results.py       NEW v7.1.1: companion summary script producing adequacy / morphology-gradient / per-galaxy ΔBIC tables from hierarchical output (partial or complete).
    │   ├── validate_v7_A.py                        A1-A8 numerical validation
    │   └── validate_v7_B.py                        B1-B7 cross-section validation
    │
    └── figures/
        ├── figure1_T4_example_grid.png         v7.0 figure naming (v7.1.x manuscript Figure 4); regenerated by figure_panels.py
        ├── figure2_ngc5055_showcase.png        v7.0 figure naming (v7.1.x manuscript Figure 3); regenerated by figure_panels.py
        ├── figure3_universal_failures.png      v7.0 figure naming (v7.1.x manuscript Figure 5); regenerated by figure_panels.py
        ├── figure4_shellfrac_vs_T.{pdf,png}    v7.0 figure naming (v7.1.x manuscript Figure 2); regenerated by figure4_300dpi.py
        ├── figure5_dbic_histogram.{pdf,png}    v7.0 figure naming (v7.1.x manuscript Figure 1); regenerated by figure5_300dpi.py
        └── T{2-9}_burkert_vs_framework.png     supplementary T-bin grids; regenerated by figure_panels.py

**Note on figure naming:** Image file names retain their v7.0 numbering (e.g., `figure5_dbic_histogram.png` is the file name for what appears as "Figure 1" in the v7.1.x PASP manuscript). This is intentional — preserving v7.0 file names maintains continuity with the frozen Zenodo v7.1.0 deposit. Paper figure numbers are assigned by `\begin{figure}` placement in the LaTeX source, not by file name. `DATA_PROVENANCE.md` documents the file-to-manuscript mapping for both v7.0 and v7.1.x.

## Reproducing the v7.1.1 results

You will need:

- Python 3.10+ with numpy, scipy, pandas, matplotlib
- The SPARC rotation-curve data files (publicly available from Lelli, McGaugh & Schombert 2016 supplementary materials)

Place SPARC rotmod files in `Rotmod_LTG/` alongside the scripts directory, then run any subset of:

    cd sparc-halo-shells/scripts/
    python3 run_canonical_fits.py                    # canonical fits, ~30 min on Apple Silicon Ultra
    python3 null_test.py                              # synthetic null test (T=2-7), ~10 min
    python3 null_test_extension.py                    # T=8/9 extension, ~5 min
    python3 run_einasto_full_sample.py                # Einasto-backbone full sample, ~30-60 min
    python3 run_einasto_robustness.py                 # Einasto 16-galaxy subsample, ~10 min
    python3 dc14_universal_failures.py                # DC14 universal failures
    python3 dc14_ngc5055.py                           # DC14 NGC 5055 showcase
    python3 run_canonical_fits_y_T.py                 # Y(T) steelman re-fit
    python3 figure4_300dpi.py                         # regenerate figure4_*
    python3 figure5_300dpi.py                         # regenerate figure5_*
    python3 figure_panels.py all                      # regenerate figure1_*, figure2_*, figure3_*, T-grids

For the v7.1.1 marginalization analyses (NEW):

    python3 ngc5055_marginalized_upsilon.py           # NGC 5055 joint-Y marginalization, ~5 min
    python3 ngc5055_marginalized_upsilon.py --fix-upsilon    # --fix-upsilon sanity check
    python3 hierarchical_upsilon_marginalization.py   # full-sample empirical-Bayes hierarchical, ~3-5 hours (or use --fast for ~1-2 hours)
    python3 summarize_hierarchical_results.py         # adequacy/morphology-gradient/ΔBIC summary tables

Total compute budget for full reproduction of v7.1.1: ~6-8 hours on a modern Apple Silicon Mac, single-process.

## Status

This package is **release v7.1.1**, validated, deposited, and ready for PASP submission:

| Component | Status |
| --- | --- |
| Manuscript text (v7.1.1, PASP-target reframe) | ✅ |
| Data CSVs | ✅ (with MD5SUMS verified) |
| Producer scripts | ✅ (paths package-relative) |
| Figures (v7.0 naming retained) | ✅ |
| Marginalization analyses (NGC 5055 single-galaxy + full-sample hierarchical) | ✅ |
| Validation pass A1-A8 (manuscript ↔ data) | ✅ from v7.0 (preserved unchanged) |
| Validation pass B1-B7 (manuscript ↔ manuscript) | ✅ from v7.0 (preserved unchanged) |
| Local seed-verification of NGC 5055 marginalization | ✅ 4 seeds, ΔBIC = +0.4815 ± 10⁻⁴ |
| Code repository | ✅ <https://github.com/RonBibb/sparc-halo-shells> |
| Zenodo concept DOI | ✅ [10.5281/zenodo.20072882](https://doi.org/10.5281/zenodo.20072882) |
| Zenodo v7.1.0 frozen DOI | ✅ [10.5281/zenodo.20100475](https://doi.org/10.5281/zenodo.20100475) (AJ-rejected snapshot, retained for citation continuity) |
| Zenodo v7.1.1 frozen DOI | ⏳ Pending (separate deposit at PASP submission) |
| PASP submission | ⏳ Pending (target: 2026-05-17/18) |

See `VALIDATION_STATUS.md` for the complete checklist.

## Citation

Until the PASP publication is final and assigned a DOI, please cite as:

> Bibb, R. (2026). *Localized Residual Structure in SPARC Rotation Curves:
> A Burkert-Backbone Model with BIC-Selected Shells.* halo_shells v7.1.1.
> Code and data: doi:10.5281/zenodo.20072882

## Companion analyses in preparation

The v7.1.1 release identifies two companion analyses in preparation as the natural follow-ups:

- **Paper II** (Bibb 2026a, in prep.): *Statistical Organization of Localized Residual Structure in SPARC Rotation Curves.* Takes the canonical-Υ shell catalog from Paper I and tests its internal statistical structure (morphology gradient, bulge correlation, radial scaling, anti-warp subsample stability, spatial-coherence destructive nulls, Einasto-backbone preservation, backbone-shift coupling).
- **Paper III** (Bibb 2026b, in prep.): *Hierarchical Bayesian Marginalization of Stellar Mass-to-Light Ratios for the SPARC Burkert-Shell Framework.* Reports the convergence and headline results of the hierarchical analysis whose data products are released here (`data/hierarchical_marginalization_*.csv`). Population-level adequacy advantage compresses from 26 to 14 galaxies but is preserved; morphology gradient survives at per-galaxy significance (Spearman ρ = −0.245, p = 0.013).

## Contact

Ron Bibb — <ronbibb@gmail.com> — ORCID 0009-0004-1153-2464

## License

Manuscript: under journal submission (license per PASP policy upon acceptance).
Data and code: MIT License (see `LICENSE`).
