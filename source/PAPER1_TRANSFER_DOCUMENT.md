# Paper 1 Transfer Document — halo_shells PASP Submission

**Date:** 2026-05-14 (end of session)
**Author:** Ron Bibb (ORCID 0009-0004-1153-2464)
**Manuscript:** *Localized Residual Structure in SPARC Rotation Curves: A Burkert-Backbone Model with BIC-Selected Shells*
**Target journal:** PASP (after AJ desk-reject; ApJ as fallback if PASP redirects)
**Soft submission deadline:** Sunday/Monday 2026-05-17 to 05-18

---

## 1. Current Manuscript State

### File on disk
- **Path:** `/mnt/user-data/outputs/PASPsample701_fixed.tex`
- **Lines:** 504
- **Size:** 91,609 bytes
- **MD5:** `d437194b8ff80c04f8301bf55e695fd8`
- **Compiles in:** AASTeX 7.0.1 (`\documentclass[linenumbers]{aastex701}`) on Overleaf
- **Last compile observed:** 22-page PDF, all 5 figures + 7 tables + 24-entry bibliography rendering correctly

### Compile status
- ✅ 0 unresolved cross-references (all 27+ `\ref{}` calls have matching `\label{}`)
- ✅ All custom macros defined (`\Yt`, `\Ydisk`, `\Ybulge`, `\Vobs`, `\Vbar`, `\Vgas`, `\Vdisk`, `\Vbulge`, `\VDM`, `\Vflat`, `\Mstar`, `\Mhalo`, `\Msun`, `\chired`)
- ✅ `\editornote{}` defined as no-op (review markup invisible in PDF)
- ✅ `\new{}` defined as passthrough
- ✅ Bibliography style `aasjournalv7`
- ✅ Float layout cleaned (tables clearpaged before bibliography, oversized figures height-constrained)
- ✅ `\phantomsection` added before `\section*` labels (fixes pdfTeX bookmark warnings)
- ✅ AAS AI-disclosure policy compliant acknowledgment text in `\begin{acknowledgments}` *(NOT YET ADDED — see Section 5 checklist item 4)*

### Recent edits applied (this session)
1. Abstract structural-qualifier sentence added: *"Because the framework reduces exactly to Burkert when BIC selects zero shells, the 26-galaxy gap against Burkert is empirical rather than structurally guaranteed."*
2. Shell terminology softened in 5 load-bearing claim sentences: "localized residual structure" → "localized residual features" or "BIC-selected localized residual components" in abstract opening/closing, §1 framework introduction, §1 MOND framing, conclusion opening
3. §1 scope paragraph added: "The scope of this paper is methodological and demonstrative..."
4. §4.4 marginalization closing extended: added clarifying sentence that structural claim survives even if adequacy gap compresses
5. §4.5 "Summary of limitations and additional channels" subsection added
6. NGC 5055 marginalization paragraphs in §3.4 confirmed complete (verified line-by-line, all paragraphs end with periods)

### Manuscript structure (current)
- §1 Introduction (with scope paragraph)
- §2 Methods (5 subsections)
- §3 Results (7 subsections including Aggregate, Win/loss, Morphology gradient, NGC 5055 showcase, Synthetic null test, Example grids, Einasto robustness)
- §4 Limitations (5 subsections including new §4.5 Summary)
- §5 Discussion (6 subsections)
- §6 Conclusion
- Acknowledgments + Software + Data and Code Availability
- 7 Tables (after `\clearpage`)
- Bibliography (`aasjournalv7`, 24 entries from `full_refs.bib`)

### Figures (in `figures/` Overleaf subdirectory)
| Paper Fig # | Section | File | Notes |
|---|---|---|---|
| 1 | §3.2 | `figure5_dbic_histogram.pdf` | ΔBIC distribution |
| 2 | §3.3 | `figure4_shellfrac_vs_T.pdf` | Morphology gradient |
| 3 | §3.4 | `figure2_ngc5055_showcase.png` | NGC 5055 single-galaxy fit |
| 4 | §3.6 | `figure1_T4_example_grid.png` | T=4 grid (height-constrained) |
| 5 | §4.1 | `figure3_universal_failures.png` | 9-galaxy panel (height-constrained) |

File numbers reflect creation order in v7.0; paper figure numbers come from `\includegraphics` placement. **Do NOT rename the image files** — preserves v7.0 Zenodo reproducibility.

---

## 2. Repo State

- **GitHub:** `github.com/RonBibb/sparc-halo-shells`
- **Current tagged release:** v7.1.0 (Zenodo version DOI: `10.5281/zenodo.20100475`)
- **Concept DOI:** `10.5281/zenodo.20072882`
- **Next tag to cut:** v7.1.1 (includes hierarchical marginalization scripts/results)

### v7.1.1 additions vs v7.1.0 (already in repo, ready to release)
- `scripts/ngc5055_marginalized_upsilon.py` (single-galaxy joint Υ marg)
- `scripts/hierarchical_upsilon_marginalization.py` (1032 lines, full-sample empirical Bayes)
- `scripts/summarize_hierarchical_results.py` (433 lines)
- `data/ngc5055_marginalized_upsilon_results.csv`
- `data/ngc5055_fixed_upsilon_recheck_results.csv` (sanity check)
- `data/hierarchical_marginalization_per_galaxy.csv` (303 rows; 101 galaxies × 3 iterations)
- `data/hierarchical_marginalization_hyperprior_history.csv`
- `data/ngc2841_marginalized_upsilon_results.csv` (supplementary; prior-escape failure mode)
- README.md and DATA_PROVENANCE.md updated to v7.1.1
- `.gitignore` committed (catches `.DS_Store` and `__pycache__/`)

### Code reproducibility
- All canonical fits reproducible with `scripts/run_canonical_pipeline.py` (or whatever your runner is named)
- Hierarchical marginalization converges in 3 iterations to bit-stable hyperparameters
- NGC 5055 marginalization is seed-stable across {1, 42, 99999, 20260513} for 0-shell and 1-shell BIC (ΔBIC = +0.4815 in every run)

---

## 3. Key Scientific Results (Canonical-Υ Numbers — Locked)

### Headline adequacy comparison (canonical Υ_disk=0.5, Υ_bulge=0.7)
- Framework: **91/102** adequate at χ²_red < 1.5 (89.2%)
- Burkert-only: **65/102** (63.7%)
- NFW (free c): **52/102** (51.0%)
- At stricter χ²_red < 1.0: 86/57/41 — gap widens

### Morphology gradient
- T=2 (n=9): 100% shell-bearing
- T=3 (n=11): 55%
- T=4 (n=17): 53%
- T=5 (n=13): 54%
- T=6 (n=16): 56%
- T=7 (n=14): 36%
- T=8 (n=6): 33%
- T=9 (n=16): 31%
- Per-galaxy T-label permutation p_two-sided = 0.002
- Per-bin Spearman ρ = -0.83, p = 0.010; Mann-Kendall p = 0.003
- **T=2-excluded sensitivity:** per-galaxy permutation p = 0.094; per-bin Spearman p = 0.052; Mann-Kendall p = 0.069 (already in §3.3)

### NGC 5055 canonical-Υ fits (Table 5)
- Burkert-only: χ² = 190.6, χ²_red = 9.53, BIC = 196.8
- NFW (free c): χ² = 604.1, χ²_red = 30.21, BIC = 610.3
- DC14 (3 free): χ² = 171.1, χ²_red = 9.01, BIC = 180.4
- Framework (1 shell at r=12.0 kpc, σ=2.91 kpc, M=3.78×10¹⁰ M_⊙): χ² = 24.9, χ²_red = 1.46, BIC = 40.3
- ΔBIC vs framework: +156.5 (Burkert), +570.0 (NFW), +140.1 (DC14)

### Synthetic null test (matched/mismatched backbone)
- Burkert-truth: 7/173 mocks select shells (4.0% aggregate FP rate)
- NFW-truth: 110/173 mocks select shells (63.6% aggregate)
- Per-bin Burkert-truth FP: ≤ 5% in 6 of 8 T-bins; uncorrelated with T-type (Spearman ρ = +0.33, p = 0.42)

### Einasto-backbone robustness (§3.7)
- Einasto framework: 89/102 adequate
- Einasto-only: 72/102 adequate
- **Adequacy gap compresses from 26 (Burkert) to 17 (Einasto)**
- Classification agreement with canonical: 90/102 (88.2%)
- Morphology gradient preserved: per-bin ρ = -0.905, p = 0.002; per-galaxy ρ = -0.347, p = 0.0004
- NGC 5055 still requires shells: Einasto-only χ²_red = 8.42; Einasto+shells χ²_red = 0.51

### DC14 universal-failure subset (Table 7)
- 10/10 fail χ²_red < 1.5 adequacy threshold
- 8/10 hit upper boundary at X = -1.5 (DC14 → NFW limit)
- Range: 2.02 (NGC 0289) to 28.58 (UGC 06787); median 4.37
- Framework achieves lowest χ²_red on 9/10 (NFW marginally beats framework on UGC 00128)

### Hierarchical Υ marginalization (v7.1.1, NOT YET IN MANUSCRIPT)
- Converged at iteration 2 (all |Δhyperparameter| < 0.01 dex)
- μ_disk = -0.315 (Υ_disk population mean = 0.485); τ_disk = 0.081 dex
- μ_bulge = -0.125 (Υ_bulge population mean = 0.749); τ_bulge = 0.081 dex
- **Framework adequate: 86/101** (vs canonical 91/102)
- **Burkert-only adequate: 72/101** (vs canonical 65/102)
- **Adequacy gap: 14** (marginalized) vs 26 (canonical) — compresses 46% but **preserved**
- Morphology gradient preserved: per-galaxy ρ = -0.245, p = 0.013
- *NGC 6674 excluded due to known fit pathology*

---

## 4. Open Vulnerabilities (Four-Reviewer Consensus)

Four independent reviewers reviewed v3.x of the manuscript. Convergence pattern was strong.

### Tier-1 vulnerabilities (worth addressing before submission)

| Vulnerability | Reviewer agreement | Mitigation status |
|---|---|---|
| NGC 5055 ΔBIC collapse under joint Υ (157 → +0.5) undermines showcase | R1, R2, R3, R4 (unanimous) | **Partially addressed.** §3.4 discloses; §4.4 has Per-galaxy joint Υ marg subsection. **Hierarchical-Υ population result (already done, not in manuscript) would fully defuse.** |
| T=2 doing disproportionate work in gradient signal | R1, R2, R4 | **Addressed in §3.3** (T=2-excluded sensitivity numbers present). Could be promoted to clearly labeled finding rather than mid-paragraph prose. |
| Base-profile shape coupling / Einasto compression (26→17) | R1, R4 | **Addressed via §3.7 "Fourth finding"** but framing is defensive ("remains substantial"). Could be inverted to strength via core-result promotion. |
| DC14 framing weakness (boundary-pegging defense reads as concession) | R1, R2 | **Acknowledged in §4.1**; calibration-domain caveat present but not fully defended. Best fix is out of scope for this submission (would require DC14 on low-mass subsample). |
| Shell terminology emotional pushback | R3, R4 | **Mostly addressed this session** — load-bearing claim sentences softened to "features/components"; interpretation context kept "structure". |
| "26 rescues, 0 reverse" structural artifact | R2 | **Addressed this session** — abstract qualifier sentence added explaining superset architecture. |
| AI disclosure (AAS policy) | R1 | **NOT YET ADDRESSED** — current acknowledgment line says "informal feedback from a reviewer who read v7.1.0" without disclosing LLM involvement. Required by AAS 2024-25 policy. |
| Long-paragraph readability in §3.3/§3.4 | R3 | **NOT YET ADDRESSED.** Pure mechanical readability fix. |
| Multi-restart convergence diagnostics thin | R2 | **NOT YET ADDRESSED.** One sentence in §2.5 closes this. |
| Universal-failure model ceiling framing | R4 | **NOT YET ADDRESSED.** Could be strengthened to "this is precisely why we flag these for 2D follow-up." |

### Tier-3 vulnerabilities (out of scope for this submission)

- **Generalized-NFW backbone test** (R1, R4): would extend Einasto robustness check. ~1-2 weeks work. Defer to Paper 2 or referee revision.
- **2D velocity-field test on a universal-failure galaxy** (R4): would convert universal-failure caveat into positive result. Months of work; Paper 3+ territory.
- **DC14 low-mass calibration-domain subsample** (R1): symmetric DC14 treatment. ~1 week.

### Four-reviewer acceptance odds estimate (consensus)

| State | PASP | ApJ |
|---|---|---|
| As-is (current state) | ~35-45% | ~25-35% |
| + Tier-1 fixes (cheap edits, no new analysis) | ~55-65% | — |
| + Hierarchical-Υ paragraph added | ~60-70% | ~40-50% |
| + Core-results promotion (Einasto + hierarchical-Υ moved out of robustness checks) | ~65-75% | ~50-60% |

---

## 5. Remaining Checklist

### Tier 1 — Editorial fixes (cheap, no new analysis)
- [x] **1. Abstract "26 vs 0 reverse" structural qualifier** ✅
- [x] **2. Shell terminology softening** ✅
- [ ] **3. §3.4 marginalization caveat forward reference** (10 min) — add one sentence at §3.4 opening: "We note upfront that the per-galaxy BIC preference reported here is sensitive to Υ marginalization; the full disclosure is in §4.4, and the load-bearing claim is the population-level adequacy comparison (§3.3), not this single-galaxy result."
- [ ] **4. AI disclosure in acknowledgments** (10 min) — required by AAS 2024-25 policy. Draft text: *"We acknowledge use of Anthropic's Claude as a sounding board for methodological feedback and presentation refinement during preparation of this manuscript; all scientific judgments, fits, and conclusions are the author's responsibility."* Consider web-searching current AAS AI disclosure language to match exactly.
- [ ] **5. Long-paragraph splits in §3.3 and §3.4** (45 min) — worst offender: §3.4 "A subtler but informative finding" paragraph (200+ words, nested clauses). Three deliberate paragraph breaks improve referee readability.
- [ ] **6. T=2-excluded sensitivity check promotion** (15 min) — already in §3.3, just reformat as clearly labeled finding rather than buried mid-paragraph.
- [ ] **7. Multi-restart convergence sentence in §2.5** (20 min) — close R2's "are you sure multi-restart found global optimum?" with one sentence quantifying restart-to-restart spread.
- [ ] **8. Universal-failure framing strengthening** (30 min) — reframe from "model ceiling" to "diagnostic flag for 2D follow-up."

### Tier 1.5 — High-leverage content addition
- [ ] **9. Hierarchical-Υ §4.4 paragraph** (2-3 hours) — write a ~200-word paragraph reporting the 86/101, 72/101, gap=14, ρ=-0.245 p=0.013 result. Numbers from `data/hierarchical_marginalization_per_galaxy.csv`. Single highest-leverage edit available. Defuses unanimous NGC 5055 critique.

### Tier 2 — Structural reframe (optional, Saturday decision)
- [ ] **10. Promote Einasto + hierarchical-Υ to core §3 results** (3-4 hours) — restructure abstract, §1 preview, §3 organization so Einasto and hierarchical-Υ are co-equal demonstrations alongside canonical Burkert, not "robustness checks." Inverts two critiques from weakness to strength. Worth doing if Saturday morning has high-quality focused time available.

### Pre-submission infrastructure tasks
- [ ] **11. ADS verification of bibliography** (~30 min) — confirm all 24 entries in `full_refs.bib` match NASA ADS bibcodes. Two flagged with `note = {VERIFY ...}` annotations: `Li2019` and `Reines2013`.
- [ ] **12. Cut v7.1.1 GitHub release tag → Zenodo deposit** (~10 min) — get frozen v7.1.1 DOI.
- [ ] **13. Update §Data and Code Availability** with v7.1.1 Zenodo DOI once issued (~5 min).
- [ ] **14. Pipeline validation pass** (~1-2 hours) — confirm canonical numbers in manuscript match latest pipeline output. **High-priority for what changed v7.1.0 → manuscript:** NGC 5055 marginalization values, hierarchical Υ numbers, §3.6 Einasto adequacy paragraph (89/72/17). Old numbers (91/65/52, morphology gradient, null test, DC14 universal-failure) should not have drifted but worth spot-checking.
- [ ] **15. Draft cover letter** (~45 min) — ~300 words. Frame as methodology paper. Lead with: BIC-selected localized residual framework, synthetic null calibration under matched/mismatched backbones, formal disclosure of per-galaxy Υ degeneracy. Mention in-prep Papers 2 and 3 briefly. Address potential PASP vs ApJ scope question proactively.
- [ ] **16. PASP ScholarOne portal walkthrough** (~1 hour) — at `https://mc04.manuscriptcentral.com/pasp`. New ScholarOne account likely needed. UAT keywords entry (6 already in manuscript). Abstract paste (277 words, under 300 limit).

### Outstanding journal clarification
- [ ] **AAS manuscript #AAS76831 journal designation** — portal shows Journal: ApJ but Desired Journal: AJ. Original clarification request was sent to `journals.manager@aas.org`. Status: pending. **For PASP, this is a new submission, so #AAS76831 is no longer relevant.** Old AJ desk-reject is closed.

---

## 6. Submission Strategy

### Three-paper program (locked)
1. **Paper 1 (this manuscript) → PASP** — methodology + canonical-Υ analysis + marginalization disclosure
2. **Paper 2 (in prep)** → AJ or PASP, 1-2 weeks after Paper 1 editorial-triage signal. Canonical-Υ shell catalog organization. Draft exists as `paper2.md` but not refined since v7.0. Topics: shell geometry (spherical vs torus rejected), two-shell galaxy candidates (16 total, NGC 2841 and UGC 11914 as best pair), SMBH proxy vs DM halo relationships.
3. **Paper 3 (TBD)** → hierarchical marginalization confirmation; full-sample Bayesian treatment. Analysis complete (v7.1.1 results), no draft yet.

**Staggered submission:** Paper 1 to PASP first. Paper 2 follows 1-2 weeks after editorial signal. Do not submit concurrently.

### Recommended timeline (3 paths, decide based on energy)

**Path A — Submit as-is Monday:** No further edits. ~35-45% PASP acceptance odds. Drag vulnerabilities into referee response.

**Path B — Tier 1 only, submit Tue/Wed:**
- Friday morning: Pipeline validation + ADS verification + Tier 1 items 3-8
- Friday afternoon: Cut v7.1.1 → Zenodo
- Saturday morning: Hierarchical-Υ paragraph (item 9)
- Saturday afternoon: Cover letter
- Sunday: Final compile review, rest
- Monday/Tuesday: Submit
- Acceptance odds: ~60-70% PASP

**Path C — Tier 1 + restructure, submit Wed:**
- Same as Path B plus:
- Saturday morning (replaces hierarchical paragraph alone): Full restructure (Einasto + hierarchical-Υ promoted to core results, abstract rewritten, §3 reorganized)
- Saturday afternoon: Cover letter
- Sunday: Final compile review, rest
- Wednesday: Submit
- Acceptance odds: ~65-75% PASP
- **Risk:** Half-done restructure is worse than no restructure. Requires sustained high-quality focus.

**Recommendation:** Path B if Saturday energy is uncertain. Path C if Saturday morning will be focused.

### PASP vs ApJ fallback

PASP fit is debatable per reviewers R1 and R3. PASP scope is methods/tools/instrumentation; this paper has substantial empirical content. **If desk-redirected to ApJ:** 2-week reset, not a rejection. Cover letter should preemptively lead with methodological contributions to defend PASP framing.

---

## 7. Key Strategic Decisions Locked Earlier

- **No M_BH framework, no two-component decomposition.** Earlier framework work (Papers A/B/C/D, k_SMBH coupling) explicitly out of scope. Paper 1 is exclusively halo_shells methodology.
- **First submission uses `linenumbers` only** (drop `trackchanges`). Use `trackchanges` only when responding to PASP referees.
- **Don't rename figure files.** Keep v7.0 naming and let LaTeX `\includegraphics{}` placement set paper figure number. Preserves v7.0 Zenodo reproducibility.
- **MNRAS not in target set** (pattern of past rejections as independent researcher). PASP primary; Open Journal of Astrophysics as long-term alternative for future work.
- **Independent researcher status:** ORCID 0009-0004-1153-2464 in author block. No institutional affiliation. PASP/AAS journals supportive of independent researchers; cover letter should not apologize for this.

---

## 8. Critical Tier Hierarchy for Claims (Internal Discipline)

Locked tier hierarchy for what the paper can claim at what confidence level:

- **Tier 1 (locked):** Shells are real — BIC-selected components capture features smooth profiles miss
- **Tier 2 (cosmetic):** Parameterization is cosmetic — Burkert vs Einasto backbone choice doesn't change the qualitative result
- **Tier 3 (substantive):** Geometry is substantive — coplanar-torus geometries decisively rejected (NGC 5055 ΔΧ² = +82)
- **Tier 4 (speculative):** Mechanism is speculative — physical origin of shells remains open

Claims should only be asserted at the tier they've earned. Paper 1 should mostly stay in Tier 1, occasionally Tier 2 (Einasto backbone substitution), with explicit "speculative" framing when discussing physical interpretation.

---

## 9. Where I Am, Where I'm Going

**Where I am:** Manuscript is structurally sound, content-locked at v3.2.6 level, compiles cleanly with all macros and cross-references working. Items 1 and 2 of Tier 1 checklist done this session. Confidence in canonical-Υ numbers high (verified pre-PASP). Confidence in v7.1.1 marginalization numbers high (seed-stable, hierarchical converged).

**Where I'm going:** Friday-Saturday execute checklist items 3-9 plus pre-submission infrastructure (ADS, Zenodo, cover letter). Submit Monday (Path B) or Wednesday (Path C). Monday/Tuesday buffer.

**Biggest leverage available:** Item 9 (hierarchical-Υ paragraph) and item 10 (core-results promotion). Both use data already in hand. Item 9 alone moves acceptance odds by ~10%. Item 10 adds another ~5-10%.

**Biggest risk:** Half-done restructure (item 10) introducing inconsistency between abstract and body. Don't attempt item 10 unless Saturday morning has sustained focus.

---

## 10. Contact / Repo Pointers

- **Manuscript:** `/mnt/user-data/outputs/PASPsample701_fixed.tex` (this session's working version)
- **Bibliography:** `/mnt/user-data/outputs/full_refs.bib`
- **README and DATA_PROVENANCE:** in `/mnt/user-data/outputs/` (already pushed to repo at v7.1.0/v7.1.1)
- **GitHub:** `github.com/RonBibb/sparc-halo-shells`
- **Zenodo concept DOI:** `10.5281/zenodo.20072882`
- **v7.1.0 version DOI:** `10.5281/zenodo.20100475` (current)
- **v7.1.1 version DOI:** TBD (to be minted Friday)
- **ORCID:** `0009-0004-1153-2464`
- **PASP portal:** `https://mc04.manuscriptcentral.com/pasp`
- **Author email (in manuscript):** `ronbibb@gmail.com`

---

*End of transfer document. Generated 2026-05-14 end-of-session.*
