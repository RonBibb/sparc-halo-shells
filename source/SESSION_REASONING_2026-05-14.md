# Session Reasoning Document — PASP Submission Preparation

**Date:** 2026-05-14 evening through late night
**Scope:** Trains of thought during PASP submission preparation; companion to PAPER1_TRANSFER_DOCUMENT.md (which captures state)
**Purpose:** Capture *why* we made each decision, not just *what* we decided. Useful for future sessions resuming this work or for understanding the strategic frame if circumstances change.

---

## 1. Where the Session Started

Came in to a manuscript that compiled but had three rendering problems: figures rendering correctly but bibliography and tables interleaving across the back-matter pages, plus some "float too large" and "stuck float" warnings in the compile log. The .tex had previously had all the macro fixes (\\Yt, \\Ydisk, \\Ybulge, \\Vobs, \\Vbar, etc.) and the structural skeleton from the migration to aastex701.

Ron's framing was practical: fix the warnings, ensure NGC 5055 marginalization paragraph was complete, add a scope paragraph and limitations paragraph.

What emerged over the course of the session was substantially bigger than that initial framing — four independent reviewer rounds came in and reshaped the strategic landscape.

---

## 2. The Compile-Layout Problem and How We Diagnosed It

The interleaving wasn't random. The pattern in the rendered PDF was: Tables 1-3 on page 20, then "REFERENCES" header appears mid-page, then bibliography entries start, then Tables 4-6 on page 21, then more bibliography, then Table 7 on page 22.

The diagnosis: LaTeX's float-deferral system. With 7 sequential `[h]` (here-please) table environments, LaTeX places what it can on the first available page, then defers the rest. When the bibliography environment is hit, the bibliography starts. The deferred table floats then get placed wherever they fit — which means around the bibliography text.

The fix was three-part:
1. Update `\graphicspath{}` to include `figures/` (Ron's actual directory) alongside `images/` and `./` for safety
2. Remove `[h]` from all table environments so LaTeX uses default `tbp` placement
3. Add `\clearpage` right before `\bibliographystyle{aasjournalv7}` to force all pending floats to flush before the bibliography begins

The key insight worth remembering: **`\clearpage` before the bibliography is the canonical fix for table-bibliography interleaving in aastex documents.** This will come up again with other AAS submissions.

For the figure warnings ("Float too large by 72.6pt"), the issue was the T=4 grid figure (891×1192 pt natively, taller than wide) requested at `width=\linewidth` which makes it 686pt tall — just over the available text height. Fix: add `height=0.85\textheight,keepaspectratio` so LaTeX caps it at 85% of page height. Same fix for the universal-failures figure as a preventative measure.

For the `pdfTeX warning (dest): name{section*.2}` warning, the issue was that `\section*{}` doesn't create a real PDF destination, so `\label{}` on it has nothing to anchor to. Fix: add `\phantomsection` before each `\section*{}` that has a label.

**Pattern worth internalizing:** when a compile warning is purely about float placement or layout, the fixes are almost always (1) remove placement specifiers, (2) add `\clearpage` or `\phantomsection` at strategic points. Cosmetic warnings about class chatter (e.g., "Repairing hyperref-unfriendly LaTeX definition") can be ignored.

---

## 3. The Three Additions Ron Requested

Three content additions came up early in the session:

**A. Scope paragraph at end of §1.** This was an editorial decision to set explicit boundaries: state what the paper is (methodological, demonstrative), what it isn't claiming (DM origin, per-galaxy event attribution), and what the principal empirical content is (population-level adequacy comparison). Three to four sentences, deliberate placement after the "paper organization" paragraph.

The drafting principle: a scope paragraph should be the first thing a skeptical referee can point to when asking "what is this paper claiming?" If it answers cleanly, the rest of the review is about evidence quality, not scope ambiguity.

**B. Marginalization section closing.** Ron flagged that the §4.4 Per-galaxy joint Υ marginalization paragraph "cuts off mid-sentence." Investigation showed the .tex actually terminated properly with a period — every paragraph in both §3.4 and §4.4 had clean sentence endings.

Two possible explanations: (1) Ron's local Overleaf compile was showing different content than my source (stale .aux cache, or he was viewing an older compile), or (2) the section *read* as cut-off because it ended on an open question ("we have not resolved... future paper should address...") without bringing closure.

Decision: assume #2 and extend the closing paragraph with two sentences that bring proper closure. The added text said: even under the conservative assumption that hierarchical-Υ compresses the adequacy gap, the framework's *structural* claim is unaffected — BIC-selected shells capture features smooth profiles cannot absorb at any Υ. The question marginalization addresses is what fraction of features admit alternative low-Υ descriptions, not whether features are present.

This was the right move for two reasons. First, even if Ron's compile was bad, fixing the content's epistemic closure was useful independently. Second, the extension does substantive work — it pre-empts the "but what if marginalization undermines everything?" critique that the four-reviewer rounds confirmed was the unanimous concern.

**C. Limitations summary in §4.** §4 already had subsections on universal-failure galaxies, per-galaxy decomposition limits, base-profile shape coupling, and stellar-population systematics. The new §4.5 "Summary of limitations and additional channels" was framed to acknowledge three additional limitation channels not covered by §§4.1-4.4: (1) BIC vs alternative model-selection criteria, (2) the architectural choice of σ/r ≤ 0.4 threshold, (3) heterogeneous SPARC uncertainties under uniform treatment.

The drafting principle for limitations sections: name the limitation, note its likely magnitude, state the load-bearing claim that survives it. Don't catastrophize, but don't elide. The §4.5 paragraph ends with the load-bearing claim restated — "the population-level adequacy comparison is the load-bearing claim, robust to the steelmans considered above; per-galaxy interpretations should be read as constrained by these channels rather than as independently verified detections."

---

## 4. The Four-Reviewer Consensus Analysis

Over the course of the session, four independent reviews came in. The convergence pattern was strong and informative.

The agreement matrix that emerged:

- **Unanimous (all 4):** NGC 5055 marginalization collapse undermines the showcase
- **3 of 4:** T=2 doing disproportionate work in the gradient signal
- **2 of 4 (different framings):** Base-profile shape coupling / Einasto compression; DC14 framing weakness; Shell terminology ambiguity
- **Single-reviewer:** AI disclosure (R1); software lacks MCMC (R1); "26 rescues, 0 reverse" structural artifact (R2); multi-restart diagnostics thin (R2); long paragraphs (R3); framework-vs-interpretation oscillation (R3); universal-failure model ceiling (R4); generalized-NFW request (R1, R4)

The first analytical move that mattered: **distinguish convergent critiques from single-reviewer concerns.** Five items showed up in two or more reviews. Those are real signal — what an actual PASP referee would also see. Single-reviewer items are real perspectives but represent individual emphasis, not collective signal.

The second analytical move: **rank by impact on acceptance odds, not by reviewer emphasis.** The most damaging critique is the one most likely to get the paper desk-rejected or sent for major revision. The unanimous NGC 5055 marginalization concern was clearly top — if a referee reads §3.4 as a load-bearing showcase that collapses in §4.4, they may not finish the paper.

The third analytical move: **separate what's expensive to fix from what's cheap to fix.** Most of the critique surface turned out to be addressable through editorial work (relocation, qualification, terminology softening) rather than new analysis. That changed the strategic frame from "we need new work" to "we need careful presentation of work we've already done."

**Pattern worth remembering:** when multiple independent reviewers converge, the convergence is the signal. When they disagree, the disagreement reveals individual emphasis. Both are useful, but they have different weights.

---

## 5. The Hierarchical-Υ Realization

This was the pivotal moment of the session.

Reviewer 2's biggest critique was about NGC 5055 marginalization collapsing the showcase and Ron not having run joint Υ marginalization on the full sample. Reviewer 4's "systematic degeneracies" concern was structurally the same critique.

Then I recalled (from memory of the project state): **Ron has already done full-sample hierarchical empirical-Bayes Υ marginalization.** The results are sitting in `data/hierarchical_marginalization_per_galaxy.csv` and the summary numbers are 86/101 framework adequate, 72/101 Burkert-only adequate, gap = 14, ρ = -0.245 p = 0.013, converged in 3 iterations. The decision had been made to keep it out of Paper 1 to keep scope clean for the PASP submission.

The realization: **the single biggest critique from the four-reviewer consensus is addressable with a 2-3 hour paragraph using data already in the repo.** This isn't new analysis. It's writing.

This reframes the whole strategic landscape. Without the hierarchical paragraph, the manuscript's strongest critique remains unrebutted. With it, the unanimous concern transforms from "your single-galaxy showcase collapses" to "we've tested this at population scale and the gap survives with magnitude 14."

This led to the conclusion that the hierarchical-Υ paragraph isn't optional — it's the single highest-leverage edit available before submission.

---

## 6. The Structural Reframing Discussion

A reviewer recommendation came in after we'd been thinking about Tier-1 fixes: "pre-emptively move the Einasto-backbone and joint-Υ marginalization tests from 'robustness checks' to core results to show you have already addressed the primary sources of systematic bias."

This wasn't another vulnerability item. It was a structural recommendation about how the paper makes its claim.

The current paper's structure makes the canonical 91/102 vs 65/102 result load-bearing. Two-thirds of the four-reviewer critique surface attacks the fragility of relying on a single canonical comparison.

The restructured paper would make a **triangulated adequacy claim** load-bearing:
- Canonical Burkert: 91/102 vs 65/102 = 26-galaxy gap
- Einasto-backbone: 89/102 vs 72/102 = 17-galaxy gap
- Hierarchical-Υ marginalized: 86/101 vs 72/101 = 14-galaxy gap

The magnitude compresses (26 → 17 → 14, a 46% reduction), but the sign doesn't change. Three different smooth-profile and stellar-population assumptions, three positive adequacy gaps. The compression itself becomes the story.

**The analytical insight worth preserving:** restructuring inverts two critiques from weakness to strength. Einasto compression currently reads as erosion of the framework's advantage; in the restructured version, it reads as "the gap survives under a more flexible backbone." Hierarchical-Υ marginalization currently reads as a §4.4 confession; in the restructured version, it reads as a co-equal core result demonstrating the gap survives marginalization.

Five of seven convergent critiques get meaningfully better under restructuring. Two get smaller improvements. None gets worse.

Trade-off analysis: 3-4 hours of focused work for a ~10% acceptance probability gain. The risk is half-doing it — a partially restructured paper with abstract promising triangulation and a body that still buries Einasto in §3.7 is worse than no restructuring at all. So the decision became conditional on Saturday energy quality.

**Pattern worth remembering:** when reviewer critiques converge on "your load-bearing claim is fragile," sometimes the fix isn't strengthening the claim but distributing the load across multiple converging claims. Triangulation defuses fragility critiques more effectively than defensive arguments do.

---

## 7. The "Largely Defused" Precision Problem

When I first laid out how restructuring interacts with the convergent critiques, I used "largely defused," "fully defused," and "inverted" without distinguishing them clearly. Ron pushed back asking what "largely defused" actually meant.

The honest answer required decomposing each critique into sub-components and tracking which sub-components the restructuring addresses.

The NGC 5055 critique decomposed into three parts:
1. **Showcase fragility** — the framing of NGC 5055 as "the showcase." Fully defused by demoting it to one of three convergent demonstrations.
2. **Numerical collapse from ΔBIC 157 → +0.5** — the actual fact of the collapse. Not defused; that fact survives any restructuring.
3. **Implicit "showcase as proof" rhetoric** — the rhetorical weight on a single galaxy. Fully defused if restructuring is committed to throughout.

So "largely defused" meant: two of three sub-components fully defused, one survives. A referee can still observe the collapse, but they can no longer make it a fatal critique because no load-bearing claim rests on the single galaxy.

This precision matters because it determines what the restructured paper still needs to say. The collapse itself still needs disclosure (§3.4 and §4.4 both still acknowledge it), but the disclosure can be matter-of-fact rather than apologetic.

**Pattern worth internalizing:** when defending a paper against critique, decompose the critique. "Your X is wrong" often has three or four sub-claims, and each may have a different remediation. Treating it as a single monolithic critique misses where the real leverage is.

---

## 8. The Cost-Benefit Triage

When Ron asked "is it worth making the common changes?" the response had to be honest about which fixes pay off and which don't.

The framework that emerged:
- **Cheap items (under 3 hours total):** Unambiguous yes — these are pure editorial work, no new analysis, no risk
- **Hierarchical-Υ paragraph (2-3 hours):** Non-negotiable — addresses the unanimous critique
- **Full restructuring (3-4 hours):** Conditional — worth it if Saturday morning has high-quality focus; risky if attempted in fragmented time
- **Out-of-scope (weeks to months):** Generalized-NFW test, 2D velocity-field follow-up, DC14 calibration-domain test — defer to Paper 2/3 or referee revision

The decision criterion that emerged for the restructuring: **not "is it worth doing in the abstract" but "do you have 4 hours of high-quality focus available."** This is the right kind of criterion because it accounts for execution risk, not just theoretical upside.

The actual recommended path: Tier 1 + 2 ("Path C") if Saturday morning is focused; Tier 1 + hierarchical paragraph alone ("Path B") otherwise. Both are defensible. Acceptance odds approximately 60% (B) vs 70% (C).

**Pattern worth remembering:** when advising on revision strategy, account for the human's actual energy state, not just the theoretical optimal. A partial restructuring at 80% energy is worse than a complete tier-1-only revision at 95% energy.

---

## 9. The Checklist Construction

The checklist that emerged from triaging four-reviewer feedback:

**Tier 1 — Editorial fixes:**
1. Abstract "26 vs 0 reverse" structural qualifier (10 min)
2. Shell terminology softening (15 min)
3. §3.4 marginalization caveat forward reference (10 min)
4. AI disclosure in acknowledgments (10 min)
5. Long-paragraph splits in §3.3 and §3.4 (45 min)
6. T=2-excluded sensitivity check promotion (15 min)
7. Multi-restart convergence sentence in §2.5 (20 min)
8. Universal-failure framing strengthening (30 min)

**Tier 1.5:**
9. Hierarchical-Υ §4.4 paragraph (2-3 hours)

**Tier 2:**
10. Promote Einasto + hierarchical-Υ to core §3 results (3-4 hours)

The ordering principle was cheapest-first. This isn't because the early items matter most — item 9 (hierarchical paragraph) is the highest-leverage by a wide margin. The cheapest-first ordering is because each completed item builds momentum, and because if energy runs out partway through, you've still made progress on the easy wins.

We executed items 1 and 2 cleanly. Item 1 was the abstract structural qualifier ("Because the framework reduces exactly to Burkert when BIC selects zero shells, the 26-galaxy gap against Burkert is empirical rather than structurally guaranteed.") — a single sentence addition with 14-word footprint, taking the abstract from 263 to 277 words, well under the 300 PASP limit.

Item 2 (shell terminology softening) was more interesting because it required principle-setting: keep "shell" as working shorthand throughout the body (the §1 generic-label paragraph already does heavy lifting explaining the terminology), but soften the load-bearing claim sentences where "structure" implies physical reality. The pattern that emerged: when the paper *claims* what's measured, say "features" or "components" (statistical descriptors). When the paper *interprets* what those might physically be, say "structure" (the hypothesis being tested).

This led to softening 5 spots: abstract opening, abstract closing, §1 framework introduction, §1 MOND framing, conclusion opening. Three "structure" usages were kept where appropriate: §3.7 Einasto methods, §5.2 Interpretation subsection (heading and content), §5.6 phenomenological-test framing.

**Pattern worth preserving:** when softening terminology, the right discipline isn't global find-replace. It's identifying which usages are load-bearing claims (where softening pays off) versus which are interpretation context (where the loaded word is appropriate). Soften the former; leave the latter alone.

---

## 10. The Item 7 Investigation

When Ron asked about item 7 (multi-restart convergence sentence), the analytical move was to recognize that this item *cannot* be written by me alone. It requires empirical data from Ron's pipeline.

The critique R2 made was sharp: the current §2.5 shows that single-restart is worse than multi-restart (NGC 2403: 728→193; NGC 2841: 1547→160). But that proves single-restart is bad, not that multi-restart found the global optimum. A skeptical referee can ask "how do I know your best fit isn't a local minimum?"

The defensive response that closes this question is some empirical statement about restart-to-restart spread — how often the best fit was hit, or how tight the BIC values cluster across the 8-24 restarts per galaxy.

Three options emerged depending on what data Ron's pipeline preserves:
- **Option A (quantitative full sample):** if per-restart BIC is logged for all galaxies
- **Option B (quantitative subsample):** if only best fit is logged, run a diagnostic on 15-20 representative galaxies (~20 min compute)
- **Option C (qualitative):** if neither is feasible, frame as "during development we observed..."

Drafted sentences for all three options so Ron could choose based on pipeline state.

**Pattern worth preserving:** some checklist items need data, not text. When the data exists, the sentence writes itself. When it doesn't, the question becomes whether the data is worth acquiring at the cost of time. For item 7, the diagnostic-run cost (~20 min) is small relative to the referee credibility it buys, so Option B is recommended if Ron's pipeline doesn't log per-restart values.

The session ended with this item still in flight — Ron hadn't confirmed which option his pipeline supports.

---

## 11. The Transfer Document Decision

When Ron asked for a transfer document, the question wasn't just "what's the state of the manuscript." It was "what does someone (Ron in 12 hours, a collaborator, a future Claude session) need to know to resume cleanly?"

The structure that emerged covered:
- Current manuscript state (file metadata, compile status, recent edits)
- Repo state (DOIs, v7.1.1 additions)
- Key scientific results (canonical numbers locked, hierarchical numbers ready)
- Open vulnerabilities (four-reviewer consensus matrix)
- Remaining checklist (16 items across editorial, content, infrastructure)
- Submission strategy (three-paper program, three timeline paths)
- Strategic decisions locked earlier (scope boundaries, formatting choices)
- Tier hierarchy for claims (internal discipline)
- Where I am, where I'm going
- Contact/repo pointers

The document came out at 304 lines / 21 KB. Long but useful — long enough to be self-contained, structured enough to navigate.

**Pattern worth preserving:** a transfer document should answer three questions for the resuming party. (1) What state am I in? (2) What's the next thing to do? (3) What's the strategic frame that determines whether the next thing is still the right next thing? Many transfer documents get (1) right but miss (2) and (3). The strategic frame especially matters because conditions change between sessions.

---

## 12. What the Session Did Not Resolve

Items 3 through 10 of the Tier 1+1.5+2 checklist are still pending. Specifically:
- Item 3 (§3.4 marginalization caveat forward reference) — planned to start when connection died
- Item 4 (AI disclosure in acknowledgments) — same
- Items 5-8 (remaining Tier 1 editorial)
- Item 9 (hierarchical-Υ paragraph) — the highest-leverage edit
- Item 10 (full restructuring) — conditional on Saturday energy

The pipeline-validation pass that Ron planned for 1-2am is also unconfirmed at session end. The data this validation would touch includes the canonical numbers (should not have drifted) and the v7.1.1 marginalization numbers (most important to verify before they appear in the manuscript).

Cover letter, ADS bibliography verification, Zenodo v7.1.1 release, and ScholarOne portal walkthrough are all still pending as pre-submission infrastructure tasks.

---

## 13. Patterns Worth Carrying Forward

Five patterns from this session that are likely to recur:

1. **`\clearpage` before `\bibliography{}` is the canonical fix for table-bibliography interleaving in AAS documents.** Add `\phantomsection` before any `\section*{}` that has a `\label{}`.

2. **When multiple independent reviewers converge, the convergence is the signal. When they disagree, the disagreement reveals individual emphasis.** Both matter, but they should be weighted differently.

3. **When defending against critique, decompose the critique.** "Your X is wrong" usually has 2-4 sub-claims with different remediations. Treating it as monolithic misses leverage.

4. **For terminology softening, identify load-bearing claim sentences (soften) versus interpretation context (leave alone).** Don't do global find-replace.

5. **For revision strategy, account for human energy state, not just theoretical optimal.** A partial restructuring at low energy is worse than complete tier-1 at high energy.

---

*End of session reasoning document.*
