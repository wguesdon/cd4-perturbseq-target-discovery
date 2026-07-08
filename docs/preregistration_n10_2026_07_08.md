# Pre-registration — N10: does the safety gate enrich for genetics-supported autoimmune targets?

**Written 2026-07-08, before `scripts/20_prior_knowledge_nomination.py` was written and before the test
was run.** Committed before the analysis. No number in the decision rule was obtained from the test it
governs. The feasibility numbers in section 3 were measured first and are labelled as such.

This is the project's target-nomination layer. It combines the within-screen gate with external prior
knowledge, which is what the user asked for and what the target-discovery literature says predicts success.

---

## 1. The question, in one sentence

Among the perturbations that clear our two-axis safety gate, does the ranking enrich for genes with human
genetic support for autoimmune disease, **relative to** the naive suppression ranking and to a
magnitude-matched background — or do genetic support and screen-native safety pull in opposite directions?

## 2. Why this is the right test, and why both answers are worth having

Human genetic support is the best-evidenced prior for drug-target success: genetically supported targets
are roughly twice as likely to reach approval (Nelson et al. 2015 *Nat Genet*; replicated King, Davis &
Degner 2019 *PLoS Genet*; refined to ~2.6x by Minikel et al. 2024 *Nature*). A target-nomination layer that
ignores it is indefensible. So we ask whether our gate, built only from the knockdown transcriptome,
aligns with or opposes that prior.

**Both outcomes are publishable and are pre-committed:**

- **ALIGN.** The gate enriches for genetics-supported targets beyond effect magnitude. The decision layer
  adds drug relevance the naive ranking lacks.
- **OPPOSE.** The gate enriches *less* than the naive ranking, or de-enriches. This would be the immunodeficiency
  result arriving from a third independent direction. Genetic association finds genes that *matter* for
  immunity; genes that matter for immunity are genes whose loss causes immunodeficiency; a safety layer must
  reject exactly those. It would say, in one figure: the standard target-discovery prior, applied to an
  inflammatory reversal signature, selects for the targets a safety layer is built to remove. That is
  "reversal is not enough" generalised past reversal, and it is the stronger result.

We commit to reporting whichever fires, with equal prominence, and to leaving this document unedited.

## 3. Data (measured before writing this, for feasibility only)

- **P1 genetic support.** Open Targets 26.06 `association_by_datatype_direct`, rows with
  `aggregationValue == "genetic_association"`, for the **14 autoimmune diseases** the source paper used
  (Zhu, Dann et al., Methods p.50; their table minus the three negative controls). MONDO ids resolved from
  the on-disk `disease.parquet`, all 14 verified present in the association table, which is fully keyed on
  `disease.parquet.id`. Fetched and filtered by `scripts/19_fetch_ot_genetic_evidence.sh` to
  `data/external/open_targets/autoimmune_genetic_assoc.parquet` (gitignored). `targetId` is Ensembl,
  mapped to symbol via the on-disk `target` parquet. **The authors' threshold is genetic evidence >= 0.1**;
  we adopt it, per disease, and define a gene as genetically supported if its `associationScore` >= 0.1 for
  **at least one** of the 14. This is a per-gene label, not the cluster-mediated one in
  `cluster_autoimmune_enrichment_results.suppl_table.csv`, which misses canonical autoimmune genes
  (PTPN22, IL2RA, TYK2, STAT4, CTLA4 all absent) and is therefore unusable as the primary label.
- **P2 loss-of-function tolerance.** `loeuf` and `prec` from `results/tables/window_score_organism_safety.csv`
  (N3, already computed). Reported. LOEUF reads dominant intolerance; `prec` reads recessive, which is what a
  potent inhibitor phenocopies. Minikel et al. 2020 *Nature* is the rationale.
- **P3 tractability.** On-disk Open Targets `target` parquet `tractability` / `safetyLiabilities` /
  `targetClass` (6,289 genes carry tractability). Reported.
- **P4 held-out function.** Schmidt & Steinhart 2022 CD4 IL-2 CRISPRi screen (`priors.schmidt_cd4_il2_screen()`).
  Validation only. It never enters the gate, the score, or any threshold (RULE #3).

**Feasibility, measured first:** the cluster-mediated authors' set gives 267 autoimmune genes in our 6,371
ranked universe but misses canonical targets; the OT per-gene fetch is the primary source and its coverage
of our universe is confirmed once the fetch completes and is reported in the results, not assumed here.

## 4. Populations and the exposure

- **Universe `U`.** The 6,371 QC-passing Stim48hr perturbations with a known `z_l2` (needed for the matched
  background). Every perturbation carries a P1 label (supported / not-supported / not-in-OT), P2, P3, P4.
- **Gate.** `safe` from `window_score.csv` (evidence floor + co-inhibitory preservation). 214 genes.
- **Naive ranking.** `-eff_mean_z`, the suppression ranking the whole project critiques.

## 5. Primary endpoint

Fixed in advance. Run once.

- **Test 1 — does the gate enrich for genetic support at all?** Among `U`, is P1-support more frequent in the
  safe set than in the evidence-passing-but-rejected set? Fisher exact, two-sided, and reported as an odds
  ratio with a 95% CI. This is descriptive of the gate.
- **Test 2 — the pre-registered primary. Gate vs naive, both top 100, against a magnitude-matched null.**
  Take the top 100 of the naive suppression ranking and the top 100 of the safe set by `window_score`. For
  each, count how many are P1-supported. Compare each observed count to the distribution from 5,000
  shortlists drawn from `U` matched on `z_l2` decile (the same machinery as `scripts/14`, `_matched_draws`,
  which must be fixed first so it does not silently truncate a top-100 draw against a thin pool — see §8).
  Report, for the naive top 100 and the safe top 100 separately: observed supported count, matched mean,
  and `P(matched >= observed)`.
- **The contrast that decides the direction:** `delta = (safe-set supported fraction) − (naive-top supported
  fraction)`, and its sign. A permutation CI on `delta` from the same 5,000 matched draws.
- **alpha = 0.05**, two-sided for the descriptive OR, one-sided is NOT used because the direction is the
  thing under test.

## 6. Decision rule — both branches pre-committed

- **ALIGN** iff the safe set is enriched for P1-support above its magnitude-matched null (`P < 0.05`) **and**
  `delta > 0` with a CI excluding 0. Report: the gate adds genetic relevance beyond effect size.
- **OPPOSE** iff `delta < 0` with a CI excluding 0, OR the safe set is *de*-enriched relative to naive.
  Report as the headline generalisation described in §2, with the immunodeficiency mechanism.
- **NULL** iff neither: the gate is orthogonal to genetic support. Report it plainly; the nomination then
  rests on the gate plus P2/P3 without a genetic-enrichment claim.

The nomination table is produced under every branch: the 214 safe genes, each carrying P1–P4, ranked by a
**transparent, pre-specified external-evidence order** (below), never by a fitted score.

## 7. The nomination ranking (external evidence, transparent, pre-specified)

The 214 safe genes are ordered by a lexicographic tier, fixed here, because a fitted composite would be
another un-validated score of the kind this project keeps retracting:

1. **Tier 1 — genetically supported AND LoF-tolerant AND tractable.** P1 `associationScore >= 0.1`;
   `prec <= 0.90` AND not `lof_intolerant` (safe to inhibit on both dominant and recessive axes); tractable
   in OT. These are the genes the literature says are most likely to succeed.
2. **Tier 2 — genetically supported, but carrying a LoF-constraint or tractability flag.** Reported with the
   flag named.
3. **Tier 3 — not genetically supported for these 14 diseases.** The exploratory novel candidates. Ranked
   within tier by `window_score` as the display default, and every one carries its selectivity annotation,
   its `prec`, and its Schmidt status, so a reader sees why to be cautious.

Within each tier, ties are broken by held-out Schmidt support, then `window_score`. No weights are fitted.
The recovered approved-drug targets (IMPDH2, PPP3R1, CD3E, IL4R, CD2, CD28) are labelled where they fall;
they validate the gate and are not counted as novel.

## 8. Falsification controls — all must fire, or the primary is void

1. **Matched-null calibration.** 500 random size-matched within-decile shortlists must give a P1-support
   count distribution whose 5% tail matches the nominal rate. If the matched-draw machinery truncates
   (the `size=min(count, pool.size)` bug in `scripts/14:_matched_draws`), the top-100 draw silently returns
   fewer than 100 rows; the fix is mandatory and the realised draw size must be asserted == 100 and logged.
2. **Positive control.** The authors' own regulator-burden core genes (Fig 6A), or any set known to be
   genetically loaded, must show enrichment where our machinery looks for it. If the pipeline cannot detect
   genetic loading where it certainly exists, it cannot report its absence.
3. **Label-shuffle.** Permuting the P1 label across `U` must destroy every enrichment. Report null mean ~0.
4. **Coverage is reported, never imputed.** Genes absent from OT (no `targetId` match, or no autoimmune
   association row at all) are a **third category**, not "not supported". The fraction of `U` and of the safe
   set that is not-in-OT is reported; a gene is never counted as unsupported because it was unmeasured by OT.

## 9. Known limitations, stated before the run

1. **Reverse causality in the ALIGN branch.** OT genetic-association scores partly derive from the same
   disease biology that motivates the screen; enrichment does not prove our gate found anything OT did not.
   The OPPOSE branch does not have this problem, which is one reason it is the stronger result.
2. **P1 is direct associations only** (`association_by_datatype_direct`), not ontology-propagated. A gene
   supported only for a child term of one of the 14 may be missed. Reported as a coverage caveat.
3. **The 0.1 threshold is the authors'**, adopted to match them, not tuned here. A sensitivity at 0.05 and
   0.2 is reported, never decisive.
4. **The screen's assay blindness carries through.** No polarising cytokines, so genetically supported
   cytokine-signalling targets (JAK2, TYK2, S1PR1) are ranked low by our directional score regardless of
   their genetics; this is a property of the ranking, documented in the source-paper briefing §1.4, and it
   will depress the ALIGN signal for a known reason.
5. **The nomination is hypothesis-generating.** No target here has a functional or in-vivo readout in this
   work. The gate is T-cell-intrinsic; organism-level safety is annotation, not proof (N3). One dataset, one
   analyst, one week. The deliverable is a ranked, evidence-annotated shortlist a wet lab could triage, not a
   validated drug target.

## 10. What is written where

`scripts/20_prior_knowledge_nomination.py`, new. Reads `results/tables/*.csv` and the on-disk OT parquets;
never touches the h5ad layers. Emits `results/tables/n10_nomination.csv` (214 genes, P1–P4, tier, rank),
`results/tables/n10_enrichment.csv` (the pre-registered tests and controls). It must `raise SystemExit(1)`
if any falsification control in §8 fails. `src/cd4_perturbseq/priors.py` gains an OT genetic-support loader
and an OT tractability loader. `scripts/14`'s `_matched_draws` truncation is fixed first, in its own commit.
The report gains a nomination section and the enrichment verdict.
