# Pre-registration — N9: is the co-inhibitory module preferentially reduced, or merely co-induced?

**Written 2026-07-08, before the analysis was run. Committed before `scripts/17_tolerance_is_special.py`
was written.** No number in this document was obtained from the test it specifies. The three quantities
quoted under "Background" were measured first and are what motivated the test; they are stated here so a
reader can see exactly how much was known in advance.

---

## 1. The question, in one sentence

Among perturbations, does a **suppression-ranked** shortlist reduce the nine-gene activation-induced
co-inhibitory module **more than it reduces a random gene module matched on both baseline expression and
stimulation-induction strength** — or is the apparent effect fully explained by the module being co-induced
with the effector program it was selected against?

## 2. What the module is, and what it is not

The nine genes are `FOXP3`, `IL10`, `IKZF2`, `CTLA4`, `PDCD1`, `LAG3`, `TIGIT`, `TGFB1`, `LRRC32`, measured
as mRNA log fold change in a 48 h anti-CD3/CD28 stimulation of **bulk primary human CD4 conventional
T cells**, knockdown versus control guide.

In that setting `CTLA4`, `PDCD1`, `LAG3` and `TIGIT` are activation-induced co-inhibitory receptors, present
on any activated Tconv. `FOXP3` is transiently induced in activated human Tconv **without conferring
suppressive function** — a human-specific phenomenon (Wang et al. 2007 *Eur J Immunol*; Tran, Ramsey &
Shevach 2007 *Blood*; Allan et al. 2007 *Int Immunol*); murine Tconv do not do this. `IKZF2` (Helios) and
`LRRC32` (GARP) are Treg-associated and also expressed on activated Tconv.

**Therefore this module is named `activation-induced co-inhibitory module`.** It is not a measurement of
tolerance, of suppressive function, or of regulatory T cell identity. No result below licenses a causal or
clinical claim, and none will be made. This limitation is larger than anything the test can settle, and it
is stated first for that reason.

## 3. Background: why the existing control is insufficient

`scripts/10_tolerance_is_real.py` already rejects one alternative and does it correctly. Its test 2 compares
the module against 200 random nine-gene modules **matched on baseline-expression decile**. In the naive
top-100 the module's z against that null is +3.87, against +1.26 in the background, p = 3.5e-13.

That null does not answer the obvious objection. The nine module genes are stimulation-induced, as are the
32 effector genes the ranking is built on. Random genes at the same *expression* level are mostly **not**
stimulation-induced, so they do not fall when TCR signalling is blocked. **Any co-induced gene set would
reproduce the result.** The objection was raised by the project's own adversarial audit, recorded as the
"tolerance tautology", and was recorded as refuted. It was refuted against an expression-matched null only.
It is not refuted.

Two further facts, measured before this document was written:

- `scripts/10`'s test 1 regresses module suppression on effector suppression and asks whether the residual is
  elevated in the naive top-100. `in_naive_top100` **is** the top-100 by effector suppression, so the residual
  is evaluated at the extreme of its own regressor, and the relationship is convex (median residual by
  efficacy decile: +0.021, +0.004, +0.005, −0.009, −0.021, −0.029, −0.025, −0.025, −0.013, +0.005 — positive
  at both ends). A linear fit through a convex curve produces a positive top-end residual for free.
  **It survives anyway**: restricted to the top efficacy decile, top-100 versus the other 537 gives module
  suppression +0.574 vs +0.283 (MWU p = 1.2e-10) and residual +0.117 vs −0.004 (p = 0.0055). Test 1 is sound.
  Its null is the problem, not its arithmetic.
- N7's magnitude-matched control (`scripts/14`) matches on `z_l2`, transcriptome-wide effect magnitude. It does
  not match on effector suppression. Our top 20 have the highest effector suppression in the pool by
  construction, and co-induction then predicts the highest module suppression. N7 therefore excludes
  *magnitude* and leaves *co-induction* untouched.

**What is already excluded, and is not re-litigated here:** effect magnitude. `rho(z_l2, module suppression)
= 0.069`. Perturbations that *induce* the effector program do not suppress the module; they induce it
(z = +9.27, `scripts/12`, test B). Global downregulation of expression-matched random genes (`scripts/10`).

## 4. Data

- Screen: `GWCD4i.DE_stats.h5ad`, `Stim48hr` arm, `qc_mask` applied → 6,371 perturbations. Layer `log_fc`
  for module suppression; layer `baseMean` for expression deciles.
- Induction: `priors.arce_stim_vs_rest()` → Arce 2024 DESeq2, AAVS1 control-guide Teff, Stimulation vs
  Resting. 16,087 genes with finite `log2FoldChange`. This is an **external, non-perturbational** definition
  of the activation program, which is why it is the right matching variable.

Measured before writing this document, and used only to establish feasibility:

| module | n present in Arce | median induction log2FC | median percentile |
|---|---|---|---|
| effector | 32/32 | +2.95 | 97.5% |
| co-inhibitory | 9/9 | +2.29 | 95.6% |

Candidate null pool: genes measured in the screen ∩ present in Arce ∩ not in either module. 344 measured
genes lie at or above the co-inhibitory module's median induction. Sufficient to draw nine-gene modules
matched on a joint decile grid.

## 5. Primary endpoint

Fixed in advance. Run once.

- **Population.** The 6,371 QC-passing `Stim48hr` perturbations.
- **Exposure.** `in_top100` = the 100 highest by `efficacy` = `−mean(log_fc)` over the 32 effector genes
  (leave-one-out on the perturbation's own gene). The remaining 6,271 are the background.
- **Null modules.** 200 random nine-gene sets. Gene *j* of each set is drawn from the pool cell matching
  module gene *j* on the **joint (baseline-expression decile × Arce induction-log2FC decile)**. The 32
  effector and 9 module genes are excluded from every pool. Cells that are empty after exclusion are
  recorded, and the fallback used is recorded; if more than 1 of the 9 positions requires a fallback the
  test is declared uninterpretable and reported as such.
- **Outcome.** For each perturbation *i*,
  `z_i = (module_suppression_i − mean_s null_suppression_{i,s}) / sd_s null_suppression_{i,s}`.
  The null mean and sd are computed **within perturbation *i***, which conditions on that perturbation's own
  effect without any regression.
- **Test.** One-sided Mann–Whitney U: `z[in_top100] > z[background]`. **α = 0.05.**
- **Co-primary, descriptive.** One-sample Wilcoxon of `z[in_top100]` against 0, and its median. If the
  module falls exactly as far as induction-matched modules, this median is ≈ 0.

**Why this test and not another.** The z is computed per perturbation against modules matched on the two
variables that could manufacture the effect: how highly the genes are expressed, and how strongly they are
switched on by the very stimulus the ranking blocks. If the effect survives, it is not co-induction. If it
does not, it is.

## 6. Decision rule — both branches are acceptable outcomes

**PASS** (primary p < 0.05 **and** median `z[in_top100]` > 0):
> The module is reduced beyond what expression- and induction-matched modules predict. The report may say
> "beyond co-induction" and must show this test beside `scripts/10`'s. It may **not** say "specifically
> destroys tolerance", "reveals", or any causal or clinical verb (RULE #8). It must still carry §2.

**FAIL** (either condition unmet):
> The effect is explained by co-induction. The report states this as a **result, not a caveat**: a
> suppression objective cannot separate the effector program from the co-inhibitory program because the two
> share an upstream cascade. Every occurrence of "specifically" is deleted from the report, README and video
> script. The word "tolerance" is replaced by "activation-induced co-inhibitory module" throughout.

**What does not change under either branch.** The practical finding is independent of N9: of the naive top
20, the gate refuses 18; nine of the 20 were independently shown to reduce IL-2 by Schmidt & Steinhart 2022;
six are refused on the co-inhibitory axis alone and would pass a homeostasis check. A triage layer scoring
only suppression is incomplete either way. **The FAIL branch is not a failure of the project.**

I commit in advance to reporting the FAIL branch as prominently as the PASS branch, and to leaving this
document in git unedited whichever fires.

## 7. Falsification controls — all must fire, or the result is void

1. **Recovery control.** Permute the induction values across genes before matching. The null degenerates to
   expression-matched-only and the primary must revert to approximately `scripts/10`'s answer (z ≈ +3.9 in the
   top-100). If it does not, the induction matching is not doing what this document claims.
2. **Positive control.** Substitute the 32-gene *effector* module as the "real" module, with its own
   expression- and induction-matched nulls. The top-100 is selected on effector suppression, so `z` must be
   strongly positive. If the pipeline cannot detect a preferentially-suppressed module where one certainly
   exists, it cannot be trusted to report its absence.
3. **Negative control.** Substitute a random nine-gene induction-matched module as the "real" module. Median
   `z[in_top100]` must be ≈ 0. Repeat 20 times; the nominal false-positive rate at α = 0.05 must be ≈ 5%.

## 8. Secondary, reported only, never decisive

- **The non-induced pair.** `TIGIT` (induction log2FC −0.04, 50th pct) and `TGFB1` (−0.54, 29th pct) are not
  stimulation-induced; the other seven are (+1.19 to +4.73). Co-induction predicts the top-100's suppression
  of the module is carried **entirely by the seven induced genes**, with no excess for `TIGIT` and `TGFB1`.
  Report per-gene suppression in top-100 versus background for all nine. **n = 2. This may corroborate the
  primary. It may never override it.**
- **The IL-2 axis, exploratory.** A candidate mechanism, if PASS fires: IL-2/STAT5 maintains `FOXP3` and Treg
  identity (Fontenot et al. 2005 *Nat Immunol*; Burchill et al. 2007 *J Immunol*; Yao et al. 2007 *Blood*) and
  does not maintain `IFNG`/`TNF`. Genes fixed a priori: `IL2`, `IL2RA`, `IL2RB`, `IL2RG`, `JAK1`, `JAK3`,
  `STAT5A`, `STAT5B`. Report `z` for perturbations of these versus other top-100 members. `JAK1` and `JAK3`
  are not perturbed in this library. n is small. **Descriptive. Hypothesis-generating. No p-value is
  promoted to a claim.**

## 9. Known limitations of this test, stated before it is run

1. §2. The module is a transcriptional proxy named by interpretation. No functional assay exists here. This
   dominates every statistical concern below.
2. Arce 2024 is bulk RNA-seq of Teff under a different stimulation protocol from the Perturb-seq atlas.
   Induction log2FC is therefore an imperfect covariate, measured in a different experiment. Matching on it
   reduces but does not eliminate the confound.
3. The joint decile grid is coarse. A residual induction gradient within a cell can survive matching.
4. `n = 9` module genes. The z has a nine-gene numerator; its variance is large and its null sd is estimated
   from 200 draws.
5. The exposure is the top-100 by effector suppression, a hard threshold, chosen to match `scripts/10` so the
   two tests are comparable. It is not tuned.
6. The assay has no polarising cytokines, so cytokine-signalling targets are near-invisible (`JAK2` rank 5392,
   `TYK2` 5618, `S1PR1` 5993, `IL4R` 6047 of 6371). This bounds what any conclusion here generalises to.
7. `qc_mask` conditions on `ontarget_significant`, a collider on target expression. The 6,371 are survivors.
   This test compares modules *within* perturbation, so the collider does not bias the z; it does bound the
   population the claim is about.

## 10. What will be written where

`scripts/17_tolerance_is_special.py`, new. It must `raise SystemExit(1)` if any falsification control in §7
fails. It writes `results/tables/tolerance_induction_matched.csv`. `scripts/10_tolerance_is_real.py` is **not
edited**; it is a load-bearing self-test and its expression-matched result stands as a separate, weaker claim.
The report gains a subsection reporting both nulls side by side, and the verdict from §6.

---

# ADDENDUM 1 — 2026-07-08, after the first run. The negative control fired.

**The original text above is unedited and stands. This addendum records that the pre-registered
falsification control §7.3 failed, what it caught, and the machinery fix. The decision rule of §6 is
NOT changed. Only the null distribution of the primary's p-value is corrected.**

## What happened

First run of `scripts/17_tolerance_is_special.py`, 200 nulls, 20 negative-control modules:

| quantity | value |
|---|---|
| primary: `z` top-100 median, induction-matched | **+2.600** |
| primary: `z` background median | +0.916 |
| primary: one-sided MWU p | 5.57e-11 |
| control §7.1 recovery, permuted induction: `z` top-100 median | +3.745 (reproduces `scripts/10`'s +3.87) |
| control §7.2 positive, effector module: `z` top-100 median | +10.093 |
| control §7.3 negative: median of 20 random modules' top-100 median `z` | +0.091 |
| **control §7.3 negative: false-positive rate at α = 0.05** | **25.0% (5 of 20)** |

§7.3 requires that rate to be "≈ 5%". It is 25%. The control **failed**.

(The script's coded gate was `fpr <= 0.25`, which is looser than the criterion registered in §7.3 and
wrongly returned PASS. The registered criterion binds, not the code. The gate has been tightened to
`fpr <= 0.10` so the code and this document agree. That the code was more permissive than its own
pre-registration is itself the kind of defect this project exists to catch, and it is recorded here
rather than quietly fixed.)

## What the control caught

The primary applies a Mann–Whitney U across 6,371 perturbations. Those 6,371 `z` values are **not
independent**. The real module is a *single fixed nine-gene set*. If that set happens to covary with
effector suppression, every top-100 `z` shifts together. **The effective unit of resampling is the
module, not the perturbation**, so the effective n is 1, not 100. The MWU p-value is anticonservative
by orders of magnitude. A 25% false-positive rate on random modules is the direct measurement of that.

This is not a defect of the induction matching. It is a defect of the inference placed on top of it,
and it was present in the same form in `scripts/10_tolerance_is_real.py`'s test 2, whose p = 3.5e-13
is anticonservative for exactly the same reason. **`scripts/10`'s verdict is not retracted — its
statistic is still the right statistic — but its p-value should not be quoted.**

## The machinery fix (does not change §6)

Same statistic. Valid null.

- **Statistic:** `T = median(z[in_top100])`, computed exactly as in §5.
- **Null:** draw `N = 200` random nine-gene "real" modules matched on the same joint
  (expression decile × induction decile) grid, each with its own 200 matched null sets, and compute
  `T` for each. Exclude each drawn module from its own nulls.
- **p:** `(1 + #{T_null >= T_real}) / (1 + N)`. This is a module-level permutation p. It is exact by
  construction and needs no independence assumption across perturbations.
- **Calibration check, new:** leave-one-out among the 200 null modules must give a uniform p, i.e.
  ≈5% of them below 0.05. This replaces §7.3 and is reported.
- **α, decision rule, both branches, and every limitation in §9: unchanged.**

The per-perturbation MWU is still computed and reported, **explicitly labelled void and
anticonservative**, so the reader can see the size of the error.

## What was already visible in the first run, and is not changed by the fix

- Induction matching reduces the top-100 `z` from **+3.745** (expression-matched, the permuted
  recovery control) to **+2.600**. Co-induction therefore accounts for roughly **31%** of the effect
  that `scripts/10` reported. It is a large part of it. It is not all of it.
- The positive control behaves: the effector module, against its own expression- and
  induction-matched nulls, gives `z` = +10.093 in the top-100. The pipeline can detect a
  preferentially-suppressed module where one certainly exists.
- **The per-gene split (§8) does not follow the co-induction prediction.** Co-induction predicts the
  two non-induced genes show no excess suppression in the naive top-100. Observed:

  | gene | induction log2FC | top-100 median | background median | MWU p |
  |---|---|---|---|---|
  | IL10 | +4.73 | 2.366 | 0.498 | 1.8e-16 |
  | LRRC32 | +3.96 | 0.020 | 0.021 | **0.43** |
  | FOXP3 | +3.04 | 0.387 | 0.040 | 3.4e-07 |
  | LAG3 | +2.78 | 0.200 | 0.007 | 9.6e-08 |
  | CTLA4 | +2.29 | 0.360 | 0.025 | 1.5e-23 |
  | PDCD1 | +1.58 | 0.562 | 0.017 | 2.4e-17 |
  | IKZF2 | +1.19 | 0.386 | 0.161 | 5.5e-03 |
  | **TIGIT** | **−0.04** | **0.930** | 0.164 | **8.1e-10** |
  | **TGFB1** | **−0.54** | 0.014 | −0.005 | 0.13 |

  `TIGIT` is not stimulation-induced and is nonetheless the second most strongly suppressed gene in
  the naive top-100. `LRRC32` is strongly induced and shows no excess at all. So the pattern is
  **1 of 2 non-induced genes suppressed, and 1 of 7 induced genes not suppressed.** These per-gene
  p-values carry the same non-independence caveat as the void primary and are descriptive only.
  n = 2 on the non-induced side. It corroborates. It cannot settle anything. Recorded because it was
  pre-registered, and because it points the other way from the hypothesis under test.

- IL-2 axis: only 4 of the 8 pre-specified genes are perturbed and QC-passing (`IL2RA`, `IL2RB`,
  `STAT5A`, `STAT5B`), and only 2 reach the naive top-100. **Reported as unavailable, not as null.**
  The exploratory mechanism in §8 cannot be assessed on this screen.
