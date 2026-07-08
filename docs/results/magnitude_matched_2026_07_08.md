# Is the risk-kill result an effect-magnitude artifact?

**2026-07-08. Tasks N1 and N2.** Scripts: `12_magnitude_matched.py`, `12a_validate_stratified.py`,
`02_risk_kill_reversal.py`. All exit 0.

**Short answer.** No. But the mechanism we asserted was wrong, the audit was half right, and the
headline should now be told through *tolerance collapse* rather than *collateral damage*.

---

## What the audit claimed

From `docs/audits/adversarial_audit_2026_07_08_findings.json`, filed critical:

> The 'collateral DE' and 'Rest DE' pillars are an effect-magnitude artifact, and the direction is
> backwards: at matched |score|, suppressors have ~66% FEWER stim DE genes than inducers.

with the supporting observation that `Spearman(n_cells_target, stim_de_genes) = -0.243`, so the
cell-count-matched background controls a confound running the opposite way. Its proposed fix:

> Replace the cell-count-matched background with an effect-magnitude-matched one (match on
> `|naive_suppression|` decile, **or** on `stim_de_genes` decile, **or** on a transcriptome-wide
> `||z||_2`) ...

The rival story would be fatal. A perturbation that moves the whole transcriptome gets a large
`|mean z|` over *any* gene set, the effector module included, and crosses the DE threshold on many
genes for the same reason. If that is all that is happening, the naive ranking selects "big effect",
not "suppresses inflammation", and every pillar is a restatement of effect size.

## Reproduce before believing (RULE #1)

The audit has been right once, wrong once, and incomplete once on this codebase.

| Claim | Result |
| --- | --- |
| `Spearman(n_cells_target, stim_de_genes) = -0.243` | **Reproduces exactly.** ρ = −0.243, p = 5.2e-86, n = 6,371 |
| "~66% fewer stim DE genes at matched \|score\|" | **Direction reproduces. Magnitude does not.** |

The direction is robust: van Elteren over all 6,371 rows, stratified on `|score|`, gives z = −7.22
(10 bins), −7.43 (20), −7.35 (50), −7.26 (100), and −7.48 restricted to the `|score|` range where
inducers actually exist. Suppressors really do carry fewer collateral DE genes at matched `|score|`.

The 66% does not. Within-decile median ratios of suppressor to inducer DE count:

| \|score\| decile | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ratio | 0.667 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.500 | 0.750 | **0.333** | 0.111 |

Pooled, the median ratio is **0.875** — 12% fewer, not 66%. The 66% is decile 8 alone. Deciles 1
through 5 give exactly 1.000. The audit reported one decile as though it were the pooled effect.

## Where the audit's conclusion comes apart

The fix text offers three controls as interchangeable: match on `|naive_suppression|` decile, on
`stim_de_genes` decile, or on transcriptome-wide `||z||_2`. They are not the same control.

> **`Spearman(z_l2, |naive_suppression|) = +0.198`.**

Score magnitude and effect magnitude share 4% of their rank variance. The audit gathered its
evidence by matching on `|score|` and stated its conclusion about *effect magnitude*. Both are
legitimate conditionings; they answer different questions, and the pillars' fate depends entirely on
which one you choose.

`z_l2` here is `sqrt(Σ z²)` over the 10,282 measured genes, excluding the perturbed gene's own
column (its knockdown z is large by construction, because `ontarget_significant` is a QC
requirement), the 32 effector genes the score is built from, and the 9 tolerance genes. So the
matching variable cannot contain the score.

Same rows, same test, only the stratifier changes:

| stratified on | collateral DE | rest DE | tolerance |
| --- | --- | --- | --- |
| `\|score\|` decile | z = −7.22 | z = −4.79 | z = +22.14 |
| `z_l2` decile | z = **−2.82** | z = **−1.22 (p = 0.22)** | z = +26.82 |

Match on effect magnitude and the collateral effect shrinks 2.6×; the resting one vanishes.

## The three tests

None draws a random number. Machinery in `src/cd4_perturbseq/stratified.py`, validated 13/13 by
`12a_validate_stratified.py` against scipy, against a 20,000-draw within-stratum permutation null,
and against 400 simulated nulls (type-I 0.037 van Elteren, 0.022 CMH).

- **A — magnitude.** Top 100 vs *every* background row, stratified on `z_l2` decile.
- **B — direction at matched magnitude.** Top 100 by suppression vs top 100 by *induction*,
  stratified on `z_l2` decile. A pillar that fires on both tails measures `|effect|`, not suppression.
- **C — direction over the whole ranking.** All 4,561 suppressors vs all 1,810 inducers, stratified on
  `|score|` decile. Exactly the test the audit's fix text asks for.

The kill rule was fixed in `SESSION_SUMMARY.md` before any number existed: **a pillar that fails both
A and B is removed from the gate.** It binds the four registered pillars only.

| Pillar | ρ with `z_l2` | A (vs magnitude-matched bg) | B (vs induction control) | C (at matched \|score\|) | Verdict |
| --- | --- | --- | --- | --- | --- |
| IEI (immunodeficiency) | — | OR 2.77, p_BH 0.0017 | OR 20.6, p_BH 0.0078 | ns | survives (annotation only) |
| collateral DE in stim | **+0.725** | z +7.04, p_BH 2.5e-12 | z +1.98, p_BH 0.040 — **fragile** | suppressors **LOWER** | survives, but entangled |
| DE at rest (raw) | +0.527 | z +3.04, p_BH 0.0017 | z −0.60, **ns** | suppressors **LOWER** | survives on A alone |
| tolerance suppression | **+0.069** | z +11.10, p_BH 3.2e-28 | z +9.27, p_BH 4.8e-20 | suppressors **HIGHER** | **survives cleanly** |

"Fragile" means the call changes with the bin count: collateral DE passes test B at 10 bins
(z = +1.98) and fails at 20 (z = +0.59, p = 0.28). A verdict that depends on an arbitrary bin count
is not a verdict. Tolerance is significant at 5, 10 and 20 bins without wavering.

**Falsification control.** 200 random shortlists of 100 reach p < 0.05 on test A at 2.5%, 7.5%, 3.0%
and 6.0% for the four pillars — the nominal 5%. The tests are not vacuous. This rules out the *noise*
null, not the *structured* nulls (rank by `z_l2` alone, by DE count alone); that is task N7.

## What this means

1. **The naive ranking really does select big-effect perturbations.** Collateral DE is 0.725
   rank-correlated with `z_l2` and is close to a restatement of effect size. Our stated mechanism —
   "reversal nominates genes that cause collateral damage" — is largely a magnitude effect.
2. **"Reversal is not enough" survives, and it survives on tolerance.** Tolerance-module suppression
   is ρ = 0.069 with effect magnitude. Inducers do not suppress tolerance; they *induce* it. It is the
   one pillar about the **direction** of a perturbation rather than its **size**, and it is the axis
   the audit attacked hardest and lost.
3. **RULE #2 again, one level deeper.** Database-membership axes all failed. Within the screen-native
   axes, the *counts of DE genes* are magnitude proxies and the *graded module score* is not. The
   contribution is a decision layer built from graded screen-native measures.

## Two schema corrections, both free

- `n_total_de_genes` counts the perturbed gene itself. The column named "collateral DE" was not
  measuring collateral DE. `.obs` carries `n_downstream`. Script 12 now uses it.
- The Rest arm is not comparable to the Stim arm. Rest rows were never QC-filtered (802 of 6,371 have
  a resting row that would fail the QC the stimulated row had to pass), and 51 perturbations have no
  resting row at all — including **IL2RB at rank 4**, the single gene `risk_kill_2026_07_08.md` uses
  to motivate the entire safety gate. A missing resting row is not a quiet zero.

## N2: the seed lottery, reproduced and retired

The old background was one 100-row sample drawn with `default_rng(0)`, discarding 6,171 of 6,271
controls. Four audit agents reported the defect with four different resample counts and four
non-significance rates (11.5%, 11.2%, 11%, 9.5%). None committed a script.

Reproduced against the code as it actually stood, 2,000 redraws:

- IEI odds ratio: seed 0 gives 3.59; median 3.59; 2.5–97.5 pct **[1.72, 14.79]**
- IEI p-value: seed 0 gives 0.0199; 2.5–97.5 pct **[0.00064, 0.178]**
- **seeds giving p ≥ 0.05: 22.4%**

The defect is twice as bad as the audit reported. Inference now runs Cochran-Mantel-Haenszel (flags)
and van Elteren (measures) over every control row, with no random numbers anywhere.

`holds = any(matched_verdicts)` is gone. Five one-sided tests at α = 0.05 carry a ~23% family-wise
false-positive rate, so `any()` made a script that advertised itself as "designed to be able to FAIL"
unfalsifiable. It is replaced by one pre-registered primary endpoint:

> **Tolerance-module suppression is higher in the naive top 100 than in the `z_l2`-decile-stratified
> background.** One-sided van Elteren, α = 0.05. **z = +11.099, p = 6.3e-29.**

Chosen structurally, not by p-value: it is the only pillar that is both a gate axis and a graded
screen-native measure, and it was validated against 200 expression-matched random modules at
`1adab65`, before the endpoint existed. One of its three justifications was known at selection time,
which the script's docstring states out loud; it is pre-registered with respect to future runs and
new data, not with respect to the run that motivated it. The other four tests are secondary and
Benjamini-Hochberg corrected, and decide nothing.

Verified falsifiable: permuting the tolerance column, and sign-flipping it, each drive the script to
exit 1.

## What this does not settle

- **N6.** `04_window_score.py` does not gate on raw rest-DE. It gates on `selectivity =
  log1p(stim_de) − log1p(rest_de)`. Conditional on effect magnitude the naive top 100 shows **no
  enrichment** for poor selectivity (test A p = 0.495; test B p = 0.908). That does not make
  selectivity a bad *filter* — NSD1, with 2,175 stimulated and 1,767 resting DE genes, really should
  be rejected — but the risk-kill argument cannot lean on homeostasis, because it never tested the
  axis it ships. A pillar test asks "is the top *enriched*"; a gate asks "*reject* this target". This
  analysis was not pre-registered, so it sets the next task and does not rewrite the gate.
  Separately, `04` median-imputes missing `rest_de_genes` straight into `selectivity`, silently.
- **N7.** The audit's `invalidates_headline: true` finding — "non-reversal null rankings pass the
  paper's own verdict function" — is only half answered. Test B is a reversal-specificity test and
  tolerance passes it. The structured nulls (rank by `z_l2`, by `stim_downstream`, by `|score|`) have
  not been run.

## Reproduce

```bash
uv run python scripts/12a_validate_stratified.py   # 13/13, the machinery is sound
uv run python scripts/12_magnitude_matched.py      # N1: the three test families
uv run python scripts/02_risk_kill_reversal.py     # N2: seed lottery, then the primary endpoint
```

All three exit 0. All three are written so they can exit 1.
