<!-- Produced 2026-07-08 by a five-lens adversarial panel (31 agents) that attacked an earlier N6
     design written by Claude Code. 25 objections raised; 17 CONFIRMED and 6 PARTIAL against the
     committed tables. The earlier design is WITHDRAWN IN FULL: its primary endpoint
     (rest_cells_ratio) was disqualified four independent ways, its `rest_not_expressed` class was
     shown not to exist, and its `indeterminate` rule was a scope error that would have voided 17
     tolerance rejections computed entirely from the stimulated arm. This document supersedes it.
     It is a PRE-SPECIFICATION WITH DISCLOSED UNBLINDING, not a blind pre-registration. See L1. -->

# N6 Pre-Registration: Validation of the Homeostasis (Selectivity) Axis

**Status:** binding. Commit this file before `scripts/15_n6_selectivity_validation.py` is written or run.
**Date:** 2026-07-08
**Supersedes:** the draft N6 design circulated to the adversarial panel. That draft is withdrawn in full. Its primary endpoint (`rest_cells_ratio` under z_l2 deciles) is disqualified. Its `indeterminate` disposition survives only in scoped form. Its `rest_not_expressed` PASS carve-out is deleted.

**This is a pre-specification with disclosed unblinding, not a blind pre-registration.** Section 7 lists every statistic already computed. The project's standing rule ("a measure computed mid-analysis must not decide an exit code") is honoured by naming, in Section 4, the one statistic that decides the exit code and that has not been computed.

---

## 1. The question

Among the perturbations the tolerance axis already keeps, does the shipped selectivity composite `log1p(stim_de_genes) - log1p(rest_de_genes) < ln(10)` carry information about harm beyond the raw resting DE burden that it reduces to inside a z_l2 stratum?

N6 does **not** re-open reversal specificity. N7 settled that: homeostasis rejects 14 of our top 20 versus 13.6 matched, P = 0.53. N6 asks only whether homeostasis is a valid **filter**.

---

## 2. FIX A: the data-provenance repair, and what survives of it

### 2.1 What is repaired (arithmetic, uncontested)

`scripts/04_window_score.py:182` imputes the screen-wide median `rest_de_genes = 3.0` into rows with no resting arm. Among evidence-passers with a measured resting row the median `rest_de_genes` is 157. All 5 `rest_absent` evidence-passers pass homeostasis on the imputed number (5/5 versus 73/258; Fisher one-sided p = 0.0021). CENPE, a mitotic kinesin, passes the entire gate on it, window_rank 131.

**FIX A1.** Delete the imputation. `selectivity` becomes `NaN` when no resting row exists. `NaN < ln(10)` evaluates `False`, which is why the 5 currently pass; the three-valued verdict in 2.3 removes that path.

### 2.2 The label is renamed, because it is not provenance

`rest_provenance == rest_measured` is a QC filter whose central term, `ontarget_significant` at rest, is a descendant of both the exposure and the outcome. Verified within z_l2 deciles over the 6,320 rows with a rest row:

* Outcome arm: `cochran_mantel_haenszel(rest_qc_pass, low_rest_cell_tertile, z_l2_decile, "less")` gives MH OR = 0.5821, z = -6.670, p = 1.28e-11.
* Exposure arm: `cochran_mantel_haenszel(rest_qc_pass, fail_homeostasis, z_l2_decile, "greater")` gives MH OR = 2.470, z = 5.435, p = 2.73e-08.

The filter depends on the exposure and on the outcome, within the stratifier, in opposite directions. Rename `rest_provenance` to `rest_qc_class`. Never call it provenance in the report.

### 2.3 Final class definition

Rebuilt from `de_stats.read_obs()` only, Rest rows deduplicated on `target_contrast_gene_name`:

```
rest_absent   : no Rest row                                              n =    51
rest_qc_fail  : has Rest row, fails distal | neighbouring | ~ontarget_significant | low_target_gex
                                                                          n =   802
rest_qc_pass  : has Rest row, passes all four                             n = 5,518
                                                                    total  = 6,371
```

The old `737 / 65` split of `rest_qc_fail` into `rest_unmeasured` and `rest_not_expressed` is **retired**. Two independent rebuilds disagree at one boundary row (736/66 versus 737/65), and `low_target_gex` mechanically entails `ontarget_significant == False`: P(ontarget_significant | low_target_gex) = 0.0115 across all 11,287 Rest rows, against 0.795 when the target is expressed. The two classes cannot be separated by the flags. D3's own wording, "fails rest QC only on low_target_gex", selects 4 rows, not 65.

**The `rest_not_expressed` PASS carve-out is deleted.** "Target absent at rest = the window itself" is an unvalidated biological inference from a scRNA-seq detection threshold, which is the exact species of assertion that got three prior mechanisms refuted. It is also inert: 1 of 286 evidence-passers (FN3K, already rejected by tolerance), 0 of the 66 safe genes, and `safe = 60` under every admission rule.

### 2.4 Does `indeterminate` survive the panel? Yes, scoped. No, as drafted.

The drafted rule ("such a row may not be called safe, and is not counted as rejected either") voids 21 rejections among the 286 evidence-passers. Seventeen of those cite **tolerance**, an axis computed entirely from the stimulated arm that never reads the resting data. That is a scope error, not a repair.

**Scoped rule, binding:**

* `homeostasis_verdict ∈ {pass, fail, indeterminate}`. It is `indeterminate` iff `rest_qc_class != rest_qc_pass`.
* `fail_evidence` and `fail_tolerance` are computed from the stimulated arm and **stand unchanged** for every row.
* `gate_outcome = rejected` if the row fails evidence or tolerance, regardless of `homeostasis_verdict`.
* `gate_outcome = safe` iff evidence passes, tolerance passes, and `homeostasis_verdict == pass`.
* `gate_outcome = indeterminate` iff evidence and tolerance both pass and `homeostasis_verdict == indeterminate`.

Only 4 rows lose a rejection under this rule, not 21. The evidence-passer pool for matched-null draws becomes 282, not 259.

### 2.5 Exactly what FIX A does

**Safe count: 66 → 60.** Removed: CENPE (`rest_absent`); EGR2, PGAP2, LETM2 (`ontarget_significant == False` at rest); HEXIM1 (`distal_offtarget_flag` at rest); QRSL1 (`neighboring_gene_KD` at rest). Six out, zero in. The imputed value 3.0 is far below the measured evidence-passer median of 173, so voiding the imputation can only raise the rejection count. FIX A is monotone by construction, and monotone toward the smaller shortlist.

Pre-declared sensitivities, reported, never decisive: `safe = 65` under FIX A1 alone; `safe = 62` if HEXIM1 and QRSL1 are admitted (Section 7, limitation L9); `safe = 60` under every rule for the retired `rest_not_expressed` class.

**The 18/20 headline: unchanged.** No naive top-20 gene is rejected only on homeostasis with an indeterminate resting arm. Zero. The three top-20 genes with an indeterminate resting arm (IL2RB, CD3G, TADA2B) are all rejected on tolerance, which stands. Re-running `scripts/14_reversal_specificity.py` on the 282-row pool, 5,000 draws, seed 0, ALPHA = 0.01:

| line | current | under scoped FIX A |
|---|---|---|
| rejected in top 20 | 18 vs 15.30 matched, P = 0.1042 | 18 vs 15.31, P = 0.1064 |
| homeostasis in top 20 | 14 vs 13.62, P = 0.5296 | 13 vs 13.43, P = 0.6912 |
| tolerance in top 20 | 12 vs 3.99, P = 0.0000 | 12 vs 4.02, P = 0.0000 |

The homeostasis count moves 14 → 13 because TADA2B's homeostasis verdict becomes indeterminate. The one surviving specificity claim, tolerance, survives at ALPHA = 0.01. Under the drafted (unscoped) rule the headline would have become 15/17 with tolerance 9 vs 3.02, P = 0.0002. That version is rejected.

### 2.6 FIX A is inert on the FAIL branch

If selectivity is demoted, the gate is `evidence + tolerance` and never reads the resting arm. `safe = 214`. D1, D2 and D3 acquire no deliverable consequence. **FIX A's 66 → 60 exists only on the PASS branch.** It is committed anyway, because it is the correct handling of the data and because the PASS branch requires it.

### 2.7 New columns written by 04

```
rest_qc_class            {rest_qc_pass, rest_qc_fail, rest_absent}
rest_qc_fail_reason      comma-joined subset of {distal_offtarget, neighboring_kd, ontarget_ns, low_target_gex}
rest_de_is_ontarget_only bool; rest_de_genes == 1 and ontarget_significant at rest.  1,590 of 5,518 = 28.8%
homeostasis_verdict      {pass, fail, indeterminate}
gate_outcome             {safe, rejected, indeterminate}
selectivity              NaN when rest_qc_class == rest_absent
```

`rest_de_is_ontarget_only` exists so that D4's `+1` is visible wherever selectivity is printed.

---

## 3. Primary endpoint

### 3.1 The incumbent is dropped. Four independent reasons.

`rest_cells_ratio` is unusable as an outcome for this exposure.

1. **It is not a ratio.** `04_window_score.py:196-197` divides `rest_cells` by the single scalar 543.0. `Spearman(rest_cells, rest_cells_ratio) = 1.000` exactly. van Elteren is a rank test, so the raw count and the "ratio" give byte-identical p. "Relative to the median" normalises nothing.
2. **It is the sample size of the exposure's denominator.** `rest_cells == obs["n_cells_target"]` of the Rest contrast for all 6,320 rows, exactly. That contrast emits `rest_de_genes`. By D4, verified at rest with 0 mismatches, `rest_de = rest_downstream + ontarget_significant`. Within a z_l2 stratum the group label **is** a threshold on `rest_de_genes`: mean within-stratum AUC = 0.996, against 0.463 for `stim_de_genes` and 0.544 for `z_l2`. Mean within-stratum variance of `log1p(stim_de)` is 0.053 against 3.768 for `log1p(rest_de)`, a ratio of 71. The endpoint is the power denominator of the exposure's only surviving input.
3. **The relation is U-shaped, so the one-sided monotone alternative is misspecified.** Among the 413 activity-matched evidence failers, `frac_rejected` by rest-cell sextile runs 0.739, 0.657, 0.544, 0.580, 0.706, 0.739. Logit quadratic term p = 5.0e-4 against a null linear term (p = 0.625). The design's own test scores that population at p = 0.496 while a permutation-calibrated both-tail dispersion test finds p = 0.0065.
4. **It is guide abundance.** `Spearman(rest_cells_ratio, stim n_cells_target) = +0.931`. The arm-to-arm coupling is +0.9353 among 2,988 null perturbations that cannot be toxic, so it is a pre-exposure common cause, not shared toxicity. Under the design's own strata, the negative-control exposure `1[n_cells_target < median]` scores z = -11.94, p = 3.8e-33. The design's machinery detects library balance at 12 sigma.

Consequently the exit code was not a property of the data. It moved with the stratifier (z_l2 -1.363, stim_de -1.021, n_cells_target -0.079, rest_de -2.407), with the bin count (p = 0.056 to 0.207 over 4 to 20 bins), and with the disposition of the filter's casualties (bracket [0.0231, 0.3081]). And a sham exposure built from `rest_de_genes` alone, with no stim term, no log1p, no ln(10) bar, beat it: z = -2.220 against -1.363.

### 3.2 No endpoint on this screen validates the selectivity composite. Say it plainly.

Two structural facts make it impossible:

* Inside a z_l2 stratum, `rejected-by-selectivity` is `high rest_de_genes` (AUC 0.996). `Spearman(z_l2, stim_de_genes) = +0.990` on the primary population, so stratifying holds the numerator fixed and 4.0% of `log1p(stim_de)` variance survives. Nothing computed on this population can distinguish the composite from a raw resting DE threshold.
* Every internal outcome is a function of `rest_cells` or `rest_de`. There is no resting-arm CD4 harm measurement in these tables that is not one of them. `viability_tier` and `log2_stim_over_rest_cells` are both derived from `rest_cells`.

**What N6 can therefore test, and all it can test, is the incremental question:** does the composite beat its own `rest_de` shadow on an outcome measured outside this screen? That is the primary below. Its expected outcome, disclosed in advance, is FAIL. It is committed anyway, because the alternative is to demote a shipped gate without a rule.

### 3.3 The primary, in full

* **Population `P204`.** `~fail_evidence & ~fail_tolerance & (rest_qc_class == "rest_qc_pass")`. n = 204: 144 homeostasis-rejected, 60 kept. `loeuf` is missing for 10 rows, which are dropped and counted; the test uses 194 rows, 136 rejected versus 58 kept.
  * Tolerance survivors only, because the decision that hangs on this axis is whether 148 genes enter the shortlist. That is the marginal contrast. The 258-row population mixes tolerance-failing rows into both arms (41/185 rejected, 13/73 kept), and tolerance failers have roughly half the resting cell count.
  * `rest_qc_pass` only, because a row FIX A calls indeterminate cannot sit in a group labelled "kept-by-selectivity".
* **Stratifier.** `deciles(z_l2, n_bins=10)`. Quantile bins (`pd.qcut`), cut on `P204`, not on the screen. Equal-width binning is out of spec.
* **Groups.** The shipped `fail_homeostasis` column, read directly off `results/tables/window_score.csv`. Verified: `selectivity == log1p(stim_de_genes) - log1p(rest_de_genes)` to 1.3e-15, and `fail_homeostasis == (selectivity < ln 10)` on all 6,371 rows. Reading the label off the shipped table closes D5 for this test.
* **Outcome.** `loeuf` from `results/tables/window_score_organism_safety.csv`. gnomAD LoF observed/expected upper bound. External to the screen. It shares no term with any DE count and cannot be caused by this screen's cell counts. `Spearman(loeuf, z_l2) = -0.047`; `Spearman(loeuf, stim_de_genes) = -0.034`.
* **Test.** `van_elteren(loeuf, in_top=fail_homeostasis, strata, alternative="less")`. Rejected genes should be **more** constrained, i.e. lower LOEUF.
* **Alpha.** 0.05, one-sided.
* **Sham comparator, the statistic that decides.** Within each z_l2 decile, relabel the top `k_s` rows by `rest_de_genes` alone as "rejected", where `k_s` is that stratum's observed number of rejected rows. Same outcome, same strata, same test. Call the deviate `z_sham`. This exposure contains no stim term, no log1p difference, and no ln(10) bar.
* **Bin stability.** The primary is recomputed at `n_bins ∈ {4, 5, 10, 20}`.

### 3.4 Why LOEUF and not something better

There is nothing better. LOEUF is external, it is on disk, it can fail, and it does not validate whatever axis you hand it (control C2). It measures organism-level LoF intolerance, not resting-CD4 damage. That gap is a permanent limitation (L13), not a defect this design can close.

---

## 4. Decision rule

Both branches are acceptable outcomes. The gate exists to be right, not to be kept.

**RETAIN** the homeostasis axis iff **all four** hold:

1. `p_sel < 0.05` at `n_bins = 10` on `P204`.
2. `p_sel < 0.05` at every `n_bins ∈ {4, 5, 10, 20}`.
3. **Dominance:** `z_sel <= z_sham - 0.5`. The margin is arbitrary and is fixed here, before the statistic exists.
4. Controls C2, C3, C4, C5 all fire as specified in Section 5.

**DEMOTE** otherwise.

### 4.1 What a RETAIN licenses, and what it does not

A RETAIN licenses exactly this sentence and no other: *the resting-arm DE burden, in the selectivity form, carries organism-level LoF-constraint information beyond a raw threshold on that burden.* It does not validate the ln(10) location, the log1p offset, or D6's implemented `(1+stim)/(1+rest) >= 10` rule. D6 remains an open documented defect on a RETAIN. The report must carry, permanently and non-optionally:

* Homeostasis has no reversal specificity. N7: 14 of ours versus 13.62 matched, P = 0.5296; 13 versus 13.43, P = 0.6912 under scoped FIX A. It is a filter, not evidence about reversal.
* On `is_iei`, the annotation closest to a CD4 therapeutic context, the axis is anti-selective: CMH one-sided less, OR = 0.202, z = -2.091, p = 0.0183, 4/148 rejected versus 6/66 kept. The axis preferentially **keeps** immunodeficiency genes.
* The NSD1 sentence is struck from the report unconditionally, on either branch. NSD1's `reject_reason` is `"homeostasis,tolerance"`, `window_rank` 3380. It is deleted by tolerance regardless of selectivity's fate and is worth nothing as evidence for the axis. If a motivating example is wanted, take one of the 148 genes homeostasis uniquely rejects.

### 4.2 What a DEMOTE costs, and where it lands

Selectivity becomes a reported per-gene annotation. The gate becomes `evidence floor + tolerance`. `homeostasis_verdict` is dropped from `gate_outcome`. FIX A becomes inert (2.6).

* **Shortlist.** `safe = 214` of 286 evidence-passers. All 66 current safe genes are retained. Net +148. The full 214 is the shortlist, ranked by `window_score`. **No cap is imposed and no post-hoc trimming is permitted.** The report displays the top 20 as it does now.
* **This is not a free simplification.** Homeostasis is near-orthogonal to tolerance. It uniquely rejects 148 of the 214 tolerance survivors (69.2%), and those 148 are, if anything, *lower* on `tolerance_loss` than the 66 it keeps: `van_elteren(tolerance_loss, fail_homeostasis, deciles(z_l2,10), "greater")` gives z = -0.895, p = 0.815; medians 0.172 versus 0.220.
* **The 148 are independently constrained.** LOEUF median 0.631 versus 0.802; `lof_intolerant` 45.3% versus 30.3%. Demotion admits 148 genes the LoF-constraint annotations dislike. That cost is reported in the abstract, not buried.
* **Held-out Schmidt screen.** Report precision at 20 for the 66-gene gate and the 214-gene gate side by side. **No threshold on that comparison is set here and none may be set later.** Setting one after seeing it would be tuning. It is descriptive. No gate change may be made on its basis.
* **Artifacts rebuilt on 214:** Figures 1 and 4, the head-to-head table in 04, the "Tier A top 20" Schmidt validation, N8's benchmark ceiling, N7's axis decomposition.

### 4.3 Honest statement about the expected outcome

On the 214-row superset the dominance contrast has already been computed: `z_sel = -2.6255`, `z_sham = -3.588`, so `Δz = +0.96`. Selectivity is worse than its own shadow by roughly one z. The 204-row value is unrun and is not expected to differ in sign. **The pre-registered outcome is DEMOTE.** We commit to the rule rather than declaring the demotion by fiat, and we do not pretend the rule was blind.

---

## 5. Falsification controls

The drafted controls are deleted. Neither is a control.

* **Within-stratum permutation of the group labels.** It reproduces the analytic p-value to Monte Carlo error and destroys causation, confounding and collider selection identically. Fed a tautology (`label = outcome < median`) it returns z = -12.5 and kills it. Fed a pure coin flip it returns the coin flip's own p. Demoted to a **software calibration check**: report null mean ≈ 0, sd ≈ 1. It may not be cited as evidence that a result is non-artefactual.
* **Sign flip.** `van_elteren(y, ~g, s) = -van_elteren(y, g, s)` exactly, to 1e-9. An algebraic identity. **Deleted entirely.**

The controls that bind. Any failure voids the endpoint and the branch is DEMOTE, which is the conservative direction.

| id | control | requirement | status |
|---|---|---|---|
| C1 | Sham dominance (Section 3.3) | `z_sel <= z_sham - 0.5` | part of the decision rule. Unrun on P204. |
| C2 | Discriminant validity: LOEUF must not validate any axis handed to it | `van_elteren(loeuf, fail_tolerance, ...)` among homeostasis survivors: `|z| < 1.0`. `van_elteren(loeuf, fail_evidence, ...)` over 6,371 must not point "less". | observed z = +0.174, p = 0.569; `p_less` = 0.995. Both fire. |
| C3 | Negative-control population | Same design, same groups, on evidence-**failers** with `rest_qc_pass & ~fail_tolerance`. Require `|z| < 2.0`. | **UNRUN.** Must be run and reported. |
| C4 | Negative-control exposure, abundance | Relabel by `1[n_cells_target < within-stratum median]`, same outcome and strata. Require `|z| < 2.0`. | observed on 214: p = 0.981. Fires. |
| C5 | Random-split calibration | 500 size-matched random within-stratum splits. Require `P(p < 0.05) ∈ [0.02, 0.09]`. | observed on 214: 0.034, median p = 0.479. Fires. |
| C6 | Bin stability | `p_sel < 0.05` at `n_bins ∈ {4,5,10,20}` | observed on 204: 0.00182 / 0.00184 / 0.00227 / 0.01899. Fires. |
| C7 | Permutation | reported, non-decisive | observed on 214: perm P = 0.0060, null mean -0.023, sd 1.004. |

Five sham axes tested on the 214 are null: reject-by-high-`stim_de_genes` p = 0.272; reject-by-low-`n_cells_target` p = 0.981; reject-by-low-`rest_cells_ratio` p = 0.896; reject-by-high-`tolerance_loss` p = 0.272; reject-by-high-`efficacy` p = 0.770. Exactly one fires: reject-by-high-`rest_de_genes`, z = -3.588. That is C1, and it is why C1 is the decision.

---

## 6. Secondary measures. Reported, BH-corrected within family, never decisive.

The family is declared here, in full, before the run. Nothing may be added or removed after. All are computed on `P204`, strata `deciles(z_l2, 10)`, `in_top = fail_homeostasis`.

| measure | test | alternative |
|---|---|---|
| `lof_intolerant` | CMH | greater |
| `systemic_risk` | CMH | greater |
| `is_ceg` | CMH | greater |
| `ubiquitous` | CMH | greater |
| `is_iei` | CMH | **less** |
| `n_nonimmune_tissues` | van Elteren | greater |
| `max_nonimmune_ntpm` | van Elteren | **greater** |
| `rest_cells_ratio` | van Elteren | less |

`benjamini_hochberg` over these 8. `max_nonimmune_ntpm` is "greater", not "less"; pinning it here prevents a post-hoc tail swap.

Also reported, not in the BH family:

* Redundancy cross-tab against tolerance on the 286 evidence-passers: homeostasis-only 148, tolerance-only 20, both 52, safe 66. Plus `van_elteren(tolerance_loss, fail_homeostasis, deciles(z_l2,10), "greater")` = 0.815. State that the axes are near-orthogonal and that demotion moves 148 genes.
* `log2_stim_over_rest_cells`, van Elteren, greater: z = +0.303, p = 0.381. Reported with the note that it cancels library representation and therefore also cancels any perturbation depleted in both arms, so it is not a resting-damage measure either.
* The incumbent primary, as an audit trail: `van_elteren(rest_cells_ratio, fail_homeostasis, deciles(z_l2,10), "less")` on the 258 gives z = -1.3629, p = 0.0865; on `P204` it gives z = -0.541, p = 0.2943. Printed with the disqualification in 3.1 attached.
* Sensitivities on the primary, reported, never decisive: `P214` (tolerance survivors, all rest classes) p = 0.0043; all 286 evidence-passers p = 0.00064.

---

## 7. KNOWN LIMITATIONS OF THIS TEST

Every objection the panel raised that was CONFIRMED or PARTIAL and is not fully closed by this design.

**L1. Unblinding. The endpoint is not blind and this document does not claim it is.** Computed in adversarial review before commit: the incumbent primary on the 258 (z = -1.3629, p = 0.0865) by five reviewers independently; the LOEUF primary on `P204` (z = -2.838, p = 0.00227) and on `P214` (z = -2.6255, p = 0.0043); its bin sweep; the entire 8-member secondary family; the 214-row sham dominance contrast (`z_sham = -3.588`). **LOEUF was chosen as primary after all eight family members were seen.** The only unrun statistics are C1 on `P204` and C3. The decision rule's thresholds were fixed before those two exist. That is the whole of the protection this document buys.

**L2. Collider selection on the population is not removed.** `P204` conditions on `rest_qc_pass`, whose central term `ontarget_significant` at rest is a descendant of both arms. Within z_l2 deciles: CMH on the outcome arm OR = 0.5821, p = 1.28e-11; on the exposure arm OR = 2.470, p = 2.73e-08. `van_elteren(rest_cells, deleted, z_l2_decile, "less")` = -8.686, p = 1.9e-18. Survival to `rest_qc_pass`, binned on the resting cell count itself, rises monotonically 69.2% → 90.6%. LOEUF being external does not make the *population* unselected. Whether the selection also depends on LOEUF is unverified. The 214 and 286 sensitivities are the only mitigation.

**L3. Within a z_l2 stratum, the group label is a threshold on `rest_de_genes`, not a property of selectivity.** Mean within-stratum AUC 0.996 versus 0.463 on the stim ingredient. `Spearman(z_l2, stim_de_genes) = +0.990` on the 258, so 4.0% of `log1p(stim_de)` variance survives stratification and 81.8% of `log1p(rest_de)` does. The primary therefore tests the denominator. C1 exists precisely because of this and is expected to fail. **A RETAIN would still not validate the composite's form.**

**L4. z_l2 is a common effect of true effect size and cell count.** `z = lfc/se`, and `se` shrinks with cells: `Spearman(z_l2, n_cells_target) = +0.341`. Conditioning on it induces a within-stratum trade-off; `Spearman(rest_de, rest_cells)` runs +0.124 marginally and -0.087 within decile. The two induced channels enter selectivity with opposite signs and largely cancil; the measured cell-count content of the adjustment on the incumbent endpoint was `dz = -0.137`, about 10% of z. z_l2 is retained because N7 established it as the magnitude confound and it is the only pre-declared matching variable. This is documented, not fixed.

**L5. FIX A2 (`rest_qc_fail` → indeterminate) is an inference, not a repair.** Its group-level control is null where it acts: among the 286 evidence-passers, `rest_unmeasured` rows pass homeostasis 7/22 (31.8%) versus 73/258 (28.3%), Fisher one-sided p = 0.4463; CMH on `stim_de` quintile z = 0.290, p = 0.386; the resting-DE depression restricted to the 280 evidence-passers is z = -1.343, p = 0.090. The screen-wide depression (z = -13.30, n = 736) does not reproduce in the population where FIX A acts. Power at n = 22 is 0.25 against a 40% pass rate. FIX A2 is kept because it is monotone conservative (6 out, 0 in), not because the harm was demonstrated. Only FIX A1 (D1, CENPE) rests on a demonstrated defect (Fisher p = 0.0021).

**L6. D2's mechanism is sign-reversed for two of the six removals.** HEXIM1 fails rest QC on `distal_offtarget_flag`, QRSL1 on `neighboring_gene_KD`. Both have `ontarget_significant == True` at rest. Off-target and neighbour contamination *adds* DE calls, inflating `rest_de` and pushing selectivity down, i.e. biased **against** passing. They passed anyway. Excluding them is not warranted by D2's stated mechanism. They are excluded under a uniform rule (the resting arm must pass the same QC the stimulated arm passed) rather than by per-gene reasoning, because per-gene reasoning after seeing the names is exactly the forbidden move. Sensitivity: `safe = 62` if they are admitted.

**L7. `rest_absent` is a censored value of a cell count, not missing data.** Median Stim48hr `n_cells_target` 104 (IQR 54 to 168) versus 563 for rows with a rest row, MWU one-sided p = 1.07e-26. Rest-row presence is 100% above the bottom stim-cell decile. The rest arm has a hard floor: min 17, p1 = 50, so `rest_cells_ratio < 17/543 = 0.0331`. The mirror set of 185 rest-only genes has median rest cells 85. For any cell-count outcome these rows must be censored at bottom rank, not dropped. This is moot for the LOEUF primary and is documented here so nobody re-derives it. The 5 rows are CENPE, DCK, ICMT, IL2RB, PAQR8, median stim cells 80 versus 267.5 for the rest of the pool. The mirror set also shows arm dropout tracks global guide abundance, so "catastrophic resting viability" overreads it.

**L8. Class counts disagree at one boundary row.** Two independent rebuilds give `rest_unmeasured` 736 or 737 and `rest_not_expressed` 65 or 66. The `rest_qc_fail` total of 802 is stable, which is why the split is retired. D2's "637 rows fail because `ontarget_significant` is False" is the count over all 802 rest-QC failures; within the 736/737 it is 575 (78.1%), and 637 double-counts the 62 `low_target_gex` rows the drafted FIX A simultaneously exempted.

**L9. The permutation control has zero power against every defect above, and the sign flip has zero power against anything.** Section 5. Do not cite either in the report as evidence of validity.

**L10. `is_iei` runs against the safety reading.** OR = 0.202, p = 0.0183. It is in the BH family and it is quoted in the report on both branches. The project has already established (commit 74da6bc) that the IUIS flag cannot separate a dangerous knockdown from a gene that is load-bearing in immunity, so it weakens the constraint secondaries rather than supplying a rival mechanism. It is not decisive.

**L11. LOEUF measures organism-level LoF intolerance, not resting-CD4 damage.** No resting-CD4 harm endpoint exists in these tables that is not a function of `rest_cells` or `rest_de`. Constructing one would require conditioning on the resting cell count while using an outcome that is not a function of it. **No such column exists.** This is stated in the report as a limitation of the screen, not of the analysis.

**L12. The 258-row incumbent endpoint's p-value sits inside a bias envelope of [0.0231, 0.3081]** generated purely by how the outcome-dependent filter's casualties are handled. Restoring the 22 `rest_qc_fail` evidence-passers gives p = 0.0484 (PASS); censoring the 5 `rest_absent` at bottom rank gives p = 0.3077 (FAIL). The exclusion set, not the data, picks the exit code. This is why the incumbent is disqualified rather than merely reported as a FAIL.

**L13. The dominance criterion is a decision rule, not a hypothesis test.** `Δz` has no analytic null. The 0.5 margin is arbitrary. It is fixed here.

**L14. Bin count and stratifier are real degrees of freedom.** On the incumbent, p moved 0.056 to 0.207 across 4 to 20 z_l2 bins, and z moved -0.079 to -2.407 across four defensible stratifiers. C6 makes bin stability a pass condition for the LOEUF primary. The stratifier is pinned to `z_l2` and may not be changed.

**L15. Demotion is expensive and is not a redundancy cleanup.** Section 4.2.

**L16. `n_cells_target` may never be used as a stratifier for a resting-cell outcome.** `Spearman(rest_cells_ratio, n_cells_target) = +0.931`; cell-count deciles absorb 71% of that outcome's variance. Under it the incumbent returns z = -0.079 whether or not the effect is present. Recorded so nobody later reads that null as reassurance.

---

## 8. Code changes, by file

Run order: `04` → `05` + `07` + `14` → `05` again → quarto.

**`scripts/04_window_score.py`**
* Delete `frame["rest_de_genes"].fillna(median)` at line 182.
* Build `rest_qc_class`, `rest_qc_fail_reason`, `rest_de_is_ontarget_only` from the deduplicated Rest `obs` rows (Section 2.7).
* `selectivity = NaN` when `rest_qc_class == "rest_absent"`.
* Replace the boolean `fail_homeostasis` with `homeostasis_verdict ∈ {pass, fail, indeterminate}`. **Keep `fail_homeostasis` as a derived boolean** (`homeostasis_verdict == "fail"`), because N6's primary reads the shipped label off `window_score.csv` and because 12 must be aligned to it.
* Replace line 233 `frame["safe"] = ~frame[fails].any(axis=1)` with the three-valued `gate_outcome` of Section 2.4. `safe` becomes `gate_outcome == "safe"`.
* `reject_reason` must not cite `homeostasis` when the verdict is indeterminate. Add the token `indeterminate` for rows whose only non-pass is an indeterminate homeostasis verdict.
* Docstring: replace "tenfold selectivity" with the implemented rule, verbatim (Section 9, D6). Keep the antiproliferation note at 191-194.
* `viability_tier` (195-208) and `log2_stim_over_rest_cells` are unchanged and remain display-only. They are **not** gate inputs. Verify that `window_score` (line 238) still reads only the three axes.

**`scripts/12_magnitude_matched.py`**
* Line 112 registers `Pillar("DE at rest (raw)", "rest_downstream", binary=False, registered=True)`. That is `n_downstream`. The shipped axis gates on `n_total_de_genes`, and 67 rows flip between them. Fix now: point the pillar at the shipped `selectivity` / `fail_homeostasis` columns, or rename the pillar to "DE at rest (downstream only)" and delete the comment claiming it tests the shipped axis. Preferred: read the shipped columns. See D5.

**`scripts/05_figures.py`**
* `REASON_LABEL`: add key `indeterminate`. A missing key raises `KeyError` on the new `reject_reason` vocabulary.
* On the DEMOTE branch: **remove** the `homeostasis` key. Any facet or legend keyed on it will render empty.
* Figures 1 and 4 rebuild on the new `gate_outcome`. The safe count in the caption is a variable, not a literal.

**`scripts/07_*.py`**
* Consumes `window_score.csv`. Rerun. Audit for any hard-coded `66`, any `df.safe.sum()` assumed to be 66, and any boolean read of `fail_homeostasis` that now has an indeterminate third state upstream.

**`scripts/14_reversal_specificity.py`**
* `AXES`: the tuple and every `reject_reason.str.contains(...)` call. `str.contains("homeostasis")` must not match `indeterminate`, and must not silently count an indeterminate homeostasis verdict as a homeostasis rejection. Use exact-token matching on the split list, not substring matching.
* The matched-null pool becomes the 282 evidence-passers whose gate outcome is determinate under the scoped rule, not 286 and not 259.
* On the DEMOTE branch: delete `homeostasis` from `AXES` and drop its matched-null column. `ALPHA = 0.01` is unchanged.
* Reissue `results/tables/reversal_specificity.csv`. Restate the headline as the exact P from the rerun. Do not carry the string `P<0.0001` forward.

**`report.qmd`**
* Every `.replace("homeostasis", ...)` must be updated in lockstep with 05 and 14, or `indeterminate` renders as a homeostasis rejection and a removed axis renders as an empty facet.
* Strike the NSD1 sentence unconditionally.
* Delete "It is a CELL-COUNT readout, not a DE-call readout" and any independence claim built on it.
* Replace "tenfold" with the implemented rule (D6).
* Print `safe` as `66 - 1 (D1, arithmetic) - 5 (FIX A2, an inferential exclusion whose group-level control is null at p = 0.446) = 60`, and carry `65` (FIX A1 only) alongside every quotation of the shortlist size.
* Add the permanent report lines from Section 4.1.

**`src/cd4_perturbseq/stratified.py`**
* No change. `van_elteren`, `cochran_mantel_haenszel`, `benjamini_hochberg`, `deciles` are used as-is. `deciles(values, n_bins=10)` is `pd.qcut`. Confirm the docstring says so.

**New: `scripts/15_n6_selectivity_validation.py`**
* Implements Sections 3, 4, 5, 6 exactly. Emits `results/tables/n6_primary.csv`, `n6_controls.csv`, `n6_secondary.csv`. Reports `n_strata_used`, `n_strata_dropped`, `n_top_used`, `n_bg_used` for every test. Any specification dropping more than 5% of the rejected rows is flagged as a protocol deviation in the output table.
* It reads only `results/tables/*.csv`. It never touches the h5ad layers.

**New: `docs/n6_prereg.md`**
* This document. Committed before `15` exists.

---

## 9. The incidental defects D4, D5, D6

Changing a threshold or a transform after seeing the data is **tuning** and is forbidden. Documenting exactly what the implemented rule is, is mandatory.

**D4. `n_total_de_genes == n_downstream + ontarget_significant`, exactly, 0 mismatches, in both conditions.**
Verified in the rest arm on all 11,287 deduplicated Rest rows and on all 6,320 matched rows.
**Verdict: DOCUMENT.** Do not change the count. Selectivity is computed from totals that include the on-target gene, and at rest that `+1` appears iff the knockdown worked. Consequences, to be stated wherever selectivity is printed:
* `rest_qc_pass` requires `ontarget_significant` at rest, so `min(rest_de_genes) = 1` by construction across all 5,518 rows, and 0.0% sit at 0. Outside it, 34.5% of the deleted rows sit at `rest_de == 0`.
* 1,590 of 5,518 `rest_qc_pass` rows (28.8%) have `rest_de_genes == 1`: their entire resting DE signal *is* the on-target gene. Surfaced as `rest_de_is_ontarget_only`.
* The 277 screen rows at `rest_de == 0` are all outside `rest_qc_pass`. 276 of them are already rejected by the evidence floor. The single evidence-passer is EGR2.

**D5. `04` gates on `n_total_de_genes`; `12_magnitude_matched.py:112` tests `homeostasis_cost` built from `n_downstream`. 67 rows flip. The pillar's comment claims it tests the shipped axis.**
**Verdict: FIX NOW**, in `12`, by pointing the pillar at the shipped columns. This is a bug in the pillar, not in the gate. Do not "fix" it by changing the gate: switching `04` from totals to downstream after learning that 67 rows flip is tuning.
N6 is insensitive to it either way. On the 258, rebuilding the group from `n_downstream` flips 1 row and moves the incumbent statistic from z = -1.3629 to -1.4295. On `P204` the primary reads the shipped `fail_homeostasis` column directly, so D5 cannot reach it.

**D6. `MIN_SELECTIVITY = np.log(10.0)` is compared against a difference of `log1p`s.**
The implemented rule is `(1 + stim_de) / (1 + rest_de) >= 10`, not `stim_de / rest_de >= 10`. The docstring and the report both say "tenfold".
**Verdict: DOCUMENT ONLY. The threshold and the transform are frozen.** Changing either now, having seen where rows fall, is tuning and is forbidden on both branches.
The mandatory documentation, verbatim, in the `04` docstring and in the report:

> The implemented homeostasis rule is `(1 + stim_de_genes) / (1 + rest_de_genes) >= 10`, applied to DE counts that include the on-target gene (D4). It is not `stim_de / rest_de >= 10`. The implied bar on `stim_de_genes` is `10 * (1 + rest_de_genes) - 1`: it is 9 at `rest_de = 0`, 19 at `rest_de = 1`, and 39 at the value 3.0 that the deleted median imputation used to supply. Across the 258 evidence-passers with a QC-passing resting arm the bar quantiles at 0/25/50/75/100 are 19 / 329 / 1,739 / 8,014 / 37,919. The word "tenfold" does not appear in this report.

On a DEMOTE, D6 becomes an annotation caveat and the rule stops gating anything. On a RETAIN, D6 remains an open, documented defect, because the primary contains no information about the threshold's location (L3).

---

## 10. DO NOT REDO

* `van_elteren(rest_cells_ratio, fail_homeostasis, deciles(z_l2, 10), "less")` on the 258 was run. z = -1.3629, p = 0.0865. It estimates the sign of a `rest_cells` / `rest_de` power-versus-viability coupling, not a homeostasis filter. **It is not evidence for or against demoting selectivity.** Do not re-run it as a decision.
* Within-stratum permutation and sign flip are not falsification controls. Do not cite them as such.
* `n_cells_target` deciles as a stratifier for any resting-cell outcome. Uninformative by construction (L16).
* Equal-width bins. Out of spec. They reach p = 0.0303 on the incumbent and 7 of 20 strata are uninformative.
* The `rest_unmeasured` / `rest_not_expressed` split. Retired (L8).
* Reading the h5ad layers. Nothing in N6 requires them.
