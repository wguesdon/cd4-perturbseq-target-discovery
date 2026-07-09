# N13 — The genetics enrichment, re-done with a continuous confidence weight

Run under RULE #9. Script `scripts/22_continuous_genetics.py`; table `results/tables/continuous_genetics.csv`.

## 1. Question and method (literature-identified)

N10 found the naive suppression ranking enriched for autoimmune genetic support (binary flag, OT score
≥ 0.1 for ≥ 1 of 14 diseases) while the safety-gated ranking was not, with a between-ranking delta whose
CI included 0. The literature review (§1.3, §2 — the "strongest single opportunity") says the fix is a
**continuous** confidence weight, because clinical success tracks *confidence in the causal gene*, not a
binary hit (Minikel et al. 2024, *Nature*, ~2.6× success lift scaling with causal-gene confidence; Nelson
et al. 2015; King, Davis & Degner 2019). A continuous weight is more statistically efficient.

This re-runs the comparison with the continuous OT association score (max across the 14 diseases) instead
of the binary flag, reports the binary result at three thresholds as a sensitivity, and confirms the
background was already correct (§1.1).

## 2. Background check (§1.1)

The matched null already draws z_l2-decile-matched shortlists from the **screened universe** (6,371
perturbations with a known z_l2), not the genome. So the "wrong background" concern the literature review
raised is not a defect here; N10 and N13 both use the correct denominator. Stated for the record.

## 3. Result

| ranking | continuous mean score | matched null | p |
|---|---|---|---|
| naive top 100 | **0.152** | 0.088 | **0.0012** |
| safe-set top 100 | 0.093 | 0.102 | 0.669 |

Binary sensitivity (naive / safe, "supported" count of 100, matched mean, p):

| threshold | naive | safe |
|---|---|---|
| 0.05 | 36, matched 24.1, **p=0.004** | 24, matched 26.5, p=0.76 |
| 0.10 | 32, matched 19.7, **p=0.002** | 21, matched 22.1, p=0.66 |
| 0.20 | 28, matched 16.7, **p=0.003** | 18, matched 19.3, p=0.68 |

Control (label-shuffle calibration, 200 shuffles): false-positive rate at α=0.05 is **0.060** — PASS.

## 4. Critical review — the marginal delta does NOT survive

The between-ranking delta (safe − naive) mean score is −0.059, and a first bootstrap CI *appeared* to
exclude 0 (upper bound near −0.002). Under RULE #9 that marginal boundary was stress-tested before being
trusted, and it fails three ways:

- **Aggregation.** With `gmax`, `gmean` and `gsum` aggregations of the per-gene score across the 14
  diseases, and with a larger bootstrap, the delta CI straddles 0 in every case (e.g. `gmax`
  [−0.119, +0.002], `gmean` [−0.081, +0.009]).
- **Bootstrap.** The CI boundary crosses 0 seed-to-seed. A result that flips with the bootstrap seed is
  not significant.
- **Concentration.** The naive-top genetic enrichment is carried by ~10 high-genetics genes (GATA3,
  CD247, STAT3, IL2RB, CD28, BHLHE40 …). Dropping the top 10 erases the delta (+0.003). These are strong
  positive regulators — the TCR signalosome and its transcription factors — which are both the strongest
  suppressors and heavily GWAS-associated. That is why the *naive* ranking is enriched, and why removing
  them removes the enrichment.

So the continuous weight **sharpened the naive-enrichment estimate** — exactly Minikel 2024's efficiency
gain — but did **not** make the safe-vs-naive difference hold up. The honest verdict is the same as N10's:
the naive ranking is robustly enriched for autoimmune genetic support; the safety-gated ranking is not; the
difference between them is not independently significant at n=100. **Do not report "the gate opposes
genetics."** The persuasive result is the consistent *direction* across three independent measures — the
IEI odds ratio, N10, and here — not any single p-value.

This was almost a fifth overclaim (a marginal CI reported as significant), caught by the stress test. It is
recorded, per RULE #0.

## 5. The tie to N11 (direction of effect)

`PTPN2` has a continuous genetic score of **0.762** (12 of the 14 diseases) — near the top of the whole
screen. The continuous method, the best genetic prior available, therefore promotes PTPN2 exactly as the
binary one did. But PTPN2 is direction-**discordant** (N11): its loss causes autoimmunity, so an inhibitor
is anti-therapeutic. **Genetics is direction-agnostic; no weighting of it fixes direction of effect.**
Genetic support is necessary, not sufficient, and the direction-of-effect axis (N11) is the gate genetics
cannot supply. This is the load-bearing conclusion of the two iterations together.

## 6. Verdict

The genetics enrichment result stands, at its true strength: naive-enriched, safe-not-enriched, delta
under-powered. The continuous weight is the correct method and is now the reported one, with the binary
thresholds as a sensitivity. It does not change the conclusion, and it does not — and cannot — resolve the
direction-of-effect problem that governs whether a genetically-supported candidate is actually the right
way round.
