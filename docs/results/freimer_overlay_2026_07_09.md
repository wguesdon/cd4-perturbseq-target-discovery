# N21 — An independent, signed, primary-human-CD4 screen replicates the efficacy axis and does not replicate the co-inhibitory axis

Pre-registered in `docs/preregistration_n21_2026_07_09.md` before the script existed and before any
Freimer sign was inspected for any gene in the nomination pool. Addendum 1 records a control failure
on the first run. Script `scripts/33_freimer_functional_overlay.py`; tables
`results/tables/freimer_overlay.csv`, `freimer_posthoc_contrast.csv`, `freimer_hypothesis_cards.csv`.

## 1. Why this instrument, and why not Schmidt

The question was whether any external dataset could still yield a defensible novel target. **Schmidt &
Steinhart 2022 cannot be used.** RULE #3 holds it out, `TRAF6` was a false lead precisely because
Schmidt was once used as a filter, and Schmidt is the **only out-of-sample validation the efficacy axis
has** (AUROC 0.702, CI [0.591, 0.814]). Spending it to nominate a gene would destroy the validation.

**Freimer 2022** (Marson lab) is a CRISPR screen in primary human CD4 T cells with three signed marker
arms: `IL2RA`, `IL2`, `CTLA4`. Two of them map directly onto the two halves of this project's thesis.
`IL2` is the activation output our **efficacy** axis scores from mRNA. `CTLA4` is one of the nine genes
in our **activation-induced co-inhibitory module**. N19 collapsed all three arms to a binary
"replicates: yes or no" and discarded the sign, so the sign was unspent.

Critically, it is **drug-label-independent and signed**. A signed functional readout does not ask
whether a drug exists, and that is the defect that killed direction inference in N17, IMPC and MGI.

> **Corrected after review.** An earlier draft called Freimer "not characterisation-bound". That is
> wrong. Its library is **targeted** — transcription factors plus candidate regulators plus controls —
> so gene *inclusion* was a human judgement about what might regulate IL2RA. That is a characterisation
> bias at the library-selection level. The defensible claim is drug-label-independent, not blind to
> prior study. It is a weaker instrument than the earlier sentence implied.

**What Freimer does not measure: autoimmune therapeutic direction.** A knockdown that lowers IL-2 could
be `CD3E`, whose loss causes immunodeficiency. Lowering IL-2 is an **efficacy** readout, not a safety
readout. Every sentence below respects that.

## 2. The first run failed its controls, and the failure found a bug in a committed script

Recorded in full in Addendum 1, not quietly fixed.

The registered positive control (the TCR-proximal signalosome) is **entirely absent** from Freimer's
focused 1,351-gene library. It was substituted with a stronger, pool-blind control: **the marker gene
inside its own arm**. Knocking out `CTLA4` must lower CTLA-4.

That control settles the sign, and it is the opposite of what this repo assumed. In every arm the
marker gene is significant on the **`pos`** side and maximally non-significant on `neg`:

| arm | marker | lfc | `neg` rank | `pos` rank | `pos` FDR |
|---|---|---|---|---|---|
| IL2RA | `IL2RA` | +2.682 | 1351 | **2** | 4.5e-4 |
| IL2 | `IL2` | +1.804 | 1347 | **5** | 2.8e-4 |
| CTLA4 | `CTLA4` | +3.787 | 1351 | **1** | 4.1e-4 |

So **`pos` = knockdown LOWERS the marker.** `scripts/27_cross_screen_concordance.py:59` had this
backwards, mislabelling **331 genes as `pos_regulator` and 141 as `neg_regulator`** in the committed
`cross_screen_concordance.csv`. **N19's statistical verdict is unaffected**, because its enrichment
test keys on `freimer_hit`, which is `min(neg|fdr, pos|fdr)` and therefore sign-agnostic; rerunning
after the fix reproduces it exactly (safe genes vs magnitude-matched null: Freimer p = 0.38, Arce
p = 0.44, either p = 0.21). Only the direction column was wrong. It is corrected.

A third control also fired correctly: MAGeCK aggregates the 593 non-targeting guides into one
pseudo-gene whose FDR is tiny by construction while its `lfc` is exactly **0.0**. An FDR-only hit rule
flags it. The rule now requires a non-zero effect.

After the corrections, all four controls pass: the marker gene lowers its own marker in all three arms;
`Non-Targeting` is silent; the `IL2` and `IL2RA` arms overlap 11 genes against 1.3 expected by chance;
no fan-out.

## 3. H2 — the efficacy axis replicates on an orthogonal, protein-level platform

> Spearman(`efficacy`, Freimer IL-2 lowering) = **+0.135**, n = **471** co-tested genes. It replicates
> under **all four** permutation nulls: stratified on z_L2 (p = 0.0115), on resting-arm disruption
> (p = 0.0065), on stimulated-arm DE burden (p = 0.0060), and on z_L2 and resting-arm disruption
> jointly (p = 0.0070). 2,000 draws each.

Stratifying the permutation on panel-wide effect magnitude means the association **cannot be bought by
effect size**. A post-hoc hit contrast agrees and is stronger: the 29 genes Freimer calls IL-2-lowering
have median `efficacy` +0.342 against +0.151 for the other 443, Mann-Whitney p = 0.0011, and it
survives magnitude matching at **p = 0.0015**.

This is an independent laboratory, a different platform, a protein readout rather than mRNA, and a
different library. It is the strongest external validation in this project.

## 4. H1 — the co-inhibitory axis does NOT replicate against CTLA-4 protein

> Spearman(`tolerance_loss`, Freimer CTLA-4 lowering) = **+0.006**, n = 471. It replicates under
> **none** of the four nulls (p = 0.89, 0.90, 0.91, 0.88). Pre-registered verdict: **DOES NOT REPLICATE.**

A post-hoc, more sensitive contrast does not rescue it. The 13 genes Freimer calls CTLA-4-lowering have
a nominally higher median `tolerance_loss` (+0.193 vs +0.112, Mann-Whitney p = 0.027), but the
difference **dies once matched on effect magnitude: p = 0.31**. The same is true of `efficacy` against
the CTLA-4 arm (matched p = 0.28) and of `tolerance_loss` against the IL-2 arm (matched p = 0.081).

**The diagnostic is calibrated, so this negative is not merely insensitivity.** The identical contrast
fires cleanly on the arm where H2 passed (29 hits, matched p = 0.0015). It has power to detect a real
effect at this sample size. It does not detect one here.

**Limitations, stated before the conclusion.** CTLA-4 is one of the nine genes in our module, not the
module. Freimer measures surface protein in its own stimulation context; we measure mRNA at 48 hours.
There are only 13 CTLA-4-lowering hits, so power is real but limited. And Freimer's library is a focused
IL2RA-regulator panel, so the co-tested set is not a random sample of our screen.

**The conclusion we draw, and the one we do not.** We do not conclude that the co-inhibitory module is
an artifact. We conclude that **it is an mRNA construct that we have validated internally and have not
been able to corroborate at protein level for its best-measured member in an independent screen**, and
that the nominal association there is explained by effect magnitude. That belongs in the Limitations,
stated before any result that rests on the axis.

There is an uncomfortable symmetry worth naming. The efficacy axis, which this project has consistently
described as the commodity half of the pipeline, is the half that replicates externally. The
co-inhibitory axis, which carries the headline, is the half that does not.

## 5. H3 — the hypothesis-target card, and why none was emitted

Of the 91-gene rebuilt nomination pool, 73 remain after excluding housekeeping and RULE #6 genes. Only
**6** of those are co-tested in the Freimer `IL2` arm: `ATF4`, `ATXN7L3`, `PTPRC`, `SMAD4`, `STAT3`,
`XBP1`.

The pre-declared coverage gate was **10**. It was not met.

> **No card was emitted, regardless of what the signs show.** The gate was fixed in advance precisely
> so that a thin, tempting result could not be promoted after the fact.

The registered sentence therefore stands:

> Even after adding a held-out, signed, primary-human-CD4 cytokine screen, no gene passed the
> pre-declared hypothesis-target standard.

The reason is structural rather than biological. Freimer's library is 1,351 genes chosen as candidate
IL2RA regulators. Our nomination pool is drawn from a genome-scale screen. The intersection is small by
construction, and a focused library cannot adjudicate a genome-scale pool. **The right instrument would
be a genome-wide, signed, functional screen in primary human CD4 T cells** — the same conclusion N19
reached from the other direction.

## 5b. The secondary, post-hoc pass — and why its two candidates fail review

**This is not the registered analysis.** Addendum 2 records that a relaxed rule was directed after the
coverage gate fired. It has no coverage gate, promotes on a significant Freimer IL-2-lowering effect,
and treats CTLA-4 loss as a liability flag rather than a disqualifier. It is post-hoc, it was run once,
and every output row is stamped with that provenance. Script `scripts/34_freimer_secondary_posthoc.py`.

It **did** return candidates. Of the six co-tested pool genes, two show a significant IL-2-lowering
knockdown effect: **`ATXN7L3`** (FDR 2.8e-4, lfc +1.263) and **`XBP1`** (FDR 2.8e-4, lfc +0.781).
Neither lowers CTLA-4. That is the **absence of a liability signal**, not evidence of preservation, and
it is labelled *no CTLA4-lowering liability in Freimer*.

**Both fail review, on grounds this pipeline had already stated.**

| gene | stim DE | rest DE | rest-DE percentile | stim:rest ratio | genetics | pocket | precedent |
|---|---|---|---|---|---|---|---|
| `ATXN7L3` | 503 | 545 | **96.2nd** | **0.92** | 0.000 | no | no |
| `XBP1` | 492 | 180 | **92.1st** | 2.73 | 0.000 | no | no |

The selectivity requirement fixed a priori in this project is a **10x** stimulated-to-resting ratio.
Neither meets it. `ATXN7L3`'s ratio is **below 1.0**: its knockdown disturbs the *unstimulated* cell
more than the stimulated one. That is inconsistent with the screen-internal therapeutic-window proxy.

`ATXN7L3` is a component of the **SAGA deubiquitinase module**, together with `USP22` — and `USP22` is
named in this very report as one of the "global transcription machinery" contaminants the naive ranking
returns. `USP22` sits at the 98.7th percentile of resting-arm disruption; `ATXN7L3` at the 96.2nd. They
are the same failure mode. `XBP1` is the master regulator of the unfolded-protein response, broadly
required, LoF-constraint unknown, with no pocket and no clinical precedent.

**And the promotion instrument itself is confounded at its threshold.** Freimer's 29 IL-2-lowering hits
carry a median of **66 resting-arm DE genes against 3 for the 442 non-hits** (Mann-Whitney
p = 9.8e-7). IL-2 is a highly induced transcript, so disabling general transcription machinery lowers
it. A Freimer IL-2 *hit call* is substantially a "this cell can no longer transcribe an induced gene"
signal.

> **Note the asymmetry, because it decides what may be claimed.** The **continuous** association
> underpinning H2 survives stratification on resting-arm disruption (p = 0.0055) and on `z_L2` and
> `rest_de_genes` jointly (p = 0.0075). **The efficacy axis genuinely replicates.** It is the
> **FDR-thresholded hit call** that is confounded — and promotion uses the hit call. So H2 stands and
> the two promoted genes do not.

**Verdict, from the three permitted in Addendum 2:**

> **(3) Freimer produces target-card candidates, but all fail review** — on housekeeping and chromatin
> biology, on resting-arm disruption far outside the a-priori selectivity requirement, and on absent
> tractability. Neither `ATXN7L3` nor `XBP1` is a Freimer-supported follow-up hypothesis worth naming
> to a reviewer.

The registered result is unchanged and governs the headline: 6 of 73 eligible genes co-tested against a
pre-declared gate of 10, and no card emitted.

## 6. What this changes

- The efficacy axis gains an orthogonal, protein-level, independent-lab replication (p = 0.0095,
  magnitude-matched). This is new and it is the strongest external result in the project.
- The co-inhibitory axis gains an honest external negative. It must now be described as internally
  validated and externally uncorroborated.
- A sign-convention bug in a committed script was found by a pre-registered control, and its blast
  radius was bounded and reported. N19's verdict survives.
- No novel target. Freimer is now spent and may not be used to validate anything else.
