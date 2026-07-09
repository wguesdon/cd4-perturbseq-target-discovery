# Pre-registration — N21: does the decision layer replicate against an independent, signed, functional T-cell screen?

Written 2026-07-09, **before** `scripts/33_freimer_functional_overlay.py` exists and before any sign in
the Freimer data has been inspected for any gene in our nomination pool. Committed before the run.
Both branches are stated in advance and both are publishable.

## 0. What this is not

It is **not** a reopening of discovery. It is one bounded overlay with a hard stop rule.

It is **not** Schmidt. **RULE #3 forbids it**, and the reason is documented: `TRAF6` was a false lead
*precisely because* Schmidt was used as a filter, and it has one module gene down. Schmidt is also the
only out-of-sample validation the efficacy axis has (AUROC 0.702, CI [0.591, 0.814], against its 33
IL-2-reducing hits). Using Schmidt to nominate would destroy that validation to buy a gene. It stays
held out.

**Freimer 2022 has never been used for direction.** N19 collapsed its three marker arms to a binary
"replicates: yes or no" and discarded the sign. That sign is unspent.

## 1. Why Freimer is the right instrument

Freimer et al. 2022 (Marson lab) is a CRISPR screen in **primary human CD4 T cells**, with three
signed marker arms: `IL2RA`, `IL2`, `CTLA4`. MAGeCK `neg|*` means knockdown **depletes** the
marker-high population, so the gene is a **positive regulator** of that marker.

This matters because it measures **both halves of our thesis at protein level, in an independent lab**:

- **`IL2`** is the activation output our **efficacy axis** scores from mRNA.
- **`CTLA4`** is one of the nine genes in our **activation-induced co-inhibitory module**.

And it is genome-agnostic to prior characterisation. That is the defect that killed direction inference
in N17 (eQTL), IMPC and MGI: every one of those instruments can only speak about genes a drug or a
syndrome has already characterised. A pooled functional screen cannot tell whether anyone has studied
the gene. It is therefore the only instrument in this project capable, in principle, of promoting an
uncharacterised gene.

**What it does NOT measure.** Autoimmune therapeutic direction. A knockdown that lowers IL-2 could be
`CD3E`, whose loss causes immunodeficiency. Lowering IL-2 is an **efficacy** readout, not a safety or
direction readout. Any claim built on this must say so in the same sentence.

## 2. The hypotheses and their decision rules, fixed now

Let the **co-tested set** be the genes present in both our rankable perturbations and a given Freimer
arm. Let `tolerance_loss` be our mRNA co-inhibitory attrition score and `efficacy` our activation
suppression score. Let `freimer_lfc` be the signed MAGeCK log fold change for that arm, where a
**negative** `neg|lfc` on the depleted side means knockdown lowers the marker.

**H1 (primary). Our mRNA co-inhibitory attrition predicts CTLA-4 protein loss in an independent screen.**
> Test: Spearman correlation between `tolerance_loss` and CTLA-4 depletion across the co-tested set,
> plus the same correlation against a background matched on our panel-wide effect magnitude `z_L2`,
> so the association cannot be bought by effect size.
>
> **Decision rule.** If rho > 0 with p < 0.05 **and** the magnitude-matched association survives, the
> co-inhibitory axis replicates at protein level, in a different lab, on a different platform. If it
> does not, the axis is an mRNA-only phenomenon and we say so. **Both outcomes are reported.**

**H2 (primary). Our efficacy axis predicts IL-2 loss in an independent screen.**
> Same test on the `IL2` arm. This is an efficacy replication, not a safety claim.

**H3 (exploratory, and labelled so). The hypothesis-target card.**
> From the **91-gene rebuilt nomination pool only** — never a genome-wide re-scan — retain genes that:
> 1. lower IL-2 in Freimer (`IL2` arm, FDR < 0.10, knockdown depletes), and
> 2. do **not** lower CTLA-4 in Freimer (`CTLA4` arm), and
> 3. are not curated approved-immunomodulator targets, and
> 4. are not housekeeping or respiratory-chain genes (by name prefix **or** function), and
> 5. are not RULE #6 competitor-overlap genes, and
> 6. are not direction-discordant.
>
> Rank by: Freimer IL-2 significance, then our internal efficacy, then weaker co-inhibitory attrition,
> then structural tractability. **Emit at most three cards.** Each card must carry an explicit
> "why this might be wrong" line. The verdict line on every card is fixed in advance and reads:
> **"Hypothesis target for follow-up, not a vetted therapeutic target."**

## 3. Controls. All must pass, or the overlay is void.

- **Positive control.** The TCR-proximal signalosome genes our naive ranking returns and our gate
  refuses (`ZAP70`, `VAV1`, `CD3G`, `CD247`, `LCK`, `ITK`, `PLCG1`), where co-tested, must show
  **IL-2 depletion** in Freimer. They are the genes an independent screen must agree lower IL-2. If
  they do not, the arms are mis-signed and nothing here is interpretable.
- **Negative control.** The `Non-Targeting` row must not register as a hit on any arm.
- **Sign control.** `IL2RA` and `IL2` arms must agree in direction for the positive-control genes.
  They measure related biology; if they disagree, the sign convention is wrong.
- **Fan-out control.** Freimer is 4,053 rows over 1,351 genes (three arms). Every merge must be
  asserted one row per gene per arm. N19 was injured by exactly this fan-out: three marker arms per
  gene inflated a spurious drug control from p = 0.17 to p = 0.001.
- **Coverage gate, pre-declared.** If fewer than **10** of the 91 nomination-pool genes are co-tested
  in the `IL2` arm, **H3 is declared underpowered and no card is emitted**, regardless of what the
  data show. H1 and H2 still run on the full co-tested set and are reported.

## 4. The stop rule

**The overlay runs once.** No threshold is re-tuned after seeing the output. No second arm is
substituted for a disappointing first. If H3 returns nothing, the result is:

> "Even after adding a held-out, signed, primary-human-CD4 cytokine screen, no gene passed the
> pre-declared hypothesis-target standard."

If H3 returns one or more genes, the result is:

> "The pipeline abstains from a vetted discovery claim, but an independent primary-human-CD4 CRISPR
> screen supports [GENE] as a hypothesis target for follow-up."

Neither sentence may be strengthened. The word "discovered" does not appear in either.

---

## ADDENDUM 1 — the first run failed two controls, and it was right to. Recorded, not quietly fixed.

Written after the first execution of `scripts/33_freimer_functional_overlay.py`, which **exited 1**.
Nothing from that run was used. Recorded here in full, in the manner of the N9 addendum.

**Failure 1: the registered positive control does not exist in the library.** Freimer is a **focused
1,351-gene IL2RA-regulator library, not genome-wide.** Of the seven TCR-proximal signalosome genes
registered as the positive control, **zero are present** (`ZAP70`, `VAV1`, `CD3G`, `CD247`, `LCK`,
`ITK`, `PLCG1` — all absent). So are `CD3E`, `PPP3R1`, `IMPDH2` and `CD28`. Only `IL4R` is in scope.
The registered control was **NOT TESTABLE**, which is a fact about the library and could have been
checked before registration. It was not. That is a registration error.

> **Substituted positive control, and why it is stronger.** In each arm, **the marker gene itself**
> must be the strongest "knockdown lowers the marker" hit: knocking out `CTLA4` must lower CTLA-4.
> This control requires no outside knowledge, is available in every arm, and is entirely independent
> of any gene in our nomination pool. It is a better control than the one registered.

**Failure 2: the sign convention in the data is the opposite of what this repo assumed.**
The substituted control settles it. In each arm the marker gene is significant on the **`pos`** side
and maximally non-significant on `neg`:

| arm | marker gene | lfc | `neg` rank | `pos` rank | `pos` FDR |
|---|---|---|---|---|---|
| IL2RA | `IL2RA` | +2.682 | 1351 | **2** | 4.5e-4 |
| IL2 | `IL2` | +1.804 | 1347 | **5** | 2.8e-4 |
| CTLA4 | `CTLA4` | +3.787 | 1351 | **1** | 4.1e-4 |

Therefore **`pos` = knockdown LOWERS the marker** (guides enriched in the marker-low bin), and
`neg` = knockdown raises it. `neg|lfc` and `pos|lfc` are identical in 100% of rows: `lfc` is a
gene-level value, not side-specific.

> **This means `scripts/27_cross_screen_concordance.py:59` has the sign backwards.** Its comment reads
> "depleted (neg) = knockdown lowers the marker", and it assigns `row_dir = "pos_regulator"` when
> `neg|fdr <= pos|fdr`. That inverts **331 genes labelled `pos_regulator` and 141 labelled
> `neg_regulator`** in the committed `cross_screen_concordance.csv`.
>
> **N19's statistical verdict is unaffected**, because its enrichment test keys on `freimer_hit`,
> which is `min(neg|fdr, pos|fdr) < 0.10` and therefore sign-agnostic. The reported non-enrichment
> (best p = 0.21) stands. Only the direction column is wrong. It is corrected, and the correction is
> reported rather than silently applied.

**Failure 3: the negative control tripped for a reason that is real but not the one it tests for.**
MAGeCK aggregates the 593 non-targeting guides into one pseudo-gene. With n = 593 its FDR is tiny by
construction, while its `lfc` is **exactly 0.0** on every arm. An FDR-only hit rule flags it. The hit
rule must therefore require a non-zero effect, not merely significance. This is a defect in the rule
as registered, and it is fixed by requiring `lfc > 0` on the lowering side.

**What is NOT changed.** No hypothesis, no threshold, no ranking rule, and no coverage gate. H1, H2
and H3 run exactly as registered in §2, and the coverage gate of 10 genes stands. Only the instrument's
sign convention, the positive control, and the zero-effect exclusion are corrected — all three from
data that is blind to our nomination pool. The overlay still runs **once** after this addendum.

---

## 5. What this costs if it fails

Nothing that matters. H1 and H2 are worth running on their own: they ask whether the two axes of the
decision layer replicate on an orthogonal platform. A negative H1 is a **substantial finding** — it
would mean the co-inhibitory axis is an mRNA-only construct — and it belongs in the Limitations
whether or not any gene is ever nominated.

Freimer is spent by this analysis. After N21 it may not be used to validate anything else.
