# Adversarial audit: the safety gate is partly tautological

Date: 2026-07-08. Method: a six-lens hostile audit of `docs/results/risk_kill_2026_07_08.md`,
run as a Claude Code multi-agent workflow. 127 agents, 1,985 tool calls, 11.8M tokens. Each of the
50 raw findings was then handed to two further agents, one instructed to refute it and one to
reproduce it empirically.

## First, a caveat about the audit itself

The workflow marked 49 of 50 findings "confirmed". That number is inflated by my own design: a
finding survived if **either** verifier failed to refute it. A one-of-two bar is not adversarial.
Reported here are only the findings I then reproduced myself, from the data, with a script that is
committed.

## The finding that matters

The audit's headline claim, verified independently in `scripts/08_verify_iei_claim.py`:

| Group | IEI rate | Fisher OR vs background |
| --- | --- | --- |
| Naive top 100 | 14.0% | 4.16 |
| **Our own approved-drug positives** | **25.0%** | **8.31** |
| Background | 3.7% | — |

The inborn-errors-of-immunity flag is **more enriched among approved drug targets than among the
perturbations we called toxic.** It therefore cannot distinguish "this knockdown is dangerous"
from "this gene is load-bearing in immunity", and a good immune drug target is, by construction, a
gene that is load-bearing in immunity.

Thirteen of our 36 curated positives are IUIS genes: `CD3D`, `CD3E`, `CD3G`, `IKZF1`, `IKZF3`,
`IL12B`, `IL17F`, `IL2RA`, `IL6R`, `JAK1`, `JAK3`, `PIK3CD`, `TYK2`.

Of the 20 rankable positives, the gate as built **rejects five approved drugs**: `CD3E` and `CD3G`
(teplizumab, muromonab), `IL2RA` (basiliximab), `PIK3CD` (leniolisib), and `TYK2`
(deucravacitinib). A 25% false-rejection rate on approved therapeutics is not a safety gate. It is
a filter that removes immunology.

## A mechanism bug underneath it

The IUIS table has a `GOF/DN` column. Sixteen entries are **gain-of-function** diseases, where
loss of function is not the pathogenic mechanism and therefore says nothing about what CRISPRi
does. `priors.iei_genes()` ignored that column and counted them all: `CXCR4`, `HCK`, `LYN`,
`NFKBIA`, `NLRC4`, `NLRP12`, `NLRP3`, `PLCG1`, `PLCG2`, `SAMD9`, `SAMD9L`, `STAT4`, `STAT6`,
`STING1`, `SYK`, `TLR8`.

`PLCG1` is one of them, and `results/figures/fig2_effective_but_rejected` labels it
"immunodeficiency gene". That label is wrong.

The symbol regex also silently drops `C4A+C4B`. Excluding the gain-of-function entries changes
nothing about the conclusion above: OR moves from 8.31 to 8.64 for the positives, and from 4.16 to
3.96 for the top 100. The tautology is not a parsing artifact.

## Three more critical findings, not yet independently reproduced

Recorded here so they are not quietly dropped. Each is scheduled for reproduction before any of it
is used.

1. **The collateral-DE and rest-DE pillars may be effect-magnitude artifacts.** The audit reports
   that at matched score magnitude, suppressors carry roughly 66% *fewer* stimulated DE genes than
   inducers, and that `Spearman(n_cells_target, stim_de_genes) = -0.243`, so the cell-count-matched
   background controls a confound that runs the opposite way. The fix is a background matched on
   transcriptome-wide effect magnitude, plus a sign-flipped control ranking by *induction*.

2. **The essentiality null is survivorship bias.** 189 of the 257 Hart essentials present in the
   sgRNA library never reach `DE_stats` at all, because their knockdown depletes cells and fails
   the authors' DE-eligibility gates. Conditioning on QC removes exactly the genes whose knockdown
   kills cells. "Zero essentials in the top 100" cannot be evidence that essentiality is the wrong
   axis when essentials were far less likely to reach the ranking at all. The correct statement is
   that **this screen cannot resolve the question**, which we independently rediscovered when only
   3 of 927 Hart nonessential control genes survived the same QC.

3. ~~**The tolerance test is close to tautological.**~~ **TESTED 2026-07-08. THE AUDIT IS WRONG
   HERE.** The audit reported that `naive_suppression` and `tolerance_suppression` correlate at
   Spearman 0.438 and that a surrogate score reproduces the effect, so the axis merely re-encodes
   suppression strength. Correlation is expected and is not the objection.
   `scripts/10_tolerance_is_real.py` runs the two tests that are:

   - **Residual.** Regress tolerance suppression on effector suppression. The residual is still
     elevated in the naive top 100, median +0.1173 versus −0.0106 in background, MWU p = 5.6e-4.
   - **Random-module null.** Draw 200 random 9-gene sets matched to the tolerance module on
     baseline-expression decile, and score every perturbation against each. The real tolerance
     module is suppressed far beyond its own null: z median +3.870 in the top 100 versus +1.258 in
     background, MWU p = 3.5e-13.

   Expression matching matters. The tolerance genes are cytokines and checkpoints, which are lowly
   expressed; an unmatched random set would be housekeeping-dominated and the module would look
   special for the wrong reason.

   The tolerance axis is real and stays in the gate. This matters, because after `is_iei` was
   demoted, six of the nine knockdowns the gate refuses are refused on tolerance alone.

   Recording this because an audit that is only ever right is not an audit. Its claim about the raw
   Mann-Whitney was fair; its conclusion did not follow.

And one that is arithmetic rather than inference: the single seed-0 matched background gives
OR 3.91, p = 0.012, but redrawing it 2,000 times gives a Monte-Carlo interval of OR 1.87–16.12 and
p 0.00032–0.129, with 11.5% of seeds non-significant. A Cochran-Mantel-Haenszel test using all
6,371 rows instead of a 100-row subsample gives OR 3.78, p = 1.1e-05.

## What survives, and what this does to the thesis

**"Reversal is not enough" survives as a descriptive claim.** The top of the naive ranking is the
TCR signalosome and the general transcription machinery. An independent lab's protein-level screen
confirms these knockdowns really do shut IL-2 down. Loss of function of `ZAP70`, `CD3G`, `CD247`,
`LCK` and `ITK` causes human immunodeficiency. None of that depends on a Fisher test.

**The statistical support for it does not survive**, and this is the second time today that the
mechanism we asserted turned out to be wrong. First the contaminants were not common-essential
genes. Now the marker we replaced essentiality with cannot separate a toxic knockdown from a drug
target.

The correct conclusion is more interesting than the one we were defending:

> On-pathway is not the same as safe, and **no binary immune-gene flag can tell them apart.**
> The genes at the top of a reversal ranking are the immune system's load-bearing nodes. Some of
> them are drugged. Muromonab was withdrawn; teplizumab causes cytokine release; the JAK
> inhibitors carry black-box warnings. The difference between a target and a hazard is not
> category membership. It is degree: how selectively, how reversibly, and at what cost to the
> resting cell.

That argument requires a **graded, measured** safety readout, not a database lookup. Selectivity
ratio, viability tier, and tolerance residual are all computed from the screen itself. The parts of
the gate that failed today are precisely the parts that were database membership tests:
DepMap/Hart essentiality, and now IUIS.

## What changes

- `is_iei` is demoted from a hard gate to a reported liability annotation, and the 25%
  false-rejection rate on approved drugs is reported alongside it.
- `priors.iei_genes()` gains the `GOF/DN` filter and recovers `C4A+C4B`.
- The gate is rebuilt on the data-native, graded axes only.
- The essentiality claim in `docs/results/risk_kill_2026_07_08.md` is retracted and replaced with
  the selection funnel.
- The tolerance claim is restated as a residual, with its effect size.
- The single-draw background is replaced with a CMH test over all rows.
- `fig2` is regenerated; `PLCG1` is no longer labelled an immunodeficiency gene.

The audit cost 11.8M tokens and found a hole in the centre of the argument four days before the
deadline. That is what it was for.
