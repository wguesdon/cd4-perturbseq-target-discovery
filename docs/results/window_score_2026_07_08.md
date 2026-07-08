# The therapeutic-window score: three failed designs and what the data forced

Date: 2026-07-08. Script: [`scripts/04_window_score.py`](../../scripts/04_window_score.py).
Output: `results/tables/window_score.csv`.

Schmidt and Steinhart 2022, a genome-wide CRISPRi screen for CD4+ IL-2 production from another
lab with a protein readout, is **held out**. It never enters the score, the gate, or any
threshold. Every number involving it below is out-of-sample.

## Three things were wrong, and the data said so

**No evidence floor.** The first score ranked by a mean of module z-scores. Of the 6,371
QC-passing perturbations, 4,131 have five or fewer significant DE genes genome-wide, and among
those the maximum number of effector-module genes significantly down is **one**. Yet 804 of them
carry a mean module z below -0.5. The score was reading sub-significant drift. `CAST` reached 1.19
on two significant DE genes. A perturbation must now significantly suppress at least three module
genes before it is called efficacious.

**Calibration was not the lever.** The hypothesis was that co-regulation among effector genes
inflated the naive mean, and that a CAMERA variance inflation factor would deflate it. Estimated
across all perturbations, `rho_bar` comes out near zero and `VIF` near 1, because the mass of null
perturbations dilutes the correlation. Estimated on the 713 perturbations with at least 100 DE
genes, `rho_bar` is 0.108 and `VIF` is 4.36. Either way the correction is nearly uniform, so it
barely moves the ranking. Held-out AUROC against Schmidt: naive 0.701, CAMERA 0.689. The
calibration made it slightly worse. It is reported because it is the honest statistic, not because
it rescued anything.

**The collateral cap was backwards.** Capping the number of DE genes in stimulated cells rejects
precisely the context-selective targets this project exists to find.

| Gene | Stim DE | Rest DE | What it is |
| --- | --- | --- | --- |
| `ITK` | 2566 | 2 | TCR kinase. A perfect window, and an immunodeficiency gene. |
| `CD3G` | 1374 | 5 | Approved target, teplizumab and muromonab. |
| `PPP3R1` | 523 | 1 | Calcineurin regulatory subunit. Ciclosporin, tacrolimus. |
| `NSD1` | 2175 | 1767 | Transcriptome collapse. |
| `STAT5B` | 1060 | 1135 | Transcriptome collapse. |

The discriminator is disruption **at rest**, never breadth in Stim. The gate now requires the
stimulated effect to exceed the resting effect at least tenfold, a ratio fixed a priori.

## Why the calcineurin benchmark looked broken

`PPP3CA` produces **one** significant DE gene and `PPP3CB` produces three. Both catalytic
subunits, both the target of ciclosporin and tacrolimus, both apparently inert. `PPP3R1`, the
obligate regulatory subunit, produces 523.

The catalytic subunits are paralogues and compensate for one another. Single-gene CRISPRi cannot
see a redundant target. The shared, non-redundant regulatory subunit shows the phenotype. This is
not a flaw in the screen or the score. It is what single-gene perturbation of a redundant enzyme
looks like, and any drug-target-recovery benchmark on this data that does not account for it will
report a false negative and blame the method.

## The result

**Efficacy is real, out of sample.** Against Schmidt's 33 significant IL-2-reducing hits:

| Score | AUROC | 95% CI |
| --- | --- | --- |
| naive, -mean(module z) | 0.701 | [0.587, 0.804] |
| efficacy, -mean(module lfc) | 0.669 | [0.557, 0.777] |
| window score | 0.688 | [0.570, 0.803] |

The window score does not beat the naive score here, and it must not. Schmidt's hits are the TCR
signalosome, which the safety gate rejects on purpose. Drug-target recovery and Schmidt validate
the **efficacy axis**. Neither can validate the **window score**. A window score that won those
comparisons would be a window score that had stopped gating.

**The gate refuses twelve knockdowns that an independent lab proved work.**

| Gene | Stim/Rest DE | Schmidt IL-2 lfc | Rejected for |
| --- | --- | --- | --- |
| `ZAP70` | 2486 / 83 | -2.78 | immune-essential, tolerance |
| `CD3G` | 1374 / 5 | -2.67 | immune-essential, tolerance |
| `CD3E` | 1586 / 4 | -1.96 | immune-essential |
| `CD247` | 828 / 5 | -1.59 | immune-essential, tolerance |
| `CD28` | 1798 / 5 | -1.73 | immune-essential |
| `ITK` | 2566 / 2 | -1.45 | immune-essential, tolerance |
| `PLCG1` | 2218 / 3 | -1.03 | immune-essential, tolerance |
| `VAV1` | 3575 / 70 | -2.92 | tolerance |
| `PTPRC` | 3226 / 55 | -0.86 | immune-essential |

They work. They would leave a patient immunodeficient. That table is the argument.

## The tiers fell out of the data

Viability is read straight off the screen. Cells carrying a guide against a gene the cell cannot
live without are simply absent. Hart core-essentials are measurably depleted at rest, 0.83 times
the median, Mann-Whitney p = 0.003. Resting cell count is therefore the context-free viability
signal, and depletion only under stimulation is antiproliferation, which is mycophenolate's
mechanism rather than a disqualification.

Five approved-drug targets clear the evidence floor. Nothing in the score or the gate was told
which drug is well tolerated.

| Gene | Drug | Resting cells | Verdict |
| --- | --- | --- | --- |
| `PPP3R1` | ciclosporin, tacrolimus | 1.96x | **passes**, non-depleting |
| `IL4R` | dupilumab | 1.85x | **passes**, non-depleting |
| `IMPDH2` | mycophenolate | 0.16x | passes, depleting at rest |
| `CD3E` | teplizumab, muromonab | 0.40x | **rejected**, immune-essential |
| `CD3G` | teplizumab, muromonab | 0.39x | **rejected**, immune-essential, tolerance |

Well-tolerated agents passing the gate: 3 of 3. Narrow-index agents passing: 0 of 2. Five points
is a direction, not a p-value, and it is reported as one.

The three strata match the clinic without being told the clinic. Non-depleting signalling blockade
is the best-tolerated class. Antiproliferatives work and require monitoring. Depleting plus
immune-essential is where muromonab was withdrawn.

## The shortlist

**Tier A**, safety-passing and not depleted at rest: `TFB1M`, `CD2`, `ACAD9`, `IBA57`, `CARS2`,
`ZNF649`, `PPP3R1`, `NDUFAF3`, `LIG3`, `POLG`, `MTHFD2`, `MTG1`, `ECSIT`, `MTCH2`, `USP24`,
`MTHFD1L`, `GRSF1`, `GAS2L1`, `VARS2`, `TMEM126B`, `RRAGC`.

The class is immunometabolism: mitochondrial biogenesis, OXPHOS assembly, one-carbon metabolism,
and `RRAGC`, the Rag GTPase through which mTORC1 senses amino acids, whose knockdown phenocopies
rapamycin. Activated T cells are metabolically demanding and resting T cells are quiescent, so
these carry huge stimulated effects and almost none at rest. `IMPDH2`, an approved drug of exactly
this class, sits in Tier B, which is the internal control for the class being real.

**Tier B**, safety-passing but depleted at rest: `ENO1`, `CCNT1`, `IMPDH2`, `TOMM70`, `COMMD5`,
`PGAP2`, `CENPE`, `RIDA`, `AP2A1`, `PTCD2`.

## `CD2`, and a rule we are not breaking

`CD2` ranks second in Tier A, is confirmed by the held-out Schmidt screen (IL-2 lfc -1.94, FDR
1.5e-4), and is the target of **alefacept**, an approved drug. It is not in our curated ground
truth. We missed it.

We are **not** adding it now. Adding a positive to the gold standard after watching it rank second
is tuning the benchmark to the result, which is the exact failure this project keeps warning about.
The curated table stays untouched. [HANDOFF #1](../handoffs/HANDOFF_01_open_targets_ground_truth.md)
asks Open Targets to find the positives we missed, independently. If it returns `CD2`, that is
evidence. If we had typed it in ourselves, it would be nothing.

## Honest limitations

- **Our safety gate is T-cell-intrinsic, not organism-level.** "Rest" here means resting CD4 T
  cells. Mitochondrial and one-carbon genes would be toxic to gut epithelium and bone marrow, and
  this screen cannot see that. The Tier A shortlist needs a tissue-breadth axis from GTEx or the
  Human Protein Atlas before anyone acts on it.
- **The shortlist is depleted overall.** Median resting cell ratio 0.80. The viability tier
  separates the worst offenders, but cell count is confounded with guide representation, which we
  cannot control without a time-zero library measurement.
- **Twenty positives, chance-level AUROC.** The drug-recovery benchmark remains underpowered. Only
  five positives even clear the evidence floor.
- **`n = 5`** for the therapeutic-index result. Directionally perfect, statistically nothing.
- Effector genes `IL13`, `IL23A`, `IL2RA` are themselves benchmark positives and module members.
  Leave-one-out removes the self-effect, but the circularity of scoring a cytokine by a module
  containing cytokines remains.
- Pseudobulk DE, single-gene library, non-polarised culture. All inherited.
