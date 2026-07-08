# Ground truth: approved immunomodulator targets

The headline benchmark asks whether our therapeutic-window ranking recovers the targets of
approved immunomodulatory drugs. That question is only meaningful if the positive set is
built correctly, and the obvious way to build it is wrong.

Source table: [`resources/ground_truth/immunomodulator_targets.csv`](../resources/ground_truth/immunomodulator_targets.csv).

## The direction-of-effect principle

CRISPRi knockdown removes protein function. A perturbation therefore mimics an **inhibitor,
an antagonist, a neutralising antibody, or a degrader**. It does not mimic an agonist, and it
does not mimic a drug whose mechanism requires the protein to be present.

So "gene X is the target of an approved immunomodulator" is not sufficient for X to be a
positive. The test is narrower: **does loss of function of X move the cell in the same
direction as the drug?** Three classes fail that test, and all three are in the table with
`include_as_positive = FALSE`.

**Agonists.** `NR3C1` is the glucocorticoid receptor, the target of prednisolone and
dexamethasone, unquestionably the most widely used immunosuppressants in this space.
Glucocorticoids are receptor *agonists*. Knocking down `NR3C1` removes the anti-inflammatory
signal. A ranking that scores `NR3C1` knockdown as a hit would be scoring a gene whose
suppression makes inflammation worse. `IL2` fails the same way: aldesleukin and low-dose IL-2
are administered as agonists.

**Drug-binding chaperones.** Tacrolimus and sirolimus bind `FKBP1A`; ciclosporin binds `PPIA`.
The drug-immunophilin complex is what inhibits calcineurin or mTOR. Knocking down `FKBP1A`
abolishes tacrolimus efficacy rather than reproducing it. The effector targets are `PPP3CA`,
`PPP3CB`, `PPP3R1` and `MTOR`, and those are the genes we score.

**Wrong cell type.** `CD80` and `CD86` are the targets of abatacept but are expressed on
antigen-presenting cells. `MS4A1` (CD20, rituximab) and `BTK` are B cell genes. They cannot be
recovered from a CD4+ T cell screen, so counting them as missed positives would understate
precision for a reason that has nothing to do with the method.

Getting this wrong in either direction corrupts the benchmark. Including `NR3C1` inflates the
apparent difficulty and then rewards a score that is anti-therapeutic. Including `MS4A1`
depresses recall for free. Both are silent failures: the numbers still come out.

## Positive set

35 genes with `include_as_positive = TRUE`, spanning the mechanisms the field actually uses:

| Mechanism class | Genes |
| --- | --- |
| Calcineurin/NFAT | `PPP3CA` `PPP3CB` `PPP3R1` |
| JAK/TYK2 | `JAK1` `JAK2` `JAK3` `TYK2` |
| mTOR | `MTOR` |
| Nucleotide synthesis | `DHODH` `IMPDH1` `IMPDH2` `DHFR` |
| Cytokines and receptors | `TNF` `IL17A` `IL17F` `IL23A` `IL12B` `IL6R` `IL4R` `IL13` `IL5` `IL2RA` |
| cAMP | `PDE4A` `PDE4B` `PDE4D` |
| Trafficking | `S1PR1` `ITGA4` `ITGAL` `CCR4` |
| TCR complex and PI3K | `CD3D` `CD3E` `CD3G` `PIK3CD` |
| Degrader archetype | `IKZF1` `IKZF3` |
| Depletion | `CD52` |

`DHFR` is deliberately retained. Methotrexate works, and `DHFR` is also core-essential. It is
the cleanest illustration of the effective-but-toxic class the safety gate is built to
separate, and the benchmark should not hide it.

## Held out, not excluded

`TNFRSF4` (OX40, phase 3) and `STAT6` (KT-621, phase 1) pass the direction test but are not
approved. They are held out as prospective validation: a good ranking should surface them
high without ever having been told about them.

## Coverage caveat

Positives only count if the gene was perturbed in the library and measured in the 10,282-gene
transcriptome. `scripts/03_build_benchmark.py` reports the intersection and the genes lost, so
precision and recall are computed against a stated, reproducible denominator.

## Provenance and open work

This table was curated by hand from the drug classes named in the strategy and literature
review. It is a v0. The planned cross-check is a handoff to Claude Science against the
official Open Targets MCP, pulling `max_phase` and mechanism-of-action per target so the
`approved` and `mechanism` columns are sourced rather than asserted. That round-trip is
logged in [`claude_tooling_log.md`](./claude_tooling_log.md).
