# HANDOFF #1 — Open Targets: can the drug-target-recovery benchmark be salvaged?

**Run in:** Claude Science web UI, with the official Open Targets MCP connector enabled.
**Log as:** CS1 in [`claude_tooling_log.md`](../claude_tooling_log.md).
**Machine:** any. Everything this handoff needs is committed under `resources/handoff_inputs/`.
A `git pull` is sufficient. The 16.8 GB h5ad is not required and must not be requested.

## Why this handoff exists

The headline benchmark was going to be "our ranking recovers the targets of approved
immunomodulators". As of 2026-07-08 that benchmark is failing, and not because the ranking is
bad. It is failing because the positive set is too small and the wrong shape.

Measured facts, from `scripts/03_export_handoff_inputs.py`:

| Stage | Positives |
| --- | --- |
| Hand-curated | 36 |
| Perturbed in the CRISPRi library | 32 |
| Also measured in the 10,282-gene transcriptome | 30 |
| Also surviving perturbation QC, therefore rankable | **20** |

`JAK1` and `JAK3` were never perturbed, so tofacitinib's targets cannot be recovered at all.
With 20 positives against 6,371 rankable genes, the naive effector-suppression score achieves
AUROC 0.542 with a bootstrap 95% CI of [0.373, 0.707]. That interval covers chance. The
benchmark currently proves nothing.

There are exactly two ways out, and this handoff tests both.

1. **Grow the positive set.** If Open Targets knows of approved immune-indication drug targets
   that are inhibitory in mechanism and that ARE in our 6,371 rankable genes, and we simply
   failed to curate them, the benchmark becomes powered.
2. **Accept that it cannot be grown, and say so.** If the intersection is genuinely tiny, then
   drug-target recovery is not a viable headline on this dataset, and we report that as a
   finding rather than torturing the number.

Either answer is publishable. The failure mode to avoid is a number that looks fine because
nobody checked its denominator.

## The direction-of-effect constraint

CRISPRi knockdown removes protein function. A perturbation therefore mimics an **inhibitor,
antagonist, neutralising antibody, or degrader** and cannot mimic an agonist. It also cannot
mimic a drug whose mechanism requires the target to be present: tacrolimus binds `FKBP1A` in
order to inhibit calcineurin, so `FKBP1A` knockdown abolishes the drug's effect rather than
reproducing it. `NR3C1` (glucocorticoid receptor) is the sharpest trap, because it is the most
prescribed immunosuppressant target in the world and glucocorticoids are **agonists**.

Any positive set that ignores this is silently wrong.

## Prompt to paste into Claude Science

> You have the official Open Targets MCP connector enabled. Cite the Open Targets release
> version in your answer.
>
> Context. I am benchmarking a genome-scale CRISPRi Perturb-seq ranking in primary human CD4+
> T cells. Knockdown removes protein function, so a perturbation can only phenocopy a drug whose
> mechanism is INHIBITOR, ANTAGONIST, NEUTRALISING, or DEGRADER. It cannot phenocopy an AGONIST,
> and it cannot phenocopy a drug that must bind its target to work through a third protein
> (tacrolimus binds FKBP1A to inhibit calcineurin, so FKBP1A knockdown abolishes efficacy rather
> than mimicking it).
>
> I attach three files:
> - `immunomodulator_targets.csv` — my hand-curated ground truth, 46 rows, 36 marked
>   `include_as_positive = TRUE` and 10 marked FALSE with a stated exclusion reason.
> - `rankable_genes.txt` — the 6,371 genes that were perturbed in the library AND survived
>   perturbation QC. Only these can appear in the ranking.
> - `measured_genes.txt` — the 10,282 genes measured in the transcriptome.
>
> TASK A, audit my curation.
> For each of the 46 gene symbols, query Open Targets and return:
> `ot_max_phase` (maximum clinical trial phase across all indications),
> `ot_max_phase_immune` (maximum phase restricted to autoimmune, inflammatory, or allergic
> indications), `ot_mechanisms` (the mechanism-of-action strings), and
> `ot_mechanism_class` (one of INHIBITOR, ANTAGONIST, NEUTRALISING, DEGRADER, AGONIST,
> MODULATOR, OTHER).
> Set `ot_lof_mimics_drug = TRUE` only for INHIBITOR, ANTAGONIST, NEUTRALISING, DEGRADER.
> Flag every row where you disagree with my `include_as_positive`, and give a one-sentence
> reason. Do not drop the 10 rows I excluded; I need to know whether you agree with the
> exclusions, especially NR3C1, IL2, FKBP1A, and PPIA.
>
> TASK B, the one that matters. Find the positives I missed.
> Query Open Targets for ALL human targets that have at least one drug with
> `maximumClinicalTrialPhase = 4` whose indication is autoimmune, inflammatory, or allergic
> disease, and whose mechanism of action is inhibitory, antagonist, neutralising, or degrading.
> Do not restrict to my list.
> Then intersect that target set with `rankable_genes.txt`.
> Return every gene that is in the intersection but NOT already marked
> `include_as_positive = TRUE` in my table, with its drug, indication, mechanism, and phase.
>
> Report three counts explicitly:
>   1. how many phase-4 immune-indication inhibitory targets Open Targets knows in total,
>   2. how many of those are in `measured_genes.txt`,
>   3. how many of those are in `rankable_genes.txt`.
> Count 3 is the maximum achievable positive-set size for this benchmark. I need that number
> more than I need anything else in this task.
>
> TASK C, sanity.
> Confirm from Open Targets that `JAK1` and `JAK3` have approved inhibitors, so that when I
> report "JAK1 and JAK3 are absent from the perturbation library" a reader understands what was
> lost.
>
> Return one CSV for Task A with my original columns plus the `ot_*` columns and
> `disagrees_with_curation`, `disagreement_reason`. Return a second CSV for Task B with columns
> `gene_symbol, drug, indication, mechanism, max_phase, in_measured, in_rankable`.
> Then have the reviewer agent verify that every count you reported is reproducible from the
> data you pulled, and give me the provenance artifact reference.

## Attachments

All three are committed. From the repo root:

- `resources/ground_truth/immunomodulator_targets.csv` (6 KB)
- `resources/handoff_inputs/rankable_genes.txt` (40 KB, 6,371 genes)
- `resources/handoff_inputs/measured_genes.txt` (65 KB, 10,282 genes)

## What we do with the answer

Save Task A to `resources/ground_truth/immunomodulator_targets.open_targets.csv` and Task B to
`resources/ground_truth/open_targets_additional_positives.csv`. Commit both next to the
hand-curated original, and commit the diff. The diff is a deliverable: it shows the positive set
was independently verified rather than tuned until the AUROC looked good.

Then:

- **If count 3 is above roughly 60**, the benchmark is powered. Rebuild it on the union set and
  report AUROC on both the hand-curated and the Open Targets positive sets, so the reader can
  see the sensitivity to the gold standard.
- **If count 3 stays near 20**, drug-target recovery is retired as the headline. It is reported
  honestly as an underpowered check, and the efficacy axis is instead validated against the
  Schmidt and Steinhart 2022 CD4+ IL-2 CRISPRi screen, where we already have AUROC 0.702 with a
  95% CI of [0.591, 0.814] over 33 independent hits.

A ranking whose apparent quality swings on the choice of gold standard is a ranking to distrust.
Reporting that sensitivity is worth more than a clean single number.
