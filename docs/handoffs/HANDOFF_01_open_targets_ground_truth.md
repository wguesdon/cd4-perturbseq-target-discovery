# HANDOFF #1 — Open Targets cross-check of the immunomodulator ground truth

**To run in:** Claude Science, with the official Open Targets MCP connector enabled.
**Log as:** CS1 in [`claude_tooling_log.md`](../claude_tooling_log.md).
**Returns:** a corrected CSV that replaces `resources/ground_truth/immunomodulator_targets.csv`.

## Why this handoff

The ground-truth table is currently hand-curated: the `approved` and `mechanism` columns are
asserted from domain knowledge, not sourced. The headline benchmark reports precision and
recall against this table, so every one of its rows is load-bearing. Open Targets carries
`maximumClinicalTrialPhase` and mechanism-of-action per target-drug pair, which turns the
assertion into a citation.

This is also the cleanest demonstration of the tooling model: Claude Code curated a candidate
set from the literature, Claude Science verifies it against an authoritative connector, and the
disagreements are the interesting part.

## Prompt to paste into Claude Science

> You have the Open Targets MCP connector enabled. I am building a ground-truth set of approved
> immunomodulator drug targets to benchmark a CRISPRi Perturb-seq target ranking in CD4+ T cells.
>
> The critical constraint is **direction of effect**. CRISPRi knockdown removes protein function,
> so a perturbation can only mimic an inhibitor, antagonist, neutralising antibody, or degrader.
> It cannot mimic an agonist, and it cannot mimic a drug whose mechanism requires the target
> protein to be present (for example tacrolimus binds FKBP1A to inhibit calcineurin, so FKBP1A
> knockdown abolishes the drug's effect rather than reproducing it).
>
> For each gene symbol in the attached table:
>
> 1. Query Open Targets for the target's known drugs. Report `maximumClinicalTrialPhase` across
>    all indications, and separately the maximum phase restricted to immune, inflammatory,
>    autoimmune, and allergic indications.
> 2. Report the mechanism of action strings Open Targets gives for each drug.
> 3. Classify the mechanism as one of: INHIBITOR, ANTAGONIST, NEUTRALISING, DEGRADER, AGONIST,
>    MODULATOR, OTHER.
> 4. Set `lof_mimics_drug = TRUE` only for INHIBITOR, ANTAGONIST, NEUTRALISING, DEGRADER.
> 5. Flag any row where your classification disagrees with the `include_as_positive` column I
>    supplied, and explain the disagreement in one sentence.
>
> Then do the reverse direction, which matters more:
>
> 6. Query Open Targets for **all** human targets with an approved drug (phase 4) whose indication
>    is autoimmune, inflammatory, or allergic disease, and whose mechanism is inhibitory,
>    antagonist, neutralising, or degrading.
> 7. Intersect that set with the 10,282 measured genes in the attached `measured_genes.txt` and
>    the perturbed gene list in `perturbed_genes.txt`.
> 8. Return the genes that pass all filters but are **missing** from my table. These are positives
>    I failed to curate, and each one is a false negative the benchmark would silently accept.
>
> Return one CSV with the same columns as the input plus: `ot_max_phase`,
> `ot_max_phase_immune`, `ot_mechanisms`, `ot_mechanism_class`, `ot_lof_mimics_drug`,
> `disagrees_with_curation`, `disagreement_reason`, `newly_added_by_open_targets`.
>
> Cite the Open Targets release version. Do not drop any of my rows, including the ones I marked
> `include_as_positive = FALSE`; I need to see whether you agree with the exclusions.

## Attachments to provide

- `resources/ground_truth/immunomodulator_targets.csv` (the curated table)
- `data/interim/measured_genes.txt` (written by `scripts/03_build_benchmark.py`)
- `data/interim/perturbed_genes.txt` (written by `scripts/03_build_benchmark.py`)

## What to do with the result

Save the returned CSV to `resources/ground_truth/immunomodulator_targets.open_targets.csv`,
commit it alongside the hand-curated original, and record the diff. The diff is a deliverable in
its own right: it shows the benchmark's positive set was independently verified rather than
tuned until the numbers looked good.

If Open Targets adds positives we missed, rerun the benchmark against both sets and report the
sensitivity of AUROC to the ground-truth definition. A ranking whose score swings on the choice
of gold standard is a ranking to distrust, and saying so is worth more than a clean number.
