# HANDOFF #2 — Reviewer agent: independently verify the risk-kill result

**Run in:** Claude Science web UI, Python sandbox plus the reviewer agent.
**Network:** not required. No connector, no `socat`. Runs offline.
**Log as:** CS2 in [`claude_tooling_log.md`](../claude_tooling_log.md).
**Machine:** any. One committed 558 KB CSV is the only input.

## Why this handoff exists

`docs/results/risk_kill_2026_07_08.md` claims the project's central empirical result. The whole
entry rests on it. It was produced by Claude Code in a single pass, and Claude Code then checked
its own work, which is not worth much.

This handoff asks Claude Science to recompute it from the raw ranked table, with no access to the
code that produced the claims, and to say plainly whether they hold.

Claude Code already suspects one defect and wants it examined from the outside rather than
excused from the inside: the cell-count-matched background is sampled **once**, with a fixed seed
of 0. The reported IEI enrichment, odds ratio 3.91 and p=0.012, rests on that single draw. If it
does not survive resampling, we need to know now and not on stage.

## Input

`results/tables/risk_kill_naive_reversal.csv`, 6,371 rows, one per rankable perturbation.

Columns: `target_contrast_gene_name`, `naive_suppression`, `tolerance_suppression`,
`stim_de_genes`, `n_cells_target`, `rest_de_genes`, `is_core_essential`, `is_iei`, `rank`,
`in_top`, `matched_background`.

`rank` is ascending by descending `naive_suppression`. `in_top` marks `rank <= 100`.
`matched_background` marks the single cell-count-decile-matched background sample that was drawn.

## The claims to verify

From `docs/results/risk_kill_2026_07_08.md`, against the matched background:

| Claim | Reported |
| --- | --- |
| IEI genes enriched in the top 100 | 14.0% vs 4.0%, Fisher OR 3.91, p=0.012 |
| Collateral DE genes in Stim, top vs background | median 356.5 vs 5.5, MWU p=2.4e-17 |
| DE genes at Rest, top vs background | median 70.0 vs 4.0, MWU p=2.6e-09 |
| Tolerance-module suppression, top vs background | median 1.09 vs 0.21, MWU p=8.1e-14 |
| Score is not driven by power | Spearman(naive_suppression, n_cells_target) = -0.109 |
| Core-essentiality is flat | 31 rankable essentials, median rank 3148 of 6371, MWU p=0.611 |

Overall conclusion drawn: a naive reversal ranking is toxic, but **not** because it favours
common-essential genes. The safety gate must therefore be re-anchored on immune essentiality
(IUIS inborn errors of immunity), resting-transcriptome disruption, and tolerance preservation,
rather than on DepMap or Hart cancer-cell essentiality.

## Prompt to paste into Claude Science

> I am attaching `risk_kill_naive_reversal.csv`, 6,371 rows, one per CRISPRi perturbation of a
> gene in primary human CD4+ T cells, ranked by how strongly the knockdown suppresses an
> inflammatory effector program in stimulated cells (`naive_suppression`, higher means stronger
> suppression; `rank` 1 is the strongest).
>
> Other columns: `is_iei` marks genes whose loss of function causes a human inborn error of
> immunity; `is_core_essential` marks Hart core-essential genes; `stim_de_genes` and
> `rest_de_genes` are the counts of significantly differentially expressed genes in the
> stimulated and resting conditions; `tolerance_suppression` is how strongly the knockdown also
> suppresses a tolerance module (FOXP3, IL10, CTLA4 and other checkpoints); `n_cells_target` is
> the number of cells carrying the guide, so it is a proxy for statistical power;
> `matched_background` marks one background sample matched to the top 100 on `n_cells_target`
> decile.
>
> The hypothesis under test is: "the top of this ranking is enriched for safety liabilities, and
> that enrichment is not an artifact of statistical power."
>
> Do the following in the sandbox. Show your code.
>
> 1. Reproduce, exactly, comparing the top 100 (`in_top`) against the supplied
>    `matched_background`: the Fisher exact odds ratio and one-sided p for `is_iei` and for
>    `is_core_essential`; and one-sided Mann-Whitney U p-values for `stim_de_genes`,
>    `rest_de_genes`, and `tolerance_suppression`. Report every number you get, including where
>    you disagree with the table above.
>
> 2. **The matched background was drawn once, with a fixed random seed.** Redraw it 1,000 times.
>    For each draw, sample from each `n_cells_target` decile as many background genes as there are
>    top-100 members in that decile, excluding the top 100 themselves. Report, for the IEI Fisher
>    test, the full distribution of the odds ratio and the p-value across those 1,000 draws:
>    median, 2.5th and 97.5th percentiles, and the fraction of draws with p < 0.05. Do the same
>    for the three Mann-Whitney tests. State clearly whether each conclusion is robust to
>    resampling or whether it was an artifact of one lucky draw.
>
> 3. Five tests were run against two backgrounds. Apply a Benjamini-Hochberg correction across
>    all of them and report which survive at FDR 0.05.
>
> 4. The claim "the naive score is indifferent to core essentiality" accepts a null hypothesis
>    (p=0.611) on only 31 rankable essential genes. Do a power analysis. With 31 positives out of
>    6,371, what effect size, expressed as a shift in median rank, would this Mann-Whitney test
>    have had 80% power to detect? Then state what can and cannot be concluded. Specifically:
>    can we rule out a modest enrichment, or only a strong one?
>
> 5. Adversarial check on confounding. `n_cells_target` is a power proxy. Both `stim_de_genes`
>    and `rest_de_genes` are counts of SIGNIFICANT genes, so they also rise with power. Test
>    whether the reported enrichment of `stim_de_genes` and `rest_de_genes` in the top 100
>    survives after conditioning on `n_cells_target`, for example by partial Spearman correlation,
>    or by logistic regression of `in_top` on the disruption measure with `log(n_cells_target)` as
>    a covariate. Note that Spearman(naive_suppression, n_cells_target) is reported as -0.109,
>    which is negative. Explain what that implies and whether it strengthens or weakens the
>    original argument.
>
> 6. Look for anything else wrong that I have not asked about. You have the raw table. Be hostile.
>
> Then hand the whole analysis to the reviewer agent. Ask it to check the arithmetic, check that
> each figure matches the code that generated it, and check that no stated conclusion outruns its
> evidence.
>
> Finish with a verdict in three lines:
> - Which claims SURVIVE independent recomputation and resampling.
> - Which claims are WEAKENED or must be restated, and how they should be restated.
> - Whether the overall conclusion, that the safety gate must be re-anchored on immune
>   essentiality rather than cancer-cell essentiality, is supported.
>
> Give me the provenance artifact reference.

## Attachment

`results/tables/risk_kill_naive_reversal.csv` (558 KB, committed)

## What we do with the answer

Save to `docs/handoffs/results/CS2_reviewer_verification.md`, including the provenance artifact
reference and the reviewer agent's findings.

If resampling weakens the IEI enrichment, `docs/results/risk_kill_2026_07_08.md` gets corrected,
not quietly patched. The correction, with the before and after numbers, is a stronger artifact
than the original claim was, and it goes in the demo. A project that catches its own error in
front of the judges is making the exact argument the Demo criterion asks for: findings you trust.

If the enrichment survives 1,000 redraws, we say so with the interval, and the headline is
sturdier than a single p-value ever made it.
