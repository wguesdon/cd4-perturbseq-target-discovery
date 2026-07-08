# Claude Science task board

The queue of work to run in the Claude Science web UI, and the protocol for working across two
machines without stepping on each other.

> **To actually run CS1 and CS2, follow
> [`claude_life_science/RUN_CS1_AND_CS2.md`](../../claude_life_science/RUN_CS1_AND_CS2.md).**
> That is the operator runbook: startup commands, attachments, prompts, and where to save the
> answers. This file is the rationale and the board. **`HANDOFF_02`'s prompt is SUPERSEDED** and
> must not be pasted; the runbook carries its replacement.

Read this first on the Mac Mini after `git pull`.

---

## Machines and who owns what

| | **um890** (Ubuntu, this repo's compute node) | **Mac Mini** |
| --- | --- | --- |
| Runs | Claude Code | Claude Science web UI |
| Has the 16.8 GB `GWCD4i.DE_stats.h5ad` | yes | **no, and never needs it** |
| Writes | `src/`, `scripts/`, `results/`, `data/`, `docs/` | **only new files** under `docs/handoffs/results/` and `resources/ground_truth/` |
| Reads | everything | the committed artifacts listed below |

**The rule that prevents merge pain: the Mac Mini only ADDS files. It never edits an existing
one.** Results land as new files in `docs/handoffs/results/`. Claude Code transcribes them into
`docs/claude_tooling_log.md` back on the um890. That way `git pull --rebase` is always clean.

Before starting work on either machine: `git pull --rebase`. After: commit and push immediately.

## Mac Mini setup, once

```bash
git clone git@github.com:wguesdon/cd4-perturbseq-target-discovery.git
cd cd4-perturbseq-target-discovery
bash scripts/fetch_priors.sh     # data/external is gitignored; this refetches the 22 prior tables
```

Claude Science needs `socat` for sandbox networking, and therefore for any connector query.
Handoffs marked **network: no** below run fine without it.

## Everything Claude Science needs is committed

No handoff requires the effect matrix. These are the inputs, all small:

| File | Size | What it is |
| --- | --- | --- |
| `resources/handoff_inputs/cs2_blind_ranking.csv` | 752 KB | **CS2 attaches THIS.** The ranked table with our analytical decisions stripped out |
| `results/tables/risk_kill_naive_reversal.csv` | 752 KB | the full ranked table. Carries `z_l2_decile`, `n_cells_decile` and `matched_background`, which **give the answer away**. Never attach it to a blind run |
| `resources/handoff_inputs/rankable_genes.txt` | 40 KB | the 6,371 genes that survived perturbation QC |
| `resources/handoff_inputs/measured_genes.txt` | 64 KB | the 10,282 measured transcriptome genes |
| `resources/handoff_inputs/perturbed_genes.txt` | 73 KB | the 11,526 genes perturbed anywhere in the library |
| `resources/handoff_inputs/ground_truth_coverage.csv` | 2 KB | per-positive coverage and naive rank |
| `resources/ground_truth/immunomodulator_targets.csv` | 8 KB | the curated benchmark ground truth |

---

## Queue

| # | Handoff | Network | Blocks | Status |
| --- | --- | --- | --- | --- |
| CS1 | [Open Targets: can the benchmark be salvaged?](./HANDOFF_01_open_targets_ground_truth.md) | **yes**, needs `socat` + Open Targets MCP | the entire benchmark design | **READY, do this first** |
| CS2 | Blind adversarial re-analysis of the risk-kill result — prompt in the [runbook](../../claude_life_science/RUN_CS1_AND_CS2.md); rationale in [HANDOFF_02](./HANDOFF_02_reviewer_verify_riskkill.md) | no | nothing, run in parallel with CS1 | **READY, prompt REWRITTEN 07-08** |
| CS3 | Evidence cards for the safety-passing shortlist, every citation verified | yes | needs the null model on the um890 first | blocked |
| CS4 | Figures rendered from the committed tables, with provenance artifacts | no | needs the final ranking | blocked |
| CS5 | **Head to head: rerun a Claude Code analysis in Claude Science and compare** | no | needs a finished analysis to rerun | candidate, see below |

### CS5 — the controlled comparison

Take one analysis that Claude Code already ran, hand Claude Science the same committed inputs and
the same question, and let it work independently. Then compare. Not "both tools are great", but a
real diff: what did each catch, what did each miss, where did they disagree, and which was right.

The risk kill analysis is the natural candidate, because CS2 already sets it up. Claude Code drew
the matched background once. If Claude Science resamples and the conclusion shifts, that is a
concrete, filmable finding about the two tools, not a marketing line.

This is a strong answer to the Claude Use criterion, which is 25% of the score and asks for
creative use that surfaces surprising capabilities. Nobody scores points for using two tools. They
score points for showing what one tool caught that the other did not.

Decide the scope over the weekend. Do not let it displace CS1, which gates the benchmark.

### CS1 — the one that matters

The headline benchmark is currently failing. Not because the ranking is bad, but because only
**20** of 36 curated positives survive into the rankable set, and `JAK1` and `JAK3` were never
perturbed. AUROC is 0.542 with a 95% CI of [0.373, 0.707], which covers chance.

CS1 asks Open Targets for the ceiling: how many phase-4, immune-indication, inhibitory-mechanism
targets exist that are also among our 6,371 rankable genes? Above roughly 60 and the benchmark is
powered. Near 20 and we retire it as the headline, honestly, and validate the efficacy axis on the
Schmidt IL-2 screen instead.

### CS2 — the trust asset

**Rewritten 2026-07-08.** The original asked Claude Science to check a defect Claude Code had
already spotted: the cell-count-matched background drawn once with `seed=0`. Claude Code has since
reproduced that defect (22.4% of seeds give p ≥ 0.05), removed it, and in the process discovered
that cell count was the wrong covariate to match on in the first place. Asking Claude Science to
confirm a correction we have already made would prove nothing.

So CS2 is now a **blind** run. Claude Science gets the measurements with our analytical decisions
stripped out — no `z_l2_decile`, no `matched_background` — and is asked to identify the confound
itself, build its own control, and reach its own verdict. Both candidate covariates are present and
neither is labelled.

That makes a genuine disagreement possible, which is the entire point. Nobody scores for using two
tools. They score for showing what one caught that the other did not, and for being able to say
which was right. An independent agent that reaches "only tolerance is direction-specific" without
being told is worth more to the Demo criterion (30%, "findings you trust") than any figure we could
draw. An independent agent that reaches something else is worth more still.

The sealed comparison, and the rule against leaking our answer mid-session, are in the runbook.

---

## Returning results

Save each result as a new file, then commit and push from the Mac:

```
docs/handoffs/results/CS1_open_targets.md          # narrative + the counts
resources/ground_truth/immunomodulator_targets.open_targets.csv    # Task A output
resources/ground_truth/open_targets_additional_positives.csv       # Task B output
docs/handoffs/results/CS2_reviewer_verification.md # narrative + revised statistics
```

In each result file, record:

1. The Claude Science **provenance artifact id or link**. This is the audit trail for the
   Claude Use criterion, which is 25% of the score.
2. The connector and its **release version** (for example, the Open Targets release).
3. Anything Claude Science **corrected or contradicted**. Log the disagreements, not just the
   agreements. A handoff that only ever confirms what we already believed is not evidence of
   anything.

Then say so here, and Claude Code will pick it up on the next `git pull`.
