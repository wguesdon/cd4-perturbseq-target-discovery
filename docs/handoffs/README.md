# Claude Science task board

The queue of work to run in the Claude Science web UI, and the protocol for working across two
machines without stepping on each other.

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
| `results/tables/risk_kill_naive_reversal.csv` | 558 KB | the full ranked table, 6,371 perturbations, all scores and flags |
| `resources/handoff_inputs/rankable_genes.txt` | 40 KB | the 6,371 genes that survived perturbation QC |
| `resources/handoff_inputs/measured_genes.txt` | 65 KB | the 10,282 measured transcriptome genes |
| `resources/handoff_inputs/perturbed_genes.txt` | 73 KB | the 11,526 genes perturbed anywhere in the library |
| `resources/handoff_inputs/ground_truth_coverage.csv` | 2 KB | per-positive coverage and naive rank |
| `resources/ground_truth/immunomodulator_targets.csv` | 6 KB | the curated benchmark ground truth |

---

## Queue

| # | Handoff | Network | Blocks | Status |
| --- | --- | --- | --- | --- |
| CS1 | [Open Targets: can the benchmark be salvaged?](./HANDOFF_01_open_targets_ground_truth.md) | **yes**, needs `socat` + Open Targets MCP | the entire benchmark design | **READY, do this first** |
| CS2 | [Reviewer agent: independently verify the risk-kill result](./HANDOFF_02_reviewer_verify_riskkill.md) | no | nothing, run in parallel with CS1 | **READY** |
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

Claude Science's reviewer agent independently recomputes the risk-kill statistics from the
committed ranked table. This is not a formality. Claude Code already knows of one defect it wants
checked from the outside: the cell-count-matched background is drawn **once**, with `seed=0`, and
the IEI enrichment p=0.012 rests on that single draw.

An independent verification that finds and quantifies that, then reports what survives, is worth
more to the Demo criterion (30%, "findings you trust") than any figure we could draw.

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
