# Claude Science runbook: CS1 and CS2

**This is the operator runbook. You are the human bridge.** Claude Code cannot drive Claude Science:
`v0.1.16-dev` has no headless mode, and its only subcommands are `serve | open | url | status | logs
| stop | update | import`. Everything below you do by hand, in a browser.

**Canonical background lives in [`../docs/handoffs/`](../docs/handoffs/).** This file is the
*procedure*; those files are the *rationale*. Where they disagree, this file wins, because
`HANDOFF_02` was written before the corrections of 2026-07-08 and is now marked SUPERSEDED.

**Two rules that matter more than the prompts.**

1. **Never paste the Anthropic org UUID, or anything from a `*.local.md`, into a Claude Science
   session.** The repo goes public.
2. **CS2 is a blind test. Do not tell it our answer.** Not in the prompt, not in a follow-up, not by
   attaching the wrong CSV. Its value is entirely in whether it reaches our conclusion without help.
   Read [Sealing the blind](#sealing-the-blind) before you start.

---

## Before you start

```bash
cd "$(git rev-parse --show-toplevel)"
git pull --rebase
git status -sb                       # must be clean

claude-science --version             # expect 0.1.16-dev or later
claude-science serve                 # or: serve --detached --no-browser, then: claude-science url
```

The web UI is at **`127.0.0.1:8000`**. Claude Science needs `bubblewrap` and `socat` installed;
both are already on this box. `socat` is what gives the sandbox network access, so it is required
for CS1 and irrelevant to CS2.

| | CS1 | CS2 |
| --- | --- | --- |
| Network | **yes**, Open Targets MCP connector | **no**, runs offline |
| Role | independent **replication** of a committed number | **blind** adversarial re-analysis |
| Blocks | nothing; the benchmark decision is already made | nothing |
| Attachments | 3 files | 1 file |
| Do it | either order | either order |

Enable the Open Targets connector under **Customize → Permissions**. That is a UI action; there is
no CLI for it.

---

## CS1 — Open Targets: an independent REPLICATION, not a first pass

**Changed 2026-07-08.** Claude Code has already answered this, offline, from the pinned Open Targets
26.06 bulk release: `scripts/16_open_targets_benchmark.py`, committed at `4f8fb4e`. So CS1 is now a
contemporaneous A/B: the same question, a different tool, a different route to the data, against a
number that is already in git and cannot be quietly adjusted afterwards.

That is a much better use of the run. Nobody scores points for querying a database. They score for
showing that two independent tools agree, or for showing exactly where they do not.

### What Claude Code found, and where it nearly went wrong

**Count 3 = 38** under the primary definition, **53** at the upper bound. The threshold, 60, was fixed
before counting. The benchmark stays retired.

The first run returned **173** and declared the benchmark powered. Three separate traps:

1. MONDO classifies multiple myeloma, AML, DLBCL and CLL as descendants of *immune system disorder*,
   because they are malignancies of immune cells. Bortezomib therefore counted as an approved immune
   drug, and `PSMB5`, the HDACs and beta-tubulin became immunomodulator targets.
2. Open Targets annotates metformin's mechanism against all forty subunits of mitochondrial complex
   I, making `NDUFA4` an approved drug target for type 1 diabetes.
3. **The real error.** `clinical_target.maxClinicalStage` is the maximum stage a *drug-target pair*
   reached across all indications. Vorinostat is approved for cutaneous T-cell lymphoma and merely
   *trialled* in Crohn's disease, so joining that approval to an immune indication made `HDAC1` a
   benchmark positive. `clinical_indication` carries the stage per indication.

A validation harness (five positive controls, four negative) refused to print the ceiling until the
filters passed. It caught trap 3 by itself.

### The question for Claude Science

**Do not tell it any of the above.** Ask it the original question and see whether it walks into the
same three traps, or into different ones. If it returns a number near 173, you have a filmable
finding about how easy it is to get this wrong. If it returns 38–53 independently, the ceiling is
solid and the benchmark is retired on two independent authorities.

The prompt below asks for `maximumClinicalTrialPhase = 4`. **Leave that wording in.** In the 26.06
schema there is no phase 4; approval is encoded as the string `APPROVAL`. Whether Claude Science
notices, or silently returns nothing, is itself the test.

### Attach

| File | Size |
| --- | --- |
| `resources/ground_truth/immunomodulator_targets.csv` | 8 KB |
| `resources/handoff_inputs/rankable_genes.txt` | 40 KB |
| `resources/handoff_inputs/measured_genes.txt` | 64 KB |

### Prompt

The prompt is in **[`../docs/handoffs/HANDOFF_01_open_targets_ground_truth.md`](../docs/handoffs/HANDOFF_01_open_targets_ground_truth.md)**,
under "Prompt to paste into Claude Science". Paste it verbatim. Do not summarise it; the
direction-of-effect constraint in it is load-bearing and easy to lose.

The trap it exists to avoid: CRISPRi knockdown removes protein function, so a perturbation can
phenocopy an **inhibitor, antagonist, neutralising antibody, or degrader**, and never an
**agonist**. `NR3C1` is the sharpest case, because glucocorticoids are the most prescribed
immunosuppressants in the world and they are agonists. `FKBP1A` is the second, because tacrolimus
must *bind* it to inhibit calcineurin, so knocking it down abolishes the drug rather than mimicking
it.

### Save the answer to

```
docs/handoffs/results/CS1_open_targets.md    narrative, the three counts, and the DIFF vs Claude Code
```

Record the **Open Targets release version** it used and the **provenance artifact id**. Then diff its
counts against `results/tables/open_targets_benchmark_ceiling.csv`, and write the disagreement up.

Do **not** overwrite `resources/ground_truth/immunomodulator_targets.open_targets.csv` or
`open_targets_additional_positives.csv`. Claude Code already wrote those, and they are the artifact
being replicated. If Claude Science produces different ones, save them beside with a `.cs1` suffix.

---

## CS2 — Blind adversarial re-analysis of the risk-kill result

**Rewritten 2026-07-08. Do NOT use the prompt in `HANDOFF_02`; it asks Claude Science to verify
claims we have since retracted** (that the gate should be re-anchored on the IUIS immunodeficiency
list, and that cancer-cell essentiality is "flat" at p = 0.611). It also asks it to resample a
seed-0 background we have already resampled and removed. Running it would produce a confident
verification of a document that no longer exists.

### Why this version

Yesterday's central result rested on a background that was matched on **cell count** and drawn
**once**, with `seed = 0`. Both were defects. We have since fixed them, and in fixing them we
reached a conclusion that reverses part of our own headline.

We are not going to tell Claude Science any of that.

Instead we hand it the same measurements with **our analytical decisions stripped out**, state the
hypothesis, and ask it to choose its own confounder and its own control. If it independently lands
where we landed, that is corroboration worth having. If it lands somewhere else, that is a
disagreement worth filming. Either way it is evidence, which a confirmation of our own summary
would not be.

`docs/handoffs/README.md` puts it plainly: *"A handoff that only ever confirms what we already
believed is not evidence of anything."*

### Sealing the blind

The committed `results/tables/risk_kill_naive_reversal.csv` carries three columns that **give the
answer away in the header**: `z_l2_decile`, `n_cells_decile`, and `matched_background`. Those are
not measurements. They are a record of which covariate we decided to stratify on and which
background we decided to draw.

So **do not attach that file.** Attach the blinded one, which is generated by
`scripts/03_export_handoff_inputs.py` and asserts on its own output that no decision column
survived:

```
resources/handoff_inputs/cs2_blind_ranking.csv     752 KB, 6,371 rows, 12 columns
```

It keeps every measurement, including **both** candidate covariates (`n_cells_target` and `z_l2`),
because without `z_l2` the reviewer cannot compute effect magnitude at all: that needs full rows of
a 16.8 GB layer it will never have. Neither covariate is labelled as the answer.

Then, in this order:

1. Run the prompt below. Let it finish.
2. **Save its verdict verbatim, before you read it critically, and before you show it anything of
   ours.** That file is the sealed record.
3. Only then, open `docs/results/magnitude_matched_2026_07_08.md` and compare.
4. Log the diff in `docs/claude_tooling_log.md`, **including the disagreements**. Especially the
   disagreements.

If you find yourself typing "actually, we used z_l2" into the session, stop. The run is over. Start
a new session for any follow-up.

### Prompt to paste into Claude Science

> I am attaching `cs2_blind_ranking.csv`: 6,371 rows, one per CRISPRi perturbation of a single gene
> in primary human CD4+ T cells. Each perturbation was assayed in resting cells and in cells
> stimulated for 48 hours.
>
> Columns:
>
> - `target_contrast_gene_name` — the knocked-down gene.
> - `naive_suppression` — how strongly the knockdown suppresses a 32-gene inflammatory effector
>   program in stimulated cells. It is the negated mean z-score over those 32 genes. Higher means
>   stronger suppression. Negative values mean the knockdown *induced* the program.
> - `rank` — 1 is the strongest suppressor. Derived from `naive_suppression`.
> - `tolerance_suppression` — the same statistic computed over a separate 9-gene regulatory
>   tolerance module (FOXP3, IL10, CTLA4, PDCD1, LAG3, TIGIT and similar checkpoints). Higher means
>   the knockdown also suppressed tolerance, which would be a liability rather than a benefit.
> - `stim_de_genes`, `rest_de_genes` — counts of significantly differentially expressed genes in
>   the stimulated and resting conditions.
> - `stim_downstream`, `rest_downstream` — the same counts excluding the perturbed gene itself.
> - `n_cells_target` — the number of cells carrying a guide against this gene.
> - `z_l2` — the L2 norm of the perturbation's z-scores across all 10,282 measured genes in the
>   stimulated condition, excluding the perturbed gene's own column and both module gene sets.
> - `is_iei` — the gene's loss of function causes a human inborn error of immunity (IUIS 2024).
> - `is_core_essential` — Hart core-essential gene.
>
> **The hypothesis under test.** "The top 100 of this ranking is enriched for safety liabilities,
> and that enrichment is a real property of ranking by suppression rather than an artifact of
> something else."
>
> The four liabilities claimed are: `is_iei`, high `stim_downstream`, high `rest_downstream`, high
> `tolerance_suppression`.
>
> Work in the sandbox and show all code.
>
> **1. Find the confound yourself.** Something other than suppression could produce all four
> apparent enrichments at once. Identify what it is, from the data. `n_cells_target` and `z_l2` are
> both candidates and I am not telling you which, if either, is the right one to control for.
> Examine how each relates to the four liabilities and to `naive_suppression`. State which you will
> control for and **justify the choice quantitatively**. If you think both, or neither, say so.
>
> **2. Build the control you just argued for**, and test each of the four liabilities in the top 100
> against it. Use every available control row rather than a subsample, and say why subsampling would
> be a mistake here. Correct for multiple comparisons and name the correction.
>
> **3. The sign-flipped control.** `naive_suppression` is signed. Rank the perturbations by
> *induction* of the effector program instead, take that top 100, and rerun all four tests against
> the same control. Any liability that fires on **both** tails of the ranking is a property of the
> magnitude of a perturbation, not of its direction, and therefore says nothing about suppression.
> Report which of the four survive this and which do not.
>
> **4. Is any of it robust?** Everything above depends on how finely you stratified. Redo tests 2
> and 3 at several stratum counts. Report any liability whose significance call changes. A verdict
> that depends on an arbitrary bin count is not a verdict.
>
> **5. Can your test fail?** Draw 200 random shortlists of 100 and run test 2 on each. Report the
> fraction reaching p < 0.05 per liability. If it is far from 0.05, your machinery is broken and
> nothing above is trustworthy.
>
> **6. Be hostile.** You have the raw table. Look for what I have not asked about: selection
> effects, columns that cannot mean what their names suggest, missing values that are being read as
> zeros, any claim that outruns its evidence. Note that `rest_de_genes` is missing for some rows.
> Find out how many, find out whether they are missing at random, and say what that does to the
> resting-disruption liability.
>
> Then hand the whole analysis to the reviewer agent. Ask it to check the arithmetic, to check that
> each number matches the code that produced it, and to check that no stated conclusion outruns its
> evidence.
>
> **Finish with a verdict in four lines.**
> - Which single covariate, if any, is the real confound, and why.
> - Which of the four liabilities survive a control for it.
> - Which of the four are properties of direction, and which merely of magnitude.
> - Whether "the top of a suppression ranking is enriched for safety liabilities" is supported, and
>   if so, by which liability specifically.
>
> Give me the provenance artifact reference.

### What we are watching for

Do not read this section until CS2 has finished and you have saved its verdict. It is here so the
comparison is fair, not so it can be anticipated.

<details>
<summary>Our answer, sealed. Open after CS2 finishes.</summary>

We found, in `scripts/12_magnitude_matched.py` and `scripts/14_reversal_specificity.py`:

- **`Spearman(n_cells_target, stim_de_genes) = -0.243`.** Cell count is a *viability* readout, not a
  statistical-power proxy. More cells means fewer DE genes. Matching on it controls the confound
  backwards. Our original design did exactly that.
- **`Spearman(z_l2, |naive_suppression|) = +0.198`.** Effect magnitude and score magnitude are not
  the same variable, and which one you stratify on decides the answer.
- Stratified on `z_l2`, collateral DE survives but is `rho = 0.725` with the stratifier and its
  direction-specificity dies when the bins go from 10 to 20. Resting DE fires equally on the
  induction tail. **Only `tolerance_suppression` survives everything** (`rho = 0.069` with `z_l2`;
  inducers *induce* tolerance rather than suppressing it).
- The headline "18 of the naive top 20 are rejected by the safety gate" is worth **2.7 rejections**
  over a magnitude-matched 20 (15.3 of 20, P = 0.104). Split by axis: homeostasis 14 vs 13.6
  (P = 0.53); tolerance 12 vs 4.0 (P < 0.0001).
- `rest_de_genes` is missing for 51 of 6,371, three of them inside the top 100, including `IL2RB`
  at rank 4 — the single gene our own results doc used to motivate the whole safety gate. It is
  missing because those perturbations have **no resting row at all**, not because the count is zero.

**The interesting questions for the log.** Does it pick `z_l2` over `n_cells_target` unprompted?
Does it notice the `-0.243` sign? Does it catch `IL2RB`? Does it reach "only tolerance is
direction-specific", or does it stop at "the pillars survive"? Does it find something we missed?

</details>

### Save the answer to

```
docs/handoffs/results/CS2_blind_reanalysis.md
```

Include the **provenance artifact id**, the sealed verdict verbatim, and then a section headed
**Disagreements** listing every point where Claude Science and Claude Code differ, with who was
right and how you decided. If there are none, say that too, and say it as a weakness.

---

## After both runs

1. `git pull --rebase`, add the result files, commit, push. The Mac Mini only ever **adds** files
   under `docs/handoffs/results/` and `resources/ground_truth/`; it never edits an existing one.
   That is what keeps the rebase clean.
2. Fill the **disagreements column** in `docs/claude_tooling_log.md`. It is a submission-checklist
   item and it is currently empty.
3. CS1's count 3 is a **replication**. Claude Code got 38 primary / 53 upper bound and the benchmark
   is already retired in `reports/report.qmd` @sec-ceiling. If Claude Science lands near 173, note
   which of the three traps it fell into. If it lands near 38–53, say so: two tools, two routes, one
   ceiling.
4. If CS2 disagrees with us anywhere that matters, **put the disagreement in the demo video**. A
   project that shows two Claude tools reaching different conclusions, and then shows which was
   right and why, is making the exact argument the Demo criterion asks for. Nobody scores points
   for using two tools. They score points for showing what one caught that the other did not.

## Reproduce the blinded attachment

```bash
uv run python scripts/03_export_handoff_inputs.py
```

It writes `resources/handoff_inputs/cs2_blind_ranking.csv` and raises `AssertionError` if any
decision column survives into it.
