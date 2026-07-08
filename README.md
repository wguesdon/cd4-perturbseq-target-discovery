# CD4+ T cell Perturb-seq target triage

Researcher Track submission for the **2026 Built with Claude: Life Sciences** hackathon.

> All work in this repository was started from scratch on 2026-07-07 and carried out during the
> event. The starting question and the public dataset predate the hackathon. The analysis does not.

## The question

Which single gene knockdowns suppress the inflammatory program of a stimulated human CD4+ T cell
**without** disabling the cell, collapsing its resting transcriptome, or destroying its tolerance
machinery?

That last clause is the whole project. Ranking perturbations by how well they reverse a disease
signature is a solved, released capability. On its own it nominates targets that work and would
leave a patient immunodeficient. **Reversal is not enough.** The contribution is the safety gate
that separates the two, and an honest account of where it fails.

## What we have found so far

Every number below is reproducible from this repository. Days 2 to 5 are still in progress.

**The naive ranking is toxic, and not for the reason we predicted.** A well-QC'd score that ranks
perturbations purely by suppression of the effector program produces a top 100 that is enriched
for genes causing human inborn errors of immunity. It is **not** enriched for cancer cell
essential genes. The 31 rankable Hart core essentials sit at median rank 3148 of 6371, p equals
0.611. A viability penalty anchored on DepMap, which is what our own plan originally specified,
would have rejected almost nothing.

The contaminants are of two kinds. Immune essential signalling machinery: `STAT5B`, `VAV1`,
`IL2RB`, `CD3G`, `CD247`, `LCK`. And global transcription machinery that collapses the resting
transcriptome: `NSD1`, `TADA2B`, `CXXC1`, `WDR82`, `USP22`. All 25 of the naive top 25 are
rejected by the safety gate. See [`docs/results/risk_kill_2026_07_08.md`](docs/results/risk_kill_2026_07_08.md).

**An independent screen confirms the toxic hits really work.** In Schmidt and Steinhart 2022, a
separate lab's genome wide CRISPRi screen with a protein level readout, `VAV1`, `CD3G` and `CD247`
knockdowns genuinely reduce IL-2 production. They are effective. They are also immunodeficiency
genes. That single table is the argument.

**The efficacy axis is real.** Our transcriptome derived suppression score reaches AUROC 0.702,
95% CI [0.591, 0.814], against the 33 significant IL-2 reducing hits in that independent screen.

**The drug target recovery benchmark is underpowered, and we say so.** Of 36 curated positives,
32 are perturbed, 30 are measured, and only 20 survive QC. `JAK1` and `JAK3` were never perturbed,
so tofacitinib's targets cannot be recovered at all. AUROC on those 20 is 0.542, 95% CI
[0.373, 0.707], which covers chance. Whether this benchmark survives is being decided by an Open
Targets query, not by us choosing a number we like.

## Honest limitations

- The screen stimulates through the TCR with no polarising cytokines. Cytokine signalling targets
  are therefore close to invisible. `JAK2` ranks 5392, `TYK2` 5618, `S1PR1` 5993, `IL4R` 6047 of
  6371. This bounds recall from above and is a property of the dataset, not of the method.
- `CD3E` and `CD3G` are simultaneously approved drug targets and inborn errors of immunity. The
  safety gate demotes the genes the benchmark counts as positives. That is correct behaviour,
  since muromonab was withdrawn, but it means drug recovery validates the **efficacy axis** and
  cannot validate the **window score**. The two are benchmarked separately.
- The released differential expression is pseudobulk, not distribution aware. We inherit that.
- The library is single gene. There are no combinatorial perturbations.
- Th2 and Th17 readouts are proxies, because the culture is not polarised.

## Data

Genome scale CD4+ T cell CRISPRi Perturb-seq, Marson lab, hosted on CZI Virtual Cell Models.

- Dataset: https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
- Reference analysis code: https://github.com/emdann/GWT_perturbseq_analysis_2025
- Preprint: https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1

We use the released per perturbation effect matrix `GWCD4i.DE_stats.h5ad`, 33,983 perturbation by
condition pairs across 10,282 measured genes. We do not reprocess the 22 million raw cells. The
matrix is 16.8 GB and is not committed. See [`data/README.md`](data/README.md).

Prior gene lists and reference screens are fetched from the authors' analysis repository by
[`scripts/fetch_priors.sh`](scripts/fetch_priors.sh), which pulls 22 small tables rather than
cloning 1.66 GB.

## How Claude Code and Claude Science are used

This is a Built with Claude entry, so how the tools are wielded is itself a graded artifact. The
division of labour is deliberate and follows what each tool is uniquely good at. The full audit
trail, one row per action, is [`docs/claude_tooling_log.md`](docs/claude_tooling_log.md). The queue
is [`docs/handoffs/README.md`](docs/handoffs/README.md).

**Claude Code orchestrates and builds.** It plans, writes and runs the pipeline, owns the
repository and git history, and runs adversarial multi agent audits against its own results. The
literature review was a six agent fan out. The risk kill result above was then attacked by a six
lens hostile audit, each finding independently refuted or reproduced by two further agents, before
any of it was believed.

**Claude Science verifies, sources, and renders.** It is a browser workbench with official
connectors, a reviewer agent, sandboxed kernels, and provenance artifacts. It is used for the four
things Claude Code cannot do as well:

| # | Handoff | What only Claude Science gives us |
| --- | --- | --- |
| CS1 | Open Targets: how many phase 4 immune inhibitory targets are in our 6,371 rankable genes | An authoritative, versioned answer to a ground truth question we would otherwise be asserting from memory |
| CS2 | Reviewer agent independently recomputes the risk kill statistics | An outside party trying to break our central claim, including resampling a background we drew only once |
| CS3 | Evidence cards for the surviving shortlist, with PubMed citation verification | Every claim on a card checked against a source by the reviewer agent, not asserted |
| CS4 | Figures rendered from the committed tables, with provenance | An artifact proving each figure matches the code that produced it |

Two principles govern this. Claude Science is given the tasks that could **break** our result, not
only the tasks that decorate it. CS2 exists specifically to attack the headline. And every handoff
runs on small committed artifacts, never on the 16.8 GB matrix, so any reader can rerun them.

Where Claude Science corrects or contradicts Claude Code, the disagreement is logged. A handoff
that only ever confirms what we already believed is not evidence of anything.

## Working across two machines

Compute and the effect matrix live on one machine. Claude Science runs in a browser on another.
Git is the only channel between them.

| | Compute node | Claude Science node |
| --- | --- | --- |
| Runs | Claude Code | Claude Science web UI |
| Has the 16.8 GB matrix | yes | no, and never needs it |
| Writes | `src/`, `scripts/`, `results/`, `data/`, `docs/` | only **new** files under `docs/handoffs/results/` and `resources/ground_truth/` |

The rule that keeps them from colliding: the Claude Science node only adds files, never edits an
existing one. Rebase therefore stays clean. Everything a handoff needs is committed and under
600 KB, listed in [`docs/handoffs/README.md`](docs/handoffs/README.md).

There is deliberately no cloud compute. One layer of the matrix is 2.8 GB and the rows we read are
0.93 GB, so the entire analysis runs locally in minutes. Nothing in this project trains a model.
Foundation and perturbation models are kept off the critical path, because the literature shows
they do not beat simple baselines here.

## Method

1. Build the activation program from two independent sources, a curated effector list and an
   external bulk RNA-seq stimulation contrast. Separate the **effector module**, which is the
   objective, from the **tolerance module**, which is a penalty. `FOXP3`, `IL10`, `IKZF2` and the
   co-inhibitory checkpoints are stimulation induced too, so a naive objective wrongly rewards
   destroying them.
2. Score effector suppression in stimulated cells, calibrated against a within condition null.
3. Gate on the three axes the data actually supports: immune essentiality from the IUIS inborn
   errors of immunity list, resting transcriptome disruption, and tolerance preservation.
4. Validate the efficacy axis against an independent CRISPRi screen with a protein readout.
   Validate the window score on therapeutic index, separately, because the two are not the same
   question.
5. Annotate survivors with druggable class and Open Targets tractability.
6. Ship a ranked table and sourced evidence cards.

Superseded planning documents are kept rather than deleted, so the reasoning is auditable. See
[`docs/strategy_2026_07_07.md`](docs/strategy_2026_07_07.md) and the corrections recorded in
[`docs/results/risk_kill_2026_07_08.md`](docs/results/risk_kill_2026_07_08.md).

## Reproducibility

Managed with [uv](https://docs.astral.sh/uv/). Python 3.13.

```bash
git clone git@github.com:wguesdon/cd4-perturbseq-target-discovery.git
cd cd4-perturbseq-target-discovery
uv sync
bash scripts/fetch_priors.sh          # 22 small prior tables, no large clone

# Optional, 16.8 GB. Only needed to regenerate the ranking from scratch.
aws s3 cp s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad \
  data/raw/ --no-sign-request

uv run python scripts/00_inspect_de_stats.py
uv run python scripts/01_build_activation_program.py
uv run python scripts/02_risk_kill_reversal.py      # exits non-zero if the headline claim fails
uv run python scripts/03_export_handoff_inputs.py
```

`scripts/02_risk_kill_reversal.py` is written so it can fail. If the naive ranking is no more
toxic than a cell count matched background, it exits non-zero and the project's headline is wrong.

## Repository layout

```
data/                 raw, interim, external (contents gitignored)
docs/                 strategy, literature, results, and the Claude Science handoff queue
resources/            committed ground truth and the small handoff input bundle
scripts/              numbered, runnable pipeline steps
src/cd4_perturbseq/   reusable package code
results/              figures and tables produced by the analysis
```

## License

MIT. See [`LICENSE`](LICENSE). The source dataset is MIT plus the CZI Acceptable Use Policy. The
preprint is CC BY.
