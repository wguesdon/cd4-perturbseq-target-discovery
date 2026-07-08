# Claude tooling log

The audit trail of how this project used Claude Code and Claude Science. One row per action.
Append a row for every Claude Code milestone and every Claude Science handoff. See
[`claude_tooling_strategy_2026_07_08.md`](./claude_tooling_strategy_2026_07_08.md) for the model.

Columns: # | date | tool | task | prompt or commit | artifact or output | result summary.

## Claude Code

| # | date | task | prompt or commit | output | result |
| --- | --- | --- | --- | --- | --- |
| CC1 | 2026-07-07 | Scaffold the private repo (uv, structure, MIT) | commit b783719 | repo skeleton | project created and pushed |
| CC2 | 2026-07-07 | Plan-mode research planning for preprint acquisition | plan file | approved plan | acquisition and comparison plan set |
| CC3 | 2026-07-07 | Six-agent literature fan-out across the idea's facets | 6 subagents | `literature/lit_review_2026_07_07.md`, `methods_toolkit.md` | cited review; toolkit; novelty white-space |
| CC4 | 2026-07-07 | Extract the source paper's future work | subagent | `literature/source_paper_future_work.md` | future work and limitations captured |
| CC5 | 2026-07-08 | Deep-read the paper figure-by-figure + reconstruct the repo | 2 subagents | strategy inputs | pert2state, released resources, zero druggability language |
| CC6 | 2026-07-07 | Competition and prior-work check | subagent | strategy inputs | space open; Weinstock 2024 is GRN prior art |
| CC7 | 2026-07-07 | Synthesize the detailed strategy | commit 4852699 | `strategy_2026_07_07.md` | benchmark-validated therapeutic-window triage |
| CC8 | 2026-07-08 | Solo adjustment | commit 96dc69e | strategy update | reversal-is-not-enough differentiation |
| CC9 | 2026-07-08 | Research Claude Science and the ecosystem | 2 subagents | `claude_tooling_strategy_2026_07_08.md` | tooling strategy and connector map |
| CC10 | 2026-07-08 | Audit the released data before trusting the plan | S3 + readme recon | — | found the suppl CSV is a Dec-2025 obs snapshot missing the confidence columns; the May-2026 h5ad is required |
| CC11 | 2026-07-08 | Locate the bundled priors without a 1.66 GB clone | `scripts/fetch_priors.sh` | 22 prior tables | found `core_essentials_hart.tsv`, the IUIS inborn-errors-of-immunity list, druggable-class lists, and the Schmidt CD4 IL-2 CRISPRi screen |
| CC12 | 2026-07-08 | Define the window objective | commit d5e6530 | `src/cd4_perturbseq/programs.py` | separated the effector module (objective) from the tolerance module (penalty); FOXP3/IL10/checkpoints are stim-induced too |
| CC13 | 2026-07-08 | Curate the benchmark ground truth | commit 4c94e42 | `resources/ground_truth/immunomodulator_targets.csv` | direction-of-effect filter: excluded NR3C1 and IL2 (agonists), FKBP1A and PPIA (drug-binding chaperones), and non-T-cell genes |
| CC14 | 2026-07-08 | Steelman and de-confound the risk-kill test | `scripts/02_risk_kill_reversal.py` | matched-background test | naive baseline gets full QC; enrichment retested against a cell-count-matched background so power cannot fake the result |

## Claude Science handoffs

Queue and two-machine protocol: [`handoffs/README.md`](./handoffs/README.md).
Claude Science runs on the Mac Mini; it never needs the 16.8 GB effect matrix.

| # | date | handoff task | handoff prompt ref | artifact id | result |
| --- | --- | --- | --- | --- | --- |
| CS1 | (ready) | Open Targets: how many phase-4 immune inhibitory targets are in our 6,371 rankable genes? Decides whether the benchmark is powered. | [HANDOFF #1](./handoffs/HANDOFF_01_open_targets_ground_truth.md) | (pending) | (pending) |
| CS2 | (ready) | Reviewer agent independently recomputes the risk-kill statistics and resamples the single-draw matched background | [HANDOFF #2](./handoffs/HANDOFF_02_reviewer_verify_riskkill.md) | (pending) | (pending) |
| CS3 | (blocked) | Evidence cards for the safety-passing shortlist | (pending null model) | | |
| CS4 | (blocked) | Tractability and degrader handle for the final shortlist | (pending CS3) | | |
