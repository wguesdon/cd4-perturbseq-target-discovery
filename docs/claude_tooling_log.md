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

## Claude Science handoffs

None yet. The first handoff will be logged here once Claude Science is installed and a task is
delegated. Template row:

| # | date | handoff task | handoff prompt ref | artifact id | result |
| --- | --- | --- | --- | --- | --- |
| CS1 | (pending) | (e.g., Open Targets druggability for top hits) | HANDOFF #1 | (artifact id) | (summary) |
