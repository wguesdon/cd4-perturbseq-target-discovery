# Judging criteria and track

The track we are on and the rubric to score every deliverable against. The full official brief is
in [`hackathon_brief.md`](./hackathon_brief.md); the detailed rubric alignment is in
[`strategy_2026_07_07.md`](./strategy_2026_07_07.md), section 5.

## Track: Researcher — "Build From the Bench"

Using Claude Science, start from a biological question, find the datasets and tools to answer it,
and submit something discrete: a finding, a trained model, or an analysis others can reproduce.
Show how Claude Science got you there.

## Criteria (Stage 1, asynchronous, 2026-07-14 to 07-15)

Standardized scoring. Weight in brackets. For each: the official ask, then how we target it.

1. Impact [25%]. Official: real-world potential; for Researcher projects, is this a finding or
   analysis others can build on, and does it fit the track problem statement.
   How we target it: a validated, reproducible target-triage method plus a novel, druggable,
   genetically supported target shortlist. Fits the CD4 Perturb-seq example directly.

2. Claude Use [25%]. Official: how creatively Claude Code was used; going beyond a basic
   application; surfacing surprising capabilities.
   How we target it: Claude Code orchestrates (multi-agent literature review, the benchmark
   harness, evidence-card generation); Claude Science analyzes over data with official connectors
   and the reviewer agent. Every handoff is logged in `claude_tooling_log.md`. See
   `claude_tooling_strategy_2026_07_08.md`.

3. Depth and Execution [20%]. Official: pushed past the first idea; sound, refined engineering;
   real craft, not a quick hack.
   How we target it: the pivot from a naive gene list to a benchmark-validated, safety-gated
   triage; the reversal-is-not-enough head-to-head; honest limitations; reproducible provenance.

4. Demo [30%, the largest weight]. Official: a working, compelling demo; findings you trust;
   genuinely cool to watch.
   How we target it: the demo spine, which is the benchmark figure recovering approved drugs, the
   naive-versus-safe head-to-head, and two or three sourced evidence cards. The reviewer agent and
   provenance artifacts are the trust signals.

## Stage 2: final live round

2026-07-16, 12:00 PM ET. The top teams' pre-recorded 3-minute demos play; winners are announced at
the 1:30 PM ET closing ceremony.

## Submission (due 2026-07-13, 9:00 PM ET)

1. Demo video, 3 minutes maximum.
2. Public open-source repository, notebook, or research write-up.
3. Written summary, 100 to 200 words.

## Scoring reminder

Demo is 30% and Claude Use is 25%, so more than half the score is the demo and the tooling story.
Budget real time for the video, and make the Claude Code and Claude Science usage visible.
