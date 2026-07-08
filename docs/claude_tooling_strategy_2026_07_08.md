# Claude Code + Claude Science: usage and showcase strategy

Date: 2026-07-08. Track: Researcher, 2026 Built with Claude: Life Sciences.

Purpose: use both Claude Code and Claude Science deliberately, and make that use visible. Claude
Use is 25% of the score, and this is the Built with Claude hackathon, so how we wield the tools
is itself a graded artifact. This document sets the operating model, the division of labor, a
task-by-task mapping, a setup checklist, and how to showcase both in the demo and writeup.

Companion: the running audit trail of every Claude Code and Claude Science action is
[`claude_tooling_log.md`](./claude_tooling_log.md).

---

## 1. The two products (accurate, brief)

- Claude Science: an AI workbench for scientists, launched 2026-06-30 (beta). A local desktop app
  that runs the same Claude models (Opus 4.8) with a Python, R, and shell sandbox and persistent
  kernels. It has a coordinating agent, specialist agents, and a reviewer agent that checks
  citations, calculations, and that figures match their code. It ships 60+ curated skills and
  connectors, produces versioned artifacts with full provenance (code, environment, and message
  history), and has scientific renderers (3D structures, genome tracks, chemical structures) that
  Claude Code lacks. It runs natively on Linux x64, so it runs directly on this workstation, no
  WSL. It needs a paid Claude account and no API key.
- Claude for Life Sciences: the connector and skill layer, plus a Claude Code marketplace
  (`anthropics/life-sciences`). The bridge between the two products is the shared Agent Skills and
  MCP standard, so a pipeline built in Claude Code can be saved as a Claude Science skill, and the
  same connectors load in both.

## 2. Operating model: Claude Code orchestrates, Claude Science executes delegated tasks

Claude Code (this session) is the orchestrator. It plans, writes and runs code, drives the repo,
spawns research subagents, and decides when a task genuinely suits Claude Science. Claude Code
cannot drive the Claude Science desktop app directly, so the integration is a human bridge, and
that bridge is also the tracking mechanism.

The loop:

1. Claude Code emits a self-contained handoff prompt block (template in section 8).
2. You paste it into Claude Science, run it, and let the reviewer agent check it.
3. You paste the result and the provenance artifact reference back into Claude Code.
4. Claude Code ingests the result, logs the round-trip in `claude_tooling_log.md`, and continues.

Every handoff is one row in the log. The log is the evidence of deliberate dual-tool use.

When to hand off to Claude Science (its unique strengths):
- Sandboxed heavy analysis over the DE matrices, with persistent kernels and provenance.
- Queries through its official connectors (Open Targets, ChEMBL, GEO, ClinVar, Reactome, UniProt,
  PDB) and single-cell skills.
- The reviewer agent, to check the benchmark math and every evidence-card citation.
- Scientific renderers, for chemical-structure and genome-track figures.
- Manuscript and figure drafting for the writeup.

What Claude Code keeps: orchestration, the repo and the reusable pipeline, the benchmark harness,
the agentic literature fan-out, the evidence-card generation logic, git and reproducibility, and
any MCP queries we choose to run here.

## 3. Division of labor

| Concern | Claude Code | Claude Science |
| --- | --- | --- |
| Role | Orchestrator and builder | Analyst over data, artifact producer |
| Best at | Repo, pipeline, benchmark harness, multi-agent research, git, evidence-card logic | Sandboxed Python/R analysis, connector queries, reviewer QC, provenance, renderers |
| Output | The shippable, reproducible tool | Auditable figures, tables, and checked results |
| Reproducibility | Git history, one-command rerun | Versioned artifacts with code, env, and history |
| Trust | Adversarial verification via subagents | Reviewer agent checks citations and numbers |

## 4. Task-by-task mapping

| Our task | Claude Code | Claude Science (connector or skill) |
| --- | --- | --- |
| Literature synthesis and novelty | Multi-agent fan-out (already done) | PubMed and Wiley connectors; reviewer agent to check citations |
| Benchmark ground truth (known drugs) | Curate and version the target list | Open Targets MCP (known drugs), ChEMBL (max_phase, mechanism) |
| Therapeutic-window scoring | Build the scorer and the pipeline | Sandboxed run over the DE matrices; single-cell skills |
| Drug-target-recovery benchmark | Write the harness, AUROC and PR | Reviewer agent verifies the metrics and the figure |
| Druggability and degrader axis | Encode the scoring logic | Open Targets tractability and PROTAC buckets; UniProt, PDB, ChEMBL; structure renderers |
| Genetics annotation | Merge into the target table | Open Targets genetic evidence; ClinVar, Reactome |
| Safety axis | DepMap and Hart logic in code | Sandboxed analysis; provenance on the safety figure |
| Evidence cards | Generate the card template and content | Provenance artifacts and reviewer-checked citations; BioRender figures |
| Reproducible pipeline | The repo, one-command rerun | Save the pipeline as a reusable Claude Science skill |

## 5. Setup checklist

1. Install Claude Science on Linux. Verify the exact command on the official get-started page,
   then run the installer and `claude-science serve`. Use a personal Claude account, since the
   hackathon requires personal-account credit eligibility.
2. Enable the connectors we actually need first: Open Targets MCP (the single most on-target
   official connector, covering druggability, genetic evidence, and known drugs), PubMed, and the
   built-in ChEMBL, GEO, ClinVar, Reactome, UniProt, and PDB databases.
3. Enable the single-cell skills: `single-cell-rna-qc` and `scvi-tools`.
4. In Claude Code, add the `anthropics/life-sciences` marketplace. Optionally add the community
   BioMCP server for a cheap variant, genetics, literature, and ChEMBL pivot layer. Prefer the
   official Open Targets MCP for anything we showcase.

## 6. How to showcase both

The demo (30%) and Claude Use (25%) reward a visible, deliberate dual-tool story. The spine:

- Show Claude Code as the conductor: it planned the work, ran a multi-agent literature review, and
  built the benchmark harness and the reusable scorer.
- Show one clean Claude Science handoff on camera: the reviewer agent verifying the benchmark
  numbers, or Open Targets returning druggability and genetics for the top hits, ending in a
  provenance artifact.
- Close on the evidence cards, each claim traceable, and the head-to-head that makes the reversal
  is not enough point.

Trust signals to name explicitly: the reviewer agent and the provenance artifacts are why a judge
should trust the findings. That directly answers the Demo criterion, findings you trust.

Reference pattern: the Claude Science launch demo reasoned over data and autonomously nominated
drug candidates, then produced an artifact. We mirror that shape: reason over the Perturb-seq
atlas, then produce validated, safety-gated target artifacts.

The tracker log is itself a showcase item. It shows exactly where each tool was used.

## 7. Honest framing and caveats

- Perturb-seq and drug-target discovery are not advertised, named Claude Science workflows. They
  are in scope through the Python and R sandbox, the single-cell skills, and the GEO and ChEMBL
  connectors, but we should describe them as our composition, not a built-in feature.
- Some connectors are community, not official. Open Targets MCP and PubMed are official. BioMCP,
  ChEMBL-MCP, and similar are community. Prefer official ones for anything on camera.
- The human bridge is a manual relay. Keep each handoff small and self-contained, and log it
  honestly, including where Claude Science corrected or improved a result.
- Do not overstate. The credibility comes from the reviewer agent and provenance, so lean on them
  rather than on claims.

## 8. Handoff prompt template (copy-ready)

```
=== CLAUDE SCIENCE HANDOFF #N — <short name> ===
Goal: <one line>
Context: <one or two lines Claude Science needs, e.g., the input file or gene list>
Steps in Claude Science:
  1. <e.g., connect Open Targets; for these genes, pull tractability, known drugs, genetic evidence>
  2. <e.g., have the reviewer agent verify the numbers and citations>
Return to Claude Code:
  - <exactly what to paste back, e.g., a TSV with columns gene, tractability, max_phase, ...>
  - the provenance artifact id or link
=== END HANDOFF #N ===
```

## 9. Tracker

Every Claude Code action and every Claude Science handoff is logged in
[`claude_tooling_log.md`](./claude_tooling_log.md), one row each, with the tool, the task, the
prompt or commit, the artifact reference, and a result summary. That file is the audit trail for
the Claude Use criterion and a demo asset.
