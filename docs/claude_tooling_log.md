# Claude tooling log

An audit trail of how this analysis was built. One row per milestone, in time order.

The column that matters is the last one. A tool that only ever agrees with you has told you nothing.
Every row that records a **disagreement** is a place where a check, a critic, or a control contradicted
the analyst and was right.

## Claude Code

| # | date | task | commit or method | output | result | disagreement |
| --- | --- | --- | --- | --- | --- | --- |
| CC1 | 07-07 | Scaffold the private repo (uv, structure, MIT) | `b783719` | repo skeleton | project created | — |
| CC2 | 07-07 | Plan-mode research planning for preprint acquisition | plan file | approved plan | acquisition plan set | — |
| CC3 | 07-07 | Six-agent literature fan-out | 6 subagents | `literature/lit_review_2026_07_07.md` | cited review; novelty white-space | — |
| CC4 | 07-07 | Extract the source paper's future work | subagent | `literature/source_paper_future_work.md` | limitations captured | — |
| CC5 | 07-08 | Deep-read the paper figure by figure | 2 subagents | strategy inputs | `pert2state`; no druggability language in the source | — |
| CC6 | 07-07 | Competition and prior-work check | subagent | strategy inputs | space open | — |
| CC7 | 07-07 | Synthesise the strategy | `4852699` | `strategy_2026_07_07.md` | benchmark-validated triage | — |
| CC8 | 07-08 | Solo adjustment | `96dc69e` | strategy update | reversal-is-not-enough framing | — |
| CC9 | 07-08 | Research the Claude tooling ecosystem | 2 subagents | `claude_tooling_strategy_2026_07_08.md` | connector map | — |
| CC10 | 07-08 | Audit the released data before trusting the plan | S3 recon | — | the Dec-2025 supplementary CSV is an obs snapshot missing the confidence columns; the May-2026 h5ad is required | **Yes.** The plan assumed the supplementary table sufficed. It does not. |
| CC11 | 07-08 | Locate bundled priors without a 1.66 GB clone | `scripts/fetch_priors.sh` | 22 prior tables | found the IUIS list and the Schmidt CD4 IL-2 screen | — |
| CC12 | 07-08 | Define the window objective | `d5e6530` | `src/.../programs.py` | effector module separated from the co-inhibitory module | — |
| CC13 | 07-08 | Curate the benchmark ground truth | `4c94e42` | `ground_truth/immunomodulator_targets.csv` | direction-of-effect filter excluded agonists and drug-binding chaperones | — |
| CC14 | 07-08 | Steelman the risk-kill test | `scripts/02` | matched-background test | enrichment retested so power cannot fake the result | — |
| CC15 | 07-08 | **127-agent adversarial audit of the analysis, against itself** | subagent fan-out | `docs/results/adversarial_audit_2026_07_08.md` | two mechanisms retracted in one day | **Yes.** Essentiality gating is a collider: only 230 of 682 library core-essentials reach the DE table at all. And the IUIS flag is *more* enriched among approved drug targets (OR 8.31) than among top-ranked perturbations (OR 4.16). Both had been proposed as safety gates. |
| CC16 | 07-08 | Test the audit's own findings rather than accept them | `scripts/08`, `scripts/11` | funnel and IEI verification | the audit was right about one collider, missed a second, and was **wrong** about the co-inhibitory axis | **Yes, in both directions.** A critic that is never checked is another assertion. |
| CC17 | 07-08 | Reproduce the seed lottery | `scripts/02::seed_lottery` | 2,000 redraws | 22.4% of seeds non-significant, against the audit's 11.5% | **Yes.** The defect was twice as bad as reported. Replaced by CMH and van Elteren over every control row. |
| CC18 | 07-08 | Read the source preprint against our own code | `20aafd1` | `source_paper_briefing_2026_07_08.md` | 147 claims checked | **Yes.** Six of our claims refuted. The platform is a probe panel, not the transcriptome; "Rest" is expanded blasts, not quiescent naive cells. |
| CC19 | 07-08 | Pre-register N9 before writing the script | `4913bfb` | `preregistration_n9_2026_07_08.md` | both branches acceptable | — |
| CC20 | 07-08 | N9 induction-matched null | `db8e6e2` | `scripts/17` | the module survives | **Yes.** The registered negative control fired on the first run: false-positive rate 25% against a nominal 5%. The effective unit is one module, not 6,371 perturbations. A p-value of 3.5e-13 was retracted. |
| CC21 | 07-08 | N6 selectivity dominance test | `00d3a47` | `scripts/18` | the axis lost to a one-variable shadow of itself | **Yes.** A gate axis was deleted on the strength of its own control. |
| CC22 | 07-08 | N8 Open Targets benchmark ceiling | `scripts/16` | `open_targets_benchmark_ceiling.csv` | 38 primary, 53 upper bound, against a threshold of 60 | **Yes.** The first run said 173 targets, "POWERED", which is the opposite of the truth. A `validate()` harness with five positive and four negative controls withheld the verdict until three inflation traps were found. |
| CC23 | 07-09 | 33-agent literature review, four streams | `374e077` | `literature/target_discovery_ideas_2026_07_08.md` | 28 ideas, prioritised | Recorded at the time: the safety classifier was down for several agents, so **its citations were not trusted** and were re-verified later (CC30). |
| CC24 | 07-09 | Evidence-card fabrication firewall | `7307b46` | `scripts/24` | cards written from a frozen record | **Yes, by design.** A deliberately fabricated card, with an invented constraint value and an invented citation, is caught on all three checks. Claude writes prose; it never sets a rank or a verdict. |
| CC25 | 07-09 | N16 synthesis, then a 7-agent adversarial critique of it | `853c604`, `f405d7f` | `final_shortlist_2026_07_09.md` | zero vetted novel nominations | **Yes.** The peer critic killed five overclaims, among them "exactly one candidate", which it identified as a filter-intersection artifact. |
| CC26 | 07-09 | N17, direction of effect via eQTL colocalisation | `76f8361` | `scripts/26` | the method self-invalidated | **Yes.** Calibration on a nine-gene known-direction panel returned 0.50, chance. The method was refused rather than shipped. A first two-gene calibration had looked as though it worked. |
| CC27 | 07-09 | 41-agent direction-carrying literature and dataset sweep | workflow | decision memo | recommended mouse knockout phenotype | **Yes.** A feasibility probe then measured its coverage and found the instrument inverts `PTPN2`, the one gene it was recruited to get right. Not run. |
| CC28 | 07-09 | N20 nomination-rule audit, pre-registered | `4db93ad`, `15a9378` | `scripts/29` | the rule re-nominates 2 of the 5 drugs it recovers | **Yes.** A self-test refused to run until it reproduced the published shortlist exactly. It caught a swapped column that had silently returned 16 genes instead of 1. |
| CC29 | 07-09 | 30-agent adversarial critique of N20 | workflow | peer-critic verdict | one real defect, five overclaims | **Yes.** The critic found that the nomination table gated on `tractable_with_precedent` while its own controls gated on structural `tractable`: the exact circularity the analysis had diagnosed one commit earlier, reintroduced by the analyst. |
| CC30 | 07-09 | Verify all 25 references by fetching them, before writing the paper | workflow | `citation_verification_2026_07_09.md`, `references.bib` | 24 of 25 pass | **Yes.** `manguso2017` is real but concerns mouse tumour-cell-intrinsic Ptpn2. It cannot support the human PTPN2 direction claim on which the analysis turns. Replaced by `wiede2011`, `wtccc2007` and `jeanpierre2024`. |
| CC31 | 07-09 | Verify the vendored comparator screens | agent and `scripts/27` | `PROVENANCE.md` | both descriptions corrected | **Yes.** `Arce 2025` is a FACS marker screen, not the fitness screen our own provenance file claimed. And it shares its 1,351-gene library with `Freimer 2022` **exactly** (Jaccard 1.0000), so "two independent screens" was true of the assay and false of the gene space. |
| CC32 | 07-09 | External peer review, acted on | `9a19ba5`, `9ac6fef` | corrected pipeline and manuscript | numbers moved; a control now fails | **Yes, comprehensively.** Two efficacy statistics were in use. The evidence floor admitted inducers. `P < 0.0001` was unsupportable from 5,000 draws. Two module p-values were floors, not estimates. The direction rules treated a cross-linking agonist as an inhibitor. Correcting them made a pre-registered control fail, and the nomination rule is now void by its own criterion. |
| CC33 | 07-09 | Trace the two reported AUROC values to their source | `scripts/35` | `auroc_validation.csv` | both recomputed on the primary statistic | **Yes.** Both had been computed on a third score. On the primary statistic the drug-recovery AUROC is 0.497 [0.339, 0.656], not 0.542. The corrected value is less favourable. |

## Claude Science handoffs

Not run. The schedule reserved the weekend for Claude Science as an adversarial reviewer against results
already committed to git, so that any disagreement would be a contemporaneous A/B rather than a first
pass. The analysis consumed that time instead, and the reviewer role was filled by the agent critiques
above (CC25, CC29) and by external peer review (CC32).

The two prepared handoffs remain in `docs/handoffs/`. `HANDOFF_02` is superseded: it asks for
verification of claims this analysis has since retracted.

## What the tooling actually bought

Nineteen of the thirty-three rows record a disagreement, and in every one of those the tool was right
and the analyst was wrong. The pattern is consistent.

**Controls written to fail caught more than critics did.** The N9 negative control, the benchmark's
`validate()` harness, and the N20 self-test each fired on their first run, before any human read the
output.

**Adversarial critique caught what controls could not, but only when the critic was given the code
rather than a summary.** The N20 peer critic found a defect the controls had passed, because the
controls and the table were computed on different columns. A critic reading a summary would have seen
neither.

**Every instrument that promised direction of effect failed for a different, principled reason.** An
expression QTL reports abundance, not activity. A systematic mouse knockout programme under-detects
autoimmunity. A full null is not a partial inhibitor. Those three failures are the paper's central
finding, and none of them is a bug.

**A citation is not evidence until it is fetched.** One reference in twenty-five was real, correctly
formatted, and unable to support the claim it was cited for.
