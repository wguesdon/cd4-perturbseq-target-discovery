# Research plan: a validated, context-selective target-triage of the CD4⁺ T cell Perturb-seq map

> **Status:** Working draft, not finalized. A candidate strategy under evaluation.
> **Track:** Researcher, 2026 Built with Claude: Life Sciences.
> **Started:** 2026-07-07. All analysis performed during the event.
> **Source:** Literature review done with Claude Science. See [`PerturbTarget-T_feasibility.md`](./PerturbTarget-T_feasibility.md).
> **Prior-work check:** how this plan sits against the authors' own stated future work is in [`strategy_vs_authors_future_work.md`](./strategy_vs_authors_future_work.md).
> **Literature review:** the methods and prior-work survey is in [`literature/lit_review_2026_07_07.md`](./literature/lit_review_2026_07_07.md); the reusable toolkit is in [`literature/methods_toolkit.md`](./literature/methods_toolkit.md).

**One-line pitch.** A reproducible pipeline that first proves it can recover the targets of approved T cell immunomodulators from the Perturb-seq data alone, then uses that validated therapeutic-window score to nominate novel, druggable, context-selective regulators of pathogenic CD4⁺ activation, each delivered as a sourced evidence card.

**Core question.** Which genes, when knocked down, selectively suppress the stimulation-induced inflammatory program in human CD4⁺ T cells while sparing the resting state, homeostasis, and viability, and of those, which are druggable and genetically supported for autoimmune disease?

This plan replaces the original "PerturbTarget-T" framing. The literature search behind it is in [PerturbTarget-T_feasibility.md](./PerturbTarget-T_feasibility.md). The change in strategy is deliberate: the original plan led with druggability annotation, which is a thin API-join layer on top of work the source paper already does. This plan leads with two things the source paper does not do and that a published-literature search found to be white space.

## Rationale

### Why this wins on novelty and rigor

The source resource ([Zhu, Dann et al. 2025](https://doi.org/10.64898/2025.12.23.696273)) already maps context-specific regulators, links them to Th1/Th2 programs, and connects them to autoimmune disease risk via GWAS and Open Targets. It never uses the word "druggable," never benchmarks its nominations against known drug targets, and never converts context-specificity into a target-selection metric. Two design choices exploit those gaps:

1. **A therapeutic-window score is the novel object, not the gene list.** The paper shows context-specificity descriptively; this plan makes it the ranking axis. A target scores well only when knockdown suppresses the pathogenic program *specifically in stimulated cells* while leaving resting-state, homeostasis, and viability modules intact. The components are established and reusable: signature reversal (Lamb 2006; Sirota 2011; Chen 2017 RGES) and selective, context-dependent dependency (DepMap, Tsherniak 2017; Freimer 2022 shows ~90% of IL2RA regulators are context-specific in T cells). But no named metric composes module-shift, stimulation-selectivity, and a viability/homeostasis penalty into one perturb-seq target score. That composition is the contribution. The literature review gives the operational definition.
2. **A positive-control benchmark converts the whole project from plausible to validated.** Before ranking novel hits, show the pipeline recovers targets of approved immunomodulators from the data: calcineurin/NFAT (tacrolimus, ciclosporin), JAK1/3 (tofacitinib), mTOR (rapamycin, sirolimus), DHODH (teriflunomide), IL2RA (basiliximab), PDE4, S1PR1. Recovering known targets to validate a screen is an established principle (Dong/Ye 2019; Frangieh 2021; Fang 2019 for immune-trait genetics), but it has not been formalized as a precision/recall benchmark for a Perturb-seq perturbation ranking against an approved-drug-target gold standard. We use ready-made positive sets (Pharos Tclin, ChEMBL max_phase=4) and standard metrics (AUROC, top-k enrichment odds ratio, precision/recall). So the contribution is the first formal drug-target-recovery evaluation of a CD4+ Perturb-seq ranking, not the idea of recovery itself. If the score ranks known targets highly, every downstream novel nomination inherits that credibility, and you get a quantitative figure instead of an unfalsifiable ranking.

The prior art this builds on, and must cite as context rather than claim to displace: [Schmidt et al. 2022, *Science*](https://doi.org/10.1126/science.abj4008) (CRISPRa/i and Perturb-seq in primary human T cells), [Freimer et al. 2020, *Nature*](https://doi.org/10.1038/s41586-020-2246-4) and [2025, *Nature*](https://doi.org/10.1038/s41586-025-08795-5) (FOXP3 circuitry), [Henriksson et al. 2019, *Cell*](https://doi.org/10.1016/j.cell.2018.11.044) (T helper screens), and [atlas-guided TF discovery, 2026, *Nature*](https://doi.org/10.1038/s41586-025-09989-7).

### What the judges want to see, and how each is delivered

1. **A causal, not correlational, use of the data.** The score is built on CRISPRi perturbation effects, so every claim is "knocking this down does X," framed correctly for loss-of-function to drug-inhibition translation.
2. **A concrete, reusable artifact.** A ranked target table plus 5 to 10 evidence cards, not a folder of UMAPs.
3. **Validation the judges can check.** The known-target benchmark is the headline figure; precision/recall against the curated immunomodulator set is a number, not an opinion.
4. **Reproducibility.** Clean repo, environment file, one command to regenerate the table and figures from released summary statistics.
5. **Visible Claude-native leverage.** The literature synthesis, ground-truth-set curation, genetics annotation, and per-target evidence cards are agent-driven synthesis tasks, which is a stronger use of the tool than writing boilerplate code.
6. **A crisp story.** One question, one validated method, a few deep target case studies.

## Approach and timeline

### Method

**Inputs.** Released gene-level, per-condition perturbation effect matrices (Rest / Stim 8h / Stim 48h) with donor and guide resolution. Do not reprocess the 22M cells. Confirmed: the released object is `GWCD4i.DE_stats.h5ad` on the public no-auth S3 bucket `s3://genome-scale-tcell-perturb-seq/marson2025_data/`, with 33,983 perturbation-by-condition rows over 10,282 genes and layers `log_fc`, `zscore`, `p_value`, `adj_p_value`, `baseMean`, `lfcSE`. Guide- and donor-resolved variants (`GWCD4i.DE_stats.by_guide.h5mu`, `.by_donors.h5mu`) feed the confidence gating. See [`../data/README.md`](../data/README.md) for the download step.

**Program module scores.** Define 6 to 8 T cell programs from canonical marker sets: activation (CD69, IL2RA, CD40LG, ICOS), inflammatory cytokines (IFNG, TNF, IL2, IL17A, CSF2), Th1 (TBX21, CXCR3), Th17 (RORC, IL23R, CCR6), proliferation (MKI67, TOP2A, PCNA), stress/apoptosis (DDIT3, JUN, FOS, BAX), tolerance (FOXP3, CTLA4, IKZF2), homeostasis/memory (IL7R, CCR7, SELL, TCF7).

**Per-perturbation program shift.** For perturbation *g* and condition *c*, ΔS(g,c,k) is the module-*k* score under perturbation minus the non-targeting-control mean in the same condition.

**Therapeutic-window score.** Reward suppression of activation, inflammation, Th1/Th17, and proliferation; penalise loss of homeostasis/tolerance and induction of stress/apoptosis. Add a context-selectivity term that rewards a strong Stim effect with a near-zero Rest effect, and a viability/essentiality penalty cross-checked against an external essentiality reference (DepMap). Keep weights transparent and equal for the MVP; show a sensitivity analysis later.

**Confidence gating.** Report donor consistency and guide consistency for every hit; treat single-donor or single-guide effects as low confidence. Rank on calibrated effect sizes with FDR control, not nominal *p*-values.

**Annotation filters, applied after the window score, not before.** Druggability and tractability from Open Targets and ChEMBL (tractability buckets, known drugs, mechanisms, chemical probes); protein-class tractability from UniProt and InterPro; surface/secreted localisation from the Human Protein Atlas; autoimmune genetics from GWAS Catalog, eQTL Catalogue, and FinnGen. This genetics layer builds on and cites the authors' own GWAS and OpenTargets disease-gene enrichment; we do not claim it as novel. Add a degrader-handle axis: Open Targets PROTAC tractability buckets; canSAR or fpocket ligandability on AlphaFold models; existing-degrader flags from PROTAC-DB and MolGlueDB; and CRBN/VHL E3 availability in CD4+ cells. This upgrades the generic druggability filter and lets an intracellular regulator score as a target even without an inhibitory pocket, the STAT6/KT-621 archetype.

**Second differentiating axis (in-silico safety liability).** Because the readout is transcriptome-wide, characterise the mechanism-based liability of inhibiting each top candidate: does knockdown collapse homeostasis/memory programs, trigger stress modules, or hit housekeeping processes? Anchor the safety readout with public priors: DepMap common-essential and Hart CEG2 core-essential lists for collateral viability disruption, and GTEx or Human Protein Atlas expression breadth as an off-tissue safety prior. The white-space check found no paper that names or benchmarks in-silico on-target safety from knockdown transcriptomes, so this is a genuine differentiator; word the safety claims humbly.

### One-week execution plan

1. **Day 1: data and ground truth.** Load the released effect matrices; confirm per-condition, per-donor, per-guide resolution; set up controls. In parallel, curate the approved-immunomodulator target set (drug, target gene, mechanism) as the benchmark. Deliverable: clean notebook, QC plots, program gene sets, ground-truth table.
2. **Day 2: program scoring and effects.** Compute module scores and ΔS(g,c,k) across Rest / Stim 8h / Stim 48h. Deliverable: program-shift matrix by perturbation and condition; a program-effect heatmap for top perturbations.
3. **Day 3: therapeutic-window score and benchmark.** Implement the window score with context-selectivity and viability penalty; run the benchmark and report precision/recall of recovering known targets. Deliverable: benchmark figure (this is the headline), first ranked target table.
4. **Day 4: mechanism and context.** Pathway enrichment on downstream genes; context-specificity plot (Rest vs Stim 8h vs Stim 48h) for top targets; perturbation-signature clustering. Deliverable: mechanism candidates, context plot.
5. **Day 5: annotation and safety axis.** Layer druggability, genetics, and the in-silico safety-liability readout onto the top 20. Deliverable: annotated shortlist, target-prioritisation scatterplot.
6. **Day 6: evidence cards and reproducibility.** Agent-generate 5 to 10 sourced evidence cards; finalise README, environment file, one-command rerun. Deliverable: evidence cards, reproducible repo.
7. **Day 7: narrative and figures.** Assemble the five figures, write limitations and a wet-lab validation plan for the top 1 to 2 candidates. Deliverable: submission report and slides.

## Outputs, risks, and success criteria

### Deliverables

1. Ranked target table (target, best context, window score, benchmark rank if known, druggability, genetics, safety-liability flag, confidence).
2. Benchmark figure: precision/recall of recovering approved immunomodulator targets.
3. Five to 10 evidence cards (window summary, context specificity, key up/down genes, pathways, druggability, genetics, safety liability, caveats, proposed validation experiment).
4. Reproducible repository with environment file and one-command regeneration.
5. Submission report and slides.

### Figures

1. Workflow schematic: Perturb-seq effects → program scoring → therapeutic-window score → benchmark validation → druggability/genetics/safety annotation → ranked targets.
2. Benchmark: recovery of known immunomodulator targets (headline).
3. Program-effect heatmap: top perturbations vs the 6 to 8 programs.
4. Context-selectivity plot: Rest vs Stim 8h vs Stim 48h for top targets.
5. Target-prioritisation scatterplot: x = inflammation suppression, y = homeostasis/stress penalty (the therapeutic window), colour = druggability, size = genetics support.

### Risks and controls

1. **Directionality.** Suppression of an inflammatory program predicts benefit from inhibition only if the gene is not required for viability in the same context; the viability penalty and essentiality cross-check guard against nominating cytotoxicity as anti-inflammatory.
2. **Overclaiming novelty.** Top hits will include known regulators; claim novelty at the level of the validated method and the shortlist, not new biology.
3. **Competition on a fresh dataset.** The preprint posted 24 December 2025 and is a CZI community resource; reanalyses will appear within months, so favour speed and a citable method over breadth.
4. **Scope creep.** Hold to one disease framing (pathogenic CD4⁺ activation in autoimmunity); do not add oncology, aging, or Treg biology.
5. **Benchmark honesty.** Report where the score fails to recover a known target as well as where it succeeds; a benchmark that only shows wins is not a benchmark.
6. **Non-polarized data.** The screen used only non-polarized culture (Rest, Stim 8h, Stim 48h). The Th1 and Th17 program scores are therefore proxies read from marker genes within stimulated cells, not measurements under polarizing conditions. The authors flag polarizing cytokine conditions as unmapped future work. State this limit and avoid strong Th-subset claims.

### What to cut from the original plan

Drop druggability as the headline (demote to a filter); drop any plan to reprocess 22M cells; drop the standalone cosine-similarity recommender as a novel model, since it overlaps with the atlas-guided-TF and RNA-fingerprinting approaches the source paper cites. Anything requiring combinatorial perturbations is infeasible: the library is single-gene.

### Literature grounding (2026-07-07)

A six-agent literature survey (see [`literature/lit_review_2026_07_07.md`](./literature/lit_review_2026_07_07.md)) confirmed the positioning and refined it:

- Reusable toolkit, no reinvention: pertpy (E-distance, Mixscape, Augur), decoupler (module scoring straight off the effect matrices), pyDESeq2 pseudobulk DE, DepMap plus Hart CEG2 for the viability floor, CMap/LINCS RGES for signature reversal, GTEx/HPA for safety breadth, Open Targets/Pharos/ChEMBL for druggability, PROTAC-DB/MolGlueDB for degraders.
- No training on the critical path. Independent benchmarks (Ahlmann-Eltze/Huber 2025; PerturBench; Arc's Virtual Cell Challenge 2025) show foundation and perturbation models mostly do not beat simple additive/linear/mean baselines for perturb-seq. If a model must appear, use Geneformer zero-shot in-silico deletion as a bounded, baseline-gated side arm.
- Type-2 option: the same machinery retargets cleanly to a type-2 / allergy frame anchored on STAT6/GATA3, with the oral STAT6 degrader KT-621 (Kymera) as a topical hook. Held pending the team-up decision; see the local collaboration notes.

### Success criteria

The submission succeeds if a judge can see, in order: a method that recovers known immunomodulator targets with reported precision/recall; a ranked shortlist of context-selective, druggable, genetically supported novel candidates; evidence cards with sourced annotations and honest caveats; and a repository that regenerates the table and figures in one command.
