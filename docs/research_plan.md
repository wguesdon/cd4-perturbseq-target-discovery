# Research plan: a validated, context-selective target-triage of the CD4⁺ T cell Perturb-seq map

> **Status:** Working draft, not finalized. A candidate strategy under evaluation.
> **Track:** Researcher, 2026 Built with Claude: Life Sciences.
> **Started:** 2026-07-07. All analysis performed during the event.
> **Source:** Literature review done with Claude Science. See [`PerturbTarget-T_feasibility.md`](./PerturbTarget-T_feasibility.md).

**One-line pitch.** A reproducible pipeline that first proves it can recover the targets of approved T cell immunomodulators from the Perturb-seq data alone, then uses that validated therapeutic-window score to nominate novel, druggable, context-selective regulators of pathogenic CD4⁺ activation, each delivered as a sourced evidence card.

**Core question.** Which genes, when knocked down, selectively suppress the stimulation-induced inflammatory program in human CD4⁺ T cells while sparing the resting state, homeostasis, and viability, and of those, which are druggable and genetically supported for autoimmune disease?

This plan replaces the original "PerturbTarget-T" framing. The literature search behind it is in [PerturbTarget-T_feasibility.md](./PerturbTarget-T_feasibility.md). The change in strategy is deliberate: the original plan led with druggability annotation, which is a thin API-join layer on top of work the source paper already does. This plan leads with two things the source paper does not do and that a published-literature search found to be white space.

## Rationale

### Why this wins on novelty and rigor

The source resource ([Zhu, Dann et al. 2025](https://doi.org/10.64898/2025.12.23.696273)) already maps context-specific regulators, links them to Th1/Th2 programs, and connects them to autoimmune disease risk via GWAS and Open Targets. It never uses the word "druggable," never benchmarks its nominations against known drug targets, and never converts context-specificity into a target-selection metric. Two design choices exploit those gaps:

1. **A therapeutic-window score is the novel object, not the gene list.** The paper shows context-specificity descriptively; this plan makes it the ranking axis. A target scores well only when knockdown suppresses the pathogenic program *specifically in stimulated cells* while leaving resting-state, homeostasis, and viability modules intact. A published search for "context-selective essential gene therapeutic window T cell activation" returns zero papers, so the scoring axis itself is the contribution.
2. **A positive-control benchmark converts the whole project from plausible to validated.** Before ranking novel hits, show the pipeline recovers targets of approved immunomodulators from the data: calcineurin/NFAT (tacrolimus, ciclosporin), JAK1/3 (tofacitinib), mTOR (rapamycin, sirolimus), DHODH (teriflunomide), IL2RA (basiliximab), PDE4, S1PR1. A search for a CRISPR-screen benchmark that recovers known immunomodulator targets returns zero. If the score ranks these highly, every downstream novel nomination inherits that credibility, and you get a quantitative figure (precision/recall on a curated ground-truth set) instead of an unfalsifiable ranking.

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

**Inputs.** Released gene-level, per-condition perturbation effect matrices (Rest / Stim 8h / Stim 48h) with donor and guide resolution. Do not reprocess the 22M cells; confirm on Day 1 that these matrices are in the released objects and in the `emdann/GWT_perturbseq_analysis_2025` repo.

**Program module scores.** Define 6 to 8 T cell programs from canonical marker sets: activation (CD69, IL2RA, CD40LG, ICOS), inflammatory cytokines (IFNG, TNF, IL2, IL17A, CSF2), Th1 (TBX21, CXCR3), Th17 (RORC, IL23R, CCR6), proliferation (MKI67, TOP2A, PCNA), stress/apoptosis (DDIT3, JUN, FOS, BAX), tolerance (FOXP3, CTLA4, IKZF2), homeostasis/memory (IL7R, CCR7, SELL, TCF7).

**Per-perturbation program shift.** For perturbation *g* and condition *c*, ΔS(g,c,k) is the module-*k* score under perturbation minus the non-targeting-control mean in the same condition.

**Therapeutic-window score.** Reward suppression of activation, inflammation, Th1/Th17, and proliferation; penalise loss of homeostasis/tolerance and induction of stress/apoptosis. Add a context-selectivity term that rewards a strong Stim effect with a near-zero Rest effect, and a viability/essentiality penalty cross-checked against an external essentiality reference (DepMap). Keep weights transparent and equal for the MVP; show a sensitivity analysis later.

**Confidence gating.** Report donor consistency and guide consistency for every hit; treat single-donor or single-guide effects as low confidence. Rank on calibrated effect sizes with FDR control, not nominal *p*-values.

**Annotation filters, applied after the window score, not before.** Druggability and tractability from Open Targets and ChEMBL (tractability buckets, known drugs, mechanisms, chemical probes); protein-class tractability from UniProt and InterPro; surface/secreted localisation from the Human Protein Atlas; autoimmune genetics from GWAS Catalog, eQTL Catalogue, and FinnGen.

**Second differentiating axis (in-silico safety liability).** Because the readout is transcriptome-wide, characterise the mechanism-based liability of inhibiting each top candidate: does knockdown collapse homeostasis/memory programs, trigger stress modules, or hit housekeeping processes? A search for in-silico on-target safety from knockdown transcriptomes returns zero papers, so this is a genuine differentiator; word the safety claims humbly.

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

### What to cut from the original plan

Drop druggability as the headline (demote to a filter); drop any plan to reprocess 22M cells; drop the standalone cosine-similarity recommender as a novel model, since it overlaps with the atlas-guided-TF and RNA-fingerprinting approaches the source paper cites. Anything requiring combinatorial perturbations is infeasible: the library is single-gene.

### Success criteria

The submission succeeds if a judge can see, in order: a method that recovers known immunomodulator targets with reported precision/recall; a ranked shortlist of context-selective, druggable, genetically supported novel candidates; evidence cards with sourced annotations and honest caveats; and a repository that regenerates the table and figures in one command.
