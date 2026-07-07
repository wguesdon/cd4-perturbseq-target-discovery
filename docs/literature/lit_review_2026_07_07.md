# Literature review: methods and prior work for a CD4+ T cell Perturb-seq target triage

Date: 2026-07-07. Track: Researcher, 2026 Built with Claude: Life Sciences.

How this was produced: six Claude research agents were fanned out in parallel, one per
facet, each returning a cited summary with an explicit white-space check. This document is
the synthesis. The parallel-agent survey is itself part of our Claude Science story.

Companion files: the reusable software and data stack is in
[`methods_toolkit.md`](./methods_toolkit.md); the source paper's own future work is in
[`source_paper_future_work.md`](./source_paper_future_work.md); the strategy is in
[`../research_plan.md`](../research_plan.md).

---

## 0. Executive summary and novelty verdict

The dataset (Zhu, Dann et al. 2025) is first of its kind: genome-scale, whole-transcriptome,
multi-context (Rest / Stim8h / Stim48h), in primary human CD4+ T cells. No prior screen
delivers all four at once. That combination is what makes our idea newly feasible.

Every building block of our idea is established and reusable. The novelty is the
composition, not the parts. Three defensible white-space contributions survive the survey:

1. A context-selective therapeutic-window composite metric. It fuses pathogenic-module
   suppression, stimulation-selectivity (Stim effect with near-zero Rest effect), and a
   viability/homeostasis penalty into one perturb-seq target score. The components are
   published; the composite is not a named metric.
2. In-silico on-target safety read from the knockdown transcriptome. No paper names or
   benchmarks it. We anchor it with DepMap and Hart core-essential lists and GTEx/HPA
   expression breadth.
3. The first formal drug-target-recovery evaluation of a Perturb-seq perturbation ranking,
   using ready-made gold-standard target sets and standard metrics (AUROC, top-k enrichment,
   precision/recall).

Two important corrections to the earlier framing:

- Recovering known targets to validate a screen is an ESTABLISHED principle (Dong/Ye 2019;
  Frangieh 2021; Fang 2019 for immune-trait genetics). Our novelty is the first formal
  precision/recall evaluation for a Perturb-seq ranking, not the idea of recovery. Soften any
  "this has never been done" language accordingly.
- Signature reversal and selective essentiality are mature. We cite them as precedent and
  claim the unification within one primary-T-cell screen.

Two strategic reads:

- Keep foundation/perturbation models off the one-week critical path. Independent benchmarks
  show they mostly do not beat simple baselines for perturb-seq. Geneformer zero-shot
  in-silico deletion is the one credible optional target-prioritization route.
- The degrader-handle idea is real, operationalizable, and topical (oral STAT6 degrader
  KT-621). It upgrades our druggability filter and anchors a possible type-2 / allergy frame.

---

## 1. Functional-genomics screens in primary human T cells

The experimental lineage our dataset extends.

- Simeonov et al. 2017, Nature. DOI 10.1038/nature23875. CRISPRa tiling of the IL2RA locus in
  CD4+ T cells; stimulation-responsive enhancer with an autoimmune-risk variant. Origin of the
  context-specific-regulator framing.
- Shifrut et al. 2018, Cell. DOI 10.1016/j.cell.2018.10.024. SLICE genome-wide KO in primary
  human T cells; separates TCR-signaling essentials from proliferation-tuning negative
  regulators. First feasibility proof.
- Henriksson et al. 2019, Cell. DOI 10.1016/j.cell.2018.11.044. Th-differentiation screens;
  activation vs differentiation crosstalk.
- Cortez et al. 2020, Nature. DOI 10.1038/s41586-020-2246-4. Treg FOXP3 screen; Usp22/Rnf20
  ubiquitin switch.
- Schmidt et al. 2022, Science. DOI 10.1126/science.abj4008. Genome-wide CRISPRa/i for IL-2 and
  IFN-gamma, with CRISPRa Perturb-seq on hits. Direct methodological precursor.
- Freimer et al. 2022, Nat Genet. DOI 10.1038/s41588-022-01106-y. CD4+ regulator networks for
  IL2RA/IL-2/CTLA4; ~90% of IL2RA regulators act in a cell-type- or stimulation-specific way,
  some with opposite effects across conditions. Strong empirical support for context windows.
- Legut et al. 2022, Nature. DOI 10.1038/s41586-022-04494-7. ORF gain-of-function screen; LTBR.
- Frangieh et al. 2021, Nat Genet. DOI 10.1038/s41588-021-00779-1. Perturb-CITE-seq; recovers
  known immune-evasion mechanisms and nominates CD58. Template for "re-derive known biology,
  therefore trust the novel hits."
- Dong/Ye et al. 2019, Cell. DOI 10.1016/j.cell.2019.07.044. In vivo CD8 screen re-identifies
  PD-1 and Tim-3 as a validation benchmark, then nominates Dhx37.
- Chen et al. 2025, Nature. DOI 10.1038/s41586-025-08795-5. Genome-wide FOXP3 screen; RBPJ-NCOR
  repressor as a context-specific negative regulator.
- Source: Zhu, Dann et al. 2025, bioRxiv. DOI 10.64898/2025.12.23.696273.

Gap: prior "drug-target" work mostly optimizes T cell function upward for cancer immunotherapy.
The inverse, context-selective suppression of a pathogenic CD4+ program while sparing rest and
viability, is shown for isolated genes but never as a systematic, benchmark-validated triage
across a genome-scale matrix.

---

## 2. Perturb-seq analysis methods and the reusable toolkit

De-facto pipeline: QC and guide assignment; Mixscape-style filtering of non-perturbed cells;
E-distance for perturbation magnitude; pseudobulk DE per donor x perturbation x condition;
program scoring; perturbation-similarity clustering; calibrated FDR against non-targeting
controls. Our released effect matrices already encode the pseudobulk-DE step.

- Peidli et al. 2024, Nat Methods (scPerturb, E-distance). DOI 10.1038/s41592-023-02144-y.
  Energy-distance is the standard scalar for perturbation magnitude; E-test must be
  cell-count-gated.
- Heumos et al. 2025, Nat Methods (pertpy). DOI 10.1038/s41592-025-02909-7. The scverse
  perturbation framework: E-distance, Mixscape, Augur, DE. Single most reusable package.
- Papalexi et al. 2021, Nat Genet (Mixscape). DOI 10.1038/s41588-021-00778-2. Per-cell
  perturbation signature and escaping-cell removal. Basis for guide-consistency gating.
- Skinnider/Squair et al. 2021, Nat Biotech (Augur). DOI 10.1038/s41587-020-0605-1. Ranks which
  states respond most, per condition.
- Badia-i-Mompel et al. 2022, Bioinformatics Advances (decoupler). DOI 10.1093/bioadv/vbac016.
  Scores a genes x contrasts statistic matrix against programs. Runs directly on our effect
  matrices. Highest-leverage tool for us.
- Aibar et al. 2017, Nat Methods (SCENIC/AUCell). DOI 10.1038/nmeth.4463; UCell (Andreatta 2021)
  as the faster per-cell scorer; Tirosh et al. 2016 score_genes / AddModuleScore as the light
  default.
- Squair et al. 2021, Nat Commun. DOI 10.1038/s41467-021-25960-2. Pseudobulk DE is required;
  naive single-cell tests inflate FDR.
- Barry et al. 2021, Genome Biol (SCEPTRE). DOI 10.1186/s13059-021-02545-2. Calibrated FDR for
  single-cell CRISPR screens.
- Zhou/Barry et al. 2023, Nat Methods (GSFA). DOI 10.1038/s41592-023-02017-4. Joint gene-module
  and perturbation-effect inference. Aligned with our module-shift framing.
- Dann et al. 2022, Nat Biotech (Milo). DOI 10.1038/s41587-021-01033-z. Differential abundance;
  compositional shift toward/away from a state. Note Emma Dann co-authored Milo and the source.

White-space: no named metric combines pathogenic-module-shift, stimulation-selectivity, and a
therapeutic-window framing into one perturb-seq target score. Operational definition to use:
`decoupler` module-score(Stim) minus module-score(Rest) per knockdown, normalized by
within-condition null variability (E-distance-style signal-to-noise), gated by donor/guide
concordance (Mixscape-style consistency).

---

## 3. Therapeutic window, selective essentiality, and signature reversal

- Signature reversal for target/drug discovery: Lamb et al. 2006, Science (Connectivity Map),
  DOI 10.1126/science.1132939; Subramanian et al. 2017, Cell (LINCS L1000, clue.io), DOI
  10.1016/j.cell.2017.10.049; Sirota et al. 2011, Sci Transl Med, DOI 10.1126/scitranslmed.3001318;
  Chen et al. 2017, Nat Commun (RGES; reversal magnitude tracks efficacy), DOI 10.1038/ncomms16022.
- Selective dependency / viability floor: Tsherniak et al. 2017, Cell (DepMap), DOI
  10.1016/j.cell.2017.06.010; Hart et al. 2017, G3 (CEG2 core-essential and nonessential lists),
  DOI 10.1534/g3.117.041277; Blomen 2015 and Wang 2015 essentialomes.
- On-target safety from expression breadth: GTEx (2013, DOI 10.1038/ng.2653) and Human Protein
  Atlas (Uhlen et al. 2015, Science, DOI 10.1126/science.1260419); broad, high tissue expression
  correlates with adverse events.
- Per-cell displacement toward/away from a state: CellOT (Bunne et al. 2023, Nat Methods, DOI
  10.1038/s41592-023-01969-x) and CINEMA-OT (Dong et al. 2023, Nat Methods, DOI
  10.1038/s41592-023-02040-5).

White-space: reversal and selective essentiality are established individually. Composing them
into a stimulation-contrast window inside one CD4+ Perturb-seq, and reading on-target safety off
the knockdown transcriptome, are the novel moves. Cite Chen 2017 and Sirota 2011 (reversal),
Tsherniak 2017 and Freimer 2022 (selective window), then claim the unification.

---

## 4. Druggability, tractability, and targeted degradation

- Druggable genome: Hopkins & Groom 2002, Nat Rev Drug Discov, DOI 10.1038/nrd892; Finan et al.
  2017, Sci Transl Med (the modern ~4,479-gene tiered list to intersect), DOI 10.1126/scitranslmed.aag1166.
- Open Targets Platform (Ochoa et al. 2021, NAR, DOI 10.1093/nar/gkaa1027). Tractability buckets
  for small molecule, antibody, and PROTAC modalities. One GraphQL call per target.
- Pharos / TCRD (Kelleher et al. 2023, NAR, DOI 10.1093/nar/gkac1033). Target Development Level:
  Tclin / Tchem / Tbio / Tdark. Tclin is a ready-made positive set for target-recovery; Tdark/Tbio
  is a novelty axis.
- ChEMBL (Zdrazil et al. 2024, NAR, DOI 10.1093/nar/gkad1004); DGIdb 5.0 (Cannon et al. 2024, NAR,
  DOI 10.1093/nar/gkad1040); structure-based ligandability via canSAR (Gingrich et al. 2025, NAR,
  DOI 10.1093/nar/gkae1050, includes AlphaFold pockets) and fpocket (2009, DOI 10.1186/1471-2105-10-168);
  DrugnomeAI ML prior (Raies et al. 2022, Commun Biol, DOI 10.1038/s42003-022-04245-4).
- Targeted degradation: Schneider et al. 2021, Nat Rev Drug Discov (the PROTACtable genome; the
  published "degrader-handle" score, now in Open Targets), DOI 10.1038/s41573-021-00245-x;
  PROTAC-DB 3.0 (Ge et al. 2025, NAR, DOI 10.1093/nar/gkae768); MolGlueDB (2025) for glue degraders,
  which reach pocketless TFs; E3 ligandability constraints (Biochemistry 2022, DOI
  10.1021/acs.biochem.1c00464): only ~2% of E3s are ligand-engaged, CRBN and VHL dominate and are
  ubiquitous.
- STAT6 degrader (the archetype): KT-621, Kymera, oral CRBN-based STAT6 degrader, Phase 1b in
  atopic dermatitis with ~94-98% degradation and ~62-63% EASI reduction; medicinal-chemistry
  precedent AK-1690 (J Med Chem 2024, DOI 10.1021/acs.jmedchem.4c01009).
- Genetic support for targets: Nelson et al. 2015, Nat Genet, DOI 10.1038/ng.3314; Minikel et al.
  2024, Nature (genetic support ~2.6x approval odds; reusable gold-standard target sets), DOI
  10.1038/s41586-024-07316-0; Fang et al. 2019, Nat Genet (immune-trait target landscape, benchmarks
  recovery of approved immune targets), DOI 10.1038/s41588-019-0456-1.

Degrader-handle score (operational): Open Targets PROTAC bucket (base) + canSAR/fpocket
ligandability + existing-degrader bonus (PROTAC-DB, MolGlueDB) + CD4-expressed E3 availability
(CRBN/VHL). Use it as a node weight that rescues undruggable intracellular regulators, the STAT6
case, not as a standalone ranker. Caveat: CRBN/VHL ubiquity means degradability does not imply
CD4-selectivity; selectivity must come from the target's own context-selective expression, which
the screen provides.

Benchmark white-space: metrics for target recovery are standard (AUROC, top-k enrichment OR,
precision/recall). Genetics-vs-target and fitness-screen-QC benchmarks exist, but a formal
recovery benchmark for a pooled/Perturb-seq perturbation ranking against an approved-drug-target
gold standard does not. That specific evaluation is our under-occupied contribution.

---

## 5. Disease biology, drug targets, and genetics (benchmark ground truth)

Known-immunomodulator target genes for the benchmark (target gene, drug, mechanism, disease;
intracellular vs surface flag matters for the degrader frame):

- Intracellular / regulatory (degrader-relevant): PPP3CA/PPP3CB/PPP3R1 with FKBP1A/PPIA and
  NFATC1/2 (tacrolimus, ciclosporin); MTOR with FKBP1A (sirolimus); JAK1/2/3, TYK2 (tofacitinib,
  baricitinib, upadacitinib); DHODH (teriflunomide); PDE4B/PDE4D (apremilast); STAT6, GATA3 (type-2).
- Surface receptors: IL2RA/CD25 (basiliximab); S1PR1/S1PR5 (fingolimod, ozanimod, etrasimod);
  ITGA4/ITGB7 (vedolizumab, natalizumab); IL4R (dupilumab); IL5RA; IL17RA; IL23R.
- Secreted cytokines: TNF (infliximab, adalimumab); IL4, IL13 (dupilumab, tralokinumab); IL5
  (mepolizumab); IL17A (secukinumab); IL23A/IL12B (ustekinumab, risankizumab).

Type-2 anchors: STAT6 is the obligate IL-4/IL-13 transducer; oral degrader KT-621 (Kymera) is the
topical proof point, further validated by human STAT6 gain-of-function genetics (JACI 2023; JEM
2023). GATA3 is the Th2 master TF; DNAzyme SB010 reached Phase 2a (Krug et al. 2015, NEJM).

Genetics stack and query recipe: GWAS Catalog and FinnGen (Kurki et al. 2023, Nature, DOI
10.1038/s41586-022-05473-8) for association; Open Targets Genetics L2G (Ghoussaini et al. 2021,
NAR, DOI 10.1093/nar/gkaa840) for causal-gene support; eQTL Catalogue and OneK1K (Yazar et al.
2022, Science, DOI 10.1126/science.abf3041) for CD4-specific colocalization. Concordance across
all three is strong genetic validation.

Risk-loci to expect: IBD/autoimmune (Th1/Th17): IL23R, IL12B, STAT3, CCR6, JAK2, TNFSF15, PTPN22,
PTPN2, IL2RA, IL6R. Allergy/AD (Th2): FLG, IL4, IL13, IL4R, TSLP, IL33, IL18R1, RORA, SMAD3,
GSDMB/ORMDL3 (17q21), STAT6.

---

## 6. Perturbation-prediction and single-cell foundation models

- Models: GEARS (Roohani et al. 2024, Nat Biotech, DOI 10.1038/s41587-023-01905-6); scGPT (Cui et
  al. 2024, Nat Methods, DOI 10.1038/s41592-024-02201-0); Geneformer (Theodoris et al. 2023, Nature;
  target-discovery protocol Nat Protoc 2026, DOI 10.1038/s41596-026-01364-8); scFoundation (Hao et
  al. 2024, Nat Methods, DOI 10.1038/s41592-024-02305-7); State (Arc, bioRxiv 2025.06.26.661135);
  CPA (Lotfollahi et al. 2023, Mol Syst Biol, DOI 10.15252/msb.202211517); Biolord; scPRAM.
- The critique: Ahlmann-Eltze, Huber & Anders 2025, Nat Methods, DOI 10.1038/s41592-025-02772-6 —
  for unseen perturbations, deep models did not beat predicting the training mean; for combinations,
  not better than an additive model. PerturBench (Wu et al. 2024, arXiv 2408.10609) — a latent-additive
  baseline ranks top. Arc's own Virtual Cell Challenge 2025 — models "not yet consistently
  outperforming naive baselines." State's own claim to beat linear baselines is a developer claim,
  not independently reproduced.
- Target-prioritization exception: Geneformer zero-shot in-silico deletion has an experimentally
  validated target-discovery precedent (cardiomyopathy, Theodoris 2023) and a 2026 protocol.

Recommendation: baselines on the critical path; models as a time-boxed, baseline-gated side arm;
Geneformer zero-shot as the most credible optional contender. A model earns a place only if it
beats the additive/linear/DE baseline on target ranking on this dataset.

---

## 7. Consolidated novelty positioning

Cite as precedent, do not claim: signature reversal (Lamb, Sirota, Chen); selective essentiality
(DepMap, Freimer); recovery-of-known-targets validation (Dong/Ye, Frangieh, Fang); PROTACtable
genome (Schneider); the source paper's context-specific-regulator mapping and GWAS/OpenTargets
enrichment (Zhu, Dann 2025).

Claim as contribution: (1) the context-selective therapeutic-window composite metric; (2)
in-silico on-target safety from the knockdown transcriptome; (3) the first formal
precision/recall drug-target-recovery evaluation of a CD4+ Perturb-seq ranking; and, as a
differentiator, (4) a degrader-handle axis fused into the target ranking.

The honest framing strengthens the pitch: we assemble established, citable primitives into a
validated, context-selective triage that no prior work delivers on this kind of dataset.
