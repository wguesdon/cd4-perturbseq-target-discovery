# Decision Memo: Built with Claude Life Sciences, Research Track

Prepared 2026-07-09 evening. Science must be done by Sunday night. Monday is video, summary, repo.
Produced by a 91-agent workflow (13 proposers, dedupe, per-candidate data verification, 5 red-team
lenses each, 3 judge panels, completeness critic). 39 proposals, 12 shortlisted, 12 survivors.

**DECISION: SatMut selected by the user 2026-07-09.**

---

## 1. The recommendation

**Build SatMut: predict, then break, a disease enhancer base by base.**

It is the one candidate whose core positive result is published fact, so it cannot end in an
abstention the way the CD4 project did, and you will know it is positive by Friday noon. Its demo is
the most legible in the field: a wet lab saturation mutagenesis heatmap and a model heatmap lighting
up in the same places, then a live database lookup that turns agreement into an on screen surprise.
It has the lowest execution risk of any option, one environment, trivial compute, tiny open data,
which for a solo builder with three days is the factor that dominates everything else.

Accept the cost with eyes open. SatMut has the lowest scientific novelty of the shortlist. The hero
heatmap is essentially a 2021 Nature Methods figure. So the framing must be surgical. Own the
reproduction as a trust signal, and make the actual contribution three things: a sub fifty dollar
fully open recipe any lab reruns in an afternoon, a live Open Targets and ClinVar cross check that
surfaces a variant the model was never told about, and one concrete new variant nomination. Ship it
as an interactive tool, not a Quarto report. That is the pattern that wins hackathons judged by
Anthropic.

---

## 2. Why it beats the runners up

Consensus across three independent judge panels. Expected score is the rubric weighted total times
calibrated P(positive by Sunday), averaged over panels.

| Project | P(positive) | Demo wattage | Expected score |
|---|---|---|---|
| SatMut (regulatory DNA) | 0.80 | High. Twin heatmaps, live lookup, base by base sweep | 48 |
| The Metric Decides (CD4) | 0.67 | Medium high. Causal search tool, jargon metric bars | 44 |
| Direct or Indirect (Krogan) | 0.57 | High. Recoloring hairball, rotating dimers | 38 |
| Beat the Book (TB) | 0.60 | High. Catalogue miss split screen, rotating Rv0678 | 36 |
| NeedleRank (neoantigen) | 0.69 | Medium. Needle enrichment strip | 36 |

SatMut is ranked first by two of three panels and second by the third. Highest P(positive), highest
expected score. Every other candidate either carries a fatal circularity, is pre empted by a named
recent paper on the exact data, or has a positive control that fires for a trivial reason. SatMut has
one real weakness, novelty, and that weakness is addressable by framing rather than by luck.

---

## 3. SatMut in detail

**Biological question.** Which single base substitutions in disease regulatory elements change their
activity, and can a sequence model redraw the measured saturation mutagenesis map of held out disease
loci well enough to nominate a causal non coding variant. This is Gladstone example (b) verbatim.

**Datasets.**

| Dataset | Role | Size | Access |
|---|---|---|---|
| Kircher 2019 satMut MPRA | Ground truth, ~17,500 SNVs, 14 elements, 9 promoters + 5 enhancers | few MB TSV | Open, no login. kircherlab.bihealth.org/satMutMPRA. Mirrors: GEO GSE126550, MaveDB urn:mavedb:00000014-a |
| Enformer weights | Frozen oracle, covers the whole critical path | ~70 MB class | Open, CC BY. HuggingFace EleutherAI/enformer-official-rough |
| gReLU | ISM plus fine tune tooling | code | Open. github.com/Genentech/gReLU |
| AlphaGenome API | Optional independent oracle | API | Free non commercial key, registration then instant. Optional only |

The Kircher release is the CAGI5 test set in open form. It exhaustively measures every SNV, so the
ground truth is unbiased. The null signal skeptic conceded that the collider that killed the CD4
project cannot operate here.

**Deliverable.** An in silico saturation mutagenesis pipeline scoring every SNV with Enformer,
predicted versus measured heatmaps with per locus Spearman, a cheap frozen or ridge readout head, one
nominated variant checked against ClinVar, GWAS Catalog, and Open Targets. Public repo plus a live
interactive front end.

**Named baseline and published number.** The floor is deltaSVM. Shigaki 2019 reports it at concordance
0.30 on single DHS, roughly half the correlation of the modern models. Two corrections before the
report freezes. Shigaki et al. were the CAGI5 organizers running post hoc integration, not the
winning team. AlphaGenome scores the CAGI5 satMut set at Pearson 0.57 cell type matched and 0.63 with
aggregation, not 0.65. Do not claim Enformer's published win as your own, since Enformer is your tool.

**Built in positive control.** Promoter saturation mutagenesis recovery. The predicted versus measured
Spearman on promoters is published fact and fires Friday regardless of the enhancer payload. The video
works even if every novel element comes up null.

**Compute plan.** One g5.2xlarge. ISM over 17,500 variants is minutes to low hours. The frozen head is
under one GPU hour. Total under forty dollars. This could run on the local NVMe workstation. No
SageMaker fan out required. Pull Enformer weights first thing and cache them so the critical path
never depends on the AlphaGenome key.

**Three minute shot list.**
1. Split screen. Left, the measured LDLR or TERT promoter heatmap, position on x, A C G T on y, red
   activating, blue silencing. Right, the model heatmap from sequence alone. They light up together.
2. Ten second caption owning it. "The field's 2021 benchmark, reproduced by any lab for under fifty
   dollars with a frozen head." Then move fast.
3. The climax. The model flags a base on a held out variant. Claude Science queries ClinVar, GWAS
   Catalog, and Open Targets live, and the hit returns as a known pathogenic or GWAS variant the model
   was never shown. Agreement becomes surprise.
4. One large legible number or a single climbing bar. Cut the deltaSVM scatter and the word Spearman.
5. Animate the ISM sweep base by base so the heatmap builds live, to compete with rotating protein
   videos.
6. Close on the one novel nomination with its Open Targets card.

**Claude Science role.** Pull Kircher records into the workbench. For each hotspot run live ClinVar,
GWAS Catalog, and Open Targets MCP lookups to auto annotate which bases already carry a drug or a
monogenic syndrome. Query AlphaGenome as an independent oracle to sanity check the map. Narrate the
reasoning, including discarding a weak association on camera.

**Claude Code role.** Orchestrate the run, write the gReLU fine tune, ISM, and motif attribution,
enforce the Friday gate by refusing to proceed if held out r is below 0.6, build the Quarto report and
the animated heatmap renderer.

**Killer risk and mitigation.** The risk is that a field literate judge reads it as rerunning
Enformer. Mitigate by leading with the visceral twin heatmap match, headlining the cost, making the
live lookup the emotional climax, and delivering one concrete new variant. Demote the enhancer beat
Shigaki sub goal. It rests on only five loci and is roughly a coin flip, so report it with a bootstrap
interval and treat it as a bonus, never the headline.

---

## 4. Day by day plan for SatMut

**Friday.** Set up the uv environment and gReLU first thing. Pull and cache Enformer weights and the
Kircher effect tables, committing the raw TSVs into the repo so the ground truth is frozen. Run zero
shot Enformer ISM on one hero promoter. Submit the AlphaGenome key request in parallel as a hedge, but
do not depend on it.

**Friday noon gate.** On a single hero promoter such as LDLR or TERT, zero shot Enformer ISM must
reach Spearman at or above 0.5 versus the Kircher measured effects. This is published fact, so a
failure means the environment, the weight load, or the coordinate mapping is broken, not the science.
First remedy: switch the oracle to the AlphaGenome API, which needs no local GPU and reaches the same
result. If Spearman is still below 0.5 by Friday 6pm, abandon to the fallback.

**Fallback project: NeedleRank.** CPU only, tiny data, positive control (immunogenic needle enrichment
on the TESLA set) fires by Saturday. Shares the recover the known needles demo shape, so a same day
pivot loses the least. Beat the Book is the alternative if you want more real world wattage, but it
needs a larger download.

**Saturday.** Run ISM across all 14 elements, produce every heatmap and per locus Spearman, lock the
positive. Fit the frozen or ridge readout head. Report the enhancer comparison honestly with a
bootstrap interval and demote it from the headline. The guaranteed deliverable is fully frozen by
Saturday night.

**Sunday.** Add the live MCP validation moment. Run ClinVar, GWAS Catalog, and Open Targets lookups on
the recovered hotspots. Nominate one variant. Run the scrambled motif control to show the map
collapses. Freeze the Quarto report with every inline number bound to a committed table. Build the
animated renderer.

**Monday.** Record the three minute video, write the 100 to 200 word summary, make the repo public.
Reserve the whole day for this.

---

## 5. What the red team could not break, and what it could

**Could not break.**
- Data access. Kircher set open, tiny, mirrored three ways, no registration on the critical path. The
  data skeptic filed no objection.
- Compute. Trivial. Minutes to low hours on one GPU, under forty dollars, feasible locally.
- The positive control. Promoter recovery is a published result, so unlike the CD4 project the demo
  cannot null out. The null signal skeptic explicitly confirmed the collider that killed CD4 cannot
  operate on this exhaustive, unbiased ground truth.
- Schedule. One environment, one dataset, one hero locus by Friday noon. The schedule skeptic could
  not land a fatal or serious verdict.

**Could break, and you must address it.**
- Novelty. The centerpiece heatmap is already an Avsec 2021 figure. Real, caps the ceiling. Fix is
  framing. Own the reproduction as a trust signal and make the contribution the open recipe, the live
  discovery, and the new variant.
- The demo climax. As first drafted it celebrates agreement, not discovery. Re cut so the live
  database lookup surfaces something the model was never told, the only genuine surprise on screen.
- The enhancer beat. Five loci, roughly a coin flip, statistically fragile. Demote it, report an
  interval, never headline it.
- Claude Science is thinner here than in a connector heavy project. Force it on screen and make the
  Open Targets lookup central, since Claude Use is a quarter of the score.
- Two factual errors to fix before freeze. AlphaGenome scores 0.57 and 0.63, not 0.65. Shigaki et al.
  were the organizers, not the winning team.

---

## 6. What was rejected and why

**Fatal on inspection.**
- Complete the Complex (Krogan/CORUM). hu.MAP edge weights are trained on CORUM and Complex Portal
  labels, so leave one subunit out recovery is a supervised leakage tautology an expert flags live.
- Blind 2022 VUS resolver. RENOVO already ran the identical prospective time split, and the ground
  truth is doubly circular. Population frequency drives the benign flips and AlphaMissense as evidence
  drives the pathogenic flips.
- PRISM precision oncology. A per drug model structurally cannot run the held out drug benchmark it
  claims to beat, and the BRAF to vemurafenib control is a lineage collider a Gladstone judge debunks
  on sight.
- Compensatory HARs. Finding already published by the cited paper, signed ISM provably cannot model
  the pairwise interaction the demo hinges on, ground truth exists for only three HARs, flagship demo
  example is from a different paper.

**Survivable but not finalists.**
- Phenocopy Cell Painting. Both proposed beat levers already baked into the published 0.277 number,
  CellProfiler ties the foundation model, cell fields not legible to the eye, Recursion license blocks
  a clean open source repo.
- Reverse the Disease (LINCS repurposing). Recovery is a tautology (an immunosuppressant reverses an
  immune signature), and Yang 2021 already did the ulcerative colitis connectivity analysis with
  clinical validation.
- NeedleRank, Beat the Book, InterfaceStack. Clean reproductions with reliable positive controls, kept
  as insurance. NeedleRank is the named Friday fallback because it is CPU only and its control fires
  fastest.

**Completeness critic additions.**
- Bloom lab immune escape time machine. High demo wattage, real gap in modality coverage. Rejected:
  never run through the gauntlet, the prepandemic forecasting framing is already EVEscape (Marks lab,
  Nature 2023), and the Omicron prediction payload is null risky in exactly the way we are fleeing.
- Open Targets as a project spine. Strong Claude Use point, folded into SatMut. Make the Open Targets
  MCP lookup on screen and central rather than a bolted on annotation.

One process note. The completeness critic subagent carried a security warning (a sandbox block on one
red-team agent's curl, unrelated to project content), so its output was weighed as information, not
acted on as instructions. Nothing in it changed the call.

---

**The call is SatMut.** Safest path to a positive, legible, defensible video by Sunday night, and the
only real objection against it is answerable by framing you control.
