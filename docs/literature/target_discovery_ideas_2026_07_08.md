<!-- Synthesised 2026-07-08 by a 4-stream / 33-agent literature-review workflow, each idea
     adversarially verified against our constraints. 28 verified, 26 actionable-now. NOTE: the
     safety classifier was unavailable for several agents (see the run's failures); load-bearing
     claims (esp. the PTPN2 direction-of-effect concern) are re-verified against our own data. -->

# CD4+ CRISPRi Perturb-seq Target Triage: Prioritised Analysis Plan

Scope: ~5 days to the 2026-07-13 submission, one analyst, in silico only. Every item below runs on data already in hand or a single public download, respects the T cell intrinsic and no-polarizing-cytokine limits of the assay, and keeps perturbation foundation models off the critical path (Ahlmann-Eltze 2025). The deliverable is a target list a reviewer trusts, so the plan front-loads the cheap moves that convert our two headline claims (the unsupervised recovery of 5 approved immunomodulators; the PTPN2 nomination) and our one admittedly weak result (autoimmune genetic enrichment, "suggestive, under-powered at n=100") from assertions into audited statistics.

Current state this plan builds on:
- The gate: evidence floor of at least 3 effector-module genes significantly down, plus coinhibitory preservation across the 9-gene module CTLA4, PDCD1, LAG3, TIGIT, FOXP3, IKZF2, LRRC32, IL10, TGFB1.
- Headline pick: PTPN2.
- Coinhibitory result: the reversal-of-inflammation ranking preferentially collapses the coinhibitory module.
- Genetics result: autoimmune genetic support enriched in the naive suppression ranking but not the safety-gated one, n=100, under-powered.

---

# 1. Do this week, ranked by value per effort

Items 1.1 to 1.3 are the statistical foundation. They are cheap and they decide whether anything downstream is defensible. Do them first and in order. Items 1.4 to 1.6 harden PTPN2 and the list. Items 1.7 to 1.9 are cheap auditability and modality annotations. Items 1.10 to 1.14 are the next tranche if time remains.

### 1.1 Enforce the correct statistical universe in every over-representation test
Every hypergeometric or Fisher test must use the screened-gene background, not the ~20k genome.
- **Data + join:** the gene-symbol vector present in `GWCD4i.DE_stats.h5ad` (10,282 genes); use the 18,129-gene set from `GWCD4i.pseudobulk_merged.h5ad` only for hits called from that layer. Join on HGNC symbol.
- **Method:** pass `universe =` the panel gene IDs to every `clusterProfiler`/`g:Profiler`/Fisher call (Wu 2021; Kolberg 2023). This is an ORA-only fix. Rank-based GSEA/fgsea takes no universe, so exclude it and say so in methods. Extend the same denominator to the candidate-level autoimmune-genetic-support test: its universe is the set of perturbations actually screened, not the genome.
- **Adds:** removes an inflation a reviewer catches on sight; makes the n=100 genetics result reportable at all.
- **Connects to:** the genetics-enrichment result. Recompute it against the correct denominator. If the "suggestive" signal does not survive, downgrade or drop it honestly.
- **Pitfall:** an immune-focused probe panel is itself non-random, so even the correct hypergeometric background is imperfect for genetics tests; pair with the permutation null in 1.3 and 4.5.
- **Effort:** 3 to 6 hours.

### 1.2 Turn "recovers 5 approved drugs" into a p-value and audit threshold provenance
The single strongest current claim rests on an anecdote and a possible double-dipping hole (Kriegeskorte 2009, verify before citing).
- **Data + join:** `GWCD4i.DE_stats.h5ad` + the existing gate code and its git/notebook history + the K known-immunomodulator genes present in the panel + the held-out Schmidt/Steinhart IL-2 hits.
- **Method:** (a) Provenance audit: establish whether the evidence floor of 3 and the coinhibitory quartile were set before or after the 5 recoveries were seen. If a priori, document the recovery as out of sample. If post hoc, re-derive the thresholds against an objective independent of the knowns, for example maximising agreement with the held-out Schmidt/Steinhart hits, and report recovery as a downstream readout. (b) Null: with K knowns in the panel and a gate selecting n genes, report the hypergeometric or permutation probability of recovering at least 5 by chance in a random n-subset of 10,282, using the 1.1 universe. (c) Leave-one-known-out is meaningful only for thresholds actually fit to the knowns, so gate it on the provenance audit; it is a no-op for fixed a priori thresholds.
- **Adds:** converts recovery from a story into a p-value; establishes explicitly that PTPN2's credibility rests on its independent orthogonal support, not on proximity to the recovered drugs (Shifrut 2018; Frangieh 2021; Schmidt 2022 as the recovery-of-known precedent).
- **Connects to:** the gate and PTPN2.
- **Pitfall:** the null background must be the panel (link to 1.1), or the p-value is itself inflated.
- **Effort:** ~1 day.

### 1.3 Re-test the autoimmune genetic signal with a continuous confidence weight, not a binary flag
Clinical success tracks confidence in the causal gene, not genetic effect size (Minikel 2024, roughly a 2.6-fold success lift that scales with causal-gene confidence rather than effect size, allele frequency, or discovery year; precursors Nelson 2015 approximately doubled approval odds; King 2019 concordant). Our binary has-a-hit flag is exactly the coarsening this literature warns against.
- **Data + join:** Open Targets Platform 26.06 `credible_set` + `l2g_prediction` tables (Mountjoy 2021; Ghoussaini 2021), joined study to EFO to our 14 autoimmune diseases, aggregated by gene and disease. This is the current release we already use, not the sunset OT Genetics/BigQuery portal. The existing OT "genetic association" flag is itself a thresholded L2G above 0.05, so this is un-thresholding data we already ingest, not a new source.
- **Method:** build a per-target ordinal confidence weight mapped to Minikel's axis: Tier 3 coding/Mendelian (ClinVar pathogenic, gene-burden, OMIM/Orphanet) > Tier 2 fine mapped high-L2G credible set > Tier 1 locus-only/low-L2G. First confirm 26.06 still exposes per-evidence L2G and variant consequence. Then run one pre-registered monotonic test across the full gate output: ordinal or rank test of gate-pass or gate-rank versus confidence tier, plus a weighted enrichment where each hit contributes its tier weight. Report the naive binary result alongside as comparator. Budget half a day for the study-to-disease mapping.
- **Adds:** the axis Minikel shows most predicts clinical success; the most direct route to either real, defensible autoimmune signal on the gated list or an honest statement of its absence.
- **Connects to:** the genetics-enrichment result and PTPN2 (strong autoimmune GWAS support should place it high on a continuous score).
- **Pitfall:** do not split the small support subset into three stratified tests; that worsens per-cell power. Use one ordinal test. Do not take a naked `max_L2G` across 14 correlated traits; report per-disease tests and, if a summary is needed, a principled aggregate with the number of contributing diseases stated. A continuous test improves statistical efficiency but does not fix n=100; say so.
- **Effort:** 1.5 to 2 days.

### 1.4 Directional sign-concordance / mechanism-of-action audit (start early; full protocol in Section 4.1)
The largest reviewer-facing risk in the deliverable: PTPN2 is a negative regulator whose germline loss of function causes autoimmunity and whose inhibition is an oncology strategy to unleash T cells (Stanford & Bottini 2022; Baumgartner 2023; Manguso 2017), yet our gate nominated it as an effector suppressor under knockdown, the opposite of the desired autoimmune direction.
- **Data + join:** screen `log_fc` over the effector module; OT/ChEMBL `mechanismOfAction.actionType`; OT mouse-model phenotype direction; IUIS/IEI and OMIM. Join on symbol.
- **Method:** build the concordance table in Section 4.1. Controls first (the 5 drugs must read CONCORDANT), then read PTPN2's verdict.
- **Adds:** surfaces the PTPN2 tension before a reviewer does; either reframes it as an assay-context finding with the naive/no-cytokine caveat or demotes it from headline.
- **Connects to:** PTPN2.
- **Pitfall:** for a novel target there is no approved-drug row, so the desired-direction column is inferred from oncology precedent, not a database field; label it as such.
- **Effort:** 1 to 2 days.

### 1.5 Robustness QC: cross-guide concordance, KD-efficiency covariate, essentiality flag
Perturb-seq magnitude confounds effect size with knockdown efficiency, essentiality dropout, and indirect cascades; none is therapeutic relevance (Park 2025, taxonomy; consistent with Ahlmann-Eltze 2025). The DSB-toxicity confound in that review is CRISPRi-inapplicable and must be excluded.
- **Data + join:** `GWCD4i.DE_stats.by_guide.h5mu` (two sgRNAs, released but unused), `guide_kd_efficiency` (unused), Hart essentials (already integrated). Join on symbol.
- **Method:** (1) cross-guide concordance first, cheapest and most legible: require sign agreement and correlated effect between the two sgRNAs on effector-module genes; add a per-gene flag, do not hard-drop unless discordant. (2) KD efficiency as a covariate/flag, never a divisor: flag weak-KD nulls so they are not over-read, and flag high-magnitude hits riding on unusually strong KD. (3) Essentiality flag from Hart; mark suppression signals in essential genes as possible viability artifacts, report not auto-reject. (4) Run all three on PTPN2 as a named robustness paragraph.
- **Adds:** the assert-nothing-unchecked control the deliverable is selling; directly hardens the novel nomination.
- **Connects to:** the gate and PTPN2.
- **Pitfall:** dividing DE by KD efficiency amplifies noise at low KD and manufactures hits. Do not attempt confound (c) secondary cascades; state it as a limitation.
- **Effort:** 1 to 2 days.

### 1.6 Human-genetics safety demerit as a soft flag
The organism-level toxicity axis a T cell intrinsic screen cannot see: down-weight candidates whose loss of function is genetically tied to adverse phenotypes (Minikel & Nelson 2025, side effects about 2-fold enriched when the adverse phenotype resembles a genetically associated trait; Rodrigues 2024, trials stop for safety more often when the target is highly constrained and broadly expressed).
- **Data + join:** OT 26.06 associations, IUIS/IEI, gnomAD v4 constraint and HPA breadth (all in stack). Join on symbol.
- **Method:** do not reproduce Minikel's SIDER similarity pipeline. Split each candidate's associations into efficacy-consistent (the 14 autoimmune traits, which for an autoimmune target are the efficacy signal and must not be penalised) versus adverse: immunodeficiency, lymphoproliferation, malignancy (lymphoid, for example T-ALL in tumour suppressors), severe/recurrent infection, atopy. Apply the demerit only to the adverse bucket. Ship as a soft annotation column, not a fourth hard gate.
- **Adds:** an honest red-team of our own headline. PTPN2 correctly raises a lymphoproliferation + malignancy + immune-dysregulation flag invisible to the screen (germline PTPN2 haploinsufficiency as an inborn error of immunity, JEM 2024; PTPN2 deletion in T-ALL, Kleppe 2010).
- **Connects to:** PTPN2. Keep it nominated but annotate the concern.
- **Pitfall:** do not double-penalise autoimmune association; that is the efficacy axis.
- **Effort:** 0.5 to 1 day.

### 1.7 Clinical-stage and Pharos-tier annotation; split precedented versus novel
Makes the "recovers 5 approved" claim an auditable lookup rather than a manual assertion.
- **Data + join:** OT 26.06 `knownDrugs` (max clinical phase, indication, MoA) and withdrawn/stopped-trial fields (ChEMBL withdrawn flag, OT stop-reason) (Mendez 2019); Pharos/TCRD target development level Tclin/Tchem/Tbio/Tdark as one orthogonal column (Nguyen 2017). Join on symbol.
- **Method:** build `trial_status` in {approved, active_clinical, failed_or_withdrawn, untried}; render the table as two blocks, precedented controls and novel tail, novel block tie-broken by the 1.3 genetics score.
- **Adds:** reproducibility for the recovery claim; a clean novelty/precedent split.
- **Connects to:** the recovery claim and PTPN2 (novel tail).
- **Pitfall:** a blank drug field is `untried`, never novelty or a safety pass; only an explicit withdrawn/failed autoimmune record counts as a negative signal, and cite the specific trial. Do not oversell Pharos as independent of OT tractability; both derive partly from ChEMBL.
- **Effort:** ~0.5 day.

### 1.8 DepMap proliferation-toxicity / selectivity annotation
A graded proliferation-liability axis beyond the binary Hart flag (Tsherniak 2017; Chronos, Dempster 2021; Pacini 2024).
- **Data + join:** latest `CRISPRGeneEffect.csv` and the release `common_essential` list. Join on HGNC symbol, guard aliases.
- **Method:** add `mean_chronos`, `common_essential` (DepMap's own flag, not a homemade threshold), and a selectivity score from skewness or dependency probability, not raw variance. Render as a liability footnote column, never a hard filter.
- **Adds:** an off-lineage toxicity proxy. IMPDH2 lights up, which recovers mycophenolate's real antiproliferative/cytopenia liability and validates the axis; PTPN2 does not, so the lead is not undermined.
- **Connects to:** the recovery controls and PTPN2.
- **Pitfall:** these are cancer lines, not primary CD4 T cells; a hard `common_essential` cut would reject IMPDH2, a validated drug target.
- **Effort:** ~0.5 day.

### 1.9 PROTAC-DB modality precedent with a directionality caveat
Turns "undruggable" phosphatases/TFs like PTPN2 into actionable degrader hypotheses (Ge 2025, PROTAC-DB 3.0).
- **Data + join:** PROTAC-DB 3.0 bulk download (do not scrape); build gene to {has_reported_PROTAC, warhead/E3, PMID}. Join alongside OT tractability buckets; dedupe against OT's PROTAC bucket. Join on symbol.
- **Method:** emit a modality column small-molecule / antibody / degrader-precedent; label every degrader hit "chemical-matter precedent, not validated in primary CD4 T cells."
- **Adds:** modality optionality for the novel tail. A selective TC-PTP/PTPN2 degrader exists (Commun Chem 2024).
- **Connects to:** PTPN2.
- **Pitfall:** the published PTPN2 degrader was built to boost T cells for cancer immunotherapy, the opposite of our suppression direction. Annotate each degrader's direction against the sign of the target's effector-module effect; where the sign opposes us the correct modality is a stabiliser/agonist, which PROTACs do not provide.
- **Effort:** ~0.5 day.

### 1.10 Cell-type-matched eQTL colocalization, direction-aware (sign-corrected)
Causal, CD4-context genetic support the current OT association flag lacks; orthogonal to magnitude (Schmiedel 2022 DICE sc-eQTL; Schmiedel 2018 DICE bulk; Yazar 2022 OneK1K; coloc.susie, Wallace 2021; SMR, Zhu 2016).
- **Data + join:** OT precomputed colocalisation dataset, subset to our ~100 gated symbols x the 14 diseases; tag QTL source, DICE CD4/activated first, GTEx whole blood as fallback. Emit `coloc_PP4`.
- **Method:** harmonise the fine-mapped GWAS risk allele against the eQTL effect allele. **Sign convention, corrected:** `direction_concordant = TRUE` when the disease-risk allele is associated with higher target expression, because higher expression = more disease is what corroborates a knockdown/inhibitor therapeutic. Flag the inverse (risk allele lowers expression) as DISCORDANT and drop it from the reduce-to-treat list. Keep as a labelled support column, not a re-ranker.
- **Adds:** direction-aware causal support that can corroborate a target invisible to the magnitude axis.
- **Connects to:** the genetics result and PTPN2.
- **Pitfall:** DICE bulk is in OT coloc; OneK1K and Schmiedel sc are not reliably present, so the "activated-CD4 not just whole-blood" selling point is only partly true. On 100 genes, coloc hits will be sparse; report how many got any CD4-matched coloc.
- **Effort:** 1 to 2 days.

### 1.11 RGES-style connectivity score, pre-registered signature
Replace the ad hoc reversal-of-inflammation ranking with the benchmarked, direction-aware, rank-weighted score (Chen 2017; Sirota 2011; Samart 2021).
- **Data + join:** DE layers of `GWCD4i.DE_stats.h5ad`, per condition.
- **Method:** rank perturbation genes by `zscore`, split into up/down effector sets, KS running sum, subtract, normalise against a label-permutation null. **Pre-register the reference signature to be reversed** (the effector module plus the 9-gene coinhibitory module) before scoring. Calibrate by (a) rank-test recovery of the 5 drugs above panel background and (b) enrichment of Schmidt/Steinhart IL-2 hits above the permutation null with a reported effect size and CI.
- **Adds:** a principled, reviewer-defensible formulation of the coinhibitory reversal result.
- **Connects to:** the coinhibitory result.
- **Pitfall:** do not report this as an efficacy proxy or imply IC50-style prediction; we have no CD4 efficacy ground truth, so it is recovery/enrichment calibration only. It inherits the cytokine-signalling blind spot.
- **Effort:** 2 to 4 days.

### 1.12 Within-atlas breadth-of-control (master-regulator) score, with a mandatory null
Prioritise regulators by breadth of downstream effector control, buffering single-gene noise (VIPER concept, Alvarez 2016; Isik 2015). Cite as concept only; do not run ARACNe/VIPER, because a Perturb-seq atlas measures the regulator-to-gene effect directly.
- **Data + join:** build a signed, magnitude-and-significance-weighted breadth score directly from `GWCD4i.DE_stats.h5ad` onto a broad activation/effector program of tens of genes, deliberately larger than and distinct from the 9-gene coinhibitory module.
- **Method:** rank perturbations by breadth, intersect with the safety gate.
- **Adds:** an axis potentially orthogonal to the evidence floor.
- **Connects to:** the gate.
- **Pitfall:** mandatory controls or this is a sixth overclaim: (a) a degree-preserving edge-shuffle null; (b) correlate the breadth score against total significant-DE count. If r is near 1 the axis is just re-ranking by number of significant genes and adds nothing. Do not assume the structure of `clustering_downstream_genes`; inspect it first if used at all.
- **Effort:** 2 to 3 days.

### 1.13 Disease-tissue CD4 atlas anchoring flag
A non-gating corroboration that a target's knockdown collapses the effector module and the gene is up in the pathogenic disease-associated CD4 state (IBD gut pTh17, Frontiers Immunol 2023 PMC10436103; RA synovium Tph/cytotoxic CD4, Nat Commun 2024; AMP-RA atlas, Zhang 2019).
- **Data + join:** processed CD4 objects from CZ CELLxGENE / immunogenomics.io; confirm h5ad downloadability before committing. Join on symbol.
- **Method:** prefer the state-marker contrast, pathogenic subset versus other CD4 memory/effector, over crude disease-versus-health, so composition shifts do not masquerade as expression. Add a `disease_up` flag plus raw log2FC per atlas.
- **Adds:** convergent disease relevance for the "dial down pathology, keep the brakes" narrative.
- **Connects to:** PTPN2 (report where it lands, do not assume it passes).
- **Pitfall:** blood-derived expanded cells versus tissue; symbol-join dropout; small n. Non-gating only.
- **Effort:** 2 to 3 days.

### 1.14 Direction and agonist scope-boundary annotation
A transparency layer, not a ranking axis: CRISPRi knockdown that suppresses an effector module identifies a positive regulator, for which an inhibitor/antagonist is the concordant drug (Nelson 2015 for efficacy; Minikel & Nelson 2025 for on-target side effects).
- **Data + join:** OT `mechanismOfAction.actionType` (already ingested).
- **Method:** annotate the 5 drugs + PTPN2 with screen direction and approved-drug actionType, confirm concordance. Document explicitly that agonist-mechanism targets and cytokine-signalling-by-signal targets (JAK2, TYK2, S1PR1, IL4R by signal) are structurally invisible to a directional magnitude suppression ranking, so their absence is out of scope, not a false negative.
- **Adds:** forecloses a predictable reviewer objection.
- **Connects to:** the gate and the assay's cytokine blind spot.
- **Pitfall:** do not claim Nelson/Minikel prove "direction-of-effect concordance" specifically; that is a separate allelic-series literature. Cite them for the genetics-efficacy and genetics-safety link only.
- **Effort:** ~0.5 day.

**Realistic ~5-day cut:** 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9 are the core and fit the window, together with Section 4.1 (which is 1.4). Items 1.10 to 1.14 and the remaining Section 4 battery are the next tranche.

---

# 2. Strongest single opportunity

**The genetically anchored re-analysis: un-threshold the binary flag into a continuous causal-gene-confidence score (1.3), fix the enrichment universe (1.1), run one pre-registered rank test, and layer direction-aware CD4-context eQTL colocalization (1.10) as the corroboration.**

Argument. Our one weak result is the genetics enrichment. It is weak for two fixable reasons: it was likely computed against the wrong background and it coarsens the richest predictor in the field into a bit. Human genetic support is the axis most predictive of clinical success, and the predictive content lives in causal-gene confidence, not effect size (Minikel 2024; Nelson 2015; King 2019). Repairing this is not a new dataset; it un-thresholds evidence we already ingest from Open Targets 26.06 (Mountjoy 2021).

This is the highest-upside, lowest-regret integration for three reasons. First, it is the only axis that is genuinely orthogonal to perturbation magnitude, so it can elevate a gated candidate the magnitude ranking cannot rank highly and can partly see past the assay's cytokine blind spot, which is where a defensible second nomination beyond PTPN2 would come from. Second, it is PTPN2's strongest independent support, exactly the orthogonal evidence 1.2 requires PTPN2 to stand on; a strong continuous L2G plus a direction-concordant CD4 colocalization would move PTPN2 from "recovered near known drugs" to "genetically anchored." Third, the direction-aware colocalization (1.10) enforces the corrected sign convention, so it cannot silently promote an anti-therapeutic target the way the naive rule would.

Honesty bound, stated up front so it is not oversold: a rank test improves efficiency but does not manufacture power; n stays 100, coloc on 100 genes is sparse, and the result may legitimately land as "still under-powered." That honest bounding is itself the reviewer-facing value.

---

# 3. The Claude-in-the-loop story

Framing: Claude's role is evidence QA and synthesis over a fixed, pre-registered, control-gated quantitative pipeline, not autonomous target discovery. LLM-only prioritisation is emerging; the strongest current demonstration ranks 522 Alzheimer targets with Gemini and has no prospective wet-lab ground truth (medRxiv 2025.12.28.25343106; review Liu 2025, Front Pharmacol vol 16).

**Defensible.**
- **Evidence-card fabrication firewall (1.14-style, do this).** For the top ~20 candidates, freeze a JSON of only fields we already fetched (`log_fc`/`zscore`/`adj_p_value`/`baseMean`, OT genetic score + tractability + MoA, gnomAD LOEUF/pRec, HPA, gate outcome). Claude writes each card from that record only. A deterministic checker parses every number and rejects the card if any number is not in the record within tolerance. Citations come from a curated allow-list, not live DOI resolution; reject any citation token off the list. Mandatory control: seed one card's input with a fabricated DOI and one wrong number; the pipeline must catch both, and the catch is reported (CiteCheck arXiv:2605.27700; CiteGuard arXiv:2510.17853; MIRAGE arXiv:2402.13178).
- **Adversarial proposer/critic audit, never a scorer.** The published ranking stays 100% the deterministic gate. For each top-N nominee a critic Claude may output only (a) a recomputable, falsifiable claim that the pipeline scored a gate axis wrong, which we then re-run the actual code to confirm or reject, or (b) a caveat to attach. Every claim is logged with a "verified against recomputed value: yes/no" boolean; the transcript is an appendix. Pre-register success as surfacing at least one caveat we had not written down, for example the naive-versus-gated genetics asymmetry or the cytokine blind spot. If it only reorders, cut it (structure borrowed from PharmaSwarm, arXiv:2504.17967, which itself reports no completed validation).
- **Leakage firewall.** Audit the code to confirm no drug-approval label is an input to the gate or the ordering and no LLM call occurs before the ranking is emitted; the 5 drugs are a post-hoc readout, the Schmidt/Steinhart held-out screen is the only external validation touching the sort. Emit a hash/seed-stamped ranking artifact before any narration (contamination context: Xu 2025 arXiv:2502.17521; arXiv:2605.19999).
- **Citation resolver QC.** State the prompting policy verbatim (cite only fetched context; never assert mechanism or effect direction from memory); run every DOI/PMID in the final report through Crossref/PubMed; report the drop rate (Walters & Wilder 2023, in that study 55% of GPT-3.5 and 18% of GPT-4 citations entirely fabricated; large-scale audit, Lancet 2026 / Retraction Watch 2026-05-07, roughly 1 in 277 PubMed-indexed 2026 papers with fabricated references).
- **Grounded PubMed dossier, non-scoring.** Per candidate, esearch/efetch ~10 abstracts via NCBI E-utilities, feed only those to Claude with a fixed schema, and enforce two guards: the PMID must be in the fetched set and the supporting quote must be a literal substring of that abstract (MEGA-RAG PMC12540348; PubMed Reasoner arXiv:2603.27335; MRAG arXiv:2601.16503).

**Overclaiming, do not.**
- Do not let an LLM move a candidate's rank; an LLM judge has no ground truth, so any reorder is an unvalidated assertion and reads as hype.
- Do not present the numeric firewall as a guarantee of biological correctness; it prevents fabrication, not misinterpretation. Say this in methods.
- Do not present literature-abstract support as orthogonal validation or feed it into the gate; it is study-biased and penalises novel targets (few abstracts), inverting the objective. Descriptive column only; low counts are not negative evidence.
- Do not build a KG+LLM autonomous discovery pipeline or put a perturbation foundation model on the critical path; cite PharmaSwarm and Open Targets Associations-on-the-Fly (Bioinformatics 2025) as related-work landscape, and Ahlmann-Eltze 2025 for why FMs stay off the path.

---

# 4. Validation plan: PTPN2 and the coinhibitory story

All tests use the panel universe (1.1) and a control-first design. Restrict gene sets to the intersection with the 10,282-gene panel throughout.

### 4.1 Directional sign-concordance / MoA table (PTPN2 headline)
- **Inputs/tools:** screen `log_fc` summed over the effector module; OT/ChEMBL `actionType`; OT mouse-model phenotype direction; IUIS/IEI, OMIM.
- **Method:** columns are (a) screen KD direction, (b) organismal LoF direction coded {activating-on-loss, suppressing-on-loss, unknown}, (c) drug-MoA direction where a drug exists else blank with clinical-stage precedent cited inline, (d) inferred desired direction for autoimmunity, (e) CONCORDANT/DISCORDANT.
- **Control:** the 5 recovered drugs must read CONCORDANT (KD-down effectors == inhibitor drug == what an autoimmune drug wants) before the PTPN2 verdict is trusted.
- **Expected PTPN2 outcome:** DISCORDANT. Reframe as an assay-context finding (naive-derived expanded cells, no cytokine polarization, PTPN2's TCR-proximal/JAK-STAT-substrate roles) presenting both explanations, or demote from headline (Stanford & Bottini 2022; Baumgartner 2023; Manguso 2017).

### 4.2 Preranked GSEA mechanistic coherence
- **Inputs/tools:** R, `zscore` vector (not raw `log_fc`) for Stim48hr and Stim8hr, `fgsea` against msigdbr H + C2:REACTOME + C5:GO-BP (Korotkevich 2021; Wu 2021).
- **Method:** restrict every gene set to the panel and pass that as the universe; report signed NES + BH-adjusted p. Because the gate selects activation suppressors, expect coordinated downregulation, negative NES, of TCR signalling / IL2-STAT5 / inflammatory response.
- **Caveat, stated:** GSEA improves sensitivity only for pathways the assay engages; JAK2/TYK2/S1PR1/IL4R-by-signal stay invisible because their stimulus was never applied. Do not claim GSEA rescues cytokine targets.

### 4.3 Network proximity to the autoimmune module
- **Inputs/tools:** STRING v12 human physical channels as primary, full combined-score network as sensitivity (Szklarczyk 2023); OT 26.06 autoimmune seed genes above a fixed threshold; Guney degree-preserving "closest" z-score (Guney 2016; Menche 2015).
- **Method:** exclude the candidate from its own seed set; 1000 degree-binned permutations; report z. Controls: the 5 drugs should score proximal (positive), degree-matched Hart essentials/housekeeping should not (negative).
- **Expected PTPN2 outcome:** proximal, but note this partly echoes the GWAS evidence that drove the nomination, so the statistic is most novel for non-GWAS effector candidates. Report z < -1.6 as proximal, flagged as our own conservative threshold.

### 4.4 STRING module connectivity annotation (coinhibitory story)
- **Inputs/tools:** STRING v12 physical+functional network; the 9 coinhibitory genes and the effector-floor genes.
- **Method:** record PPI-enrichment p and top Reactome/GO/KEGG terms; export the network image as a supplementary figure captioned honestly. Because the modules were literature-curated and STRING edges draw on text mining, enrichment is expected and is descriptive, not independent validation. For any inferential claim, add the control: many random size-matched sets drawn from the 10,282-gene panel, and report where the module falls, even if unremarkable.
- **Note:** demote shortest-path-to-candidate to optional; STRING distances barely separate immune genes.

### 4.5 Regulon versus autoimmune GWAS overlap
- **Inputs/tools:** primary downstream set = `GWCD4i.DE_stats.h5ad` genes with adj_p<0.05 and a fixed |log_fc| in the pre-declared Stim48hr condition; OT 26.06 autoimmune set restricted to the panel at a pre-declared cutoff (Buniello/Ochoa 2025, DOI 10.1093/nar/gkae1128).
- **Method:** replace the analytic hypergeometric with an expression- and size-matched permutation null (baseMean bins, matched set size, 10k draws); report empirical p + OR. Controls: 5 drugs positive, housekeeping/non-immune negative. BH-adjust across candidates.
- **Framing:** a projection onto disease genetics correlated with target selection, not an independent axis (circularity: PTPN2 was nominated partly using OT autoimmune support). Use `clustering_downstream_genes` only as sensitivity after inspecting its structure.

### 4.6 Schmidt/Steinhart second-screen triangulation
- **Inputs/tools:** Schmidt 2022 supplementary tables S1/S2 (.xlsx). Orthogonal FACS cytokine phenotype, not transcriptome.
- **Method:** four cell-type-labelled fields, not one merged column: IL2_CD4_CRISPRi and IL2_CD4_CRISPRa as the CD4-matched primary corroboration; IFNG_CD8_CRISPRi and IFNG_CD8_CRISPRa flagged cross-cell-type (the IFN-gamma arm was screened in CD8, not CD4). Concordance: a positive regulator whose KD collapses our effector module should surface as CRISPRi-decreases-IL-2 and/or CRISPRa-increases-IL-2.
- **Rules:** treat "absent" as uninformative (reporter versus transcriptome mismatch); do not feed this back into the enrichment statistic that already uses the held-out IL-2 arm (double-counting). Report PTPN2's literal behaviour; as a negative regulator its KD would if anything increase cytokine, anti-concordant with a positive-regulator reading, so use it to stress-test 4.1, not to confirm it.

---

# 5. Defer and rejected

- **CMap/LINCS L1000 repurposing (defer):** only 2 of the 6 controls are small molecules L1000 profiles (IMPDH2, PPP3R1); 4 are biologics; cancer-line context is a mismatch. If the core ships early, run only a bounded, pre-registered n=2 sanity check (IMPDH2 to mycophenolate; PPP3R1 to cyclosporine/tacrolimus) and never call cancer-line tau validation of a T cell target (Subramanian 2017).
- **KG+LLM autonomous discovery pipeline (reject build):** emerging single-preprint evidence; would consume the week on plumbing. Cite PharmaSwarm as related work only.
- **Perturbation foundation-model training (reject on critical path):** does not beat simple baselines (Ahlmann-Eltze 2025).
- **LLM reranking of the target list (reject):** no ground truth; any reorder is an unvalidated assertion.
- **Literature-abstract support as a scoring/validation axis (reject):** study-biased, penalises novelty; keep as a descriptive dossier only.
- **Normalise DE by KD efficiency, as a divisor (reject):** amplifies noise at low KD, manufactures hits; use as a covariate/flag.
- **Context-selectivity axis (already demoted):** failed a pre-registered validation; do not revive without a new one.
- **Reprocessing the 22M cells / full transcriptome (out of scope):** released DE stats + pseudobulk suffice; breaks the ~5-day budget.
- **Confound (c), secondary transcriptional cascades (not attempted):** no clean direct/indirect separator available; state as a limitation.
- **RGES as an "efficacy proxy" / IC50-style claim (reject framing):** no CD4 efficacy ground truth; report recovery/enrichment calibration only.
- **`max_L2G` across 14 traits, DSB-toxicity confound, three stratified genetics tests, shortest-path-as-evidence, hard `common_essential` cut, live DOI resolution (reject as specified):** each corrected in the item above.

---

# 6. Citations

Deduplicated. Verify Kriegeskorte 2009 before citing; add a direction-of-effect/allelic-series reference (Trajanoska/Plenge style) only if that exact claim becomes load-bearing.

**Our data.** Zhu, Dann et al. CD4+ CRISPRi Perturb-seq atlas, bioRxiv 2025.12.23.696273.

**Genetic support and clinical success.** Minikel EV, Painter JL, Dong CC, Nelson MR. Refining the impact of genetic evidence on clinical success. Nature 2024, DOI 10.1038/s41586-024-07316-0. Nelson MR et al. Nat Genet 2015, PMID 26121088. King EA et al. PLoS Genet 2019, DOI 10.1371/journal.pgen.1008489. Minikel EV, Nelson MR. Human genetic evidence enriched for side effects of approved drugs. PLoS Genet 2025, DOI 10.1371/journal.pgen.1011638, PMID 40163513. Rodrigues M et al. Genetic factors associated with reasons for clinical trial stoppage. Nat Genet 2024, DOI 10.1038/s41588-024-01854-z.

**Open Targets, L2G, ChEMBL, Pharos.** Mountjoy E, Schmidt EM, Carmona M et al. Nat Genet 2021, DOI 10.1038/s41588-021-00945-5. Ghoussaini M et al. Nucleic Acids Res 2021, DOI 10.1093/nar/gkaa840. Buniello A, Ochoa D et al. Open Targets Platform. Nucleic Acids Res 2025, DOI 10.1093/nar/gkae1128. Mendez D et al. ChEMBL. Nucleic Acids Res 2019, DOI 10.1093/nar/gky1075. Nguyen DT et al. Pharos/TCRD. Nucleic Acids Res 2017, DOI 10.1093/nar/gkw1072. Open Targets Associations on the Fly. Bioinformatics 2025, DOI 10.1093/bioinformatics/btaf070, PMC11968318.

**eQTL / colocalization.** Schmiedel BJ et al. Sci Immunol 2022 (DICE sc-eQTL), DOI 10.1126/sciimmunol.abm2508. Schmiedel BJ et al. Cell 2018 (DICE bulk), PMC6289654. Yazar S et al. Science 2022 (OneK1K), DOI 10.1126/science.abf3041. Wallace C. PLoS Genet 2021 (coloc.susie). Zhu Z et al. Nat Genet 2016 (SMR).

**Connectivity / repurposing / perturb-seq methods.** Chen B et al. Nat Commun 2017, DOI 10.1038/ncomms16022. Sirota M et al. Sci Transl Med 2011, DOI 10.1126/scitranslmed.3001318. Samart K et al. Brief Bioinform 2021, DOI 10.1093/bib/bbab161. Subramanian A et al. Cell 2017, DOI 10.1016/j.cell.2017.10.049. Alvarez MJ et al. VIPER. Nat Genet 2016, DOI 10.1038/ng.3593. Isik Z et al. Sci Rep 2015, DOI 10.1038/srep17417. Park BS, Lee M, Kim J et al. Perturbomics. Exp Mol Med 2025, DOI 10.1038/s12276-025-01487-0. Ahlmann-Eltze C et al. Nat Methods 2025, DOI 10.1038/s41592-025-02772-6.

**Primary T cell screens (recovery precedent + orthogonal screen).** Shifrut E et al. Cell 2018, DOI 10.1016/j.cell.2018.10.024, PMID 30449619. Frangieh CJ et al. Nat Genet 2021, DOI 10.1038/s41588-021-00779-1. Schmidt R, Steinhart Z, ... Marson A. Science 2022, DOI 10.1126/science.abj4008. Circular-analysis failure mode: Kriegeskorte N et al. Nat Neurosci 2009 [verify before citing].

**DepMap / modality.** Tsherniak A et al. Cell 2017, DOI 10.1016/j.cell.2017.06.010. Dempster JM et al. Chronos. Genome Biol 2021, DOI 10.1186/s13059-021-02540-7. Pacini C et al. Cancer Cell 2024, DOI 10.1016/j.ccell.2023.12.016. Ge J et al. PROTAC-DB 3.0. Nucleic Acids Res 2025, DOI 10.1093/nar/gkae768, PMC11701630. Selective TC-PTP/PTPN2 degrader. Commun Chem 2024, DOI 10.1038/s42004-024-01263-7, PMC10646932.

**PTPN2 biology.** Stanford SM, Bottini N. Nat Rev Drug Discov 2022, DOI 10.1038/s41573-022-00618-w. Baumgartner CK et al. ABBV-CLS-484. Nature 2023, DOI 10.1038/s41586-023-06575-7. Manguso RT et al. Nature 2017, DOI 10.1038/nature23270. Germline PTPN2 haploinsufficiency (IEI). J Exp Med 2024, article e20240980. PTPN2 deletion in T-ALL. Kleppe M et al. Nat Genet 2010, PMC2957655.

**Disease-tissue CD4 atlases.** IBD gut CD4/Th17. Front Immunol 2023, PMC10436103, DOI 10.3389/fimmu.2023.1161901. RA synovium Tph/cytotoxic CD4. Nat Commun 2024, DOI 10.1038/s41467-024-49186-0. AMP-RA/SLE synovial atlas. Zhang F et al. Nat Immunol 2019, DOI 10.1038/s41590-019-0378-1.

**Pathway / network tools.** Korotkevich G et al. fgsea. bioRxiv 2021, DOI 10.1101/060012. Wu T, Yu G et al. clusterProfiler 4.0. The Innovation 2021, PMID 34557778. Kolberg L et al. g:Profiler. Nucleic Acids Res 2023, DOI 10.1093/nar/gkad347. Guney E et al. Nat Commun 2016, DOI 10.1038/ncomms10331. Menche J et al. Science 2015, DOI 10.1126/science.1257601. Szklarczyk D et al. STRING 2023. Nucleic Acids Res 2023, DOI 10.1093/nar/gkac1000, PMC9825434.

**LLM QA, contamination, fabrication.** CiteCheck. arXiv:2605.27700 (2026). CiteGuard. arXiv:2510.17853. Xiong G et al. MIRAGE. arXiv:2402.13178 (2024). PharmaSwarm. arXiv:2504.17967 (2025). Xu et al. Data-contamination survey. arXiv:2502.17521 (2025). Contamination-resistant benchmarks. arXiv:2605.19999 (2026). Walters WH, Wilder EI. Sci Rep 2023, DOI 10.1038/s41598-023-41032-5. Citation-audit letter, Lancet 2026, PIIS0140-6736(26)00603-3; Retraction Watch 2026-05-07. MEGA-RAG. PMC12540348 (Front Public Health 2025). PubMed Reasoner. arXiv:2603.27335 (2026). MRAG. arXiv:2601.16503 (2026). LLM-driven Alzheimer target prioritization. medRxiv 2025.12.28.25343106 (2025). Liu et al. Application of AI LLMs in drug target discovery. Front Pharmacol 2025, 16:1597351, DOI 10.3389/fphar.2025.1597351.
