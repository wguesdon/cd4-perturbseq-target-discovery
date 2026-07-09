# N16 — The honest final shortlist: what new drug targets did we find?

Run under RULE #9 (literature → implement → critical review → report), as a synthesis of every committed
axis rather than a new measurement. Script `scripts/25_final_shortlist.py`; table
`results/tables/final_shortlist.csv`. The critical review used a 7-agent adversarial workflow
(`wf_9807ec74-7b2`): four independent lenses on the lead candidate, a triage of the wider set, and a
hostile peer-critic. Every number below was verified against a committed table before it was written
(`scripts/_verify_icam2_claims.py`).

## 1. The question and the credibility bar (literature)

The competition's literal question is: given the screen, the safety gate, direction of effect, human
genetics, LoF tolerance, tractability, and the recovery audit, **what new drug targets did we find?**

The target-discovery literature sets the bar, and it is not a single screen hit. Human genetic support
roughly doubles the odds a target reaches approval (Nelson 2015; King 2019; Minikel 2024). For a
knockdown-nominated gene the therapeutic direction must be loss-of-function, because CRISPRi knockdown
mimics an inhibitor: a positive regulator of activation is concordant, a negative regulator (a brake) is
the wrong way round (Manguso 2017, the PTPN2 lesson). The gene must tolerate loss of function in humans
(gnomAD LOEUF / pRec; Minikel 2020) and be tractable. A candidate that clears the screen but fails any of
these is not a credible novel nomination. We encode this bar as explicit, auditable buckets, never a
fitted score (RULE #2).

## 2. The buckets (implement)

The 214 gate-passing safe genes, bucketed by the pre-stated rules. Four honesty controls pass: no known
drug and no discordant gene reaches the shortlist; every recovered drug lands in validation; PTPN2 and
RC3H1 are demoted.

| bucket | n | definition |
|---|---|---|
| VALIDATION | 6 | recovered approved-drug targets (IMPDH2, CD2, CD3E, CD28, PPP3R1, IL4R) |
| DEMOTED | 2 | direction-DISCORDANT negative regulators (PTPN2, RC3H1) |
| SHORTLIST | 1 | genetically supported + LoF-tolerant + tractable + not a known drug + not discordant |
| WIDER | 44 | genetically supported, novel, concordant-or-unknown, but LoF-constrained or not tractable |
| GATE_ONLY | 161 | safe and gate-passing, no autoimmune genetic flag, direction-UNKNOWN |

**Direction of effect is adjudicable for only 23 of 214 safe genes** (21 concordant — most are recovered
drugs; 2 discordant). The remaining **191 are direction-UNKNOWN** and cannot be vetted in silico in this
budget (N11/N12: eQTL-colocalisation is 21 GB and out of budget; a GO fallback could not even flag PTPN2).

## 3. Recovery, stated with its denominator (not an unqualified "6")

The pipeline recovers approved-drug targets, but the honest statistic is the **pre-registered N14
permutation test on the assay-visible positive set (20 of 36 curated positives are rankable in this
screen)**, not the count of drug targets in the safe set:

- **Efficacy floor:** 5 of the 20 pass, vs 0.90 expected under a magnitude-matched null (5.6×,
  hypergeometric p = 0.0016, permutation p = 0.0017).
- **Full gate:** 4 of 20 pass, vs 0.67 expected (perm p = 0.0037).
- **The safety axis adds nothing beyond the efficacy floor:** 4 observed vs 3.74 expected, p = 0.63.
  Recovery is driven by the efficacy axis, out of sample (Schmidt held out, RULE #3).
- **Caveat, stated next to the number:** the drug-target-recovery AUROC on the 20-target set is
  0.542 [0.373, 0.707], which spans chance. Recall is capped by assay blindness — a single TCR-only
  stimulation with no polarising cytokines pushes cytokine-signalling targets (JAK2, TYK2, S1PR1, IL4R)
  to the bottom **by construction** — so 15 of 20 curated positives are assay-invisible here, not screen
  failures. The six approved-drug targets in the VALIDATION bucket include CD2 and CD28, which were held
  **out** of the tested positive set (RULE #4); "6" is an annotation, "5 of 20 above chance" is the test.

## 4. The single intersection candidate: ICAM2 (and why it is not a nomination)

ICAM2 is the one gene at the intersection of all filters. Its full, verified profile:

| axis | value | reading |
|---|---|---|
| screen efficacy | z 0.705, window rank 44/6371, 4 module genes down | genuine on-target suppression (skeptic lens: OK, not a magnitude artifact) |
| direction of effect | **UNRESOLVED** | the pivotal gap; see below |
| autoimmune genetics | 0.283 across 4 diseases | real but the **weakest** in the nominated set (recovered drugs are 0.68–0.92); direction-agnostic |
| LoF tolerance | LOEUF 1.167, pRec 0.861 | tolerant (consistent with ICAM1 redundancy) |
| tractability | tractable, but `sm_pocket=False`, `clinical_precedent=False` | no direct ICAM2 drug or clinical precedent |
| homeostasis | **`fail_homeostasis=True`**, viability tier **depleting-at-rest**, broadly expressed (max non-immune nTPM 96.5, 45 non-immune tissues) | fails the project's own thesis — the six recovered drugs are all `fail_homeostasis=False` |

**Direction of effect — UNRESOLVED, with active counter-evidence (both direction lenses agreed).** This is
not merely "no data". The concordant evidence (ICAM-2 as an LFA-1 costimulatory ligand; the approved
LFA-1/ICAM axis drugs efalizumab and lifitegrast) sits in the **APC and endothelial compartments** that a
CD4-T-cell-intrinsic knockdown does not model. The T-cell-intrinsic and whole-organism evidence points the
other way or to redundancy: ICAM-2 is the lower-affinity, constitutively expressed ligand, largely buffered
by ICAM-1; a documented mouse ICAM-2 knockout shows *prolonged* eosinophilic airway inflammation (loss
*worsens* inflammation, discordant-leaning); and ICAM-2 ligation has been linked to activation-induced
T-cell death (blockade could *raise* activity). So ICAM2 cannot be distinguished from a PTPN2-like brake,
whose inhibition would worsen disease.

**Mechanism-class safety precedent.** ICAM2 is novel only at the exact-gene level. The LFA-1/ICAM
adhesion-blockade class carries a serious history: efalizumab (anti-LFA-1) was withdrawn for progressive
multifocal leukoencephalopathy (PML), a class-level signal. An ICAM2 programme would inherit that caution.

**Conclusion for ICAM2:** it is direction-unresolved, antiproliferative at rest, broadly expressed, and
mechanism-class-flagged. It is a **hypothesis for per-gene review, not a vetted nomination.**

> Note on "exactly one". The count of one is a definitional artefact of requiring genetics AND
> tractability AND LoF-tolerance AND non-discordant direction at once. The 161 GATE_ONLY genes — including
> the top window ranks PRKAR1A (1), KRR1 (2), ENO1 (5), all housekeeping/proliferation genes — are
> excluded only for lacking a genetics flag and are equally direction-UNKNOWN. "One" measures a filter
> intersection, not a robust discovery.

## 5. The wider watchlist (44 genes), triaged honestly

Every gene here is genetically supported but misses strict LoF-tolerant AND tractable together, so each is
**inhibitor-only and safety-flagged** (an inhibitor with a therapeutic window, never a knockout). The
strongest concordant secondaries, with their honest caveats:

- **STAT3** (genetics 0.78 / 8 diseases, CONCORDANT, tractable). Human STAT3 gain-of-function causes
  early-onset autoimmunity — direct genetic evidence that lowering STAT3 activity is protective. But it is
  LoF-intolerant (biallelic/dominant-negative LoF → Hyper-IgE/Job syndrome); inhibitor only, narrow window.
- **HDAC7** (0.74 / 6, CONCORDANT). Best chemical handle (approved HDAC-inhibitor precedent, class-IIa
  selectivity plausible). LoF-intolerant; pan-HDAC liabilities demand isoform selectivity.
- **PTPRC / CD45** (0.83 / 7, the highest genetics in the whole set, CONCORDANT). But CD45 deficiency
  causes SCID, and it is a proximal TCR-signalling node expected to collapse the resting transcriptome —
  the same failure mode as the LCK/CD3/VAV1 contaminants. Genetically top, mechanistically the riskiest.
- **MALT1** (0.56 / 1, CONCORDANT). The clearest tractability exemplar (clinical-stage allosteric protease
  inhibitors). Thin genetics; chronic inhibition can reduce Tregs and drive IPEX-like autoimmunity.
- **TBX21 / T-bet** (0.62 / 5, CONCORDANT Th1). On-target for type-1 inflammation, but a transcription
  factor with limited tractability and a type-2-skew liability.

**Modality mismatch, flagged so they are not mis-sold as inhibitor targets.** EGR2 (0.63 / 10, the broadest
genetics in the set) and FOXO1 (0.51 / 7) are tolerance/quiescence promoters. Their loss is *discordant*
for an inhibitor — knockdown would worsen autoimmunity — so the modality would be an agonist/stabiliser.
Attractive genetics, wrong direction for our screen.

**Do not headline (RULE #6 and hazard):** STAT6, GATA3 (type-2 competitor turf), STAT5A (IL-2/Treg
adjacent, discordant, weak genetics), PTEN (tumour suppressor; inhibition is oncogenic).

## 6. Critical review — does it survive literature comparison?

- **The direction verdict for ICAM2 agrees with the adhesion-molecule literature**, which is genuinely
  two-sided: costimulatory on APCs, redundant/suppressive or dispensable T-cell-intrinsically. The honest
  reading is UNRESOLVED, and the compartment mismatch is why a T-cell-intrinsic screen cannot settle it.
- **The recovery result survives, but only as stated:** efficacy-axis enrichment is significant and out of
  sample; the safety axis adds nothing (p = 0.63); the AUROC spans chance. We do not claim the safety gate
  improves drug recovery.
- **"Correctly demotes" was too strong.** Direction was adjudicated for 23 of 214 genes via a hand-curated
  list that was wrong twice (STAT3, MLST8), caught only by the drug control. PTPN2 and RC3H1 are
  literature-defensible demotions, but the method is two curated judgements, not a validated pipeline.
- **"Audited decision layer" was too strong.** The adversarial audit's actual outcome was to *retract* two
  central safety axes (DepMap/Hart essentiality as survivorship-biased; IUIS/IEI as tautological, 25%
  false-rejection on approved drugs). The honest phrase is **adversarially stress-tested and substantially
  revised** — a decision layer whose failure modes are mapped, not one that passed audit.
- **Scope limits, stated plainly.** All evidence is in silico on one released pseudobulk DE matrix
  (single-gene, non-combinatorial, 4 donors, one TCR-only stimulation). CRISPRi models partial constitutive
  LoF, not reversible pharmacology. No wet-lab or orthogonal-cohort confirmation of any novel call. The
  earlier "gate opposes genetics" idea was walked back (N10: the safe-vs-naive delta CI straddles 0).

## 7. Verdict

**The pipeline surfaces zero *vetted* novel drug-target nominations from this dataset in one week.** That is
the honest answer to "what new targets did we find", and it is not a failure — it is the correct output of
a decision layer built to reject what it cannot support:

1. It **recovers** approved immunomodulator targets above chance on the efficacy axis (5 of 20 assay-visible,
   5.6×, permutation p = 0.0017), out of sample.
2. It **demotes** two genetically-supported genes whose therapeutic direction is discordant (PTPN2, RC3H1),
   a worked example of why genetics alone is not enough.
3. It yields **one filter-intersection candidate (ICAM2)** that is direction-unresolved, antiproliferative
   at rest, broadly expressed, and mechanism-class-flagged — a hypothesis for review, not a nomination — and
   a **44-gene, direction-UNKNOWN watchlist** whose strongest concordant members (STAT3, HDAC7, PTPRC,
   MALT1, TBX21) are LoF-intolerant, so inhibitor-only and safety-flagged.

**The deliverable is the adversarially stress-tested decision layer, not a discovery.** Its contribution is
that it says no — to reversal alone, to genetics alone, to a knockdown screen's fixed direction — and shows
its work each time. Direction of effect for the 191 UNKNOWN genes is the gating question a T-cell-intrinsic
knockdown screen cannot answer on its own, and resolving it (eQTL-coloc or per-gene curation) is the first
step before any of this watchlist becomes a nomination.

---

*Provenance note (RULE #9 §1 discipline).* The ICAM2 direction evidence and the wider-set triage were
assembled by an automated 7-agent adversarial literature review. Specific primary citations (mouse ICAM-2
knockout phenotypes, ICAM-2/AICD) are recorded in the workflow journal and require manual confirmation
before any external use; the textbook-level facts cited here (ICAM-2 as an LFA-1 ligand, ICAM-1 redundancy,
efalizumab's PML withdrawal) are well established. Numbers are from committed tables, verified by
`scripts/_verify_icam2_claims.py`.
