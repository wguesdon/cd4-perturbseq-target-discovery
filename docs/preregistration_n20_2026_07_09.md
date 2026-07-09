# Pre-registration — N20: is "zero vetted novel targets" a property of the biology or of the rule?

Written 2026-07-09, before `scripts/29_nomination_recalibration.py` exists. Committed before the
rebuilt nomination layer is written. Both branches are acceptable and are stated here in advance.

## 0. Provenance, stated first

This is a **post-hoc re-analysis.** It was prompted by a challenge from the user, not by a planned
step: *given the size of the dataset, if we do not identify a single new potential target, something
must be wrong in our method.*

Three measurements were made in scratch scripts **before** this document was written, and are recorded
here rather than presented later as if they had been predicted. They are the motivation, not the
result:

1. Applying the N16 SHORTLIST criteria to the six approved-drug targets the pipeline recovers, only
   **3 of 6** would have been nominated.
2. **21 of 21** CONCORDANT direction verdicts rest on a basis that requires the gene to be already
   characterised: 15 read `approved LoF-mimicking drug (OT actionType)`, 6 read
   `LoF causes immunodeficiency (IUIS)`.
3. `src/cd4_perturbseq/priors.py:344` defines `tractable = sm_pocket or sm_family or ab`, omitting the
   `clinical_precedent` field computed on the line above.

The point of pre-registering now is that the **decision rules below were fixed before the rebuilt
nomination layer was run**, and both outcomes are publishable. RULE #8 governs every sentence that
leaves this analysis.

## 1. The three hypotheses

**H1 — the rule fails its own positive class.**
The N16 SHORTLIST rule is `supported & lof_tolerant & tractable & ~is_known_drug & ~discordant`.
Define its *sensitivity* as the fraction of the recovered approved-drug targets (n = 6: IMPDH2, CD2,
CD3E, CD28, PPP3R1, IL4R) that satisfy `supported & lof_tolerant & tractable`, i.e. that would have
been nominated had we not already known they were drugs.

> **Decision rule, fixed in advance.** If sensitivity < 5/6, the rule is rejected as a *nomination
> gate*, and the "zero vetted novel targets" conclusion is reported as an artifact of the rule rather
> than as a statement about the screen. If sensitivity ≥ 5/6, the rule stands and the zero result
> stands with it.

**H2 — direction-adjudicability is a collider on prior characterisation.**
Let `p` be the proportion of CONCORDANT verdicts whose `direction_basis` requires either an existing
drug or a described monogenic syndrome.

> **Decision rule, fixed in advance.** If `p > 0.80`, then conditioning a nomination on direction
> concordance selects for genes that are already drug targets, "novel AND concordant" is close to a
> contradiction in terms, and direction of effect may be used **only to demote** (rule out discordant
> genes) and **never to promote**. If `p ≤ 0.80`, direction may remain a promoting criterion.

**H3 — `tractable` omits the strongest tractability evidence.**
This is a code defect, not a hypothesis. It is fixed, and every gene whose flag changes is named.
The fix is independent of which gene it promotes: `clinical_precedent` is computed and discarded,
which cannot have been intended.

## 2. The rebuilt nomination gate, and how it is chosen

The replacement rule is selected on **positive-class sensitivity alone.** That criterion never looks at
which novel gene a rule promotes, so it cannot be tuned to manufacture a hit. Among rules formed as an
AND over the available conjuncts, we take the one with the most conjuncts whose sensitivity is 6/6.

> **NOMINATE = passes the screen gate (evidence floor + tolerance axis) AND `tractable` AND NOT
> direction-discordant.**

Genetic support and LoF tolerance are **demoted from conjuncts to reported annotations**, for reasons
that are internal to this dataset and do not depend on the outcome:

- **Genetic support.** IMPDH2, CD3E and PPP3R1 each have an Open Targets autoimmune genetic score of
  exactly 0.000. Requiring genetic support would reject mycophenolate, anti-CD3 and ciclosporin. This
  agrees with N13, which found the safe ranking is not enriched for genetic support (p = 0.67). The
  literature treats genetic support as a prior on approval odds, roughly a doubling (Nelson 2015;
  King 2019; Minikel 2024), not as a necessary condition.
- **LoF tolerance.** IMPDH2 has `prec` = 0.99883, so it is recessive-intolerant and fails
  `lof_tolerant`. Mycophenolate is nonetheless approved against it. An approved drug against a
  recessive-intolerant gene demonstrates that LoF intolerance cannot gate a nomination. It constrains
  the **modality** — an inhibitor with a therapeutic window rather than a knockout or a degrader — and
  is reported as such.

Output is a **ranked, tiered candidate table with per-gene liabilities**, not a binary verdict. The
ranking is transparent and lexicographic, never a fitted score (RULE #2).

## 3. Controls. All must pass, or the analysis is void.

- **Positive control (recovery is preserved).** The rebuilt gate must still recover the curated
  approved-drug targets above chance, tested against the screened-gene background with the N14
  permutation machinery. If loosening the rule destroys the recovery signal, the loosening is wrong.
- **Negative control (discordance still demotes).** `PTPN2` and `RC3H1` must not appear anywhere in the
  nominated set. If either does, the direction demotion has been broken.
- **Negative control (the flag carries information).** Permuting the `tractable` column across genes
  must destroy the recovery enrichment. If a shuffled flag performs as well, the flag is decorative.
- **Ground-truth guard (RULE #4).** No member of the 36-gene curated positive set may be presented as a
  novel nomination. Recovered drugs are validation, and are labelled so.
- **Housekeeping guard, pre-committed.** The gate alone promotes proliferation and housekeeping genes;
  N16's peer critic already named `PRKAR1A`, `KRR1` and `ENO1` at the top window ranks. We will **count
  and report** how many nominations are core housekeeping or mitochondrial-complex genes. This number is
  published whatever it is. If the top of the ranking is housekeeping, the ranking is stated to be unfit
  for purpose.

## 4. What will be claimed, and what will not

**Will be claimed, if the controls pass.** That "zero vetted novel targets" was an artifact of a
nomination rule with 50% sensitivity on its own positive class, and of a collider that made direction
concordance nearly equivalent to being an existing drug target. That a ranked candidate table with
explicit liabilities is the correct deliverable for a triage layer.

**Will not be claimed.** That any gene here is a validated novel drug target. That a candidate is safe.
That direction of effect has been resolved for the 191 direction-UNKNOWN genes; it has not, and N17
showed a lightweight resolution fails its own calibration at chance.

**Stated wherever the leads are named.** The strongest candidates are concordant *because inhibitors
against them already exist*. They are therefore **repurposing candidates, not virgin targets.** Genes
with no prior precedent remain direction-unresolved. That limitation precedes the result, not follows it.
