# Three-minute demo video: shot list and narration

Target 3:00. Read at roughly 145 words per minute. The narration below is 420 words, about 2 min 54 s
at that pace, which leaves room to breathe on the two visual beats.

Every number spoken is in `results/tables/` and in the paper. **Nothing in this script may be spoken
that is not on screen.**

Register: a seminar, not a pitch. No adjectives about our own work.

---

## 0:00 – 0:20 · The hook

**On screen:** the paper's title block, then `scripts/29_nomination_recalibration.py` running in a
terminal and exiting with a red non-zero status.

> I asked a genome-scale Perturb-seq atlas of human CD4 T cells for drug targets. The interesting result
> is not a target. It is that the pipeline learned when not to nominate one, and it can show you why.

---

## 0:20 – 0:45 · What the data is, and what it is not

**On screen:** Figure 1A, the top twenty perturbations by efficacy, coloured by triage outcome.

> Genome-scale CRISPRi, primary human CD4 T cells, four donors. A probe panel of ten thousand genes,
> not the transcriptome. Messenger RNA only. Rank perturbations by how well they suppress the activation
> program and you get the T cell receptor signalosome: `VAV1`, `CD3G`, `CD247`. An independent screen
> confirms these knockdowns really do reduce interleukin-2. They work. That is the problem.

---

## 0:45 – 1:15 · The layer, stated at its true strength

**On screen:** Figure 1B, observed 14 against a magnitude-matched 4.6.

> The layer refuses fourteen of those twenty, against four point six for a shortlist matched on
> effect magnitude. The rejection is specific to one axis: loss of an activation-induced co-inhibitory
> transcript module. Not tolerance. Not Treg identity. Not safety. `FOXP3` is transiently induced in
> activated conventional T cells without conferring suppressive function. We measured messenger RNA.

---

## 1:15 – 1:50 · The audit, on camera

**On screen:** the p-value genealogy table scrolling; then `results/tables/pvalue_genealogy.csv`.

> Then it was audited. A per-perturbation p-value of three times ten to the minus thirteen treated six
> thousand values as independent, when the module is one fixed nine-gene set. Its own negative control
> fired at a fifteen percent false-positive rate. Retracted. A benchmark first said one hundred and
> seventy-three approved targets. The honest count is thirty-eight. Retracted. Peer review found two
> different efficacy statistics in one pipeline, and an evidence floor that admitted perturbations which
> *induce* the program it was supposed to suppress. Both corrected. This table tracks nineteen statistics.
> Five are retracted. Every one carries the reason it changed.

---

## 1:50 – 2:20 · What survived

**On screen:** Figure 4, the GSEA sign inversion; then Figure 5, the Freimer comparison.

> What survives is narrow and it replicates. Against a control matched on effect magnitude, the efficacy
> axis inverts interleukin-2 and TNF signatures while suppressing proliferation identically. And on an
> independent, protein-level CD4 screen, the efficacy axis is concordant under four stratified nulls.
>
> The co-inhibitory axis is not. Zero of four. The commodity half of the pipeline is the half that
> replicates. The half carrying the headline is the half that does not.

---

## 2:20 – 2:45 · The temptation, and the refusal

**On screen:** `ATXN7L3` and `XBP1` appearing as candidates, then red flags firing beside them.

> So I relaxed the rule, to see whether an independent screen could rescue a named target. It returned
> two. Both failed criteria this pipeline had already written down. One disturbs the unstimulated cell
> more than the stimulated one. And the comparator library contains only transcription factors and
> chromatin modifiers, so there was nothing else it could have returned.

---

## 2:45 – 3:00 · The result

**On screen:** the paper's abstract, the sentence about abstention highlighted.

> No vetted novel target. Direction of effect can veto a candidate; it cannot promote one, because every
> annotation that carries direction requires a drug to already exist. That is not a gap a better ranking
> closes. Abstention is the finding.

---

## Recording notes

- **Do not** show a candidate gene without its liabilities on the same frame.
- The terminal beat at 0:00 is real: `scripts/29` exits non-zero because a pre-registered control fails.
  Do not fake it, and do not cut away before the exit status is visible.
- Figures are in the rendered report; screenshot them from `reports/report.html` so what is on screen is
  what a reader would see.
- If a number in the paper changes before recording, this script is wrong. Re-read it against
  `results/tables/pvalue_genealogy.csv` first.
