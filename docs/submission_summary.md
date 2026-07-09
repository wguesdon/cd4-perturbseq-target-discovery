# Written summary

Required length: 100 to 200 words. The text below is **181 words**. Every number appears in
`results/tables/` and in the paper.

---

Genome-scale Perturb-seq in primary human CD4 T cells is powerful, but not a drug-target engine by
default. I built an audit-ready triage and abstention layer for CRISPRi perturbations that suppress a
polyclonal activation program, while penalising loss of an activation-induced co-inhibitory messenger
RNA module. Fourteen of the twenty highest-ranked perturbations are refused, against 4.6 for a
shortlist matched on effect magnitude. The layer recovers approved immunomodulator targets above chance
at the evidence floor, and gene-set enrichment shows the efficacy axis inverts IL-2/STAT5 and TNF
signatures where an equal-magnitude control does not. The same controls prevented overclaiming.
Essentiality, immune-disease annotation, genetics, colocalisation, mouse phenotype and cross-screen
evidence each failed calibration, exposed assay blindness, or proved characterisation-bound. On an
independent protein-level CD4 screen the efficacy axis is concordant; the co-inhibitory axis is not.
Several headlines were retracted after negative controls, seed-sensitivity tests, self-tests and peer
review exposed errors. The result is not a vetted novel target. It is a reproducible decision layer
that recovers known pharmacology, treats missing direction as missing rather than favourable, and
refuses to nominate when the evidence is uncalibrated.

---

## Word count check

```bash
sed -n '/^Genome-scale/,/uncalibrated\.$/p' docs/submission_summary.md | wc -w
```
