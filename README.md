# A screen-native triage layer for a CD4+ T cell Perturb-seq atlas

An exploratory computational reanalysis of the genome-scale CRISPRi Perturb-seq atlas of primary human
CD4+ T cells released by [Zhu, Dann et al. 2025](https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1)
(preprint). We ask whether a target-triage layer built only from quantities the screen measures can
nominate a drug target that is not already known.

**It cannot, and the reasons are the result.**

> **Read the paper:** [`reports/report.html`](reports/report.html) (self-contained), source at
> [`reports/report.qmd`](reports/report.qmd).

All work in this repository was started on 2026-07-07 and carried out during the
[Built with Claude: Life Sciences](docs/hackathon_brief.md) event. The question and the dataset predate
it; the analysis does not.

---

## What this is, and what it is not

It **is** an audit-ready decision layer that recovers known pharmacology above chance, exposes what the
assay cannot see, and declines to nominate when the evidence is uncalibrated.

It is **not** a validated method for target discovery, safety assessment, or therapeutic-window
estimation. The central co-inhibitory criterion is a **messenger RNA annotation**. Nothing here was
validated experimentally.

## The result

Of 12,779 library genes, 6,371 perturbations survive **our own** quality-control mask, not the source
authors'. Of these, 263 clear an evidence floor requiring both a minimum number of significantly
decreased effector-module genes and net suppression, and 197 additionally fall below the 75th percentile
of co-inhibitory suppression.

| finding | value |
| --- | --- |
| Top 20 by efficacy, rejected by the layer | **14 of 20**, against 4.6 for a magnitude-matched shortlist (P = 2.0e-4, the 1/5001 floor) |
| Context-selectivity axis | 14 versus 14.2 matched, P = 0.65. No information beyond effect magnitude. Demoted to an annotation. |
| Co-inhibitory module versus 200 matched random modules | exceeds every one, so the permutation P equals its floor of 1/201. The co-regulation confound is open. |
| Recovery of approved immunomodulator targets | 4 of 20 eligible positives at the evidence floor (p = 0.0081) |
| Do the safety axes add recovery? | **No.** 3 observed, 3.00 expected, p = 0.74 |
| Rank discrimination on drug targets | AUROC 0.497 [0.339, 0.656]. Chance. |
| Held-out validation of the efficacy axis | AUROC 0.701 [0.583, 0.803] against Schmidt 2022 IL-2 hits |
| Independent CD4 screen (Freimer 2022) | efficacy axis concordant (Spearman +0.115, 4 of 4 stratified nulls); **co-inhibitory axis not** (+0.025, 0 of 4) |
| Direction of effect | annotated for **17 of 197** genes; every concordant annotation requires a pre-existing approved drug |
| Novel targets nominated | **none** |

The efficacy axis, the commodity half of the pipeline, is the half that shows concordance on an
orthogonal protein-level screen. The co-inhibitory axis, which carries the headline, is the half that
does not.

## What the screen cannot see

- The stimulation carries no polarising cytokines, so cytokine-signalling targets rank near the bottom,
  and two targets of an approved JAK inhibitor were never perturbed.
- Single-gene CRISPRi is blind to redundant targets. The two catalytic calcineurin paralogues yield
  fewer than four differentially expressed genes each; the non-redundant regulatory subunit yields 523.
- Essentiality cannot be evaluated as a safety axis in a screen where essentiality governs entry into
  the analysable set. A comparison among the survivors of that selection is not causal in either
  direction.
- A drug-recovery benchmark cannot be powered here. Against a threshold of 60 fixed before counting,
  the primary count is 38 and the most permissive upper bound is 53.

## Claims withdrawn during this work

Recorded in full in the paper's supplement and in `results/tables/pvalue_genealogy.csv`, which is
validated at build time by `scripts/30_pvalue_genealogy.py`.

- An **inborn-errors-of-immunity enrichment** among top-ranked perturbations was used to argue that the
  ranking is unsafe. The same flag is *more* enriched among approved drug targets than among the
  top-ranked perturbations. It is an annotation, never a gate.
- A per-perturbation p-value of 3.5e-13 treated 6,371 z statistics as independent when the module is a
  single fixed nine-gene set. Its false-positive rate on random modules was 15% against a nominal 5%.
- A drug-recovery benchmark count of 173 "powered" targets was inflated by three annotation artifacts.
  The honest count is 38.
- Two efficacy statistics were in use, a z score for ranking and a mean log fold change for everything
  else. They share 9 of their top 20. Unified to one signed statistic.
- The evidence floor admitted perturbations that *induce* the effector module, one of them a curated
  positive. A direction requirement was added.
- `P < 0.0001` was reported from 5,000 draws without the `+1` correction.
- `BINDING AGENT` and `CROSS-LINKING AGENT` were treated as loss-of-function-mimicking drug mechanisms.
  An anti-CD3 cross-linking antibody is the agonist that stimulates these very cells.
- A rule inferring therapeutic concordance from an inborn error of immunity called `PTEN` concordant,
  although PTEN loss causes lymphoid hyperplasia and autoimmunity.
- Two areas under the ROC curve were computed on a third score, distinct from the primary statistic.

## Reproducing

```bash
uv sync                                   # Python 3.13
bash scripts/fetch_priors.sh              # prior gene lists and reference screens
bash scripts/06_fetch_safety_priors.sh    # gnomAD, Human Protein Atlas
bash scripts/15_fetch_open_targets.sh     # Open Targets, release pinned to 26.06

# the effect matrix (16.8 GB), not vendored
aws s3 cp --no-sign-request \
  s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad data/raw/

uv run python scripts/01_build_activation_program.py
uv run python scripts/04_window_score.py
# see the Methods section of reports/report.qmd for the full order

QUARTO_PYTHON="$PWD/.venv/bin/python" quarto render reports/report.qmd --to html
```

Several scripts are written to be able to fail. `scripts/02_risk_kill_reversal.py` exits non-zero if its
pre-registered endpoint breaks. `scripts/29_nomination_recalibration.py` **currently exits non-zero**,
because a pre-registered control on the rebuilt nomination rule fails for want of power; the rule is void
by its own criterion, and the paper says so. `scripts/30` and `scripts/31` fail if any quoted statistic
or denominator has drifted from its source table, and the paper's abstract fails the render on drift.

## Layout

| path | contents |
| --- | --- |
| `reports/report.qmd` | the manuscript; every number read from a committed table at render time |
| `scripts/` | numbered pipeline, one concern per script |
| `src/cd4_perturbseq/` | shared readers for the effect matrix and the external priors |
| `results/tables/` | every committed intermediate; the paper reads these, never the raw data |
| `docs/preregistration_*.md` | decision rules fixed before the analyses that used them |
| `docs/results/*.md` | per-analysis write-ups, including the negative results |
| `docs/citation_verification_2026_07_09.md` | every reference fetched and checked before it was cited |
| `docs/claude_tooling_log.md` | how the analysis was built, and where the tooling disagreed with itself |

## Limitations

Donor and guide uncertainty is not propagated: results pool four donors and generally two guides per
target, and `PPP3R1`, used here as a recovered positive, has two guides whose signatures anticorrelate.
The co-inhibitory cut is a quantile, so it rejects a quarter of the evidence-passing set by construction
and is calibrated against no external phenotype. The enrichment analysis is partly circular. The whole
analysis is exploratory rather than prospectively specified. The paper's Limitations section states these
before any result that rests on them.

## Licence

MIT. The source atlas is CC BY 4.0. Vendored reference screens carry their original attribution in
`resources/external_screens/PROVENANCE.md`.
