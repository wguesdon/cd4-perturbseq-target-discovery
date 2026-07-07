# CD4+ T cell Perturb-seq Target Discovery

Research Track submission for the **2026 Built with Claude: Life Sciences** hackathon.

> All work in this repository was started from scratch on 2026-07-07 and carried out
> during the event. The starting question and the public dataset predate the hackathon.
> The analysis does not.

## The question

Which single gene knockouts reshape the state of primary human CD4+ T cells in ways
that point to new immune drug targets? This project mines a genome scale CRISPR
Perturb-seq screen to rank perturbations by how strongly they push T cells toward, or
away from, activation and differentiation programs. The output is a prioritized,
reproducible target list with the evidence behind each call.

## Data

Primary dataset: genome scale CD4+ T cell Perturb-seq from the Marson lab, hosted on
CZI Virtual Cell Models.

- Dataset: https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
- Reference analysis code: https://github.com/emdann/GWT_perturbseq_analysis_2025
- Preprint: https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1

Raw and processed data files are **not** committed. See [`data/README.md`](data/README.md)
for the download steps and the expected layout under `data/`.

## Approach

The full plan lives in [`docs/research_plan.md`](docs/research_plan.md). It is a working
draft, not yet finalized. In short:

1. Use the released per-condition perturbation effect matrices. Do not reprocess raw cells.
2. Score 6 to 8 canonical T cell programs for each perturbation and condition.
3. Compute a therapeutic-window score. Reward context-selective suppression of the
   inflammatory program. Penalise loss of homeostasis and viability.
4. Benchmark the score by recovering targets of approved immunomodulators. This is the
   headline validation, reported as precision and recall.
5. Layer druggability and autoimmune genetics onto the top hits.
6. Deliver a ranked target table and 5 to 10 sourced evidence cards others can reproduce.

## Reproducibility

The project is managed with [uv](https://docs.astral.sh/uv/). Python 3.13.

```bash
git clone git@github.com:wguesdon/cd4-perturbseq-target-discovery.git
cd cd4-perturbseq-target-discovery
uv sync
# then follow data/README.md to fetch the dataset
```

Analysis dependencies are added on demand with `uv add`. See the research plan for the
intended stack.

## Repository layout

```
data/        raw, interim, and processed data (contents gitignored)
notebooks/   exploratory and analysis notebooks
src/         reusable package code (cd4_perturbseq)
results/     figures and tables produced by the analysis
docs/        research plan and method notes
reports/     final write up and submission materials
```

## How Claude Science was used

Documented as the work proceeds in [`docs/research_plan.md`](docs/research_plan.md) and
the final write up under `reports/`.

## License

MIT. See [`LICENSE`](LICENSE).
