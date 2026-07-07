# Data

Data files are not committed. This folder documents where the data comes from and how
to fetch it, so the analysis stays reproducible.

## Primary dataset

Genome scale CD4+ T cell Perturb-seq from the Marson lab.

- Portal: https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
- Reference analysis code: https://github.com/emdann/GWT_perturbseq_analysis_2025
- Preprint: https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1

## Expected layout

```
data/
  raw/         # exactly as downloaded, never edited by hand
  interim/     # intermediate objects (filtered, normalized)
  processed/   # analysis-ready tables and per-perturbation results
```

## Download

The processed data sits on a public, no-auth S3 bucket:
`s3://genome-scale-tcell-perturb-seq/marson2025_data/`. Two objects cover our needs.

```bash
# Per-perturbation x per-condition effect matrix, ~16.8 GB. This is the main input.
aws s3 cp --no-sign-request \
  s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad \
  data/raw/GWCD4i.DE_stats.h5ad

# Lightweight obs-level summary, ~4.8 MB: per-perturbation metadata and DE counts.
aws s3 cp --no-sign-request \
  s3://genome-scale-tcell-perturb-seq/marson2025_data/suppl_tables/DE_stats.suppl_table.csv \
  data/raw/DE_stats.suppl_table.csv

# Browse the full bucket to confirm keys and find the guide/donor-resolved variants:
aws s3 ls --no-sign-request s3://genome-scale-tcell-perturb-seq/marson2025_data/
```

Guide- and donor-resolved variants (`GWCD4i.DE_stats.by_guide.h5mu`,
`GWCD4i.DE_stats.by_donors.h5mu`) live in the same bucket for the confidence checks. Raw
reads are under SRA SRP643211 / GEO GSE314342.

Record the file name, size, and access date below after downloading, so the provenance
is clear.

| file | size | source | date accessed |
| ---- | ---- | ------ | ------------- |
|      |      |        |               |
