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

Fill in the exact command once the download source is confirmed from the portal.

```bash
# Example placeholder — replace with the confirmed download step:
# curl -L "<dataset-url>" -o data/raw/cd4_perturbseq.h5ad
```

Record the file name, size, and access date here after downloading, so the provenance
is clear.

| file | size | source | date accessed |
| ---- | ---- | ------ | ------------- |
|      |      |        |               |
