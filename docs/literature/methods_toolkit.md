# Methods toolkit: packages, APIs, and reference data to reuse

Actionable companion to [`lit_review_2026_07_07.md`](./lit_review_2026_07_07.md). Everything
here is established and reusable, so we compose rather than reinvent. Python-first, since the
project uses uv and Python 3.13.

## Python packages (add with `uv add`)

| Package | Use | Notes |
| --- | --- | --- |
| `scanpy`, `anndata`, `mudata` | Backbone for AnnData / MuData objects | The released files are `.h5ad` and `.h5mu`. |
| `pertpy` | E-distance + E-test, Mixscape, Augur, DE | The single most reusable perturbation package. |
| `decoupler` | Score the DE effect matrices against gene programs | Runs directly on a genes-by-contrasts statistic matrix. Highest leverage for us. |
| `pydeseq2` | Pseudobulk differential expression | Only if we recompute; the released matrices already have DE. |
| `ucell` / `pyucell` or `sc.tl.score_genes` | Per-cell module scores | UCell scales better than AUCell. |
| `gseapy` | Enrichment (GO, pathways) on downstream genes | For the mechanism step. |
| `scikit-learn`, `statsmodels` | Benchmark metrics, calibration, sensitivity | AUROC, precision/recall, enrichment ORs. |

Optional / if time allows: `sceptre` (R) for calibrated CRISPR-screen FDR; `GSFA` (R) for
de-novo module inference; `miloR` (R) for compositional shifts; `fpocket` (CLI) for
ligandability on AlphaFold models.

## APIs to query (per candidate gene)

- Open Targets Platform GraphQL: `https://api.platform.opentargets.org/api/v4/graphql`.
  One call returns `tractability` for small-molecule, antibody, and PROTAC modalities.
- Pharos / TCRD GraphQL: `https://pharos-api.ncats.io/graphql`. Returns Target Development
  Level. Tclin is a ready-made positive set for the benchmark.
- ChEMBL: `chembl_webresource_client` Python package. Existing bioactives, potency,
  mechanism, `max_phase` (max_phase=4 is an approved-target positive set).
- DGIdb GraphQL: `https://dgidb.org/api`. Drug-gene interactions and druggable categories.
- Open Targets Genetics L2G, GWAS Catalog, FinnGen, eQTL Catalogue, OneK1K: genetic support
  and CD4-specific eQTL colocalization.

## Reference gene lists and data to download

- Finan et al. 2017 druggable genome (tiered gene list) for a hard druggability filter.
- Hart CEG2 core-essential and nonessential gene lists (hart-lab.org) for the viability floor
  and hit-calling calibration.
- DepMap Chronos/CERES gene-effect and selective-dependency scores (depmap.org) for the
  common-essential penalty and the "selective dependency" statistical shape.
- GTEx and Human Protein Atlas expression breadth for the off-tissue safety prior.
- PROTAC-DB 3.0 and MolGlueDB for existing-degrader flags and E3 ligand precedent.
- The source data on public no-auth S3: `s3://genome-scale-tcell-perturb-seq/marson2025_data/`
  (see [`../../data/README.md`](../../data/README.md)). Main input: `GWCD4i.DE_stats.h5ad`.

## Benchmark positive sets (pick one, report all)

- Pharos Tclin (approved-drug mechanism-of-action targets).
- ChEMBL `max_phase = 4`.
- Open Targets known-drugs; DrugBank approved targets.
- Curated immunomodulator target list in the lit review, section 5.

Construct lenient/moderate gold-standard sets in the style of Fang 2019 and Minikel 2024, and
report AUROC, top-k enrichment odds ratio, and precision/recall, so reviewers recognize the
framework.

## Do not put on the critical path

Foundation / perturbation models (scGPT, scFoundation, State, GEARS, CPA) mostly do not beat
simple additive/linear/mean baselines for perturb-seq. If a model must appear, use Geneformer
zero-shot in-silico deletion as a bounded, baseline-gated side arm, and accept it only if it
beats the additive/linear/DE baseline on target ranking on this dataset.
