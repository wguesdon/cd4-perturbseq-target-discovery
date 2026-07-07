# Source paper: future work and limitations

Notes on the paper our project builds on, compiled to check our novelty claims against what
the authors themselves say is next. Quoted excerpts are reproduced under the paper's CC BY
license with attribution.

## Citation

Zhu R., Dann E., Yan J., Reyes Retana J., Goto R., Guitche R. C., Petersen L. K., Ota M.,
Pritchard J. K., Marson A. (2025). Genome-scale perturb-seq in primary human CD4+ T cells
maps context-specific regulators of T cell programs and human immune traits. bioRxiv.
DOI: 10.64898/2025.12.23.696273. Posted 2025-12-24. License: CC BY.

- Landing page: https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1
- Data and code, MIT plus CZI Acceptable Use Policy: https://github.com/emdann/GWT_perturbseq_analysis_2025
- Dataset portal: https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq

Note: the seventh author reads "Lillian K. Petersen" on bioRxiv and "Lillian Brixi" on one
mirror. Confirm before formal citation.

## Abstract (verbatim)

> "Gene regulatory networks encode the fundamental logic of cellular functions, but
> systematic network mapping remains challenging, especially in cell states relevant to
> human biology and disease. Here, we perturbed all expressed genes across 22 million
> primary human CD4+ T cells from four donors and developed a probe-based perturb-seq
> platform to measure the transcriptome effects in cells at rest and after stimulation.
> These data allow us to map genes that regulate known and novel pathways, including novel
> regulators of cytokine production. Importantly, active regulators and the gene programs
> they control change dramatically across stimulation conditions. Perturbation signatures
> enabled us to model T cell states observed in population-scale transcriptomic atlases,
> nominating regulators of Th1 and Th2 polarization and of age-related T cell phenotypes.
> Finally, we leveraged perturb-seq to implicate context-specific gene regulatory pathways
> in autoimmune disease risk. Our data provide a foundational resource to decode human
> immune function and genetic variation and for new approaches to study gene regulatory
> networks."

## What the paper already does

- Maps context-specific regulators, including novel regulators of cytokine production.
- Shows active regulators and their gene programs change across stimulation conditions.
- Models atlas T cell states from perturbation signatures. Nominates Th1 and Th2
  polarization regulators and age-related phenotypes.
- Implicates context-specific pathways in autoimmune disease risk, using GWAS, OpenTargets,
  the OneK1K single-cell eQTL cohort, and UK Biobank loss-of-function burden.

## Authors' stated future work (Discussion, verbatim)

> "Over the past two decades, human population genetics have associated genetic variation
> to human traits. More recently, single-cell genomic technologies have mapped diverse
> cell states in human health and disease. Now, perturb-seq in human primary cells has
> potential to map how genetic variation controls cell states, offering new hope that we
> can systematically link genome sequences to cell programs to human health outcomes."

> "As advanced perturbation and measurement technologies converge, we envision this
> framework extending to new perturbation modalities and cellular contexts, scaling from
> tens of millions to hundreds of millions, or even billions, of cells."

> "We make the data publicly available to facilitate these efforts and accelerate
> community-driven discovery."

## Authors' stated limitations (Discussion, verbatim)

> "Despite performing large-scale perturb-seq with 2 gRNAs per target across multiple human
> blood donors, future work would still be required to fully assess potential off-target
> guide effects with additional guides per target gene and to differentiate donor-specific
> biology from experimental batch effects."

> "This enhances statistical robustness, but potentially masks heterogeneity of response to
> perturbation across cells, which could be resolved by future distribution-aware
> analyses."

> "Our study focuses on a non-polarized culture condition; consequently, the regulatory
> rewiring driven by polarizing cytokines – central to CD4+ T cell plasticity – remains
> unmapped."

## Future work, distilled

1. More guides per target. Separate off-target effects and donor biology from batch.
2. Single-cell distribution-aware analyses beyond pseudobulk.
3. Extend to polarizing cytokine conditions and other contexts and modalities.
4. Scale cell numbers by orders of magnitude.
5. Tighter integration of perturb-seq with population genetics and cell atlases.

## Data availability, distilled

- Public no-auth S3 bucket: `s3://genome-scale-tcell-perturb-seq/marson2025_data/`.
- Effect matrix: `GWCD4i.DE_stats.h5ad`, ~16.8 GB, 33,983 perturbation-by-condition rows,
  10,282 genes, layers `log_fc`, `zscore`, `p_value`, `adj_p_value`, `baseMean`, `lfcSE`.
- Obs-level summary CSV: `DE_stats.suppl_table.csv`, ~4.8 MB.
- Guide and donor resolved: `GWCD4i.DE_stats.by_guide.h5mu`, `GWCD4i.DE_stats.by_donors.h5mu`.
- Design: 4 donors; Rest / Stim8hr / Stim48hr; genome-wide CRISPRi; single-gene; ~2 guides.
- Raw reads: SRA SRP643211 / GEO GSE314342, to be released.

## Access note

bioRxiv and SSRN block automated fetching with HTTP 403. The quotes above were retrieved
via the bioRxiv API, Europe PMC (record PPR1230028), and a reader render of the PDF. The
exact "Data Availability" statement could not be pulled verbatim; the details above come
from the authors' `data_sharing_readme.md` and the CZI page, which are authoritative for
the released artifacts.
