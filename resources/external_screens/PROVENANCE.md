# External T-cell CRISPR screens (independent replication data)

Vendored for N19 cross-screen concordance. Both files are the MAGeCK gene-level outputs bundled in the
source paper's analysis repository `emdann/GWT_perturbseq_analysis_2025@master` (MIT licence), redistributed
here with attribution for reproducibility.

- **Freimer2022_Screen.csv** — Freimer et al. 2022 (Marson lab). CRISPR screen for regulators of IL2RA
  (CD25) surface expression in primary human CD4 T cells. MAGeCK `neg|*` = knockdown depletes the selected
  (IL2RA-high) population (a positive regulator of IL2RA); `pos|*` = enriches (a negative regulator).
- **Arce2025_Screen.csv** — Arce M.M., Umhoefer J.M. et al. (2025) "Central control of dynamic gene circuits
  governs T cell rest and activation", *Nature* 637(8047):930–939, doi:10.1038/s41586-024-08314-y.
  Citation confirmed 2026-07-09 against the PMC full text (PMC11754113).
  **It is NOT a fitness or proliferation screen.** An earlier version of this file said so, and that was
  wrong. It is a **FACS marker screen**: cells were gated on the sgRNA-library marker and the top and
  bottom 20% of **IL-2Rα (CD25)** expressing cells were sorted. Guides are ranked by their shift in
  IL-2Rα, not by dropout or abundance. Three arms: resting Teff (IL-2Rα-low), restimulated Teff, resting
  Treg (IL-2Rα-high). MAGeCK `neg|*` therefore flags genes whose loss **lowers IL-2Rα**, i.e. positive
  regulators of the marker. It does not mean "depleted from a proliferating population"; no such
  population exists in this design.

> **★ The two vendored screens share one library.** Freimer 2022 and Arce 2025 target the **identical
> 1,351 genes** (Jaccard = 1.0000, verified). Both are Marson-lab panels of transcription factors and
> chromatin modifiers, approximately 1,350 trans-factor genes, **not genome-wide**. They are independent
> in *readout* and *cell state*, and not at all independent in *gene space*.
>
> Two consequences that must accompany any use of them:
> 1. Describing them as "two independent screens" is accurate only about the assay. Any statement about
>    coverage, overlap, or replication rate refers to the same 1,351 genes twice.
> 2. The 472 genes co-tested against our genome-scale screen are a **transcription-factor and
>    chromatin-modifier panel**, not a random sample of the perturbed genome. Any candidate that emerges
>    from a comparison restricted to this set is enriched for transcription machinery by construction.

**Held out, NOT vendored:** `SchmidtSteinhart2022_CRISPRi_screen_gene_phenotypes.csv` and
`Schmidt2022_hits_Supplementary_table_2.xlsx` remain the validation set (RULE #3) and must never enter a
discovery score. They are read only where the pipeline already uses Schmidt for held-out validation.

These screens measure T-cell *function/fitness*, not autoimmune therapeutic direction. Concordance here is
evidence a gene is a genuine, replicated functional regulator (not a screen artifact); it does not resolve
direction of effect (N17 showed that remains unresolved).
