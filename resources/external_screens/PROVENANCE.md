# External T-cell CRISPR screens (independent replication data)

Vendored for N19 cross-screen concordance. Both files are the MAGeCK gene-level outputs bundled in the
source paper's analysis repository `emdann/GWT_perturbseq_analysis_2025@master` (MIT licence), redistributed
here with attribution for reproducibility.

- **Freimer2022_Screen.csv** — Freimer et al. 2022 (Marson lab). CRISPR screen for regulators of IL2RA
  (CD25) surface expression in primary human CD4 T cells. MAGeCK `neg|*` = knockdown depletes the selected
  (IL2RA-high) population (a positive regulator of IL2RA); `pos|*` = enriches (a negative regulator).
- **Arce2025_Screen.csv** — Arce et al. 2025 (from the source repo's comparison set; primary citation to be
  confirmed before external use). Fitness/proliferation screen across Resting_Teff, Stimulated_Teff and
  Resting_Treg. MAGeCK `neg|*` = knockdown depletes cells (gene required for fitness/expansion); `pos|*` =
  enriches (a suppressor of fitness).

**Held out, NOT vendored:** `SchmidtSteinhart2022_CRISPRi_screen_gene_phenotypes.csv` and
`Schmidt2022_hits_Supplementary_table_2.xlsx` remain the validation set (RULE #3) and must never enter a
discovery score. They are read only where the pipeline already uses Schmidt for held-out validation.

These screens measure T-cell *function/fitness*, not autoimmune therapeutic direction. Concordance here is
evidence a gene is a genuine, replicated functional regulator (not a screen artifact); it does not resolve
direction of effect (N17 showed that remains unresolved).
