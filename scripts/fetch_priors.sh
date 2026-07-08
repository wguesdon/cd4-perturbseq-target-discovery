#!/usr/bin/env bash
# Fetch the prior gene lists and reference screens bundled by the source paper's
# analysis repo (emdann/GWT_perturbseq_analysis_2025, branch master).
#
# We pull individual files over raw.githubusercontent rather than cloning, because
# the repo is 1.66 GB and we need ~20 small tables.
#
# Usage: bash scripts/fetch_priors.sh
set -euo pipefail

REPO_RAW="https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/master"
DEST="data/external/gwt_priors"
mkdir -p "$DEST"

FILES=(
  # --- safety axis anchors ---
  "metadata/gene_lists/core_essentials_hart.tsv"          # Hart CEG2 core-essential genes
  "metadata/IUIS-IEI-list-July-2024V2.csv"                # inborn errors of immunity (human LoF -> immunodeficiency)
  "metadata/gene_lists/clinvar_path_likelypath.tsv"       # ClinVar pathogenic / likely pathogenic

  # --- druggable-class annotation ---
  "metadata/gene_lists/kinases.tsv"
  "metadata/gene_lists/gpcr_union.tsv"
  "metadata/gene_lists/ion_channels.tsv"
  "metadata/gene_lists/enzymes.tsv"
  "metadata/gene_lists/transporters.tsv"
  "metadata/gene_lists/nuclear_receptors.tsv"
  "metadata/gene_lists/catalytic_receptors.tsv"

  # --- program definition ---
  "metadata/immune_effector_genes.csv"
  "metadata/gene_lists/cytokines.tsv"
  "metadata/gene_lists/cytokine_receptors.tsv"
  "metadata/Lambert_2018_HumanTF.csv"

  # --- genetics ---
  "metadata/gene_lists/gwascatalog.tsv"

  # --- orthogonal external screens / signatures ---
  "metadata/Arce2024_20230130_DESeq2_output_AAVS1_Teff_Stimulation_vs_Resting.csv"
  "metadata/SchmidtSteinhart2022_CRISPRi_screen_gene_phenotypes.csv"
  "metadata/Freimer2022_Screen.csv"
  "metadata/Umhoefer2025_FOXP3_Teff.csv"

  # --- authors' own descriptive context-specificity + cluster annotation ---
  "metadata/nde75ntotal50_cluster_specificity.csv"
  "metadata/suppl_tables/clustering_results_and_annotations.csv"
  "metadata/suppl_tables/cluster_autoimmune_enrichment_results.suppl_table.csv"
)

echo "Fetching ${#FILES[@]} prior files -> $DEST"
for f in "${FILES[@]}"; do
  out="$DEST/$(basename "$f")"
  if [[ -s "$out" ]]; then
    echo "  skip (exists): $(basename "$f")"
    continue
  fi
  if curl -fsSL --retry 3 --retry-delay 2 "$REPO_RAW/$f" -o "$out"; then
    printf '  ok  %8s  %s\n' "$(du -h "$out" | cut -f1)" "$(basename "$f")"
  else
    echo "  FAIL             $f" >&2
    rm -f "$out"
  fi
done

echo
echo "Done. Provenance: emdann/GWT_perturbseq_analysis_2025@master (MIT)."
