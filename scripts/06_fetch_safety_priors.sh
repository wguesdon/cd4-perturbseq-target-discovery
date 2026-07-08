#!/usr/bin/env bash
# Fetch the organism-level safety priors.
#
# Our safety gate, as first built, was T-cell-intrinsic: "rest" meant resting CD4 T cells. It
# could not see that a mitochondrial or one-carbon gene would be toxic to gut epithelium and
# bone marrow. These four files close that gap.
#
#   gnomAD constraint  Is this gene loss-of-function intolerant in living humans? LOEUF is the
#                      upper bound of the observed/expected LoF ratio. Below ~0.35 means strong
#                      selection against heterozygous LoF, which is the clearest organism-level
#                      warning that exists.
#   HPA tissue RNA     Is this gene expressed everywhere, or only in immune cells? A ubiquitous
#                      housekeeping gene has no therapeutic window outside the target tissue.
#   Hart CEGv2         The FULL 684-gene core-essential set. The copy bundled by the source
#                      paper's repo is a 283-symbol subset with only 68 genes in this library,
#                      so it could never adjudicate essentiality.
#   Hart NEGv1         927 nonessential genes. The negative control our viability axis has never
#                      had: core essentials should deplete at rest, nonessentials should not.
#
# DepMap is deliberately absent. Its download API sits behind a browser bot-check and returns an
# HTML verification page to curl, so it cannot be fetched reproducibly from a script. Saying so
# is better than shipping a pipeline that silently parses HTML as CSV.
#
# Usage: bash scripts/06_fetch_safety_priors.sh
set -euo pipefail

DEST="data/external/safety_priors"
mkdir -p "$DEST"

fetch() {
  local url="$1" out="$DEST/$2"
  if [[ -s "$out" ]]; then
    printf '  skip (exists) %8s  %s\n' "$(du -h "$out" | cut -f1)" "$2"
    return
  fi
  if curl -fsSL --retry 3 --retry-delay 2 "$url" -o "$out"; then
    printf '  ok            %8s  %s\n' "$(du -h "$out" | cut -f1)" "$2"
  else
    echo "  FAIL  $2  <- $url" >&2
    rm -f "$out"
    return 1
  fi
}

echo "Fetching organism-level safety priors -> $DEST"

# gnomAD v4.1 constraint does not exist at the v4.1 prefix; v4.0 is the current published table.
fetch "https://storage.googleapis.com/gcp-public-data--gnomad/release/v4.0/constraint/gnomad.v4.0.constraint_metrics.tsv" \
      "gnomad.v4.0.constraint_metrics.tsv"

fetch "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip" \
      "rna_tissue_consensus.tsv.zip"

fetch "https://raw.githubusercontent.com/hart-lab/bagel/master/CEGv2.txt" "hart_CEGv2.txt"
fetch "https://raw.githubusercontent.com/hart-lab/bagel/master/NEGv1.txt" "hart_NEGv1.txt"

if [[ -s "$DEST/rna_tissue_consensus.tsv.zip" && ! -s "$DEST/rna_tissue_consensus.tsv" ]]; then
  echo "  unzipping HPA consensus table"
  unzip -o -q "$DEST/rna_tissue_consensus.tsv.zip" -d "$DEST"
  printf '  ok            %8s  %s\n' "$(du -h "$DEST/rna_tissue_consensus.tsv" | cut -f1)" "rna_tissue_consensus.tsv"
fi

echo
echo "Done. Provenance:"
echo "  gnomAD v4.0 constraint metrics (Broad Institute, public GCS bucket)"
echo "  Human Protein Atlas consensus tissue RNA (proteinatlas.org, CC BY-SA 3.0)"
echo "  Hart CEGv2 / NEGv1 (hart-lab/bagel, MIT)"
echo "  DepMap: NOT fetched. Its API is behind a browser bot-check."
