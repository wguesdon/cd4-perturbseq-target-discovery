#!/usr/bin/env bash
# Fetch the Open Targets genetic-association scores for the 14 autoimmune diseases the source paper
# used (Zhu, Dann et al., Methods p.50: "genetic evidence score >= 0.1"), so N10 can test whether
# the safety-gated ranking enriches for genetics-supported autoimmune targets.
#
# The full association_by_datatype_direct dataset is ~1 GB across 15 parquet parts. We stream each
# part, keep only rows whose diseaseId is one of the 14 and whose datatypeId is genetic_association,
# append to one small parquet, and delete the part. Peak disk stays ~150 MB.
#
# Public EBI FTP, no authentication, pinned to 26.06 to match scripts/15. Written to
# data/external/open_targets/ (gitignored). The derived per-gene label is committed under results/.
#
# Usage:
#     bash scripts/19_fetch_ot_genetic_evidence.sh

set -euo pipefail

OT_VERSION="26.06"
BASE="https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/${OT_VERSION}/output"
DATASET="association_by_datatype_direct"
DEST="$(cd "$(dirname "$0")/.." && pwd)/data/external/open_targets"
OUT="${DEST}/autoimmune_genetic_assoc.parquet"

mkdir -p "${DEST}"
if [[ -s "${OUT}" ]]; then
    echo "have ${OUT}; delete it to refetch"
    exit 0
fi

parts=$(curl -fsSL --max-time 60 "${BASE}/${DATASET}/" \
        | grep -oE 'href="[^"]+\.parquet"' | sed 's/href="//; s/"//' | sort -u)
n=$(echo "${parts}" | grep -c parquet)
echo "association_by_datatype_direct: ${n} parts, filtering to 14 autoimmune diseases x genetic_association"

tmp="${DEST}/.assoc_part.parquet"
i=0
for part in ${parts}; do
    echo "  [$((i + 1))/${n}] ${part}"
    curl -fsSL --retry 3 --max-time 600 -o "${tmp}" "${BASE}/${DATASET}/${part}"
    FILTER_PART="${tmp}" ACCUM="${OUT}" uv run python "$(dirname "$0")/_filter_ot_assoc.py"
    rm -f "${tmp}"
    i=$((i + 1))
done

echo
uv run python -c "import pandas as pd; d=pd.read_parquet('${OUT}'); print('rows:', len(d), 'unique targets:', d.targetId.nunique(), 'diseases:', d.diseaseId.nunique())"
echo "wrote ${OUT}"
