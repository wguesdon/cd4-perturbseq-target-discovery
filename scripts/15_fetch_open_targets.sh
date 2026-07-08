#!/usr/bin/env bash
# Fetch the Open Targets tables needed to bound the drug-target-recovery benchmark.
#
# Public, no authentication, no API key, ~22 MB. This is the same evidence the Open Targets MCP
# connector serves, taken from the bulk release so the result is reproducible from a pinned
# version rather than from whatever the API returned on the day.
#
# The version is PINNED. Open Targets reissues data quarterly and target-disease-drug links move
# between releases. A benchmark ceiling quoted without a release version is not a number.
#
# Written to data/external/open_targets/, which is gitignored. Derived tables under
# results/tables/ and resources/ground_truth/ are committed instead.
#
# Usage:
#     bash scripts/15_fetch_open_targets.sh

set -euo pipefail

OT_VERSION="26.06"
BASE="https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/${OT_VERSION}/output"
DEST="$(cd "$(dirname "$0")/.." && pwd)/data/external/open_targets"

mkdir -p "$DEST"
echo "Open Targets ${OT_VERSION} -> ${DEST}"
echo "${OT_VERSION}" > "${DEST}/VERSION"

fetch_single() {
    local dataset="$1" filename="$2"
    if [[ -s "${DEST}/${filename}" ]]; then
        echo "  have  ${filename}"
        return
    fi
    echo "  get   ${filename}"
    curl -fsSL --retry 3 --max-time 300 -o "${DEST}/${filename}" "${BASE}/${dataset}/${filename}"
}

# One parquet each.
#
# clinical_INDICATION, not clinical_TARGET. The latter reports the maximum clinical stage a
# drug-target pair reached across ALL indications, so joining it to an immune disease conjoins an
# approval earned in oncology with a trial that never finished. Vorinostat is approved for cutaneous
# T-cell lymphoma and merely trialled in Crohn's; read that way it becomes an approved
# immunomodulator and HDAC1 becomes a benchmark positive. clinical_indication carries the stage per
# indication, so the approval and the disease are the same fact.
fetch_single clinical_indication  clinical_indication.parquet
fetch_single disease              disease.parquet

# drug_molecule and drug_mechanism_of_action are partitioned; discover the parts.
fetch_parts() {
    local dataset="$1" prefix="$2"
    local parts
    parts=$(curl -fsSL --max-time 60 "${BASE}/${dataset}/" \
            | grep -oE 'href="[^"]+\.parquet"' | sed 's/href="//; s/"//' | sort -u)
    if [[ -z "${parts}" ]]; then
        echo "  ERROR: no parquet parts found under ${dataset}/" >&2
        exit 1
    fi
    local i=0
    for part in ${parts}; do
        local out="${prefix}_${i}.parquet"
        if [[ -s "${DEST}/${out}" ]]; then
            echo "  have  ${out}"
        else
            echo "  get   ${out}  (${part})"
            curl -fsSL --retry 3 --max-time 300 -o "${DEST}/${out}" "${BASE}/${dataset}/${part}"
        fi
        i=$((i + 1))
    done
}

fetch_parts drug_mechanism_of_action moa
fetch_parts drug_molecule           drug_molecule

# 85 MB. Needed only to map Ensembl ids to gene symbols for targets OUTSIDE our 10,282-gene
# measured panel. Without it, "Open Targets does not know this drug target" is indistinguishable
# from "our transcriptome does not measure it", and IL17A, DHODH, JAK1 and JAK3 all look like
# curation errors when they are nothing of the kind.
fetch_parts target target

echo
du -sh "${DEST}"
ls -1 "${DEST}"
echo
echo "Next: uv run python scripts/16_open_targets_benchmark.py"
