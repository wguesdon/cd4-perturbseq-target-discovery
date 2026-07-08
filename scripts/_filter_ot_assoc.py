"""Filter one Open Targets association part to the 14 autoimmune diseases and the genetic datatype.

Helper for scripts/19_fetch_ot_genetic_evidence.sh. Reads the part named by $FILTER_PART, keeps rows
whose diseaseId is one of the 14 autoimmune MONDO ids and whose aggregationValue is
``genetic_association``, and appends the kept rows to the parquet at $ACCUM.
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

# The 14 autoimmune diseases from Zhu, Dann et al. (their table minus the three negative controls:
# age-related macular degeneration, chronic kidney disease, coronary artery disease). MONDO ids
# resolved from data/external/open_targets/disease.parquet, all 14 verified present in the
# association table.
AUTOIMMUNE = {
    "MONDO_0005011": "Crohn disease",
    "MONDO_0007699": "Hashimoto thyroiditis",
    "MONDO_0004979": "asthma",
    "MONDO_0004980": "atopic eczema",
    "MONDO_0005083": "psoriasis",
    "MONDO_0005101": "ulcerative colitis",
    "MONDO_0005130": "celiac disease",
    "MONDO_0005147": "type 1 diabetes mellitus",
    "MONDO_0005265": "inflammatory bowel disease",
    "MONDO_0005301": "multiple sclerosis",
    "MONDO_0005306": "ankylosing spondylitis",
    "MONDO_0007179": "autoimmune disease",
    "MONDO_0007915": "systemic lupus erythematosus",
    "MONDO_0008383": "rheumatoid arthritis",
}


def main() -> None:
    """Filter the part and append to the accumulator."""
    part = Path(os.environ["FILTER_PART"])
    accum = Path(os.environ["ACCUM"])

    df = pd.read_parquet(part, columns=["diseaseId", "targetId", "aggregationValue", "associationScore"])
    keep = df[
        (df["aggregationValue"] == "genetic_association") & (df["diseaseId"].isin(AUTOIMMUNE))
    ][["diseaseId", "targetId", "associationScore"]].copy()

    if accum.exists() and accum.stat().st_size:
        keep = pd.concat([pd.read_parquet(accum), keep], ignore_index=True)
    keep.to_parquet(accum, index=False)
    print(f"    kept {len(keep):,} genetic-association rows so far "
          f"({keep['targetId'].nunique():,} targets, {keep['diseaseId'].nunique()} diseases)")


if __name__ == "__main__":
    main()
