"""Confirm the IUIS parsing fix does what it claims and nothing else."""

from __future__ import annotations

import pathlib

import sys

# Resolve the checkout from this file, so the script runs from any clone.
REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, f"{REPO}/src")

from cd4_perturbseq import priors  # noqa: E402

lof = priors.iei_genes()
allmech = priors.iei_genes(lof_only=False)

print(f"iei_genes(lof_only=True) : {len(lof)}")
print(f"iei_genes(lof_only=False): {len(allmech)}")
print(f"dropped as gain-of-function: {len(allmech - lof)}")
print("  " + ", ".join(sorted(allmech - lof)))

print("\nC4A / C4B recovered?")
for g in ("C4A", "C4B"):
    print(f"  {g}: {'YES' if g in lof else 'NO'}")

print("\nGOF genes must NOT be flagged (they were, before):")
for g in ("PLCG1", "STAT6", "SYK", "NLRP3", "CXCR4"):
    print(f"  {g}: {'STILL FLAGGED (bug)' if g in lof else 'correctly excluded'}")

print("\nreal LoF immunodeficiency genes must STILL be flagged:")
for g in ("ZAP70", "CD3E", "CD3G", "CD247", "ITK", "JAK3", "IL2RA", "RAG1", "LCK"):
    print(f"  {g}: {'flagged' if g in lof else '!! LOST !!'}")

print("\nentries that still fail the symbol regex (should be cytogenetic, not genes):")
import pandas as pd  # noqa: E402

raw = pd.read_csv(f"{REPO}/data/external/gwt_priors/IUIS-IEI-list-July-2024V2.csv")
raw.columns = [c.strip() for c in raw.columns]
defects = raw["Genetic defect"].astype(str).str.strip()
exploded = defects.str.split("+").explode().str.strip()
unmatched = sorted(set(exploded[~exploded.str.fullmatch(r"[A-Z][A-Z0-9\-]{1,14}").fillna(False)]))
for u in unmatched:
    print(f"  {u!r}")
