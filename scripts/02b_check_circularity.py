"""Quantify the self-effect artifact and the ground-truth circularity.

Two suspicions:
1. A perturbation of gene G is scored by the mean z-score over module M. If G is in M, its
   own on-target knockdown z-score (large, negative, and REQUIRED significant by QC) enters
   its own score and inflates it.
2. Ground-truth positives (IL17A, IL2, IFNG, TNF...) may themselves be members of the
   effector module, making their recovery self-fulfilling.

Uses only small CSVs. Never touches the h5ad.
"""

from __future__ import annotations

import pathlib

import sys

import pandas as pd

# Resolve the checkout from this file, so the script runs from any clone.
REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, f"{REPO}/src")

from cd4_perturbseq import programs  # noqa: E402

program = pd.read_csv(f"{REPO}/data/interim/activation_program.csv")
effector = set(programs.effector_core(program))
tolerance = set(programs.tolerance_module(program))

ranked = pd.read_csv(f"{REPO}/results/tables/risk_kill_naive_reversal.csv")
measured = set(open(f"{REPO}/data/interim/measured_genes.txt").read().split())
perturbed = set(open(f"{REPO}/data/interim/perturbed_genes.txt").read().split())

eff_measured = sorted(effector & measured)
tol_measured = sorted(tolerance & measured)
print(f"effector core (measured): {len(eff_measured)}")
print(f"tolerance     (measured): {len(tol_measured)}")

# --- 1. self-effect: which module members are themselves perturbed & ranked? ---
ranked_genes = set(ranked["target_contrast_gene_name"])
eff_self = sorted(set(eff_measured) & ranked_genes)
tol_self = sorted(set(tol_measured) & ranked_genes)
print(f"\n[1] effector-core genes that are ALSO perturbed and rankable: {len(eff_self)}")
print("   ", ", ".join(eff_self))
print(f"    tolerance genes that are ALSO perturbed and rankable: {len(tol_self)}")
print("   ", ", ".join(tol_self))

sub = ranked[ranked["target_contrast_gene_name"].isin(eff_self)][
    ["rank", "target_contrast_gene_name", "naive_suppression"]
].sort_values("rank")
print("\n    where do effector-core self-members rank?")
print(sub.to_string(index=False))

n = len(eff_measured)
print(f"\n    module size n={n}. An on-target z of -10 shifts the mean by {10/n:+.3f};")
print(f"    z of -20 -> {20/n:+.3f}; z of -30 -> {30/n:+.3f} (added to naive_suppression).")

top100 = ranked.nsmallest(100, "rank")
n_self_top100 = top100["target_contrast_gene_name"].isin(eff_self).sum()
print(f"    effector-core self-members inside the top 100: {n_self_top100}")

# --- 2. ground-truth circularity ---
gt = pd.read_csv(f"{REPO}/resources/ground_truth/immunomodulator_targets.csv")
pos = gt[gt["include_as_positive"]]["gene_symbol"].tolist()
print(f"\n[2] ground-truth positives: {len(pos)}")

pos_perturbed = sorted(set(pos) & perturbed)
pos_measured = sorted(set(pos) & measured)
pos_rankable = sorted(set(pos) & ranked_genes)
print(f"    ... perturbed in the library: {len(pos_perturbed)}")
print(f"    ... measured in the transcriptome: {len(pos_measured)}")
print(f"    ... perturbed AND surviving QC (rankable): {len(pos_rankable)}")
print("    rankable positives:", ", ".join(pos_rankable))
print("    LOST (perturbed but not rankable):", ", ".join(sorted(set(pos_perturbed) - set(pos_rankable))))
print("    NEVER PERTURBED:", ", ".join(sorted(set(pos) - set(pos_perturbed))))

circular = sorted(set(pos) & effector)
print(f"\n    positives that are THEMSELVES effector-module members: {len(circular)}")
print("   ", ", ".join(circular))
circular_rankable = sorted(set(pos_rankable) & effector)
print(f"    of the rankable positives, circular ones: {len(circular_rankable)} -> {', '.join(circular_rankable)}")

if pos_rankable:
    r = ranked[ranked["target_contrast_gene_name"].isin(pos_rankable)][
        ["rank", "target_contrast_gene_name", "naive_suppression"]
    ].sort_values("rank")
    r["in_effector_module"] = r["target_contrast_gene_name"].isin(effector)
    print("\n    where do the rankable positives sit in the NAIVE ranking?")
    print(r.to_string(index=False))
    print(f"\n    total rankable perturbations: {len(ranked)}")
