"""Does the cell dropout happen only in stimulated cells, or at rest too?

Resting CD4 T cells are quiescent. Stimulated ones proliferate hard. So:

  depleted in Stim, normal at Rest  ->  antiproliferative. Mycophenolate and methotrexate work
                                        exactly this way. This IS a therapeutic window.
  depleted in BOTH                  ->  core-essential. The cell dies regardless of context.
                                        That is chemotherapy, not an immunomodulator.

Discriminator: log2(n_cells_target in Stim48hr / n_cells_target at Rest), per perturbed gene.
Negative means stimulation-selective dropout.

Guides are assigned before stimulation and each condition is a separate culture from the same
donors, so the resting arm is a within-experiment control for guide representation.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd
from scipy import stats

REPO = "/home/will/Documents/Github/cd4-perturbseq-target-discovery"
sys.path.insert(0, f"{REPO}/src")

from cd4_perturbseq import de_stats, priors  # noqa: E402

pd.set_option("display.width", 200)

obs = de_stats.read_obs()
obs["gene_name"] = obs["target_contrast_gene_name"].astype(str)

cells = obs.pivot_table(index="gene_name", columns="culture_condition",
                        values="n_cells_target", aggfunc="first", observed=True)
cells = cells.dropna(subset=["Rest", "Stim48hr"])
cells["log2_stim_over_rest"] = np.log2(cells["Stim48hr"] / cells["Rest"])
print(f"genes with both conditions: {len(cells)}")
print(f"global median log2(Stim/Rest) = {cells['log2_stim_over_rest'].median():+.3f}")

essential = priors.core_essential_genes()
cells["is_essential"] = cells.index.isin(essential)

ess = cells[cells["is_essential"]]
non = cells[~cells["is_essential"]]
print(f"\nHart core-essentials present: {len(ess)}")
print(f"  essential   median log2(Stim/Rest) = {ess['log2_stim_over_rest'].median():+.3f}")
print(f"  other       median log2(Stim/Rest) = {non['log2_stim_over_rest'].median():+.3f}")
_, p = stats.mannwhitneyu(ess["log2_stim_over_rest"], non["log2_stim_over_rest"], alternative="less")
print(f"  MWU essentials MORE stim-depleted: p = {p:.3g}")

# Absolute resting depletion: essentials should also be low at rest relative to the median.
rest_med = cells["Rest"].median()
print(f"\n  median resting cell count overall: {rest_med:.0f}")
print(f"  essential   median resting cells = {ess['Rest'].median():.0f}  ({ess['Rest'].median()/rest_med:.2f}x)")
print(f"  other       median resting cells = {non['Rest'].median():.0f}  ({non['Rest'].median()/rest_med:.2f}x)")
_, p2 = stats.mannwhitneyu(ess["Rest"], non["Rest"], alternative="less")
print(f"  MWU essentials depleted AT REST: p = {p2:.3g}")
print("  -> resting depletion is the context-free viability signal we want to gate on.")

frame = pd.read_csv(f"{REPO}/results/tables/window_score.csv")
merged = frame.merge(cells.reset_index()[["gene_name", "Rest", "Stim48hr", "log2_stim_over_rest"]],
                     on="gene_name", how="left")
merged["rest_cells_ratio"] = merged["Rest"] / rest_med

print("\n=== the discriminator, applied to the genes that matter ===")
probe = ["PPP3R1", "IL4R", "CD2", "IMPDH2", "CD3E", "CD3G", "ENO1", "MTHFD2", "TFB1M",
         "CARS2", "VARS2", "TOMM70", "CCNT1", "POLG", "RRAGC", "NDUFAF3", "ACAD9"]
cols = ["gene_name", "Rest", "Stim48hr", "rest_cells_ratio", "log2_stim_over_rest", "safe"]
view = merged[merged["gene_name"].isin(probe)][cols].sort_values("rest_cells_ratio")
print(view.to_string(index=False))

print("\n  interpretation:")
print("   rest_cells_ratio near or above 1 + negative log2(Stim/Rest) = stimulation-selective")
print("     dropout. The cell is fine until you activate it. Mycophenolate's mechanism.")
print("   rest_cells_ratio well below 1 = the cell dies at rest too. Core-essential toxicity.")

safe25 = merged[merged["safe"]].nlargest(25, "window_score")
print(f"\n=== current top-25 shortlist, viability audit ===")
print(f"  median resting cell ratio: {safe25['rest_cells_ratio'].median():.2f}")
print(f"  how many are depleted AT REST (ratio < 0.5)? {(safe25['rest_cells_ratio'] < 0.5).sum()} of 25")
depleted = safe25[safe25["rest_cells_ratio"] < 0.5]["gene_name"].tolist()
print(f"  they are: {', '.join(depleted) if depleted else 'none'}")
survivors = safe25[safe25["rest_cells_ratio"] >= 0.5]["gene_name"].tolist()
print(f"\n  shortlist genes that are NOT depleted at rest: {', '.join(survivors) if survivors else 'none'}")
