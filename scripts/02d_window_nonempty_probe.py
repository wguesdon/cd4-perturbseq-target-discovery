"""Is the therapeutic window EMPTY?

The crude window proxy dropped AUROC-vs-Schmidt from 0.702 to 0.502. Two readings:

  (a) The safety gate is correctly demoting toxic-but-effective genes (CD3, VAV1, LCK), and
      Schmidt IL-2 hits are mostly those genes, so a low AUROC against them is EXPECTED and
      not a failure.
  (b) The safety gate deletes everything that works, the window is empty, and the project
      has no novel hits to show.

These are distinguished by one question: after gating, does any gene retain BOTH strong
effector suppression AND independent Schmidt IL-2 support?

Small CSVs only. Never touches the h5ad.
"""

from __future__ import annotations

import pathlib

import sys

import numpy as np
import pandas as pd

# Resolve the checkout from this file, so the script runs from any clone.
REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, f"{REPO}/src")

from cd4_perturbseq import priors  # noqa: E402

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 30)

ranked = pd.read_csv(f"{REPO}/results/tables/risk_kill_naive_reversal.csv").rename(
    columns={"target_contrast_gene_name": "gene_name"}
)
schmidt = priors.schmidt_cd4_il2_screen()
iei = priors.iei_genes()
essential = priors.core_essential_genes()
classes = priors.druggable_classes()

df = ranked.merge(schmidt, on="gene_name", how="left")
df["rest_de_genes"] = df["rest_de_genes"].fillna(df["rest_de_genes"].median())

print("distributions over the 6371 rankable perturbations:")
for col in ("naive_suppression", "stim_de_genes", "rest_de_genes", "tolerance_suppression"):
    q = df[col].quantile([0.5, 0.75, 0.9, 0.95, 0.99]).round(2)
    print(f"  {col:22s} p50={q[0.5]:8.2f} p75={q[0.75]:8.2f} p90={q[0.90]:8.2f} p95={q[0.95]:8.2f} p99={q[0.99]:8.2f}")

# ---------------------------------------------------------------- the safety gate
REST_MAX = df["rest_de_genes"].quantile(0.75)
STIM_MAX = df["stim_de_genes"].quantile(0.90)
TOL_MAX = df["tolerance_suppression"].quantile(0.75)

df["fail_iei"] = df["gene_name"].isin(iei)
df["fail_essential"] = df["gene_name"].isin(essential)
df["fail_rest"] = df["rest_de_genes"] > REST_MAX
df["fail_collateral"] = df["stim_de_genes"] > STIM_MAX
df["fail_tolerance"] = df["tolerance_suppression"] > TOL_MAX

fail_cols = ["fail_iei", "fail_essential", "fail_rest", "fail_collateral", "fail_tolerance"]
df["safe"] = ~df[fail_cols].any(axis=1)

print(f"\ngate thresholds: rest_de<={REST_MAX:.0f}  stim_de<={STIM_MAX:.0f}  tol_supp<={TOL_MAX:.2f}")
print(f"safety-passing perturbations: {int(df['safe'].sum())} / {len(df)}")

top100 = df.nsmallest(100, "rank")
print(f"\nHEAD-TO-HEAD: of the naive top 100, how many pass the safety gate? {int(top100['safe'].sum())}")
print("rejection reasons among the naive top 100:")
for c in fail_cols:
    print(f"  {c:18s} {int(top100[c].sum()):3d}")

# ---------------------------------------------------------------- do safe hits have efficacy?
df["il2_support"] = (df["il2_neg_fdr"] < 0.05) & (df["il2_lfc"] < 0)

safe = df[df["safe"]].sort_values("naive_suppression", ascending=False)
print(f"\nsafe genes with independent Schmidt IL-2 support: {int(safe['il2_support'].sum())}")
print(f"unsafe genes with Schmidt IL-2 support:            {int(df.loc[~df['safe'],'il2_support'].sum())}")

def annotate(g: str) -> str:
    """Return the druggable classes a gene belongs to."""
    hits = [name for name, members in classes.items() if g in members]
    return ",".join(hits) if hits else "-"

top_safe = safe.head(30).copy()
top_safe["class"] = top_safe["gene_name"].map(annotate)
cols = ["rank", "gene_name", "naive_suppression", "stim_de_genes", "rest_de_genes",
        "tolerance_suppression", "il2_lfc", "il2_neg_fdr", "il2_support", "class"]
print("\n=== TOP 30 SAFETY-PASSING BY EFFECTOR SUPPRESSION ===")
print(top_safe[cols].to_string(index=False))

print("\n=== safe AND independently IL-2-supported (the money list) ===")
money = safe[safe["il2_support"]].copy()
money["class"] = money["gene_name"].map(annotate)
if len(money):
    print(money[cols].to_string(index=False))
else:
    print("EMPTY. The therapeutic window contains nothing with independent efficacy support.")

# ---------------------------------------------------------------- what did the gate reject?
print("\n=== naive top 25, with gate verdict (the demo figure) ===")
t25 = df.nsmallest(25, "rank").copy()
t25["verdict"] = np.where(t25["safe"], "PASS", "REJECT")
t25["reason"] = t25.apply(
    lambda r: ",".join(c.replace("fail_", "") for c in fail_cols if r[c]) or "-", axis=1
)
print(t25[["rank", "gene_name", "naive_suppression", "verdict", "reason", "il2_support"]].to_string(index=False))
