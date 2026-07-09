"""Is the effector-suppression score measuring anything real?

Three independent checks, all on small CSVs. Never touches the h5ad.

1. AUROC of the naive suppression score against the 20 rankable approved-drug positives.
   Tiny positive set, so bootstrap the CI rather than trusting the point estimate.
2. Concordance with Schmidt & Steinhart 2022, a genome-wide CRISPRi screen for CD4+ IL-2
   PRODUCTION. This is a protein-level readout from a different lab, different assay,
   different modality. If our transcriptome-derived score is real, knockdowns that suppress
   our effector program should also reduce IL-2 protein: naive_suppression should correlate
   NEGATIVELY with the IL-2 log fold change.
3. AUROC against Schmidt's own significant IL-2-reducing hits. Large positive set, tight CI.
   This is the efficacy validation that the 20-positive drug benchmark cannot deliver.
"""

from __future__ import annotations

import pathlib

import sys

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

# Resolve the checkout from this file, so the script runs from any clone.
REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, f"{REPO}/src")

from cd4_perturbseq import priors  # noqa: E402

RNG = np.random.default_rng(0)


def auroc_ci(y: np.ndarray, score: np.ndarray, n_boot: int = 2000) -> tuple[float, float, float]:
    """Bootstrap a percentile CI for AUROC.

    Args:
        y: Binary labels.
        score: Ranking score, higher is more positive.
        n_boot: Bootstrap resamples.

    Returns:
        Tuple of (auroc, lower 2.5th percentile, upper 97.5th percentile).
    """
    point = roc_auc_score(y, score)
    n = len(y)
    boots = []
    for _ in range(n_boot):
        idx = RNG.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            continue
        boots.append(roc_auc_score(y[idx], score[idx]))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return point, lo, hi


ranked = pd.read_csv(f"{REPO}/results/tables/risk_kill_naive_reversal.csv")
ranked = ranked.rename(columns={"target_contrast_gene_name": "gene_name"})
print(f"rankable perturbations: {len(ranked)}")

# ---------------------------------------------------------------- 1. drug positives
gt = pd.read_csv(f"{REPO}/resources/ground_truth/immunomodulator_targets.csv")
positives = set(gt.loc[gt["include_as_positive"] == True, "gene_symbol"])  # noqa: E712

ranked["is_positive"] = ranked["gene_name"].isin(positives)
n_pos = int(ranked["is_positive"].sum())
print(f"\n[1] approved-drug positives rankable: {n_pos}")

auc, lo, hi = auroc_ci(ranked["is_positive"].to_numpy(), ranked["naive_suppression"].to_numpy())
print(f"    AUROC naive_suppression = {auc:.3f}   95% CI [{lo:.3f}, {hi:.3f}]")
print(f"    {'ABOVE chance' if lo > 0.5 else 'NOT distinguishable from chance (CI covers 0.5)'}")

# precision at k
for k in (50, 100, 250, 500):
    hits = int(ranked.nsmallest(k, "rank")["is_positive"].sum())
    expected = n_pos * k / len(ranked)
    print(f"    top-{k:<4d} positives: {hits:2d}   expected by chance {expected:4.1f}   enrichment {hits/max(expected,1e-9):5.2f}x")

# ---------------------------------------------------------------- 2. Schmidt concordance
schmidt = priors.schmidt_cd4_il2_screen()
print(f"\n[2] Schmidt CD4+ IL-2 screen genes: {len(schmidt)}")

merged = ranked.merge(schmidt, on="gene_name", how="inner")
merged = merged.dropna(subset=["naive_suppression", "il2_lfc"])
print(f"    overlap with our rankable perturbations: {len(merged)}")

rho, p = stats.spearmanr(merged["naive_suppression"], merged["il2_lfc"])
print(f"    spearman(naive_suppression, il2_lfc) = {rho:+.3f}  p={p:.3g}")
print("    expected NEGATIVE: suppressing the effector program should reduce IL-2 protein")
print(f"    -> {'CONCORDANT' if rho < 0 and p < 0.05 else 'NOT concordant'}")

# ---------------------------------------------------------------- 3. Schmidt hits as positives
merged["il2_hit"] = (merged["il2_neg_fdr"] < 0.05) & (merged["il2_lfc"] < 0)
n_hit = int(merged["il2_hit"].sum())
print(f"\n[3] Schmidt significant IL-2-reducing hits among our rankable set: {n_hit}")

if n_hit >= 10:
    auc2, lo2, hi2 = auroc_ci(merged["il2_hit"].to_numpy(), merged["naive_suppression"].to_numpy())
    print(f"    AUROC naive_suppression vs Schmidt IL-2 hits = {auc2:.3f}  95% CI [{lo2:.3f}, {hi2:.3f}]")
    print(f"    {'ABOVE chance' if lo2 > 0.5 else 'not above chance'}")

    # Does the safety axis destroy efficacy? Compare AUROC of a naive-minus-penalty proxy.
    z = lambda s: (s - s.mean()) / s.std()  # noqa: E731
    merged["rest_pen"] = z(np.log1p(merged["rest_de_genes"].fillna(0)))
    merged["tol_pen"] = z(merged["tolerance_suppression"].fillna(0))
    merged["window_proxy"] = z(merged["naive_suppression"]) - merged["rest_pen"] - merged["tol_pen"]
    auc3, lo3, hi3 = auroc_ci(merged["il2_hit"].to_numpy(), merged["window_proxy"].to_numpy())
    print(f"\n    crude window proxy (efficacy - rest_penalty - tolerance_penalty)")
    print(f"    AUROC vs Schmidt IL-2 hits = {auc3:.3f}  95% CI [{lo3:.3f}, {hi3:.3f}]")
    print(f"    delta vs naive = {auc3 - auc2:+.3f}")
    print("    (a big drop means the safety gate is throwing away real efficacy;")
    print("     a small drop means it is removing toxicity at little efficacy cost)")
