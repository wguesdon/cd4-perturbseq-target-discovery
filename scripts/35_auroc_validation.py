"""Rank-discrimination of the primary efficacy statistic, on two labelled sets.

Both areas under the receiver operating characteristic curve quoted in the manuscript were computed
on `naive_suppression`, a score written by `scripts/02_risk_kill_reversal.py` that is neither the
primary efficacy statistic nor the quantity `scripts/04_window_score.py` validates. The three scores
order perturbations differently. This script computes both curves on the single primary statistic and
writes them to a committed table, so the manuscript stops quoting literals.

Two labelled sets, with different standing:

- **Curated approved-immunomodulator targets.** The evaluation set. Discrimination here is expected to
  be weak, because most positives are not eligible under the quality-control mask.
- **Schmidt & Steinhart 2022 IL-2-reducing hits.** HELD OUT (RULE #3). Used only to validate the
  efficacy axis out of sample. It is never a score, a filter, or a threshold. Its hits are the
  receptor-proximal signalosome, which the triage layer rejects on purpose, so a triage score that
  beat the efficacy axis here would be suspect rather than reassuring.

Usage:
    uv run python scripts/35_auroc_validation.py
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from cd4_perturbseq import paths, priors  # noqa: E402

SEED = 0
N_BOOT = 2_000
RNG = np.random.default_rng(SEED)


def auroc_ci(labels: np.ndarray, scores: np.ndarray) -> tuple[float, float, float]:
    """Area under the ROC curve with a bootstrap percentile interval.

    Args:
        labels: Boolean positives.
        scores: Ranking score; higher means more positive.

    Returns:
        Tuple of (point estimate, lower bound, upper bound) at the 95% level.
    """
    point = float(roc_auc_score(labels, scores))
    boots = []
    n = labels.size
    for _ in range(N_BOOT):
        idx = RNG.integers(0, n, n)
        if np.unique(labels[idx]).size < 2:
            continue
        boots.append(roc_auc_score(labels[idx], scores[idx]))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return point, float(lo), float(hi)


def main() -> None:
    """Compute both curves on the primary statistic and write the table."""
    window = pd.read_csv(paths.TABLES / "window_score.csv")
    truth = pd.read_csv(paths.GROUND_TRUTH / "immunomodulator_targets.csv")
    positives = set(truth.loc[truth["include_as_positive"], "gene_symbol"])

    rows = []

    labels = window["gene_name"].isin(positives).to_numpy()
    auc, lo, hi = auroc_ci(labels, window["efficacy"].to_numpy())
    rows.append({
        "labelled_set": "curated approved-immunomodulator targets",
        "n_positives": int(labels.sum()), "n_total": len(window),
        "auroc": round(auc, 3), "ci_low": round(lo, 3), "ci_high": round(hi, 3),
        "spans_chance": bool(lo <= 0.5 <= hi),
        "role": "evaluation set; discrimination expected to be weak because most positives are "
                "ineligible under the quality-control mask",
    })

    schmidt = priors.schmidt_cd4_il2_screen()
    merged = window.merge(schmidt, on="gene_name", how="inner")
    merged["il2_hit"] = (merged["il2_neg_fdr"] < 0.05) & (merged["il2_lfc"] < 0)
    labels = merged["il2_hit"].to_numpy()
    auc, lo, hi = auroc_ci(labels, merged["efficacy"].to_numpy())
    rows.append({
        "labelled_set": "Schmidt 2022 IL-2-reducing hits (HELD OUT)",
        "n_positives": int(labels.sum()), "n_total": len(merged),
        "auroc": round(auc, 3), "ci_low": round(lo, 3), "ci_high": round(hi, 3),
        "spans_chance": bool(lo <= 0.5 <= hi),
        "role": "out-of-sample validation of the efficacy axis only; never a score, filter or threshold",
    })

    frame = pd.DataFrame(rows)
    out = paths.TABLES / "auroc_validation.csv"
    frame.to_csv(out, index=False)

    print("Rank discrimination of the primary efficacy statistic\n")
    print(frame[["labelled_set", "n_positives", "n_total", "auroc", "ci_low", "ci_high",
                 "spans_chance"]].to_string(index=False))
    print(f"\nwrote {out}")
    print("\nSuperseded: AUROC 0.542 [0.373, 0.707] and 0.702 [0.591, 0.814] were computed on")
    print("`naive_suppression`, a different score. Both are recomputed above on the primary statistic.")


if __name__ == "__main__":
    main()
