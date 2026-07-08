"""RISK-KILL: does a naive reversal ranking nominate targets a safety gate rejects?

The project's headline claim is "reversal is not enough": ranking perturbations purely
by how strongly they suppress an inflammatory program will surface genes that are
common-essential or that collapse the resting transcriptome, and those are
phenotype-but-toxic false positives rather than drug targets.

This script tests that claim directly, and it is designed to be able to FAIL. If the
top of the naive ranking is no more essential, no more immunodeficiency-linked, and no
more transcriptome-disrupting than background, the headline is wrong and the strategy
needs rewriting.

Naive score, per perturbed gene, in stimulated cells:

    suppression = -mean(zscore over effector program genes)   [Stim48hr rows]

Two design choices keep the result honest.

First, the naive baseline is STEELMANNED, not a straw man. We apply the QC any competent
analyst would apply before ranking: drop flagged off-target perturbations and require a
significant on-target knockdown. The claim is that even a well-QC'd reversal ranking is
toxic, which is a much stronger claim than "an unfiltered one is".

Second, we check the obvious confound. Perturbations assayed in more cells have more
statistical power, so they get more extreme z-scores AND more significant DE genes. That
alone could manufacture the enrichment we are looking for. We report the cell-count
distribution in the top K versus background, and repeat every test against a background
matched on cell-count decile.

Enrichment of liabilities in the top K is tested against background with Fisher's exact
test (binary flags) and the Mann-Whitney U test (continuous disruption measures).

Usage:
    uv run python scripts/02_risk_kill_reversal.py --top-k 100
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy import stats

from cd4_perturbseq import de_stats, paths, priors, programs

CONDITION = "Stim48hr"


def _fisher(flag: pd.Series, in_top: pd.Series) -> tuple[float, float, int, int]:
    """Fisher's exact test for a binary liability flag enriched in the top of a ranking.

    Args:
        flag: Boolean series, True where the perturbed gene carries the liability.
        in_top: Boolean series, True where the perturbed gene is in the top K.

    Returns:
        Tuple of (odds ratio, p-value, count in top, count in background).
    """
    a = int((flag & in_top).sum())
    b = int((~flag & in_top).sum())
    c = int((flag & ~in_top).sum())
    d = int((~flag & ~in_top).sum())
    odds, pval = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    return odds, pval, a, c


def apply_qc(stim: pd.DataFrame) -> pd.DataFrame:
    """Steelman the naive baseline by applying routine perturbation QC.

    Drops perturbations with an off-target flag and keeps only those with a significant
    on-target knockdown, since ranking perturbations that never knocked their target down
    would be ranking noise. Column names differ between the December 2025 supplementary
    CSV and the May 2026 h5ad, so each filter is applied only if its column exists.

    Args:
        stim: Rows of the obs frame for the stimulated condition.

    Returns:
        The QC-passing subset.
    """
    before = len(stim)
    offtarget_cols = [c for c in ("offtarget_flag", "distal_offtarget_flag", "neighboring_gene_KD") if c in stim.columns]
    for col in offtarget_cols:
        stim = stim[~stim[col].astype(bool)]
        print(f"  QC drop {col:24s} -> {len(stim):6d} remain")

    if "ontarget_significant" in stim.columns:
        stim = stim[stim["ontarget_significant"].astype(bool)]
        print(f"  QC keep ontarget_significant     -> {len(stim):6d} remain")
    if "low_target_gex" in stim.columns:
        stim = stim[~stim["low_target_gex"].astype(bool)]
        print(f"  QC drop low_target_gex           -> {len(stim):6d} remain")

    print(f"  QC total: {before} -> {len(stim)}")
    return stim


def _matched_background(frame: pd.DataFrame, n_bins: int = 10) -> pd.Series:
    """Flag a background matched to the top K on cell-count decile.

    For each decile of ``n_cells_target``, sample as many background perturbations as
    there are top-K members in that decile, so power cannot drive the enrichment.

    Args:
        frame: Ranked frame carrying ``n_cells_target`` and ``in_top``.
        n_bins: Number of quantile bins.

    Returns:
        Boolean series, True for the matched background rows.
    """
    rng = np.random.default_rng(0)
    bins = pd.qcut(frame["n_cells_target"], n_bins, labels=False, duplicates="drop")
    matched = pd.Series(False, index=frame.index)
    for b in sorted(pd.unique(bins.dropna())):
        in_bin = bins == b
        n_top = int((in_bin & frame["in_top"]).sum())
        if n_top == 0:
            continue
        pool = frame.index[in_bin & ~frame["in_top"]]
        take = min(n_top, len(pool))
        if take:
            matched.loc[rng.choice(pool, size=take, replace=False)] = True
    return matched


def build_scores(top_k: int) -> pd.DataFrame:
    """Compute the naive suppression ranking and attach safety annotations.

    Args:
        top_k: Number of top-ranked perturbations treated as the naive shortlist.

    Returns:
        DataFrame indexed by perturbed gene with the naive score, liability flags,
        and disruption measures.
    """
    obs = de_stats.read_obs()
    var = de_stats.read_var()

    gene_names = var["gene_name"].astype(str).to_numpy()
    name_to_idx = {g: i for i, g in enumerate(gene_names)}

    program = programs.load_activation_program()
    effector = [g for g in programs.effector_core(program) if g in name_to_idx]
    tolerance = [g for g in programs.tolerance_module(program) if g in name_to_idx]

    print(f"effector core genes measured: {len(effector)}")
    print(f"tolerance genes measured:     {len(tolerance)}")

    wanted = effector + tolerance
    columns = de_stats.read_layer_columns([name_to_idx[g] for g in wanted], layer="zscore")
    effector_block = columns[:, : len(effector)]
    tolerance_block = columns[:, len(effector) :]

    obs = obs.reset_index(drop=True)
    obs["effector_z"] = np.nanmean(effector_block, axis=1)
    obs["tolerance_z"] = np.nanmean(tolerance_block, axis=1)

    stim = obs[obs["culture_condition"] == CONDITION].copy()
    rest = obs[obs["culture_condition"] == "Rest"].copy()

    print(f"\nQC on the naive baseline ({CONDITION}):")
    stim = apply_qc(stim)

    stim["naive_suppression"] = -stim["effector_z"]
    stim["tolerance_suppression"] = -stim["tolerance_z"]

    gene_col = "target_contrast_gene_name"
    keep = [gene_col, "naive_suppression", "tolerance_suppression", "n_total_de_genes", "n_cells_target"]
    frame = stim[keep].copy()
    frame = frame.rename(columns={"n_total_de_genes": "stim_de_genes"})

    rest_de = rest.drop_duplicates(gene_col).set_index(gene_col)["n_total_de_genes"]
    frame["rest_de_genes"] = frame[gene_col].map(rest_de)

    essential = priors.core_essential_genes()
    iei = priors.iei_genes()
    frame["is_core_essential"] = frame[gene_col].isin(essential)
    frame["is_iei"] = frame[gene_col].isin(iei)

    frame = frame.dropna(subset=["naive_suppression", "n_cells_target"]).copy()
    frame = frame.sort_values("naive_suppression", ascending=False).reset_index(drop=True)
    frame["rank"] = np.arange(1, len(frame) + 1)
    frame["in_top"] = frame["rank"] <= top_k
    frame["matched_background"] = _matched_background(frame)
    return frame


def report(frame: pd.DataFrame, top_k: int) -> bool:
    """Print the enrichment tests and return whether the empirical bet holds.

    Args:
        frame: Output of :func:`build_scores`.
        top_k: Size of the naive shortlist.

    Returns:
        True if the top of the naive ranking is significantly enriched for at least
        one liability, meaning the "reversal is not enough" claim is supported.
    """
    in_top = frame["in_top"]
    n_bg = int((~in_top).sum())
    print(f"\nranked perturbations: {len(frame)}  top-K: {top_k}  background: {n_bg}")

    print(f"\n--- top {min(25, top_k)} by naive suppression ---")
    cols = [
        "rank",
        "target_contrast_gene_name",
        "naive_suppression",
        "is_core_essential",
        "is_iei",
        "stim_de_genes",
        "rest_de_genes",
    ]
    with pd.option_context("display.width", 140, "display.max_columns", 20):
        print(frame.head(min(25, top_k))[cols].to_string(index=False))

    print("\n--- confound check: does statistical power drive the ranking? ---")
    rho, rho_p = stats.spearmanr(frame["naive_suppression"], frame["n_cells_target"])
    print(f"  spearman(naive_suppression, n_cells_target) = {rho:+.3f}  p={rho_p:.3g}")
    print(
        f"  n_cells_target  top median {frame.loc[in_top, 'n_cells_target'].median():9.1f}   "
        f"bg median {frame.loc[~in_top, 'n_cells_target'].median():9.1f}"
    )
    if abs(rho) > 0.2:
        print("  WARNING: cell count correlates with the score. The matched test below is the")
        print("           one that counts; the unmatched test may be inflated by power alone.")

    def run_tests(bg_mask: pd.Series, label: str) -> list[bool]:
        """Run every enrichment test of the top K against a given background."""
        out: list[bool] = []
        n_background = int(bg_mask.sum())
        print(f"\n--- {label} (background n={n_background}) ---")
        print("  binary liabilities, Fisher one-sided greater:")
        for name, col in (("core-essential", "is_core_essential"), ("IEI (immunodeficiency)", "is_iei")):
            subset = frame[in_top | bg_mask]
            odds, pval, n_top, n_back = _fisher(subset[col], subset["in_top"])
            top_rate = n_top / max(int(in_top.sum()), 1)
            bg_rate = n_back / max(n_background, 1)
            flag = pval < 0.05 and odds > 1
            out.append(flag)
            print(
                f"    {name:24s} top {n_top:3d} ({top_rate:6.1%})  bg {n_back:4d} ({bg_rate:6.1%})  "
                f"OR={odds:6.2f}  p={pval:.3g}  {'ENRICHED' if flag else 'ns'}"
            )

        print("  continuous disruption, Mann-Whitney U one-sided greater:")
        for name, col in (
            ("stim DE genes (collateral)", "stim_de_genes"),
            ("rest DE genes (homeostasis)", "rest_de_genes"),
            ("tolerance suppression", "tolerance_suppression"),
        ):
            top_vals = frame.loc[in_top, col].dropna()
            bg_vals = frame.loc[bg_mask, col].dropna()
            if len(top_vals) < 5 or len(bg_vals) < 5:
                print(f"    {name:28s} insufficient data")
                continue
            _, pval = stats.mannwhitneyu(top_vals, bg_vals, alternative="greater")
            flag = pval < 0.05
            out.append(flag)
            print(
                f"    {name:28s} top median {top_vals.median():9.2f}  "
                f"bg median {bg_vals.median():9.2f}  p={pval:.3g}  {'HIGHER' if flag else 'ns'}"
            )
        return out

    run_tests(~in_top, "unmatched background")
    matched_verdicts = run_tests(frame["matched_background"], "cell-count-matched background")

    # The matched background is the test that decides it. Power cannot explain it away.
    holds = any(matched_verdicts)
    print("\n" + "=" * 72)
    if holds:
        print("VERDICT: the empirical bet HOLDS (on the cell-count-matched background).")
        print("The naive reversal ranking is enriched for liabilities a safety gate rejects.")
        print("'Reversal is not enough' is supported by the data.")
    else:
        print("VERDICT: the empirical bet FAILS on the matched background.")
        print("Power, not toxicity, explains any unmatched enrichment. The headline needs a rewrite.")
    print("=" * 72)
    return holds


def main() -> None:
    """Run the risk-kill test and persist the annotated ranking."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top-k", type=int, default=100)
    args = parser.parse_args()

    paths.ensure_dirs()
    frame = build_scores(args.top_k)
    holds = report(frame, args.top_k)

    out = paths.TABLES / "risk_kill_naive_reversal.csv"
    frame.to_csv(out, index=False)
    print(f"\nwrote {out}")

    raise SystemExit(0 if holds else 1)


if __name__ == "__main__":
    main()
