"""Is the tolerance axis real, or is it just re-encoding suppression strength?

Six of the nine knockdowns the gate now refuses are refused on the tolerance axis alone. So the
gate stands or falls on this question.

The worry, raised by the adversarial audit: a knockdown that broadly downregulates the
transcriptome will show a low mean log fold change over ANY gene set, including the tolerance
module. If so, "collapses tolerance" means nothing more than "suppresses strongly", the axis
re-encodes the efficacy score, and the gate rejects the very hits it is meant to find.

Two tests, in increasing severity.

**Test 1, the residual.** Regress tolerance suppression on effector suppression. Is the residual
still elevated in the perturbations we reject? A raw correlation proves nothing; the residual is
the claim.

**Test 2, the random-module null.** This is the one that matters. Sample many random gene sets of
the same size as the tolerance module, matched on baseline expression. Compute each perturbation's
mean log fold change over each random set. If the real tolerance module is suppressed no more than
a random module of the same size and expression, then the axis is measuring global downregulation
and nothing else, and it must be removed from the gate.

Usage:
    uv run python scripts/10_tolerance_is_real.py --n-null 200
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy import stats

from cd4_perturbseq import de_stats, paths, programs

STIM = "Stim48hr"
FDR = 0.10
RNG = np.random.default_rng(0)


qc_mask = de_stats.qc_mask
"""Routine perturbation QC. Single definition in :mod:`cd4_perturbseq.de_stats`."""


def _expression_matched_sets(
    baseline: np.ndarray, module_idx: np.ndarray, n_sets: int, exclude: set[int]
) -> np.ndarray:
    """Draw random gene sets matched to the module on baseline-expression decile.

    Matching on expression matters: the tolerance genes are cytokines and checkpoints, which are
    lowly expressed. An unmatched random set would be dominated by housekeeping genes whose log
    fold changes are far less variable, and the module would look special for the wrong reason.

    Args:
        baseline: Per-gene baseline expression, used only for decile matching.
        module_idx: Positional indices of the real module's genes.
        n_sets: How many null sets to draw.
        exclude: Gene indices never to draw (the real module, and the effector module).

    Returns:
        Array of shape ``(n_sets, len(module_idx))`` of gene indices.
    """
    finite = np.isfinite(baseline)
    deciles = np.full(baseline.size, -1, dtype=np.int64)
    deciles[finite] = pd.qcut(baseline[finite], 10, labels=False, duplicates="drop")

    pools: dict[int, np.ndarray] = {}
    for d in range(deciles.max() + 1):
        candidates = np.flatnonzero((deciles == d))
        pools[d] = np.array([c for c in candidates if c not in exclude])

    wanted = [deciles[i] for i in module_idx]
    sets = np.empty((n_sets, module_idx.size), dtype=np.int64)
    for s in range(n_sets):
        for j, d in enumerate(wanted):
            pool = pools.get(int(d), np.array([], dtype=np.int64))
            if pool.size == 0:  # decile empty after exclusion; fall back to any gene
                pool = np.array([c for c in range(baseline.size) if c not in exclude])
            sets[s, j] = RNG.choice(pool)
    return sets


def main() -> None:
    """Run both tests and print a verdict on whether tolerance may stay in the gate."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-null", type=int, default=200)
    args = parser.parse_args()
    paths.ensure_dirs()

    obs = de_stats.read_obs().reset_index(drop=True)
    var = de_stats.read_var()
    names = var["gene_name"].astype(str).to_numpy()
    idx_of = {g: i for i, g in enumerate(names)}

    program = programs.load_activation_program()
    effector = [g for g in programs.effector_core(program) if g in idx_of]
    tolerance = [g for g in programs.tolerance_module(program) if g in idx_of]
    eff_idx = np.array([idx_of[g] for g in effector])
    tol_idx = np.array([idx_of[g] for g in tolerance])
    print(f"effector module {len(effector)} genes, tolerance module {len(tolerance)} genes")

    rows = obs.index[(obs["culture_condition"] == STIM) & qc_mask(obs)].to_numpy()
    sub = obs.loc[rows].reset_index(drop=True)
    print(f"QC-passing {STIM} perturbations: {rows.size}")

    # Baseline expression per gene, for decile matching. Median across perturbations is fine.
    base_block = de_stats.read_layer_columns(np.arange(len(names)), layer="baseMean")[rows[:200]]
    baseline = np.nanmedian(base_block, axis=0)
    del base_block

    lfc_eff = de_stats.read_layer_columns(eff_idx, layer="log_fc")[rows]
    lfc_tol = de_stats.read_layer_columns(tol_idx, layer="log_fc")[rows]

    efficacy = -np.nanmean(lfc_eff, axis=1)
    tolerance_loss = -np.nanmean(lfc_tol, axis=1)

    # ------------------------------------------------------------------ test 1: the residual
    rho, p = stats.spearmanr(efficacy, tolerance_loss)
    print(f"\n[1] Spearman(efficacy, tolerance_loss) = {rho:+.3f}  p={p:.3g}")
    print("    A high correlation is expected and is NOT itself the objection.")

    slope, intercept, *_ = stats.linregress(efficacy, tolerance_loss)
    residual = tolerance_loss - (slope * efficacy + intercept)

    naive_rank = pd.Series(-np.nanmean(
        np.where(np.isnan(lfc_eff), np.nan, lfc_eff), axis=1)).rank(ascending=False)
    in_top = (naive_rank <= 100).to_numpy()
    _, p_res = stats.mannwhitneyu(residual[in_top], residual[~in_top], alternative="greater")
    print(f"    residual tolerance_loss, top-100 median {np.median(residual[in_top]):+.4f} "
          f"vs background {np.median(residual[~in_top]):+.4f}   MWU p={p_res:.3g}")
    print(f"    -> {'residual survives' if p_res < 0.05 else 'RESIDUAL DOES NOT SURVIVE'}")

    # ------------------------------------------------------------------ test 2: random-module null
    print(f"\n[2] random-module null, {args.n_null} expression-matched sets of {len(tolerance)} genes")
    exclude = set(eff_idx.tolist()) | set(tol_idx.tolist())
    null_sets = _expression_matched_sets(baseline, tol_idx, args.n_null, exclude)

    needed = np.unique(null_sets)
    print(f"    reading {needed.size} distinct null genes ...")
    lfc_null_all = de_stats.read_layer_columns(needed, layer="log_fc")[rows]
    position = {g: k for k, g in enumerate(needed)}

    null_means = np.empty((rows.size, args.n_null), dtype=np.float64)
    for s in range(args.n_null):
        cols = [position[g] for g in null_sets[s]]
        null_means[:, s] = -np.nanmean(lfc_null_all[:, cols], axis=1)
    del lfc_null_all

    null_mu = null_means.mean(axis=1)
    null_sd = null_means.std(axis=1, ddof=1)
    tolerance_z = (tolerance_loss - null_mu) / np.maximum(null_sd, 1e-9)

    print(f"    tolerance_loss median            {np.median(tolerance_loss):+.4f}")
    print(f"    random-module median (per-pert)  {np.median(null_mu):+.4f}")
    print(f"    tolerance z vs its own null: median {np.median(tolerance_z):+.3f}")

    frac_exceeds = float(np.mean(tolerance_z > 1.645))
    print(f"    perturbations where tolerance is suppressed beyond the 95th pct of random: {frac_exceeds:.1%}")

    _, p_top = stats.mannwhitneyu(tolerance_z[in_top], tolerance_z[~in_top], alternative="greater")
    print(f"    tolerance_z, top-100 median {np.median(tolerance_z[in_top]):+.3f} "
          f"vs background {np.median(tolerance_z[~in_top]):+.3f}   MWU p={p_top:.3g}")

    real = p_top < 0.05 and np.median(tolerance_z[in_top]) > 0
    print("\n" + "=" * 76)
    if real:
        print("VERDICT: the tolerance axis is REAL. The tolerance module is suppressed more than")
        print("an expression-matched random module of the same size, and more so in the top 100.")
        print("It may stay in the gate.")
    else:
        print("VERDICT: the tolerance axis is NOT distinguishable from global downregulation.")
        print("It re-encodes suppression strength. REMOVE IT FROM THE GATE: it would reject the")
        print("strongest hits for the crime of being strong.")
    print("=" * 76)

    out = paths.TABLES / "tolerance_null_test.csv"
    pd.DataFrame({
        "gene_name": sub["target_contrast_gene_name"].astype(str),
        "efficacy": efficacy,
        "tolerance_loss": tolerance_loss,
        "tolerance_residual": residual,
        "tolerance_z_vs_random_module": tolerance_z,
        "in_naive_top100": in_top,
    }).to_csv(out, index=False)
    print(f"\nwrote {out}")
    raise SystemExit(0 if real else 1)


if __name__ == "__main__":
    main()
