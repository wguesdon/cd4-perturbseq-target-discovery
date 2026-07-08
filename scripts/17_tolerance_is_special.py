"""N9. Is the co-inhibitory module preferentially reduced, or merely co-induced?

Pre-registered in ``docs/preregistration_n9_2026_07_08.md``, committed before this file was written.
Read it first, including **ADDENDUM 1**. Nothing here may be tuned against its own output.

``scripts/10_tolerance_is_real.py`` shows the module falls further than 200 random nine-gene modules
**matched on baseline expression**. That null cannot answer the obvious objection. The nine module
genes are stimulation-induced, and so are the 32 effector genes the ranking is built on. Random genes
at the same expression level are mostly not stimulation-induced, so they do not fall when TCR
signalling is blocked. Any co-induced gene set would reproduce the result. So the null is redrawn,
matched on baseline expression **and** on how strongly the gene is switched on by stimulation, taken
from an external non-perturbational experiment (Arce 2024, AAVS1 control-guide Teff, Stim vs Rest).

**The first run failed its own negative control (25% false positives at nominal 5%), and the control
was right.** A Mann-Whitney across 6,371 perturbations treats them as independent. They are not: the
real module is ONE fixed nine-gene set, so if it covaries with effector suppression every top-100 z
shifts together. The effective unit of resampling is the MODULE, not the perturbation. The MWU is
retained below, computed and printed, and **labelled void**, so the size of the error is visible.
The decisive test is a module-level permutation: draw many random induction-matched nine-gene modules,
compute the same statistic for each, and rank the real one among them. That p is exact by construction.

Both outcomes were written down in advance. The practical finding does not depend on which fires: of
the naive top 20, six of the nine knockdowns Schmidt & Steinhart 2022 independently showed reduce IL-2
are refused on this axis alone.

The module is named for what it is. In a 48 h stimulation of bulk primary human CD4 *conventional*
T cells, CTLA4/PDCD1/LAG3/TIGIT are activation-induced co-inhibitory receptors, and FOXP3 is
transiently induced without conferring suppressive function (Wang 2007; Tran, Ramsey & Shevach 2007;
Allan 2007). This is mRNA. It is not tolerance.

Usage:
    uv run python scripts/17_tolerance_is_special.py [--n-null 200] [--n-modules 200]
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy import stats

from cd4_perturbseq import de_stats, paths, priors, programs

STIM = "Stim48hr"
TOP_N = 100
"""Matches ``scripts/10``'s exposure exactly, so the two nulls are comparable. Not tuned."""

ALPHA = 0.05
N_BINS = 10
BASELINE_ROWS = 200
MAX_FALLBACKS = 1
"""Pre-registered: if more than one of the nine positions cannot be matched, the test is void."""

MWU_FPR_MAX = 0.10
"""ADDENDUM 1: the coded gate was 0.25, looser than the registered "approximately 5%". Tightened so
the code and the pre-registration agree. The per-perturbation MWU is expected to FAIL this."""

N_POSITIVE = 50
NULLS_POSITIVE = 100
"""The positive control is a sanity check, so it runs at lower resolution than the primary."""

IL2_AXIS = ("IL2", "IL2RA", "IL2RB", "IL2RG", "JAK1", "JAK3", "STAT5A", "STAT5B")
"""Fixed a priori. Exploratory only; no p-value from this set is promoted to a claim."""

SEED = 0


def _decile(values: np.ndarray, n_bins: int = N_BINS) -> np.ndarray:
    """Quantile-bin, with -1 wherever the value is not finite.

    Args:
        values: Per-gene quantity to bin.
        n_bins: Number of quantile bins.

    Returns:
        Integer bin label per gene; -1 where ``values`` is not finite.
    """
    out = np.full(values.size, -1, dtype=np.int64)
    finite = np.isfinite(values)
    if finite.sum() >= n_bins:
        out[finite] = pd.qcut(values[finite], n_bins, labels=False, duplicates="drop")
    return out


def _pools(
    exp_dec: np.ndarray, ind_dec: np.ndarray, exclude: set[int]
) -> dict[tuple[int, int], np.ndarray]:
    """Build the joint (expression decile x induction decile) candidate pools.

    Args:
        exp_dec: Per-gene baseline-expression decile.
        ind_dec: Per-gene stimulation-induction decile.
        exclude: Gene indices never eligible.

    Returns:
        Mapping from (expression decile, induction decile) to eligible gene indices.
    """
    eligible = (exp_dec >= 0) & (ind_dec >= 0)
    if exclude:
        eligible[np.fromiter(exclude, dtype=np.int64, count=len(exclude))] = False
    pools: dict[tuple[int, int], list[int]] = {}
    for g in np.flatnonzero(eligible):
        pools.setdefault((int(exp_dec[g]), int(ind_dec[g])), []).append(int(g))
    return {k: np.array(v, dtype=np.int64) for k, v in pools.items()}


def _draw_matched(
    module_idx: np.ndarray,
    exp_dec: np.ndarray,
    ind_dec: np.ndarray,
    pools: dict[tuple[int, int], np.ndarray],
    n_sets: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, int, list[str]]:
    """Draw ``n_sets`` gene sets matched to ``module_idx`` position by position.

    Each position is matched on the joint (expression decile, induction decile) cell of the
    corresponding real module gene. Fallback ladder: same induction decile with expression decile
    within +/-1, then same induction decile with any expression decile. Induction is the variable
    under test, so it is the last thing relaxed.

    Args:
        module_idx: Positional gene indices of the real module.
        exp_dec: Per-gene baseline-expression decile.
        ind_dec: Per-gene stimulation-induction decile.
        pools: Output of :func:`_pools`.
        n_sets: Number of null sets to draw.
        rng: Random generator.

    Returns:
        Tuple of (sets of shape ``(n_sets, len(module_idx))``, number of fallback positions, notes).
    """
    candidates: list[np.ndarray] = []
    fallbacks = 0
    notes: list[str] = []

    for g in module_idx:
        e, i = int(exp_dec[g]), int(ind_dec[g])
        pool = pools.get((e, i), np.array([], dtype=np.int64))
        if pool.size == 0:
            near = [p for p in (pools.get((e - 1, i)), pools.get((e + 1, i))) if p is not None and p.size]
            if near:
                pool = np.concatenate(near)
                fallbacks += 1
                notes.append(f"gene idx {int(g)}: cell ({e},{i}) empty, used expression decile {e}+/-1")
        if pool.size == 0:
            same_ind = [v for (_, ii), v in pools.items() if ii == i]
            if same_ind:
                pool = np.concatenate(same_ind)
                fallbacks += 1
                notes.append(f"gene idx {int(g)}: relaxed to induction decile {i}, any expression")
        if pool.size == 0:
            raise RuntimeError(f"no eligible null gene for module gene {int(g)} (cell {e},{i})")
        candidates.append(pool)

    sets = np.empty((n_sets, module_idx.size), dtype=np.int64)
    for s in range(n_sets):
        for j, pool in enumerate(candidates):
            sets[s, j] = rng.choice(pool)
    return sets, fallbacks, notes


def _suppression(lfc: np.ndarray, sets: np.ndarray, self_idx: np.ndarray) -> np.ndarray:
    """Mean module suppression per perturbation, for many modules at once.

    Each perturbation's own gene is masked wherever it appears in a module, because its on-target
    knockdown is large, negative, and required to be significant by QC.

    Args:
        lfc: Log fold changes, shape ``(n_pert, n_genes)``, gene axis in ``var`` order.
        sets: Gene indices, shape ``(n_sets, module_size)``.
        self_idx: Per-perturbation gene index of its own target, or -1.

    Returns:
        ``-mean(log_fc)`` over each module, shape ``(n_pert, n_sets)``. Larger means more suppressed.
    """
    n_pert = lfc.shape[0]
    n_sets, m = sets.shape
    block = lfc[:, sets.ravel()].astype(np.float64).reshape(n_pert, n_sets, m)
    block[sets[None, :, :] == self_idx[:, None, None]] = np.nan
    with np.errstate(invalid="ignore"):
        return -np.nanmean(block, axis=2)


def _z(real: np.ndarray, nulls: np.ndarray) -> np.ndarray:
    """Standardise a real module against its own nulls, within perturbation.

    Args:
        real: Real-module suppression, shape ``(n_pert,)``.
        nulls: Null-module suppression, shape ``(n_pert, n_sets)``.

    Returns:
        z per perturbation. The null mean and sd are computed per perturbation, which conditions on
        that perturbation's own effect size without fitting any regression.
    """
    mu = nulls.mean(axis=1)
    sd = nulls.std(axis=1, ddof=1)
    return (real - mu) / np.maximum(sd, 1e-9)


def _stat(z: np.ndarray, in_top: np.ndarray) -> float:
    """The module-level statistic: median z among the naive top.

    Args:
        z: Per-perturbation z.
        in_top: Boolean exposure.

    Returns:
        Median of the finite z values inside the exposure.
    """
    v = z[in_top]
    return float(np.median(v[np.isfinite(v)]))


def _mwu(z: np.ndarray, in_top: np.ndarray) -> tuple[float, float, float]:
    """One-sided Mann-Whitney U, top versus background. VOID as an inference; see ADDENDUM 1.

    Args:
        z: Per-perturbation statistic.
        in_top: Boolean exposure.

    Returns:
        Tuple of (median top, median background, one-sided p).
    """
    a, b = z[in_top], z[~in_top]
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    _, p = stats.mannwhitneyu(a, b, alternative="greater")
    return float(np.median(a)), float(np.median(b)), float(p)


def _module_null(
    module_idx: np.ndarray,
    exp_dec: np.ndarray,
    ind_dec: np.ndarray,
    base_exclude: set[int],
    lfc: np.ndarray,
    self_idx: np.ndarray,
    in_top: np.ndarray,
    n_modules: int,
    n_null: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Null distribution of the module-level statistic.

    Draws ``n_modules`` random modules matched to ``module_idx`` on the joint decile grid. Each is
    treated as though it were the real module: it gets its own ``n_null`` matched null sets, from
    which its z and then its statistic are computed. Each drawn module is excluded from its own nulls.

    Args:
        module_idx: The real module, used only for its decile cells.
        exp_dec: Per-gene expression decile.
        ind_dec: Per-gene induction decile.
        base_exclude: Genes never eligible (both real modules).
        lfc: Log fold changes.
        self_idx: Per-perturbation own-gene index.
        in_top: Boolean exposure.
        n_modules: Size of the null distribution.
        n_null: Null sets per drawn module.
        rng: Random generator.

    Returns:
        Tuple of (statistics of shape ``(n_modules,)``, the drawn modules of shape
        ``(n_modules, module_size)``).
    """
    pools = _pools(exp_dec, ind_dec, base_exclude)
    fakes, _, _ = _draw_matched(module_idx, exp_dec, ind_dec, pools, n_modules, rng)
    fake_real = _suppression(lfc, fakes, self_idx)

    stats_out = np.empty(n_modules, dtype=np.float64)
    for k in range(n_modules):
        pools_k = _pools(exp_dec, ind_dec, base_exclude | set(fakes[k].tolist()))
        nulls_k, _, _ = _draw_matched(fakes[k], exp_dec, ind_dec, pools_k, n_null, rng)
        stats_out[k] = _stat(_z(fake_real[:, k], _suppression(lfc, nulls_k, self_idx)), in_top)
    return stats_out, fakes


def _perm_p(observed: float, null: np.ndarray) -> float:
    """Exact permutation p-value, one-sided greater, with the +1 correction.

    Args:
        observed: The real module's statistic.
        null: Statistics of the null modules.

    Returns:
        ``(1 + #{null >= observed}) / (1 + n)``.
    """
    return float((1 + int((null >= observed).sum())) / (1 + null.size))


def _rho_bar(lfc: np.ndarray, sets: np.ndarray) -> np.ndarray:
    """Mean pairwise correlation among a module's genes, across perturbations.

    The real module is a co-regulated biological program. The null modules are arbitrary gene sets
    matched on expression and induction. If co-regulation alone inflated the statistic, this is where
    it would show. Median-based statistics are robust to variance inflation, but the objection is
    worth measuring rather than asserting away.

    Args:
        lfc: Log fold changes, shape ``(n_pert, n_genes)``.
        sets: Gene indices, shape ``(n_sets, module_size)``.

    Returns:
        Mean off-diagonal Pearson correlation per set, shape ``(n_sets,)``.
    """
    out = np.empty(sets.shape[0], dtype=np.float64)
    iu = np.triu_indices(sets.shape[1], 1)
    for k, s in enumerate(sets):
        block = lfc[:, s].astype(np.float64)
        block = block[np.isfinite(block).all(axis=1)]
        with np.errstate(invalid="ignore", divide="ignore"):
            c = np.corrcoef(block.T)
        out[k] = float(np.nanmean(c[iu]))
    return out


def main() -> None:
    """Run the corrected primary, its controls, and the pre-registered secondaries."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-null", type=int, default=200)
    parser.add_argument("--n-modules", type=int, default=200)
    args = parser.parse_args()
    paths.ensure_dirs()
    rng = np.random.default_rng(SEED)

    obs = de_stats.read_obs().reset_index(drop=True)
    var = de_stats.read_var()
    names = var["gene_name"].astype(str).to_numpy()
    idx_of = {g: i for i, g in enumerate(names)}

    program = programs.load_activation_program()
    effector = [g for g in programs.effector_core(program) if g in idx_of]
    tolerance = [g for g in programs.tolerance_module(program) if g in idx_of]
    eff_idx = np.array([idx_of[g] for g in effector])
    tol_idx = np.array([idx_of[g] for g in tolerance])

    gene_col = "target_contrast_gene_name"
    obs[gene_col] = obs[gene_col].astype(str)
    rows = obs.index[(obs["culture_condition"] == STIM) & de_stats.qc_mask(obs)].to_numpy()
    sub = obs.loc[rows].reset_index(drop=True)
    print(f"effector {len(effector)} genes, co-inhibitory module {len(tolerance)} genes")
    print(f"QC-passing {STIM} perturbations: {rows.size}")

    base_block = de_stats.read_layer_rows(rows[:BASELINE_ROWS], layer="baseMean")
    baseline = np.nanmedian(base_block.astype(np.float64), axis=0)
    del base_block

    arce = priors.arce_stim_vs_rest().drop_duplicates("gene_name").set_index("gene_name")
    induction = np.array([arce["log2FoldChange"].get(g, np.nan) for g in names], dtype=np.float64)
    print(f"induction available for {int(np.isfinite(induction).sum())} of {names.size} measured genes")

    exp_dec, ind_dec = _decile(baseline), _decile(induction)
    both = set(eff_idx.tolist()) | set(tol_idx.tolist())
    pools = _pools(exp_dec, ind_dec, both)
    print(f"joint (expression x induction) cells with >=1 eligible gene: {len(pools)}")

    print(f"reading log_fc, {rows.size} perturbations x {names.size} genes ...")
    lfc = de_stats.read_layer_rows(rows, layer="log_fc")
    self_idx = np.array([idx_of.get(g, -1) for g in sub[gene_col]], dtype=np.int64)
    print(f"perturbations targeting one of their own module genes (masked): "
          f"{int(np.isin(self_idx, tol_idx).sum())}")

    efficacy = _suppression(lfc, eff_idx[None, :], self_idx)[:, 0]
    in_top = pd.Series(-efficacy).rank(ascending=True, method="first").to_numpy() <= TOP_N
    assert in_top.sum() == TOP_N, "exposure is not the top-100 by effector suppression"
    # This project has written the same rank-direction bug twice, both times flattering its own
    # result. The assertion is cheap and it fires before anything is reported.
    assert np.nanmin(efficacy[in_top]) >= np.nanmax(efficacy[~in_top]) - 1e-9, "rank direction inverted"

    # ------------------------------------------------------------------ the real module
    null_ind, fb, notes = _draw_matched(tol_idx, exp_dec, ind_dec, pools, args.n_null, rng)
    for n in notes:
        print(f"  FALLBACK {n}")
    if fb > MAX_FALLBACKS:
        print(f"\nVOID: {fb} of {tol_idx.size} positions needed a fallback (max {MAX_FALLBACKS}).")
        raise SystemExit(1)

    module_supp = _suppression(lfc, tol_idx[None, :], self_idx)[:, 0]
    z_ind = _z(module_supp, _suppression(lfc, null_ind, self_idx))
    T_real = _stat(z_ind, in_top)

    # ------------------------------------------------------------------ the VOID per-perturbation MWU
    mwu_top, mwu_bg, p_mwu = _mwu(z_ind, in_top)

    # ------------------------------------------------------------------ CORRECTED PRIMARY
    print(f"\nbuilding the module-level null: {args.n_modules} random induction-matched modules, "
          f"{args.n_null} nulls each ...")
    T_null, fakes = _module_null(tol_idx, exp_dec, ind_dec, both, lfc, self_idx, in_top,
                                 args.n_modules, args.n_null, rng)
    p_module = _perm_p(T_real, T_null)

    print("\n" + "=" * 78)
    print("PRIMARY, corrected per ADDENDUM 1 (module is the unit of resampling)")
    print(f"  statistic T = median z among the naive top-{TOP_N}")
    print(f"  real co-inhibitory module         T = {T_real:+.3f}")
    print(f"  {args.n_modules} random induction-matched modules  "
          f"T: median {np.median(T_null):+.3f}, 95th pct {np.quantile(T_null, 0.95):+.3f}, "
          f"max {T_null.max():+.3f}")
    print(f"  module-level permutation p = {p_module:.4g}   (alpha {ALPHA})")
    print(f"\n  [VOID, shown so the size of the error is visible] per-perturbation MWU")
    print(f"    z top-{TOP_N} median {mwu_top:+.3f} vs background {mwu_bg:+.3f}, p = {p_mwu:.3g}")
    print(f"    This p treats {rows.size:,} correlated z values as independent. See ADDENDUM 1.")

    # ------------------------------------------------------------------ CONTROLS
    ind_perm = ind_dec.copy()
    movable = ind_perm >= 0
    ind_perm[movable] = rng.permutation(ind_perm[movable])
    pools_perm = _pools(exp_dec, ind_perm, both)
    null_perm, _, _ = _draw_matched(tol_idx, exp_dec, ind_perm, pools_perm, args.n_null, rng)
    z_perm = _z(module_supp, _suppression(lfc, null_perm, self_idx))
    T_perm = _stat(z_perm, in_top)

    eff_supp = efficacy
    null_eff, _, _ = _draw_matched(eff_idx, exp_dec, ind_dec, pools, args.n_null, rng)
    z_eff = _z(eff_supp, _suppression(lfc, null_eff, self_idx))
    T_eff = _stat(z_eff, in_top)
    T_eff_null, _ = _module_null(eff_idx, exp_dec, ind_dec, both, lfc, self_idx, in_top,
                                 N_POSITIVE, NULLS_POSITIVE, rng)
    p_eff = _perm_p(T_eff, T_eff_null)

    # Calibration: leave-one-out among the null modules. Uniform by construction; verify it.
    loo_p = np.array([
        (1 + int((np.delete(T_null, k) >= T_null[k]).sum())) / T_null.size for k in range(T_null.size)
    ])
    loo_fpr = float((loo_p < ALPHA).mean())
    mwu_fpr_hits = 0
    for k in range(min(20, args.n_modules)):
        pools_k = _pools(exp_dec, ind_dec, both | set(fakes[k].tolist()))
        nulls_k, _, _ = _draw_matched(fakes[k], exp_dec, ind_dec, pools_k, args.n_null, rng)
        zk = _z(_suppression(lfc, fakes[k][None, :], self_idx)[:, 0],
                _suppression(lfc, nulls_k, self_idx))
        mwu_fpr_hits += int(_mwu(zk, in_top)[2] < ALPHA)
    mwu_fpr = mwu_fpr_hits / min(20, args.n_modules)

    print("\nFALSIFICATION CONTROLS")
    print(f"  [1] recovery: permute induction -> the grid degenerates to expression-matched")
    print(f"      T = {T_perm:+.3f}   (expect ~ +3.9, reproducing scripts/10)")
    print(f"      => co-induction accounts for {(T_perm - T_real) / T_perm:.0%} of scripts/10's effect")
    print(f"  [2] positive: the effector module, against its own matched nulls")
    print(f"      T = {T_eff:+.3f}   module-level p = {p_eff:.4g} (vs {N_POSITIVE} random 32-gene modules)")
    print(f"  [3a] negative, MODULE-LEVEL: leave-one-out among the {args.n_modules} null modules")
    print(f"      false-positive rate at alpha={ALPHA}: {loo_fpr:.1%}   (expect ~5%, exact by construction)")
    print(f"  [3b] negative, PER-PERTURBATION MWU on the first {min(20, args.n_modules)} null modules")
    print(f"      false-positive rate at alpha={ALPHA}: {mwu_fpr:.1%}   "
          f"(nominal 5%; gate {MWU_FPR_MAX:.0%}) -> the MWU is {'OK' if mwu_fpr <= MWU_FPR_MAX else 'VOID'}")

    ctrl_recovery = T_perm > 2.0
    ctrl_positive = T_eff > 2.0 and p_eff < ALPHA
    ctrl_negative = 0.0 <= loo_fpr <= 0.15
    ok = ctrl_recovery and ctrl_positive and ctrl_negative
    print(f"\n  recovery {'PASS' if ctrl_recovery else 'FAIL'}   "
          f"positive {'PASS' if ctrl_positive else 'FAIL'}   "
          f"negative(module-level) {'PASS' if ctrl_negative else 'FAIL'}")

    # ------------------------------------------------------------------ POST-HOC diagnostic
    # NOT pre-registered. Added after the primary was run, because a reviewer will raise it and it is
    # better found here than there. It characterises a limitation. It CANNOT change the verdict: the
    # primary was registered in advance and this was not. If it contradicts the primary, BOTH are
    # reported and the contradiction is the finding.
    rho_real = float(_rho_bar(lfc, tol_idx[None, :])[0])
    rho_null = _rho_bar(lfc, fakes)
    rs = stats.spearmanr(rho_null, T_null)
    hi = rho_null >= np.quantile(rho_null, 0.75)
    p_cond = _perm_p(T_real, T_null[hi])

    print("\nPOST-HOC DIAGNOSTIC (not pre-registered; cannot change the verdict)")
    print("  The real module is a co-regulated program; the nulls are arbitrary matched gene sets.")
    print("  Does within-module co-regulation, rather than targeting, produce the statistic?")
    print(f"    mean pairwise r among the 9 real genes' log_fc : {rho_real:+.3f}")
    print(f"    same, across the {args.n_modules} null modules  : median {np.median(rho_null):+.3f}, "
          f"95th pct {np.quantile(rho_null, 0.95):+.3f}, max {rho_null.max():+.3f}")
    print(f"    Spearman(rho_bar, T) across null modules       : {rs.statistic:+.3f}  p={rs.pvalue:.3g}")
    print(f"    -> co-regulation {'DOES' if rs.pvalue < 0.05 and rs.statistic > 0.2 else 'does not'}"
          f" predict the statistic among the nulls")
    print(f"    p restricted to the most co-regulated quartile of nulls (n={int(hi.sum())}): {p_cond:.4g}")

    fit = stats.linregress(rho_null, T_null)
    pred = float(fit.intercept + fit.slope * rho_real)
    print(f"    OLS  T_null ~ rho_bar: slope {fit.slope:+.2f}, R^2 {fit.rvalue ** 2:.3f}, p {fit.pvalue:.3g}")
    print(f"    T predicted at the real module's rho_bar ({rho_real:.3f}) : {pred:+.3f}")
    print(f"    T observed                                        : {T_real:+.3f}")
    print(f"    excess over the co-regulation trend               : {T_real - pred:+.3f} "
          f"({1 - pred / T_real:.0%} of the observed statistic is unexplained by it)")
    print(f"    CAVEAT: rho_bar_real ({rho_real:.3f}) exceeds EVERY null ({rho_null.max():.3f}), so this")
    print("    is an extrapolation beyond the support of the null. It is reported as one. Co-regulation")
    print("    cannot be matched on this screen; it is a limitation, not a solved confound.")

    # ------------------------------------------------------------------ SECONDARY, reported only
    print("\nSECONDARY (reported, never decisive)")
    print(f"  Per-gene suppression, naive top-{TOP_N} vs background. Co-induction predicts the two")
    print("  NON-induced genes (TIGIT, TGFB1) show no excess. n=2: corroborates, never overrides.")
    print("  These p-values carry the same non-independence caveat as the void MWU. Descriptive only.")
    per_gene = []
    for j, g in enumerate(tolerance):
        gi = tol_idx[j]
        v = -lfc[:, gi].astype(np.float64)
        v[self_idx == gi] = np.nan
        a, b = v[in_top], v[~in_top]
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        _, p = stats.mannwhitneyu(a, b, alternative="greater")
        per_gene.append({
            "gene_name": g, "induction_log2fc": float(induction[gi]),
            "induced": bool(induction[gi] > 1.0),
            "top_median": float(np.median(a)), "bg_median": float(np.median(b)), "mwu_p": float(p),
        })
    pg = pd.DataFrame(per_gene).sort_values("induction_log2fc", ascending=False)
    print(pg.to_string(index=False))
    non_induced = pg[~pg["induced"]]
    print(f"\n  non-induced genes suppressed in the top-{TOP_N} at p<{ALPHA}: "
          f"{int((non_induced['mwu_p'] < ALPHA).sum())} of {len(non_induced)}   "
          f"(co-induction predicts 0; n=2 cannot settle it)")
    induced = pg[pg["induced"]]
    print(f"  induced genes NOT suppressed at p<{ALPHA}: "
          f"{int((induced['mwu_p'] >= ALPHA).sum())} of {len(induced)}")

    il2 = [g for g in IL2_AXIS if g in set(sub[gene_col])]
    is_il2 = sub[gene_col].isin(il2).to_numpy()
    n_il2_top = int((is_il2 & in_top).sum())
    print(f"\n  IL-2 axis perturbed and QC-passing: {len(il2)} of {len(IL2_AXIS)} ({', '.join(il2)})")
    print(f"  of which in the naive top-{TOP_N}: {n_il2_top}")
    print("  too few to compare. Reported as UNAVAILABLE, not as null. The exploratory mechanism"
          if n_il2_top < 3 else "  n is small; exploratory only")
    print("  in section 8 of the pre-registration cannot be assessed on this screen.")

    # ------------------------------------------------------------------ persist and rule
    pd.DataFrame({
        "gene_name": sub[gene_col],
        "efficacy": efficacy,
        "module_suppression": module_supp,
        "z_induction_matched": z_ind,
        "z_expression_matched_recovery_control": z_perm,
        "z_effector_positive_control": z_eff,
        "in_naive_top100": in_top,
        "is_il2_axis": is_il2,
    }).to_csv(paths.TABLES / "tolerance_induction_matched.csv", index=False)
    pg.to_csv(paths.TABLES / "tolerance_per_gene_induction.csv", index=False)
    pd.DataFrame({"module_index": np.arange(T_null.size), "T_null": T_null}).to_csv(
        paths.TABLES / "tolerance_module_level_null.csv", index=False)
    pd.DataFrame([{
        "T_real": T_real, "p_module_level": p_module, "n_modules": args.n_modules,
        "n_null_per_module": args.n_null, "T_null_median": float(np.median(T_null)),
        "T_null_p95": float(np.quantile(T_null, 0.95)), "T_null_max": float(T_null.max()),
        "T_recovery_expression_matched": T_perm, "T_effector_positive_control": T_eff,
        "p_effector_positive_control": p_eff, "loo_fpr_module_level": loo_fpr,
        "mwu_fpr_per_perturbation": mwu_fpr, "mwu_p_void": p_mwu,
        "coinduction_share_of_scripts10_effect": (T_perm - T_real) / T_perm,
        "rho_bar_real": rho_real, "rho_bar_null_median": float(np.median(rho_null)),
        "rho_bar_null_max": float(rho_null.max()),
        "spearman_rhobar_vs_T_null": float(rs.statistic), "spearman_rhobar_vs_T_null_p": float(rs.pvalue),
        "p_module_level_coregulated_quartile": p_cond,
    }]).to_csv(paths.TABLES / "tolerance_induction_verdict.csv", index=False)
    print(f"\nwrote {paths.TABLES / 'tolerance_induction_verdict.csv'}")

    if not ok:
        print("\n" + "=" * 78)
        print("VOID: a falsification control did not fire. The primary is uninterpretable.")
        print("Do not report the primary. Fix the machinery, not the conclusion.")
        print("=" * 78)
        raise SystemExit(1)

    passed = p_module < ALPHA and T_real > 0
    print("\n" + "=" * 78)
    if passed:
        print("VERDICT: PASS. The co-inhibitory module is reduced in the naive top beyond what")
        print(f"expression- AND induction-matched modules predict (p = {p_module:.4g}, module-level).")
        print(f"Co-induction accounts for {(T_perm - T_real) / T_perm:.0%} of the effect scripts/10")
        print("reported, and it is not all of it. The report may say 'beyond co-induction'.")
        print("It may NOT say 'specifically destroys tolerance', nor any causal or clinical verb.")
        print("This is mRNA in a 48h in-vitro stimulation of bulk primary human CD4 Tconv.")
    else:
        print("VERDICT: FAIL. The effect is explained by co-induction. A suppression objective")
        print("cannot separate the effector program from the co-inhibitory program, because the")
        print("two share an upstream cascade. This is a RESULT, not a caveat: no ranking that")
        print("maximises effector suppression can avoid it, so an explicit gate is required.")
        print("Delete the word 'specifically' everywhere. The practical finding is untouched:")
        print("six of nine Schmidt-validated IL-2 reducers in the naive top 20 are refused on")
        print("this axis alone.")
    print("=" * 78)
    print("\nThe largest limitation is not statistical. See section 2 of the pre-registration:")
    print("this module is named by interpretation, and no functional assay exists in this screen.")
    raise SystemExit(0)


if __name__ == "__main__":
    main()
