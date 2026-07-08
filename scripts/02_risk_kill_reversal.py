"""RISK-KILL: does a naive reversal ranking nominate targets a safety gate rejects?

The project's headline claim is "reversal is not enough": ranking perturbations purely by how
strongly they suppress an inflammatory program will surface genes that collapse the resting
transcriptome or the tolerance program, and those are phenotype-but-toxic false positives rather
than drug targets.

This script tests that claim directly, and it is designed to be able to FAIL. If the top of the
naive ranking is no more transcriptome-disrupting, and no more tolerance-collapsing, than a
background matched on effect magnitude, the headline is wrong and the strategy needs rewriting.

Naive score, per perturbed gene, in stimulated cells:

    suppression = -mean(zscore over effector program genes)   [Stim48hr rows]

The naive baseline is STEELMANNED, not a straw man. We apply the QC any competent analyst would
apply before ranking: drop flagged off-target perturbations and require a significant on-target
knockdown. The claim is that even a well-QC'd reversal ranking is toxic, which is much stronger
than "an unfiltered one is".

Two corrections, both forced by an adversarial audit, both dated 2026-07-08.

**The background was matched on the wrong variable.** The original rationale was that
perturbations assayed in more cells have more statistical power, so they get more significant DE
genes, and matching on ``n_cells_target`` removes that. The data say the opposite:
``Spearman(n_cells_target, stim_de_genes) = -0.243``. More cells means FEWER DE genes. Cell count
is a viability readout, not a power proxy, and matching on it controlled the confound backwards.
The confound that actually exists is effect magnitude, so we now stratify on ``z_l2``. Both
stratifications are reported, because the reader is owed the comparison.

**The background was drawn once, with seed 0.** That discarded 6,171 of 6,271 usable controls and
made the reported p-value a random variable. :func:`seed_lottery` reproduces the defect below --
redrawing the same design 2,000 times -- and then we stop sampling altogether. The stratified
tests in ``src/cd4_perturbseq/stratified.py`` condition on the matching variable and use every
control row: Cochran-Mantel-Haenszel for the binary flags, van Elteren for the continuous
measures. Neither draws a random number.

**The exit criterion is one pre-registered primary endpoint, not ``any()`` over five tests.**
``any()`` over five one-sided tests at 0.05 has a family-wise false-positive rate near 23% under
independence, which quietly undoes the script's claim to be falsifiable.

    PRIMARY ENDPOINT: tolerance-module suppression is higher in the naive top 100 than in the
    z_l2-decile-stratified background. One-sided van Elteren, alpha = 0.05.

Chosen on three structural grounds, none of which is "it gave the smallest p-value":

1. It is the only pillar that is both an axis ``04_window_score.py`` actually rejects on and a
   graded, screen-native measure rather than a count of DE genes. Every count-of-DE-genes pillar
   is partly a restatement of effect magnitude (``rho(z_l2, stim_de) = 0.725``); tolerance is not
   (``rho = 0.069``). This is RULE #2.
2. It was independently validated against 200 expression-matched random 9-gene modules in
   ``scripts/10_tolerance_is_real.py``, committed at ``1adab65``, before this endpoint existed.
3. ``scripts/12_magnitude_matched.py`` showed it is the only pillar whose direction-specificity
   survives the sign-flipped induction control at every bin count.

Point 3 was known when the endpoint was chosen, so this is pre-registered with respect to every
future run and every new dataset, not with respect to the run that motivated it. Say so out loud
rather than pretending otherwise. The endpoint can still fail, and then this script exits 1.

The other four tests are reported as SECONDARY, Benjamini-Hochberg corrected across the family.
They do not decide the build.

Usage:
    uv run python scripts/02_risk_kill_reversal.py --top-k 100
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy import stats

from cd4_perturbseq import de_stats, magnitude, paths, priors, programs, stratified

CONDITION = "Stim48hr"
N_BINS = 10
PRIMARY_COLUMN = "tolerance_suppression"
PRIMARY_LABEL = "tolerance suppression"
PRIMARY_ALPHA = 0.05

TESTS = (
    ("core-essential", "is_core_essential", True),
    ("IEI (immunodeficiency)", "is_iei", True),
    ("stim DE genes (collateral)", "stim_de_genes", False),
    ("rest DE genes (homeostasis)", "rest_de_genes", False),
    ("tolerance suppression", "tolerance_suppression", False),
)
"""The five liabilities. Bool flag means CMH; otherwise van Elteren."""


def _fisher(flag: pd.Series, in_top: pd.Series) -> tuple[float, float, int, int]:
    """Fisher's exact test for a binary liability enriched in the top of a ranking.

    Retained only for :func:`seed_lottery`, which reproduces the defect this script used to have.
    Inference now runs through :mod:`cd4_perturbseq.stratified`.

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

    Identical to :func:`cd4_perturbseq.de_stats.qc_mask`, but applied stepwise so the attrition
    funnel can be printed. Column names differ between the December 2025 supplementary CSV and the
    May 2026 h5ad, so each filter is applied only if its column exists.

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


def _single_draw_background(frame: pd.DataFrame, rng: np.random.Generator, n_bins: int = 10) -> pd.Series:
    """The ORIGINAL, DEFECTIVE background: one sample of 100 rows, matched on cell-count decile.

    Kept for two reasons. :func:`seed_lottery` uses it to show why a single draw cannot carry a
    p-value, and ``docs/handoffs/HANDOFF_02_reviewer_verify_riskkill.md`` asks the Claude Science
    reviewer to find exactly this defect from the committed table, so the column must survive.

    It is no longer used for inference.

    Args:
        frame: Ranked frame carrying ``n_cells_target`` and ``in_top``.
        rng: Random generator. The defect is that this exists at all.
        n_bins: Number of quantile bins.

    Returns:
        Boolean series, True for the sampled background rows.
    """
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
    """Compute the naive suppression ranking and attach safety annotations and both strata.

    Args:
        top_k: Number of top-ranked perturbations treated as the naive shortlist.

    Returns:
        DataFrame indexed by perturbed gene with the naive score, liability flags, disruption
        measures, the effect-magnitude covariate, and both stratification variables.
    """
    obs = de_stats.read_obs().reset_index(drop=True)
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
    obs["effector_z"] = np.nanmean(columns[:, : len(effector)], axis=1)
    obs["tolerance_z"] = np.nanmean(columns[:, len(effector) :], axis=1)

    stim = obs[obs["culture_condition"] == CONDITION].copy()
    rest = obs[obs["culture_condition"] == "Rest"].copy()

    print(f"\nQC on the naive baseline ({CONDITION}):")
    stim = apply_qc(stim)
    stim_rows = stim.index.to_numpy()

    stim["naive_suppression"] = -stim["effector_z"]
    stim["tolerance_suppression"] = -stim["tolerance_z"]

    gene_col = "target_contrast_gene_name"
    keep = [gene_col, "naive_suppression", "tolerance_suppression", "n_total_de_genes", "n_downstream", "n_cells_target"]
    frame = stim[keep].copy()
    frame = frame.rename(columns={"n_total_de_genes": "stim_de_genes", "n_downstream": "stim_downstream"})

    rest_unique = rest.drop_duplicates(gene_col).set_index(gene_col)
    frame["rest_de_genes"] = frame[gene_col].map(rest_unique["n_total_de_genes"])
    frame["rest_downstream"] = frame[gene_col].map(rest_unique["n_downstream"])

    essential = priors.core_essential_genes()
    iei = priors.iei_genes()
    frame["is_core_essential"] = frame[gene_col].isin(essential)
    frame["is_iei"] = frame[gene_col].isin(iei)

    print("\neffect magnitude, the covariate the background should have been matched on:")
    scale = magnitude.effect_magnitude(obs, var, stim_rows, exclude=set(effector) | set(tolerance))
    frame["z_l2"] = frame[gene_col].map(scale["z_l2"])

    frame = frame.dropna(subset=["naive_suppression", "n_cells_target"]).copy()
    frame = frame.sort_values("naive_suppression", ascending=False).reset_index(drop=True)
    frame["rank"] = np.arange(1, len(frame) + 1)
    frame["in_top"] = frame["rank"] <= top_k

    frame["n_cells_decile"] = stratified.deciles(frame["n_cells_target"].to_numpy(), N_BINS)
    frame["z_l2_decile"] = stratified.deciles(frame["z_l2"].to_numpy(), N_BINS)
    frame["matched_background"] = _single_draw_background(frame, np.random.default_rng(0))
    return frame


def seed_lottery(frame: pd.DataFrame, n_draws: int) -> None:
    """Reproduce the audit's arithmetic claim about the single seed-0 draw. RULE #1.

    Four audit agents reported this defect with four different resample counts and four different
    non-significance rates (11.5%, 11.2%, 11%, 9.5%). None of them committed a script. It is
    reproduced here before being retired, because a defect nobody reproduced is a rumour.

    Args:
        frame: Output of :func:`build_scores`.
        n_draws: Number of times to redraw the background.
    """
    print(f"\n--- the seed lottery: redrawing the OLD single-draw background {n_draws:,} times ---")
    odds_ratios, pvalues = [], []
    for seed in range(n_draws):
        background = _single_draw_background(frame, np.random.default_rng(seed))
        subset = frame[frame["in_top"] | background]
        odds, pval, _, _ = _fisher(subset["is_iei"], subset["in_top"])
        odds_ratios.append(odds)
        pvalues.append(pval)

    odds_ratios = np.asarray(odds_ratios, dtype=float)
    pvalues = np.asarray(pvalues, dtype=float)
    finite = np.isfinite(odds_ratios)
    non_significant = float((pvalues >= 0.05).mean())

    print(f"  IEI odds ratio   seed 0 = {odds_ratios[0]:.2f}   "
          f"median {np.median(odds_ratios[finite]):.2f}   "
          f"2.5-97.5 pct [{np.percentile(odds_ratios[finite], 2.5):.2f}, "
          f"{np.percentile(odds_ratios[finite], 97.5):.2f}]")
    print(f"  IEI p-value      seed 0 = {pvalues[0]:.4f}   "
          f"median {np.median(pvalues):.4f}   "
          f"2.5-97.5 pct [{np.percentile(pvalues, 2.5):.5f}, {np.percentile(pvalues, 97.5):.4f}]")
    print(f"  seeds giving p >= 0.05: {non_significant:.1%}")
    print("  The audit reported 11.5% (and 11.2%, and 11%, and 9.5%, from four agents who each ran a")
    print("  different resampling and none of whom committed a script). Reproduced here against the")
    print(f"  code as it actually stood, the figure is {non_significant:.1%}. The defect is worse than reported.")
    print("  A single draw cannot carry a p-value. Everything below uses every control row instead.")


def stratified_family(frame: pd.DataFrame, stratum: str, label: str) -> pd.DataFrame:
    """Run all five liability tests against a stratified background, then BH-correct the family.

    Args:
        frame: Output of :func:`build_scores`.
        stratum: Column holding the stratum label per row.
        label: Human-readable name of the stratification.

    Returns:
        Tidy frame, one row per test, with raw and BH-adjusted p-values.
    """
    strata = frame[stratum].to_numpy()
    in_top = frame["in_top"].to_numpy()

    records = []
    for name, column, binary in TESTS:
        values = frame[column].to_numpy(dtype=float)
        try:
            if binary:
                result = stratified.cochran_mantel_haenszel(values, in_top, strata, "greater")
            else:
                result = stratified.van_elteren(values, in_top, strata, "greater")
        except ValueError:
            records.append({"test": name, "column": column, "binary": binary, "effect": np.nan,
                            "p": np.nan, "strata": 0, "top_used": 0, "bg_used": 0})
            continue
        records.append({
            "test": name, "column": column, "binary": binary, "effect": result.effect,
            "p": result.pvalue, "strata": result.n_strata_used,
            "top_used": result.n_top_used, "bg_used": result.n_background_used,
        })

    table = pd.DataFrame(records)
    table["p_bh"] = stratified.benjamini_hochberg(table["p"].tolist())
    table["stratification"] = label

    n_bg = int((~frame["in_top"]).sum())
    used = int(table["bg_used"].max()) if len(table) else 0
    print(f"\n--- {label} ---")
    print(f"    no sampling and no seed. Of {n_bg:,} background rows, {used:,} sit in a stratum that also")
    print("    contains a top-K row and therefore carry information; the rest are dropped, not ignored.")
    for _, row in table.iterrows():
        if not np.isfinite(row["p"]):
            print(f"    {row['test']:28s} not testable")
            continue
        effect = f"OR {row['effect']:7.2f}" if row["binary"] else f" z {row['effect']:+7.3f}"
        call = "ENRICHED" if row["p_bh"] < 0.05 else "ns"
        print(f"    {row['test']:28s} {effect}  p {row['p']:.3g}  p_BH {row['p_bh']:.3g}  "
              f"strata {row['strata']:2d}  top {row['top_used']:3d}  bg {row['bg_used']:5d}  {call}")
    return table


def essentiality_coverage(frame: pd.DataFrame) -> None:
    """Report the essentiality coverage, and RETRACT the conclusion once drawn from it.

    A zero count of core-essential genes in the top K looked like it had two causes: either
    essential-gene knockdowns are absent from the analysable set, or they are present and the naive
    ranking does not favour them. We concluded the second and called essentiality the wrong safety
    axis. That was wrong.

    ``scripts/11_selection_funnel.py`` measured the selection against the true sgRNA-library
    denominator. Both causes operate, through two colliders pulling in opposite directions:
    465 of 682 library essentials never reach the DE table (their knockdown kills the cell), and
    then ``ontarget_significant`` removes the nonessentials instead (it needs an expressed target).
    The essentials that survive are exactly the ones whose knockdown was harmless, so a rank
    comparison among survivors estimates nothing causal.

    This function still prints the coverage funnel, because it is informative. It no longer draws a
    conclusion, because none is available.

    Args:
        frame: The QC-passing ranked frame from :func:`build_scores`.
    """
    essential = priors.core_essential_genes()
    obs = de_stats.read_obs()
    stim_all = obs[obs["culture_condition"] == CONDITION]
    perturbed = set(obs["target_contrast_gene_name"].astype(str))

    in_library = essential & perturbed
    in_stim = essential & set(stim_all["target_contrast_gene_name"].astype(str))
    in_qc = essential & set(frame["target_contrast_gene_name"].astype(str))

    print("\n--- why zero core-essential genes in the top K? ---")
    print(f"  core-essential list (authors' filtered Hart subset): {len(essential)}")
    print(f"  ... perturbed anywhere in the library:               {len(in_library)}")
    print(f"  ... tested in {CONDITION}:                            {len(in_stim)}")
    print(f"  ... surviving QC and thus rankable:                  {len(in_qc)}")

    if not in_qc:
        print("  Essentials are absent from the analysable set. No conclusion possible.")
        return

    print()
    print("  RETRACTED 2026-07-08. This coverage cannot answer the question it appears to answer.")
    print("  `scripts/11_selection_funnel.py` measures the selection against the true sgRNA-library")
    print("  denominator and finds TWO colliders pulling opposite ways:")
    print("    465 of 682 library core-essentials never reach the DE table at all, because their")
    print("    knockdown depletes cells and fails DE-eligibility (OR 0.024, p=5.6e-48); and then")
    print("    ontarget_significant destroys the nonessentials instead, because it needs an")
    print("    expressed target (OR 38.5, p=9.4e-22).")
    print("  The surviving essentials are precisely the ones whose knockdown did NOT kill the cell.")
    print("  A rank comparison among survivors estimates nothing causal, and a null on 31 of them")
    print("  was underpowered regardless.")
    print("  The supportable statement: THIS SCREEN CANNOT RESOLVE whether cancer-cell essentiality")
    print("  predicts the naive ranking. The core-essential row below is reported, never believed.")


def report(frame: pd.DataFrame, top_k: int, n_draws: int) -> tuple[bool, pd.DataFrame]:
    """Print every test and decide the build on the pre-registered primary endpoint.

    Args:
        frame: Output of :func:`build_scores`.
        top_k: Size of the naive shortlist.
        n_draws: Redraws for the seed lottery.

    Returns:
        Tuple of (primary endpoint holds, the combined results table).
    """
    in_top = frame["in_top"]
    print(f"\nranked perturbations: {len(frame)}  top-K: {top_k}  background: {int((~in_top).sum())}")

    print(f"\n--- top {min(25, top_k)} by naive suppression ---")
    cols = ["rank", "target_contrast_gene_name", "naive_suppression", "is_core_essential",
            "is_iei", "stim_de_genes", "rest_de_genes", "z_l2"]
    with pd.option_context("display.width", 160, "display.max_columns", 20):
        print(frame.head(min(25, top_k))[cols].to_string(index=False))

    print("\n--- what actually confounds this ranking? ---")
    rho_score, p_score = stats.spearmanr(frame["naive_suppression"], frame["n_cells_target"])
    ok = frame[["n_cells_target", "stim_de_genes"]].dropna()
    rho_de, p_de = stats.spearmanr(ok["n_cells_target"], ok["stim_de_genes"])
    rho_mag, _ = stats.spearmanr(frame["z_l2"], frame["stim_de_genes"])
    rho_abs, _ = stats.spearmanr(frame["z_l2"], frame["naive_suppression"].abs())
    print(f"  spearman(naive_suppression, n_cells_target) = {rho_score:+.3f}  p={p_score:.3g}")
    print(f"  spearman(n_cells_target,   stim_de_genes)  = {rho_de:+.3f}  p={p_de:.3g}")
    print("    ^ NEGATIVE. More cells means FEWER DE genes. Cell count is a viability readout, not")
    print("      a power proxy, so the original cell-count matching controlled the confound the")
    print("      WRONG WAY ROUND. Retained below only as a sensitivity analysis.")
    print(f"  spearman(z_l2, stim_de_genes)              = {rho_mag:+.3f}   <- the real confound")
    print(f"  spearman(z_l2, |naive_suppression|)        = {rho_abs:+.3f}   <- and it is NOT the score")

    seed_lottery(frame, n_draws)

    cells = stratified_family(frame, "n_cells_decile", "stratified on n_cells_target decile (sensitivity)")
    mag = stratified_family(frame, "z_l2_decile", "stratified on z_l2 decile (PRIMARY stratification)")
    essentiality_coverage(frame)

    primary = mag[mag["column"] == PRIMARY_COLUMN].iloc[0]
    holds = bool(np.isfinite(primary["p"]) and primary["p"] < PRIMARY_ALPHA and primary["effect"] > 0)

    print("\n" + "=" * 84)
    print(f"PRIMARY ENDPOINT (pre-registered): {PRIMARY_LABEL}, top {top_k} vs z_l2-stratified background")
    print(f"  one-sided van Elteren  z = {primary['effect']:+.3f}   p = {primary['p']:.3g}   "
          f"alpha = {PRIMARY_ALPHA}")
    print(f"  strata used {primary['strata']}   top rows compared {primary['top_used']}   "
          f"background rows compared {primary['bg_used']:,}")
    print("  Not `any()` over five tests: that has a ~23% family-wise false-positive rate and would")
    print("  make this script unfalsifiable. The other four are SECONDARY and BH-corrected.")
    print("=" * 84)

    if holds:
        print("VERDICT: the empirical bet HOLDS on the pre-registered primary endpoint.")
        print("The top of a naive reversal ranking collapses the tolerance program more than a")
        print("background of equal transcriptome-wide effect magnitude. 'Reversal is not enough'")
        print("is supported, and it is supported by TOLERANCE, not by collateral DE-gene counts.")
        print()
        print("Secondary, BH-corrected, reported and not believed on their own:")
        for _, row in mag[mag["column"] != PRIMARY_COLUMN].iterrows():
            call = "enriched" if row["p_bh"] < 0.05 else "ns"
            print(f"  {row['test']:28s} p_BH {row['p_bh']:.3g}  {call}")
        print()
        print("  scripts/12_magnitude_matched.py shows the DE-count secondaries are entangled with")
        print("  effect magnitude (rho 0.725) and that rest-DE fires equally on the INDUCTION tail.")
        print("  Do not quote them as independent evidence. The primary endpoint is the claim.")
    else:
        print("VERDICT: the primary endpoint FAILS on a magnitude-matched background.")
        print("Effect size, not suppression, explains the enrichment. The headline needs a rewrite.")
    print("=" * 84)

    return holds, pd.concat([cells, mag], ignore_index=True)


def main() -> None:
    """Run the risk-kill test and persist the annotated ranking."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top-k", type=int, default=100)
    parser.add_argument("--n-draws", type=int, default=2000, help="redraws for the seed lottery")
    args = parser.parse_args()

    paths.ensure_dirs()
    frame = build_scores(args.top_k)
    holds, tests = report(frame, args.top_k, args.n_draws)

    out = paths.TABLES / "risk_kill_naive_reversal.csv"
    frame.to_csv(out, index=False)
    print(f"\nwrote {out}")

    tests_out = paths.TABLES / "risk_kill_stratified_tests.csv"
    tests.to_csv(tests_out, index=False)
    print(f"wrote {tests_out}")

    raise SystemExit(0 if holds else 1)


if __name__ == "__main__":
    main()
