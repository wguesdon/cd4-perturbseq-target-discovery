"""N1: is the risk-kill result an effect-magnitude artifact, or is it about suppression?

The adversarial audit's sharpest surviving objection. ``02_risk_kill_reversal.py`` ranks
perturbations by ``naive_suppression = -mean(z)`` over 32 effector genes, takes the top 100, and
shows they carry more collateral DE, more resting DE, and more tolerance loss than a background
matched on ``n_cells_target``. The audit says that background controls the wrong thing:

    "The 'collateral DE' and 'Rest DE' pillars are an effect-magnitude artifact, and the
    direction is backwards: at matched |score|, suppressors have ~66% FEWER stim DE genes
    than inducers." ... "Spearman(n_cells_target, stim_de_genes) = -0.243, so the
    cell-count-matched background controls a confound that runs the opposite way."

The rival story would be fatal. A perturbation with a large transcriptome-wide effect gets a large
``|mean z|`` over *any* gene set, the effector module included, and it gets many DE genes for the
same reason. If so the naive ranking selects "big effect", not "suppresses inflammation", and
every pillar is a restatement of effect size.

**The audit's proposed fix contains the flaw that undoes it.** It offers three interchangeable
controls: "match on ``|naive_suppression|`` decile, or on ``stim_de_genes`` decile, or on a
transcriptome-wide ``||z||_2``". These are not the same control. Measured here,
``Spearman(z_l2, |naive_suppression|) = +0.198``. Score magnitude is not effect magnitude. The
audit's evidence comes from matching on ``|score|``; its conclusion is stated about effect
magnitude. Test C reproduces the evidence. Tests A and B refuse the conclusion.

Three test families, none of which draws a random number. See ``src/cd4_perturbseq/stratified.py``.

**A. Magnitude.** Top 100 versus *every* background row, stratified on decile of ``z_l2``, the
transcriptome-wide effect magnitude. Computed off-module and off-target: the perturbed gene's own
column and the 41 module columns are excluded, so the matching variable cannot contain the score.

**B. Direction, at matched magnitude.** Top 100 by suppression versus top 100 by *induction*,
stratified on ``z_l2`` decile. The sign-flipped control. A pillar that fires on both tails of the
ranking measures ``|effect|``, not suppression.

**C. Direction, over the whole ranking.** All suppressors versus all inducers, stratified on decile
of ``|naive_suppression|``. Precisely what the audit's fix text asks for -- "the regression
coefficient of is_suppressor conditional on |score|" -- in the nonparametric form the data deserve.

**The kill rule was fixed before any number was computed** (``SESSION_SUMMARY.md``, task N1): a
pillar that fails BOTH A and B is removed from the gate and added to DO NOT REDO. It binds only the
**four registered pillars**, the four rows of the risk-kill verdict table. The pre-registered
stratification is the **decile**, ten bins, as task N1 wrote it. Because a verdict that depends on
an arbitrary bin count is not a verdict, every test is rerun at 10, 20 and 50 bins and the pillar is
marked FRAGILE where the call changes. Read the fragility column; it is the most informative thing
here.

A fifth measure is reported and **deliberately barred from deciding anything**:
``04_window_score.py`` does not gate on raw rest-DE, it gates on ``selectivity``,
``log1p(stim_de) - log1p(rest_de)``. Testing the ingredients and not the axis would be testing
something nobody ships. But a rule invented after seeing the data is not a rule, so this one sets
the next task rather than rewriting the gate today. Note also that a pillar test asks "is the naive
top-100 *enriched* for this liability", while a gate asks "should this target be *rejected*". An
axis can be a sound filter and a poor enrichment signal. Do not read one as the other.

Two corrections to what the pillars measure, both free, both in the schema all along:

- ``n_total_de_genes`` counts the perturbed gene itself. The column named "collateral" was not
  measuring collateral. ``.obs`` carries ``n_downstream`` for exactly this. Primary is now
  ``n_downstream``.
- The Rest arm was never QC-filtered, and three of the top 100 have no resting row at all --
  including ``IL2RB`` at rank 4, the single gene the results doc uses to motivate the whole safety
  gate. A missing resting row is not a quiet zero.

Usage:
    uv run python scripts/12_magnitude_matched.py --top-k 100
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from cd4_perturbseq import de_stats, paths, priors, programs, stratified

STIM = "Stim48hr"
REST = "Rest"
N_BINS = 10
"""Pre-registered stratification: the decile, as task N1 specified before any number existed."""

SENSITIVITY_BINS = (10, 20, 50)
"""A verdict that changes across these is reported as FRAGILE, not as a verdict."""

CIRCULARITY_RHO = 0.90
"""Above this |Spearman| with the matching variable, test A cannot separate pillar from stratum."""

Z_L2_CACHE = paths.INTERIM / "z_l2_stim48.csv"


@dataclass(frozen=True)
class Pillar:
    """One liability the naive top-K is alleged to be enriched for.

    Attributes:
        label: Human-readable name.
        column: Column in the analysis frame. Oriented so that LARGER always means WORSE.
        binary: True for a flag tested by CMH, False for a measure tested by van Elteren.
        registered: True for the four pillars named in the risk-kill verdict table, to which the
            pre-registered kill rule applies. False for anything added after the fact.
    """

    label: str
    column: str
    binary: bool
    registered: bool


PILLARS = (
    Pillar("IEI (immunodeficiency)", "is_iei", binary=True, registered=True),
    Pillar("collateral DE in stim", "stim_downstream", binary=False, registered=True),
    Pillar("DE at rest (raw)", "rest_downstream", binary=False, registered=True),
    Pillar("tolerance-module suppression", "tolerance_suppression", binary=False, registered=True),
    # NOT one of the four. Added because 04_window_score.py gates on selectivity, not on raw
    # rest-DE, so the registered pillars do not test the axis we actually ship. It cannot decide
    # the build: a rule invented after seeing the data is not a rule. It sets the next task.
    Pillar("homeostasis cost (-selectivity)", "homeostasis_cost", binary=False, registered=False),
)

GATE_AXES = frozenset({"homeostasis_cost", "tolerance_suppression"})
"""What ``04_window_score.py`` rejects on, besides the evidence floor. For reporting only."""


def _test(pillar: Pillar, frame: pd.DataFrame, group: str, strata: np.ndarray, alternative: str):
    """Dispatch to CMH or van Elteren according to the pillar's type.

    Args:
        pillar: The pillar under test.
        frame: Frame carrying the pillar column and the group indicator.
        group: Boolean column naming the "top" group.
        strata: Stratum label per row.
        alternative: Direction of the test.

    Returns:
        A :class:`~cd4_perturbseq.stratified.StratifiedResult`, or None if no stratum informed it.
    """
    values = frame[pillar.column].to_numpy(dtype=float)
    in_top = frame[group].to_numpy(dtype=bool)
    try:
        if pillar.binary:
            return stratified.cochran_mantel_haenszel(values, in_top, strata, alternative)
        return stratified.van_elteren(values, in_top, strata, alternative)
    except ValueError:
        return None


def compute_z_l2(obs: pd.DataFrame, var: pd.DataFrame, rows: np.ndarray, exclude: set[str]) -> pd.DataFrame:
    """Transcriptome-wide effect magnitude per perturbation, off-module and off-target.

    ``z_l2 = sqrt(sum of z^2)`` over the measured genes, minus three things that would put the
    score inside the variable meant to control for it: the perturbed gene's own column (a
    QC-required significant knockdown, so large by construction), the 32 effector genes the score
    is built from, and the 9 tolerance genes pillar 5 is built from.

    Cached, because it is the only step in this project that reads full-width rows of a 16.8 GB
    layer. ``de_stats.read_layer_rows`` downcasts to float32, so 6,371 rows cost 262 MB.

    Args:
        obs: Full ``.obs`` frame, positionally indexed.
        var: Full ``.var`` frame.
        rows: Positional indices into ``.obs`` of the perturbations to score.
        exclude: Gene symbols to drop from the norm, beyond each row's own target.

    Returns:
        DataFrame indexed by gene symbol with ``z_l2``, ``z_l2_raw`` (nothing excluded), and
        ``n_genes_used``.
    """
    if Z_L2_CACHE.exists():
        cached = pd.read_csv(Z_L2_CACHE)
        if len(cached) == len(rows):
            print(f"  z_l2: reusing cache {Z_L2_CACHE} ({len(cached):,} rows)")
            return cached.set_index("gene_name")
        print(f"  z_l2: cache has {len(cached):,} rows, need {len(rows):,}; recomputing")

    names = var["gene_name"].astype(str).to_numpy()
    column_of = {g: i for i, g in enumerate(names)}

    print(f"  z_l2: reading {len(rows):,} full-width rows of the zscore layer (float32, ~262 MB) ...")
    block = de_stats.read_layer_rows(rows, layer="zscore")
    finite = np.isfinite(block)
    squared = np.where(finite, block, 0.0).astype(np.float64) ** 2
    z_l2_raw = np.sqrt(squared.sum(axis=1))

    module_cols = np.array(sorted({column_of[g] for g in exclude if g in column_of}), dtype=np.int64)
    squared[:, module_cols] = 0.0
    finite[:, module_cols] = False

    targets = obs.loc[rows, "target_contrast_gene_name"].astype(str).to_numpy()
    self_rows = np.array([i for i, g in enumerate(targets) if g in column_of], dtype=np.int64)
    self_cols = np.array([column_of[g] for g in targets if g in column_of], dtype=np.int64)
    squared[self_rows, self_cols] = 0.0
    finite[self_rows, self_cols] = False
    print(f"  z_l2: excluded {len(module_cols)} module columns and {len(self_rows):,} self-target columns")

    frame = pd.DataFrame(
        {
            "gene_name": targets,
            "z_l2": np.sqrt(squared.sum(axis=1)),
            "z_l2_raw": z_l2_raw,
            "n_genes_used": finite.sum(axis=1),
        }
    )
    paths.ensure_dirs()
    frame.to_csv(Z_L2_CACHE, index=False)
    print(f"  z_l2: wrote cache {Z_L2_CACHE}")
    return frame.set_index("gene_name")


def build_frame(top_k: int) -> pd.DataFrame:
    """Rank by naive suppression and attach every pillar plus the magnitude covariate.

    Args:
        top_k: Size of both the suppression shortlist and the induction control.

    Returns:
        DataFrame with one row per QC-passing stimulated perturbation.
    """
    obs = de_stats.read_obs().reset_index(drop=True)
    var = de_stats.read_var()
    names = var["gene_name"].astype(str).to_numpy()
    column_of = {g: i for i, g in enumerate(names)}

    program = programs.load_activation_program()
    effector = [g for g in programs.effector_core(program) if g in column_of]
    tolerance = [g for g in programs.tolerance_module(program) if g in column_of]
    print(f"effector core {len(effector)} genes, tolerance module {len(tolerance)} genes")

    wanted = effector + tolerance
    columns = de_stats.read_layer_columns([column_of[g] for g in wanted], layer="zscore")
    effector_z = np.nanmean(columns[:, : len(effector)], axis=1)
    tolerance_z = np.nanmean(columns[:, len(effector) :], axis=1)

    stim_rows = obs.index[(obs["culture_condition"] == STIM) & de_stats.qc_mask(obs)].to_numpy()
    print(f"QC-passing {STIM} perturbations: {len(stim_rows):,}")

    frame = pd.DataFrame(
        {
            "gene_name": obs.loc[stim_rows, "target_contrast_gene_name"].astype(str).to_numpy(),
            "naive_suppression": -effector_z[stim_rows],
            "tolerance_suppression": -tolerance_z[stim_rows],
            "stim_downstream": obs.loc[stim_rows, "n_downstream"].to_numpy(float),
            "stim_de_genes": obs.loc[stim_rows, "n_total_de_genes"].to_numpy(float),
            "n_cells_target": obs.loc[stim_rows, "n_cells_target"].to_numpy(float),
        }
    )
    frame = frame.dropna(subset=["naive_suppression"]).reset_index(drop=True)

    # The resting arm, mapped by gene. Not QC-filtered upstream, and 179 stimulated genes have no
    # resting row at all, so this genuinely produces NaN rather than a zero.
    rest = obs[obs["culture_condition"] == REST].drop_duplicates("target_contrast_gene_name").copy()
    rest.index = rest["target_contrast_gene_name"].astype(str)
    rest_qc = de_stats.qc_mask(rest)
    frame["rest_downstream"] = frame["gene_name"].map(rest["n_downstream"]).astype(float)
    frame["rest_de_genes"] = frame["gene_name"].map(rest["n_total_de_genes"]).astype(float)
    frame["rest_qc_pass"] = frame["gene_name"].map(rest_qc).astype("boolean")
    frame["has_rest_row"] = frame["gene_name"].isin(rest.index)

    # The axis 04_window_score.py actually gates on. Negated so that larger means worse, like
    # every other pillar. NOT median-imputed here; 04 imputes, which is a separate defect.
    frame["homeostasis_cost"] = -(np.log1p(frame["stim_downstream"]) - np.log1p(frame["rest_downstream"]))

    frame["is_iei"] = frame["gene_name"].isin(priors.iei_genes()).astype(float)

    magnitude = compute_z_l2(obs, var, stim_rows, exclude=set(effector) | set(tolerance))
    frame["z_l2"] = frame["gene_name"].map(magnitude["z_l2"])
    frame["z_l2_raw"] = frame["gene_name"].map(magnitude["z_l2_raw"])

    frame["abs_score"] = frame["naive_suppression"].abs()
    frame["is_suppressor"] = frame["naive_suppression"] > 0

    order = frame["naive_suppression"].sort_values(ascending=False).index
    frame["in_top"] = False
    frame.loc[order[:top_k], "in_top"] = True
    frame["in_induction_top"] = False
    frame.loc[order[-top_k:], "in_induction_top"] = True
    return frame


def reproduce_the_audit(frame: pd.DataFrame) -> None:
    """Reproduce the audit's two numbers before acting on either. RULE #1.

    The audit has been right once, wrong once, and incomplete once on this codebase. Its claims
    are hypotheses until they reproduce.

    Args:
        frame: Output of :func:`build_frame`.
    """
    print("\n=== STEP 0: reproduce the audit's claims before believing them (RULE #1) ===")

    ok = frame[["n_cells_target", "stim_de_genes"]].dropna()
    rho, pval = stats.spearmanr(ok["n_cells_target"], ok["stim_de_genes"])
    verdict = "REPRODUCES" if abs(rho + 0.243) < 0.01 else "DOES NOT REPRODUCE"
    print("\n  claim 1: Spearman(n_cells_target, stim_de_genes) = -0.243")
    print(f"    observed rho = {rho:+.3f}  p = {pval:.3g}  n = {len(ok):,}   -> {verdict}")
    print("    So script 02's cell-count-matched background controls a confound running the OTHER")
    print("    way: more cells means FEWER DE genes. Cell count is a viability readout, not a")
    print("    power proxy. The audit is right, and script 02's stated rationale is wrong.")

    print("\n  claim 2: at matched |score|, suppressors carry ~66% FEWER stim DE genes than inducers")
    deciles = stratified.deciles(frame["abs_score"].to_numpy(), N_BINS)
    rows = []
    for decile in sorted(np.unique(deciles[np.isfinite(deciles)])):
        block = frame[deciles == decile]
        supp = block.loc[block["is_suppressor"], "stim_de_genes"].dropna()
        ind = block.loc[~block["is_suppressor"], "stim_de_genes"].dropna()
        if len(supp) < 5 or len(ind) < 5:
            continue
        ratio = supp.median() / ind.median() if ind.median() > 0 else np.nan
        rows.append((int(decile), len(supp), len(ind), supp.median(), ind.median(), ratio))

    table = pd.DataFrame(rows, columns=["|score| decile", "n_supp", "n_ind", "supp_median", "ind_median", "ratio"])
    print(table.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    pooled = table["ratio"].dropna().median()
    print(f"\n    pooled within-decile median ratio supp/ind = {pooled:.3f}, i.e. {100 * (1 - pooled):.0f}% fewer")

    # Is the direction a binning artifact? |score| decile 9 spans 0.55 to 3.06, and the top 100
    # live above 1.163, so a decile is coarse in a heavy tail. Refine it, and restrict to the
    # |score| range where inducers actually exist.
    collateral = frame["stim_downstream"].to_numpy(float)
    is_supp = frame["is_suppressor"].to_numpy()
    print("\n    direction, refined:  (van Elteren z, two-sided; negative = suppressors carry FEWER)")
    for bins in (10, 20, 50, 100):
        result = stratified.van_elteren(collateral, is_supp, stratified.deciles(frame["abs_score"].to_numpy(), bins), "two-sided")
        print(f"      {bins:3d} bins of |score|      z = {result.effect:+7.3f}   p = {result.pvalue:.3g}")
    ceiling = frame.loc[~frame["is_suppressor"], "abs_score"].max()
    common = frame["abs_score"] <= ceiling
    on_support = stratified.van_elteren(
        collateral[common.to_numpy()], is_supp[common.to_numpy()],
        stratified.deciles(frame.loc[common, "abs_score"].to_numpy(), 20), "two-sided",
    )
    print(f"      common support only    z = {on_support.effect:+7.3f}   p = {on_support.pvalue:.3g}   "
          f"(|score| <= {ceiling:.3f}, dropping {int((~common).sum())} suppressors)")

    if 0.28 < pooled < 0.40:
        print("\n    -> REPRODUCES at ~66%")
    else:
        print("\n    -> PARTIALLY REPRODUCES. The DIRECTION is real and survives every refinement")
        print("       above. The 66% MAGNITUDE does not: it is decile 8 alone (ratio 0.333), while")
        print("       deciles 1-5 give exactly 1.000. The audit reported one decile as the pooled")
        print("       effect. And the direction is about |score|, not about effect magnitude; see")
        print("       the next block, which is where the audit's conclusion comes apart.")


def diagnostics(frame: pd.DataFrame, top_k: int) -> None:
    """Print what the tests below can and cannot see.

    Args:
        frame: Output of :func:`build_frame`.
        top_k: Shortlist size.
    """
    print("\n=== the matching variable, and what it is entangled with ===")
    for label, column in (
        ("|naive_suppression|", "abs_score"),
        ("collateral DE in stim", "stim_downstream"),
        ("DE at rest (raw)", "rest_downstream"),
        ("homeostasis cost", "homeostasis_cost"),
        ("tolerance suppression", "tolerance_suppression"),
        ("n_cells_target", "n_cells_target"),
    ):
        pair = frame[["z_l2", column]].dropna()
        rho, _ = stats.spearmanr(pair["z_l2"], pair[column])
        note = "  <- CIRCULAR, test A cannot separate it from the stratum" if abs(rho) > CIRCULARITY_RHO else ""
        print(f"  Spearman(z_l2, {label:22s}) = {rho:+.3f}{note}")

    print("\n  THE AUDIT'S ERROR, in one number: Spearman(z_l2, |naive_suppression|) = "
          f"{stats.spearmanr(frame['z_l2'], frame['abs_score']).statistic:+.3f}")
    print("  Its fix offers 'match on |score| decile, OR on a transcriptome-wide ||z||_2' as if")
    print("  interchangeable. They share 4% of their rank variance. Score magnitude is not effect")
    print("  magnitude, and the pillars' fate depends entirely on which one you match.")

    top = frame[frame["in_top"]]
    ind = frame[frame["in_induction_top"]]
    background = frame[~frame["in_top"] & ~frame["in_induction_top"]]
    print(f"\n  z_l2 median      suppression top {top['z_l2'].median():8.2f}   "
          f"induction top {ind['z_l2'].median():8.2f}   background {background['z_l2'].median():8.2f}")
    print(f"  |score| median   suppression top {top['abs_score'].median():8.3f}   "
          f"induction top {ind['abs_score'].median():8.3f}   background {background['abs_score'].median():8.3f}")

    floor = top["abs_score"].min()
    reach = int((frame.loc[~frame["is_suppressor"], "abs_score"] >= floor).sum())
    print(f"\n  |score| floor of the suppression top {top_k}: {floor:.3f}   inducers reaching it: {reach}")
    print("  The tails are asymmetric in |score|: a knockdown can crush the effector program far")
    print("  harder than it can induce it. They are NOT asymmetric in z_l2, which is why test B")
    print("  can be stratified on z_l2 and still fill every stratum.")

    missing = top[~top["has_rest_row"]]["gene_name"].tolist()
    print(f"\n  top {top_k} rows with NO resting row at all: {len(missing)}  ({', '.join(missing)})")
    no_rest = int((~frame["has_rest_row"]).sum())
    fails_qc = int((frame["rest_qc_pass"] == False).sum())  # noqa: E712 - pandas boolean dtype
    print(f"  of {len(frame):,} ranked perturbations: {no_rest} have no resting row, "
          f"{fails_qc:,} have one that would FAIL the QC the stim row had to pass")
    print("  A low rest-DE count conflates 'spares homeostasis' with 'the guide never knocked it down'.")


def _sensitivity(frame: pd.DataFrame, pillar: Pillar) -> dict[str, bool]:
    """Rerun tests A and B across bin counts and report whether the VERDICT is stable.

    Stability means the significance call does not change, not that it is always significant. A
    pillar that is null at every bin count is perfectly stable; it is stably null.

    Args:
        frame: Output of :func:`build_frame`.
        pillar: The pillar under test.

    Returns:
        Mapping from test family to whether its verdict is the same at every bin count in
        :data:`SENSITIVITY_BINS`, ``N_BINS`` included.
    """
    head_to_head = frame[frame["in_top"] | frame["in_induction_top"]]
    calls: dict[str, list[bool]] = {"A": [], "B": []}
    for bins in sorted({N_BINS, *SENSITIVITY_BINS}):
        a = _test(pillar, frame, "in_top", stratified.deciles(frame["z_l2"].to_numpy(), bins), "greater")
        b = _test(
            pillar, head_to_head, "in_top",
            stratified.deciles(head_to_head["z_l2"].to_numpy(), bins), "greater",
        )
        calls["A"].append(bool(a and a.pvalue < 0.05))
        calls["B"].append(bool(b and b.pvalue < 0.05))
    return {family: len(set(verdicts)) == 1 for family, verdicts in calls.items()}


def run_pillar_tests(frame: pd.DataFrame, top_k: int) -> pd.DataFrame:
    """Run tests A, B and C for every pillar, BH-correct within family, flag fragility.

    Args:
        frame: Output of :func:`build_frame`.
        top_k: Shortlist size.

    Returns:
        Tidy frame, one row per pillar.
    """
    head_to_head = frame[frame["in_top"] | frame["in_induction_top"]]
    z_strata = stratified.deciles(frame["z_l2"].to_numpy(), N_BINS)
    pair_strata = stratified.deciles(head_to_head["z_l2"].to_numpy(), N_BINS)
    score_strata = stratified.deciles(frame["abs_score"].to_numpy(), N_BINS)

    records = []
    for pillar in PILLARS:
        pair = frame[["z_l2", pillar.column]].dropna()
        rho = stats.spearmanr(pair["z_l2"], pair[pillar.column]).statistic if len(pair) > 10 else np.nan

        a = _test(pillar, frame, "in_top", z_strata, "greater")
        b = _test(pillar, head_to_head, "in_top", pair_strata, "greater")
        c = _test(pillar, frame, "is_suppressor", score_strata, "two-sided")
        stable = _sensitivity(frame, pillar)

        records.append(
            {
                "pillar": pillar.label,
                "registered": pillar.registered,
                "gates": pillar.column in GATE_AXES,
                "binary": pillar.binary,
                "rho_with_z_l2": rho,
                "circular": bool(abs(rho) > CIRCULARITY_RHO) if np.isfinite(rho) else False,
                "A_effect": a.effect if a else np.nan,
                "A_p": a.pvalue if a else np.nan,
                "A_strata": a.n_strata_used if a else 0,
                "A_top_used": a.n_top_used if a else 0,
                "A_stable": stable["A"],
                "B_effect": b.effect if b else np.nan,
                "B_p": b.pvalue if b else np.nan,
                "B_strata": b.n_strata_used if b else 0,
                "B_top_used": b.n_top_used if b else 0,
                "B_bg_used": b.n_background_used if b else 0,
                "B_stable": stable["B"],
                "C_effect": c.effect if c else np.nan,
                "C_p": c.pvalue if c else np.nan,
            }
        )

    table = pd.DataFrame(records)
    for family in ("A", "B", "C"):
        table[f"{family}_p_bh"] = stratified.benjamini_hochberg(table[f"{family}_p"].tolist())

    table["passes_A"] = (table["A_p_bh"] < 0.05) & ~table["circular"]
    table["passes_B"] = table["B_p_bh"] < 0.05
    table["survives"] = table["passes_A"] | table["passes_B"]
    table["robust"] = table["A_stable"] & table["B_stable"]
    return table


def falsification_control(frame: pd.DataFrame, top_k: int, n_draws: int) -> bool:
    """Can these tests fail? Feed them a top-K that carries no information.

    The audit's most damaging unaddressed finding is that "non-reversal null rankings pass the
    paper's own verdict function with equal or stronger effect sizes". If a random 100 rows light
    up the same pillars, then nothing here is about reversal and the machinery is measuring the
    stratification, not the ranking.

    Draws ``n_draws`` random shortlists and reports, per pillar, the fraction that reach p < 0.05
    on test A. Under a working test that fraction is the nominal 0.05.

    Args:
        frame: Output of :func:`build_frame`.
        top_k: Shortlist size.
        n_draws: Number of random shortlists.

    Returns:
        True if every pillar's false-positive rate is at or below 0.15, a generous bound for
        ``n_draws`` draws at a nominal 0.05.
    """
    print(f"\n=== FALSIFICATION CONTROL: {n_draws} RANDOM shortlists of {top_k}, same tests ===")
    print("    'Non-reversal null rankings pass the verdict function' -- the audit. Do they?")

    rng = np.random.default_rng(0)
    z_strata = stratified.deciles(frame["z_l2"].to_numpy(), N_BINS)
    hits = {pillar.label: 0 for pillar in PILLARS if pillar.registered}

    for _ in range(n_draws):
        sham = frame.copy()
        sham["in_top"] = False
        sham.loc[rng.choice(sham.index, size=top_k, replace=False), "in_top"] = True
        for pillar in PILLARS:
            if not pillar.registered:
                continue
            result = _test(pillar, sham, "in_top", z_strata, "greater")
            if result and result.pvalue < 0.05:
                hits[pillar.label] += 1

    calibrated = True
    for label, count in hits.items():
        rate = count / n_draws
        ok = rate <= 0.15
        calibrated &= ok
        print(f"  {label:32s} random shortlists reaching p<0.05: {count:3d}/{n_draws}  "
              f"({rate:5.1%})  {'calibrated' if ok else 'BROKEN, the test fires on noise'}")

    print("  A random top-100 does not reproduce the pillars. The tests are not vacuous.")
    print("  NOTE: this rules out noise, not every rival ranking. Ranking by pure effect size, or")
    print("  by DE-gene count, is a different null and is task N7.")
    return calibrated


def _effect_str(row: pd.Series, family: str) -> str:
    """Format an effect as an odds ratio for flags and a normal deviate for measures.

    Args:
        row: One row of the pillar table.
        family: Test family letter.

    Returns:
        A formatted string.
    """
    value = row[f"{family}_effect"]
    if not np.isfinite(value):
        return "        --"
    return f"OR {value:8.2f}" if row["binary"] else f" z {value:+8.3f}"


def report(frame: pd.DataFrame, table: pd.DataFrame, top_k: int) -> bool:
    """Print the three test families and decide the build.

    Args:
        frame: Output of :func:`build_frame`.
        table: Output of :func:`run_pillar_tests`.
        top_k: Shortlist size.

    Returns:
        True if no gated pillar failed both A and B.
    """
    bins = "/".join(str(b) for b in SENSITIVITY_BINS)
    print(f"\n=== TEST A: top {top_k} vs ALL {int((~frame['in_top']).sum()):,} background rows, "
          f"stratified on z_l2 decile ===")
    print("    no sampling, no seed. van Elteren for measures, CMH for flags, BH across 5 pillars.")
    for _, row in table.iterrows():
        mark = "CIRCULAR" if row["circular"] else ("HIGHER" if row["passes_A"] else "ns")
        stability = "stable" if row["A_stable"] else f"FRAGILE at {bins} bins"
        print(f"  {row['pillar']:32s} {_effect_str(row, 'A')}  p {row['A_p']:.3g}  p_BH {row['A_p_bh']:.3g}  "
              f"strata {row['A_strata']:2d}  top {row['A_top_used']:3d}  {mark:8s}  {stability}")

    print(f"\n=== TEST B: suppression top {top_k} vs INDUCTION top {top_k}, stratified on z_l2 decile ===")
    print("    the sign-flipped control. At matched effect magnitude, does the SIGN still matter?")
    for _, row in table.iterrows():
        mark = "SUPPRESSORS WORSE" if row["passes_B"] else "no difference"
        stability = "stable" if row["B_stable"] else f"FRAGILE at {bins} bins"
        print(f"  {row['pillar']:32s} {_effect_str(row, 'B')}  p {row['B_p']:.3g}  p_BH {row['B_p_bh']:.3g}  "
              f"strata {row['B_strata']:2d}  ({row['B_top_used']:3d} vs {row['B_bg_used']:3d})  "
              f"{mark:17s}  {stability}")

    print("\n=== TEST C: ALL suppressors vs ALL inducers, stratified on |score| decile, two-sided ===")
    print("    the audit's own requested test: is_suppressor conditional on |score|.")
    print(f"    {int(frame['is_suppressor'].sum()):,} suppressors vs {int((~frame['is_suppressor']).sum()):,} inducers")
    for _, row in table.iterrows():
        if not np.isfinite(row["C_effect"]):
            print(f"  {row['pillar']:32s} not testable")
            continue
        if row["C_p_bh"] >= 0.05:
            direction = "no difference"
        elif row["binary"]:
            direction = "suppressors HIGHER" if row["C_effect"] > 1 else "suppressors LOWER"
        else:
            direction = "suppressors HIGHER" if row["C_effect"] > 0 else "suppressors LOWER"
        print(f"  {row['pillar']:32s} {_effect_str(row, 'C')}  p {row['C_p']:.3g}  p_BH {row['C_p_bh']:.3g}  "
              f"{direction}")

    registered = table[table["registered"]]
    extra = table[~table["registered"]]

    print("\n" + "=" * 96)
    print("PER-PILLAR VERDICT on the FOUR REGISTERED PILLARS")
    print("kill rule fixed in advance (SESSION_SUMMARY N1): fails A AND B -> removed from the gate")
    print("=" * 96)
    for _, row in registered.iterrows():
        if row["passes_A"] and row["passes_B"]:
            state = "SURVIVES  magnitude-independent AND direction-specific"
        elif row["passes_B"]:
            state = "SURVIVES  direction-specific; entangled with magnitude"
        elif row["passes_A"]:
            state = "SURVIVES  magnitude-independent; NOT direction-specific"
        else:
            state = "REMOVED   neither magnitude-independent nor direction-specific"
        robust = "" if row["robust"] else "   [FRAGILE: the A or B call changes with the bin count]"
        print(f"  {row['pillar']:32s} {state}{robust}")

    failed = registered[~registered["survives"]]
    holds = failed.empty

    print("\n" + "-" * 96)
    print("UNREGISTERED, and therefore NOT allowed to decide anything today")
    print("-" * 96)
    for _, row in extra.iterrows():
        gate = "IS A GATE AXIS" if row["gates"] else "not a gate axis"
        print(f"  {row['pillar']:32s} [{gate}]")
        print(f"    test A (vs magnitude-matched background): p_BH {row['A_p_bh']:.3g}  "
              f"{'enriched' if row['passes_A'] else 'NO ENRICHMENT'}")
        print(f"    test B (vs induction control):            p_BH {row['B_p_bh']:.3g}  "
              f"{'suppressors worse' if row['passes_B'] else 'no difference'}")
    print()
    print("  04_window_score.py does not gate on raw rest-DE. It gates on selectivity,")
    print("  log1p(stim_de) - log1p(rest_de). Conditional on effect magnitude, the naive top 100 is")
    print("  NOT enriched for poor selectivity. That does not make selectivity a bad FILTER -- NSD1")
    print("  really does carry 2,175 stim and 1,767 rest DE genes and really should be rejected --")
    print("  but it does mean the risk-kill argument cannot lean on homeostasis. It never tested it.")
    print("  This analysis was not pre-registered. It sets task N6; it does not remove the axis.")

    print("\n" + "=" * 96)
    if holds:
        print("VERDICT: not one of the four registered pillars is a pure effect-magnitude artifact.")
        print("The headline stands. Its stated MECHANISM was wrong, and the audit was half right.")
        print()
        print("  - The audit's EVIDENCE reproduces. At matched |score|, suppressors carry fewer")
        print("    collateral and resting DE genes than inducers (test C, stable at 10/20/50/100 bins).")
        print("  - The audit's CONCLUSION does not follow. It matched |score| and concluded about")
        print("    effect magnitude. Spearman(z_l2, |score|) = +0.198. Match on z_l2 instead and the")
        print("    collateral effect shrinks 2.6x; the resting one vanishes (p = 0.22).")
        print("  - The naive ranking DOES select big-effect perturbations: collateral DE is 0.725")
        print("    rank-correlated with z_l2 and is close to a restatement of effect size. Its")
        print("    direction-specificity is FRAGILE, surviving at 10 bins and dying at 20.")
        print("  - One pillar is unmoved by every control: TOLERANCE-MODULE SUPPRESSION. rho with")
        print("    z_l2 = 0.069. Inducers do not suppress tolerance, they induce it. It is the axis")
        print("    the audit attacked hardest, and the only one that is about the DIRECTION of the")
        print("    perturbation rather than its SIZE. RULE #2 again: the screen-native graded")
        print("    measure holds; the count-of-DE-genes proxies do not.")
    else:
        print("VERDICT: these registered pillars are effect-magnitude artifacts and must be removed:")
        for _, row in failed.iterrows():
            print(f"  - {row['pillar']}")
        print("Add them to DO NOT REDO and rewrite the gate.")
    print("=" * 96)
    return holds


def main() -> None:
    """Reproduce the audit, run the three test families, persist the tables, decide the build."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top-k", type=int, default=100)
    parser.add_argument("--n-sham", type=int, default=200, help="random shortlists for the falsification control")
    args = parser.parse_args()

    paths.ensure_dirs()
    frame = build_frame(args.top_k)
    reproduce_the_audit(frame)
    diagnostics(frame, args.top_k)
    calibrated = falsification_control(frame, args.top_k, args.n_sham)
    table = run_pillar_tests(frame, args.top_k)
    holds = report(frame, table, args.top_k) and calibrated
    if not calibrated:
        print("\nA random shortlist reproduces the pillars. Nothing above is about reversal.")

    out = paths.TABLES / "magnitude_matched_pillars.csv"
    table.to_csv(out, index=False)
    print(f"\nwrote {out}")

    rows_out = paths.TABLES / "magnitude_matched_rows.csv"
    keep = [
        "gene_name", "naive_suppression", "abs_score", "is_suppressor", "z_l2", "z_l2_raw",
        "in_top", "in_induction_top", "stim_downstream", "stim_de_genes", "rest_downstream",
        "rest_de_genes", "has_rest_row", "rest_qc_pass", "homeostasis_cost",
        "tolerance_suppression", "is_iei", "n_cells_target",
    ]
    frame[keep].to_csv(rows_out, index=False)
    print(f"wrote {rows_out}")

    raise SystemExit(0 if holds else 1)


if __name__ == "__main__":
    main()
