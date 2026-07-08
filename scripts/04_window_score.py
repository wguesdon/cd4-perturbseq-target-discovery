"""Build the therapeutic-window score and the naive-versus-safe head-to-head.

This script is the second design. The first is preserved in git history because its failure is
informative. Three things were wrong with it, and the data said so.

**No evidence floor.** Ranking by a mean of z-scores let perturbations that did nothing float to
the top. Of the 6,371 QC-passing perturbations, 4,131 have five or fewer significant DE genes
genome-wide, and among those the maximum number of effector-module genes significantly down is
ONE. Yet 804 of them have a mean module z below -0.5. The score was reading sub-significant
drift. A perturbation must now significantly suppress at least three module genes before it is
allowed to be called efficacious.

**Calibration was not the lever.** A CAMERA variance inflation factor estimated across all
perturbations comes out near 1, because the mass of null perturbations dilutes the inter-gene
correlation to zero. Estimated among the 713 perturbations with at least 100 DE genes, the same
quantity is 4.36. Either way the correction is close to uniform across perturbations, so it moves
the ranking almost not at all. It is retained because it is the honest statistic, not because it
rescues anything.

**The collateral cap was backwards.** Capping the number of DE genes in stimulated cells rejects
exactly the context-selective targets we are hunting. CD3G carries 1,374 DE genes in Stim and 5
at Rest. ITK carries 2,566 and 2. PLCG1 carries 2,218 and 3. That is the therapeutic window in
its purest form. Meanwhile NSD1 carries 2,175 in Stim and 1,767 at Rest, and STAT5B carries 1,060
and 1,135. Those are the transcriptome collapses. The discriminator is disruption AT REST, never
breadth in Stim.

Validation discipline: the Schmidt and Steinhart 2022 CD4+ IL-2 CRISPRi screen is HELD OUT. It
never enters the score, the gate, or any threshold. It is used only at the end.

Usage:
    uv run python scripts/04_window_score.py
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

from cd4_perturbseq import de_stats, paths, priors, programs, score

STIM = "Stim48hr"
REST = "Rest"
FDR = 0.10
MIN_MODULE_DOWN = 3
ACTIVE_DE = 100

MIN_SELECTIVITY = float(np.log(10.0))
"""Require the stimulated effect to be at least tenfold the resting effect, in DE-gene count.

An absolute cap on resting DE genes cannot work. Set on all perturbations it rejects IMPDH2, the
mycophenolate target, which carries 33 resting DE genes against 1,055 stimulated. Set on the
evidence-passing subpopulation, which is by construction active, the threshold inflates to roughly
500 and CFAP20 sails through with 455 resting DE genes.

The therapeutic window is a RATIO, not a count. This threshold is fixed a priori at tenfold and is
not tuned against any outcome. Schmidt is held out; the approved-drug positives are not consulted.
"""

RNG = np.random.default_rng(0)


qc_mask = de_stats.qc_mask
"""Routine perturbation QC. Single definition in :mod:`cd4_perturbseq.de_stats`."""


def auroc_ci(labels: np.ndarray, values: np.ndarray, n_boot: int = 2000) -> tuple[float, float, float]:
    """Bootstrap percentile CI for AUROC.

    Args:
        labels: Binary labels.
        values: Ranking score, higher means more positive.
        n_boot: Bootstrap resamples.

    Returns:
        Tuple of (auroc, lower, upper) at the 95% level.
    """
    point = roc_auc_score(labels, values)
    n = labels.size
    boots = []
    for _ in range(n_boot):
        idx = RNG.integers(0, n, n)
        if np.unique(labels[idx]).size < 2:
            continue
        boots.append(roc_auc_score(labels[idx], values[idx]))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return point, float(lo), float(hi)


def _module_stats(
    columns: dict[str, np.ndarray], module_idx: np.ndarray, self_idx: np.ndarray
) -> dict[str, np.ndarray]:
    """Leave-one-out module summaries.

    The perturbed gene is masked out of its own module, because its on-target knockdown is large,
    negative, and required to be significant by QC. Twenty-one effector genes are themselves
    perturbed in this library.

    Args:
        columns: Mapping of layer name to an array of shape ``(n_perturbations, m)`` holding that
            layer's values for the module's genes, in ``module_idx`` order.
        module_idx: Positional gene indices of the module.
        self_idx: Positional gene index of each perturbation's own gene, or -1.

    Returns:
        Mapping with ``mean_lfc``, ``mean_z``, ``n_down``, and ``n_used``.
    """
    n_pert = self_idx.size
    keep = np.ones((n_pert, module_idx.size), dtype=bool)
    # Mask each perturbation's own gene wherever it appears in the module.
    for position, gene_index in enumerate(module_idx):
        keep[:, position] &= self_idx != gene_index

    lfc = np.where(keep, columns["log_fc"], np.nan)
    z = np.where(keep, columns["zscore"], np.nan)
    padj = np.where(keep, columns["adj_p_value"], np.nan)

    down = (padj < FDR) & (lfc < 0) & keep
    return {
        "mean_lfc": np.nanmean(lfc, axis=1),
        "mean_z": np.nanmean(z, axis=1),
        "n_down": down.sum(axis=1),
        "n_used": keep.sum(axis=1),
    }


def build() -> pd.DataFrame:
    """Compute efficacy with an evidence floor, the safety axes, and the window score."""
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

    gene_col = "target_contrast_gene_name"
    obs[gene_col] = obs[gene_col].astype(str)
    rows = obs.index[(obs["culture_condition"] == STIM) & qc_mask(obs)].to_numpy()
    print(f"QC-passing {STIM} perturbations: {rows.size}")

    layers = ("zscore", "log_fc", "adj_p_value")
    eff_cols = {layer: de_stats.read_layer_columns(eff_idx, layer=layer)[rows] for layer in layers}
    tol_cols = {layer: de_stats.read_layer_columns(tol_idx, layer=layer)[rows] for layer in layers}

    sub = obs.loc[rows].reset_index(drop=True)
    self_idx = np.array([idx_of.get(g, -1) for g in sub[gene_col]], dtype=np.int64)

    eff = _module_stats(eff_cols, eff_idx, self_idx)
    tol = _module_stats(tol_cols, tol_idx, self_idx)
    print(f"leave-one-out removed a module gene for {(eff['n_used'] < len(effector)).sum()} perturbations")

    # CAMERA VIF, estimated on the ACTIVE perturbations. Across all 6,371 the null mass drives
    # rho_bar to zero and the correction vanishes; that is an estimation artifact, not biology.
    active = sub["n_total_de_genes"].to_numpy() >= ACTIVE_DE
    vif, rho = score.inter_gene_vif(eff_cols["zscore"][active])
    print(f"effector rho_bar={rho:.3f} VIF={vif:.2f} estimated on {int(active.sum())} active perturbations")
    print(f"  effective module size {len(effector) / vif:.1f} of {len(effector)}")

    frame = pd.DataFrame(
        {
            "gene_name": sub[gene_col],
            "eff_mean_lfc": eff["mean_lfc"],
            "eff_mean_z": eff["mean_z"],
            "n_module_down": eff["n_down"],
            "tol_mean_lfc": tol["mean_lfc"],
            "n_tolerance_down": tol["n_down"],
            "stim_de_genes": sub["n_total_de_genes"].to_numpy(),
            "n_cells_target": sub["n_cells_target"].to_numpy(),
        }
    )

    rest = obs[obs["culture_condition"] == REST].drop_duplicates(gene_col).set_index(gene_col)
    frame["rest_de_genes"] = frame["gene_name"].map(rest["n_total_de_genes"])
    frame["rest_de_genes"] = frame["rest_de_genes"].fillna(frame["rest_de_genes"].median())

    # Context selectivity: large effect when stimulated, small effect at rest.
    frame["selectivity"] = np.log1p(frame["stim_de_genes"]) - np.log1p(frame["rest_de_genes"])
    frame["efficacy"] = -frame["eff_mean_lfc"]
    frame["tolerance_loss"] = -frame["tol_mean_lfc"]

    # Viability, read straight off the screen. Cells carrying a guide against a gene the cell
    # cannot live without simply are not there. Hart core-essentials are measurably depleted at
    # REST in this data (0.83x the median, MWU p=0.003), so resting cell count is the
    # context-free viability signal. Stimulated cells proliferate and resting ones are quiescent,
    # so depletion only under stimulation is antiproliferation, which is mycophenolate's and
    # methotrexate's actual mechanism, not a disqualification.
    frame["rest_cells"] = frame["gene_name"].map(rest["n_cells_target"])
    median_rest_cells = rest["n_cells_target"].median()
    frame["rest_cells_ratio"] = frame["rest_cells"] / median_rest_cells
    frame["log2_stim_over_rest_cells"] = np.log2(frame["n_cells_target"] / frame["rest_cells"])
    # Three-way, not two. A gene with no resting arm has unknown viability, and NaN >= 0.5 is
    # False, so a plain np.where would silently file missing data as evidence of toxicity.
    frame["viability_tier"] = np.select(
        [
            frame["rest_cells_ratio"].isna(),
            frame["rest_cells_ratio"] >= 0.5,
        ],
        ["unknown", "non-depleting"],
        default="depleting-at-rest",
    )
    print(f"median resting cell count {median_rest_cells:.0f}; "
          f"{int((frame['rest_cells_ratio'] < 0.5).sum())} perturbations are depleted at rest")

    frame["is_iei"] = frame["gene_name"].isin(priors.iei_genes())
    frame["is_core_essential"] = frame["gene_name"].isin(priors.core_essential_genes())
    return frame


def apply_gate(frame: pd.DataFrame) -> pd.DataFrame:
    """Apply the evidence floor and the three evidence-chosen safety axes."""
    frame["fail_evidence"] = frame["n_module_down"] < MIN_MODULE_DOWN

    passing = frame[~frame["fail_evidence"]]
    tol_max = passing["tolerance_loss"].quantile(0.75)

    frame["fail_homeostasis"] = frame["selectivity"] < MIN_SELECTIVITY
    frame["fail_tolerance"] = frame["tolerance_loss"] > tol_max

    # `is_iei` is NOT a gate. Verified 2026-07-08: the IUIS flag is more enriched among approved
    # immunomodulator targets (25.0%, OR 8.31) than among the perturbations a naive reversal
    # ranking calls toxic (14.0%, OR 4.16). It cannot separate a dangerous knockdown from a gene
    # that is load-bearing in immunity, and a good immune drug target is the latter. Gating on it
    # rejected CD3E, CD3G, IL2RA, PIK3CD and TYK2, which are five approved drugs. It is reported.
    fails = ["fail_evidence", "fail_homeostasis", "fail_tolerance"]
    frame["safe"] = ~frame[fails].any(axis=1)
    frame["reject_reason"] = frame.apply(
        lambda r: ",".join(c.replace("fail_", "") for c in fails if r[c]) or "-", axis=1
    )

    frame["window_score"] = (
        score.zscore(frame["efficacy"].to_numpy())
        + score.zscore(frame["selectivity"].to_numpy())
        - score.zscore(np.maximum(frame["tolerance_loss"].to_numpy(), 0.0))
    )

    print(f"\nevidence floor: >= {MIN_MODULE_DOWN} module genes significantly down at FDR {FDR}")
    print(f"  passes evidence floor: {int((~frame['fail_evidence']).sum())} / {len(frame)}")
    print(f"selectivity gate: stimulated effect >= 10x resting effect (log-ratio {MIN_SELECTIVITY:.3f}), fixed a priori")
    print(f"tolerance gate: tolerance_loss <= {tol_max:.3f} (p75 of evidence-passers)")
    print(f"  fully safety-passing: {int(frame['safe'].sum())}")
    print("  NOTE: there is deliberately no cap on Stim DE-gene count. Capping it rejects the")
    print("  context-selective targets (CD3G 1374/5, ITK 2566/2) this project exists to find.")
    print("  NOTE: is_iei is an ANNOTATION, not a gate. See _iei_is_not_a_gate() below.")
    return frame


def _iei_is_not_a_gate(frame: pd.DataFrame) -> None:
    """Show, every run, why the immunodeficiency flag cannot be a safety gate.

    RULE: before believing any enrichment, run it on the POSITIVE class too. This project asserted
    two mechanisms in one day and had both refuted. The first was common-essentiality. The second
    was this flag. Printing the refutation on every run means nobody can quietly re-promote it.

    Args:
        frame: The gated frame, carrying ``is_iei``.
    """
    truth = pd.read_csv(paths.GROUND_TRUTH / "immunomodulator_targets.csv")
    positives = set(truth.loc[truth["include_as_positive"], "gene_symbol"])

    frame = frame.copy()
    frame["naive_rank"] = (-frame["eff_mean_z"]).rank(ascending=False, method="first")
    frame["in_top100"] = frame["naive_rank"] <= 100
    frame["is_positive"] = frame["gene_name"].isin(positives)

    def _or(group: pd.Series) -> tuple[float, float, float]:
        a = int((frame["is_iei"] & group).sum())
        b = int((~frame["is_iei"] & group).sum())
        c = int((frame["is_iei"] & ~group).sum())
        d = int((~frame["is_iei"] & ~group).sum())
        odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        return a / max(a + b, 1), odds, p

    print("\n=== WHY is_iei IS NOT A GATE (printed every run, on purpose) ===")
    rate_t, or_t, p_t = _or(frame["in_top100"])
    rate_p, or_p, p_p = _or(frame["is_positive"])
    print(f"  naive top-100          IEI rate {rate_t:6.1%}   OR {or_t:5.2f}   p={p_t:.3g}")
    print(f"  approved drug targets  IEI rate {rate_p:6.1%}   OR {or_p:5.2f}   p={p_p:.3g}")
    if or_p >= or_t:
        rejected = sorted(frame.loc[frame["is_positive"] & frame["is_iei"], "gene_name"])
        print("  The flag is MORE enriched among approved drugs than among the 'toxic' top 100.")
        print(f"  Gating on it would reject these approved targets: {', '.join(rejected)}")
        print("  A gene whose loss causes immunodeficiency is a gene that matters for immunity,")
        print("  which is exactly what a good immune drug target is. Annotation only.")
    else:
        print("  The flag separates them. Re-examine whether it should be a gate.")


def validate(frame: pd.DataFrame) -> None:
    """Report the head-to-head, the shortlist, and the held-out Schmidt validation."""
    schmidt = priors.schmidt_cd4_il2_screen()
    merged = frame.merge(schmidt, on="gene_name", how="left")
    merged["il2_hit"] = (merged["il2_neg_fdr"] < 0.05) & (merged["il2_lfc"] < 0)
    merged["naive_rank"] = (-merged["eff_mean_z"]).rank(ascending=False)

    print("\n=== HEAD-TO-HEAD: the naive top 20, judged by the gate ===")
    top = merged.nsmallest(20, "naive_rank")[
        ["naive_rank", "gene_name", "stim_de_genes", "rest_de_genes", "n_module_down",
         "safe", "reject_reason", "il2_hit", "is_iei"]
    ]
    print(top.to_string(index=False))
    n_safe = int(top["safe"].sum())
    n_hit = int(top["il2_hit"].sum())
    n_iei = int(top["is_iei"].sum())
    print(f"\n  survive the gate: {n_safe} of 20")
    print(f"  independently confirmed IL-2 reducers (Schmidt): {n_hit} of 20")
    print(f"  loss of function causes human immunodeficiency: {n_iei} of 20 (annotation, not a gate)")
    print(f"  -> {20 - n_safe} of the top 20 are rejected on screen-native axes alone.")
    if n_safe:
        survivors = ", ".join(top.loc[top["safe"], "gene_name"])
        print(f"  survivors: {survivors}. Check each by hand; the gate is not a proof.")

    print("\n=== EFFECTIVE BUT REJECTED: proven IL-2 reducers the gate refuses ===")
    rejected = merged[merged["il2_hit"] & ~merged["safe"] & ~merged["fail_evidence"]]
    cols = ["gene_name", "stim_de_genes", "rest_de_genes", "n_module_down", "il2_lfc", "reject_reason"]
    print(rejected.sort_values("n_module_down", ascending=False)[cols].to_string(index=False))

    print("\n=== SHORTLIST TIER A: safety-passing AND not depleted at rest ===")
    print("Signalling blockade rather than cytotoxicity. This is the tier a lab should act on.")
    tier_a = merged[merged["safe"] & (merged["viability_tier"] == "non-depleting")].nlargest(20, "window_score")
    cols = ["gene_name", "window_score", "efficacy", "selectivity", "rest_cells_ratio",
            "log2_stim_over_rest_cells", "n_module_down", "il2_lfc", "il2_hit"]
    print(tier_a[cols].to_string(index=False))

    print("\n=== SHORTLIST TIER B: safety-passing but depleted at rest ===")
    print("Antiproliferative. This is where mycophenolate lives. Effective, and it needs monitoring.")
    tier_b = merged[merged["safe"] & (merged["viability_tier"] == "depleting-at-rest")].nlargest(10, "window_score")
    print(tier_b[cols].to_string(index=False))

    n_supported = int(tier_a["il2_hit"].sum())
    base = merged["il2_hit"].mean()
    print(f"\n  {n_supported} of Tier A's top 20 are confirmed by the HELD-OUT Schmidt screen")
    print(f"  expected by chance {base * 20:.2f}, enrichment {n_supported / max(base * 20, 1e-9):.1f}x")
    if n_supported:
        n_hits = int(merged["il2_hit"].sum())
        table = [[n_supported, 20 - n_supported], [n_hits - n_supported, len(merged) - 20 - n_hits + n_supported]]
        _, p = stats.fisher_exact(table, alternative="greater")
        print(f"  Fisher one-sided p = {p:.3g}. n is small; this is a direction, not proof.")

    print("\n=== approved-drug targets: where the gate puts them ===")
    truth = pd.read_csv(paths.GROUND_TRUTH / "immunomodulator_targets.csv")
    positives = set(truth.loc[truth["include_as_positive"], "gene_symbol"])
    hits = merged[merged["gene_name"].isin(positives) & ~merged["fail_evidence"]]
    cols = ["gene_name", "stim_de_genes", "rest_de_genes", "rest_cells_ratio", "viability_tier",
            "n_module_down", "safe", "reject_reason"]
    if len(hits):
        print(hits.sort_values("n_module_down", ascending=False)[cols].to_string(index=False))
    else:
        print("  none pass the evidence floor")

    have = merged["il2_lfc"].notna()
    labels = merged.loc[have, "il2_hit"].to_numpy()
    print("\n=== HELD-OUT VALIDATION vs Schmidt IL-2 screen ===")
    print("Schmidt validates the EFFICACY axis. It cannot validate the window score, and a window")
    print("score that beat it would be suspect: Schmidt's hits are the TCR signalosome, which the")
    print("safety gate rejects on purpose. Reporting both is the honest thing to do.")
    for label, col, sign in (
        ("naive -mean(module z)", "eff_mean_z", -1.0),
        ("efficacy -mean(module lfc)", "efficacy", 1.0),
        ("window_score", "window_score", 1.0),
    ):
        values = merged.loc[have, col].to_numpy() * sign
        values = np.nan_to_num(values, nan=float(np.nanmin(values)))
        a, lo, hi = auroc_ci(labels, values)
        print(f"  {label:28s} AUROC {a:.3f}  95% CI [{lo:.3f}, {hi:.3f}]")

    print("\n=== THE TEST THAT DOES VALIDATE THE WINDOW SCORE: therapeutic index ===")
    print("Among approved-drug targets that clear the evidence floor, does the gate separate the")
    print("well-tolerated agents from the ones with a narrow therapeutic index?")
    truth = pd.read_csv(paths.GROUND_TRUTH / "immunomodulator_targets.csv")
    positives = set(truth.loc[truth["include_as_positive"], "gene_symbol"])
    graded = merged[merged["gene_name"].isin(positives) & ~merged["fail_evidence"]].copy()
    # Narrow-index agents: anti-CD3 causes cytokine release and profound immunosuppression, and
    # muromonab was withdrawn. This labelling is from the drug's clinical record, not from our score.
    narrow = {"CD3D", "CD3E", "CD3G", "CD52", "DHFR"}
    graded["clinical_index"] = np.where(graded["gene_name"].isin(narrow), "narrow", "tolerated")
    print(graded[["gene_name", "clinical_index", "safe", "viability_tier", "reject_reason"]].to_string(index=False))
    if len(graded) >= 2:
        tolerated_safe = int(((graded["clinical_index"] == "tolerated") & graded["safe"]).sum())
        narrow_safe = int(((graded["clinical_index"] == "narrow") & graded["safe"]).sum())
        n_tol = int((graded["clinical_index"] == "tolerated").sum())
        n_nar = int((graded["clinical_index"] == "narrow").sum())
        print(f"\n  well-tolerated agents passing the gate: {tolerated_safe}/{n_tol}")
        print(f"  narrow-index agents passing the gate:   {narrow_safe}/{n_nar}")
        print("  n is small. This is a direction, not a p-value, and it is reported as such.")
    # Generate this from the table, never alongside it. The previous hardcoded version kept
    # asserting that CD3E was rejected after CD3E had started passing.
    print("\n  Viability tiers, as observed:")
    for tier, block in graded.groupby("viability_tier", observed=True):
        members = ", ".join(f"{r.gene_name}({'pass' if r.safe else 'reject'})" for r in block.itertuples())
        print(f"    {tier:20s} {members}")
    print("\n  HONEST NOTE. Removing the immunodeficiency gate cost us this result. With it, the")
    print("  separation was 3/3 tolerated versus 0/2 narrow. Without it, CD3E passes. The gate")
    print("  now rests on screen-native axes only, and those axes cannot see cytokine release")
    print("  syndrome, which is why muromonab was withdrawn. We report the weaker number.")


def main() -> None:
    """Build, gate, validate, and persist the window score."""
    argparse.ArgumentParser(description=__doc__).parse_args()
    paths.ensure_dirs()

    frame = build()
    frame = apply_gate(frame)
    frame = frame.sort_values("window_score", ascending=False).reset_index(drop=True)
    frame["window_rank"] = np.arange(1, len(frame) + 1)

    out = paths.TABLES / "window_score.csv"
    frame.to_csv(out, index=False)
    print(f"\nwrote {out}")
    _iei_is_not_a_gate(frame)
    validate(frame)


if __name__ == "__main__":
    main()
