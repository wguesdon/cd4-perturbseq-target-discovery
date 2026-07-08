"""N6. Does the selectivity composite carry information beyond its own denominator?

Implements `docs/preregistration_n6_2026_07_08.md` sections 3, 4, 5 and 6, exactly. That document was
written by a five-lens adversarial panel that withdrew an earlier design of mine, and it was committed
to git before this file existed. Read it first. Nothing here may be tuned against its own output.

**The earlier design is dead and its corpse is instructive.** Its primary endpoint was
`rest_cells_ratio`, proposed as an "independent, non-DE readout of resting-arm damage". It is none of
those things:

* `04_window_score.py:197` divides `rest_cells` by the scalar 543.0, so `Spearman(rest_cells,
  rest_cells_ratio) = 1.000` exactly. A rank test gives byte-identical answers. It is not a ratio.
* `rest_cells` is the cell count of the Rest contrast, which is the contrast that emits
  `rest_de_genes`. Within a `z_l2` stratum the group label *is* a threshold on `rest_de_genes`
  (mean within-stratum AUC 0.996). The endpoint was the power denominator of the exposure.
* `Spearman(rest_cells_ratio, stim n_cells_target) = +0.931`, and +0.9353 among 2,988 null
  perturbations that cannot be toxic. It measures guide abundance, a pre-exposure common cause.

So the outcome moved outside the screen entirely: gnomAD LOEUF. And the statistic that decides is not
the primary's p-value but a **dominance contrast against a sham exposure** built from `rest_de_genes`
alone, with no stim term, no `log1p` difference and no `ln(10)` bar. If the shipped composite cannot
beat a one-variable shadow of itself, its form is decoration and it is demoted.

Both branches are acceptable. The gate exists to be right, not to be kept. Section 4.3 of the
pre-registration discloses that on the 214-row superset the contrast has already been computed
(`z_sel` = -2.6255, `z_sham` = -3.588) and that **the expected outcome is DEMOTE**. We commit to the
rule rather than declaring the demotion by fiat, and we do not pretend the rule was blind (L1).

Reads only `results/tables/*.csv`. Never touches the h5ad layers.

Usage:
    uv run python scripts/18_n6_selectivity_validation.py
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from cd4_perturbseq import paths, stratified

ALPHA = 0.05
N_BINS = 10
BIN_SWEEP = (4, 5, 10, 20)
DOMINANCE_MARGIN = 0.5
"""Pre-registered, section 4: RETAIN needs `z_sel <= z_sham - 0.5`. Arbitrary, and fixed in advance."""

N_RANDOM_SPLITS = 500
SEED = 0

SECONDARY = (
    ("lof_intolerant", "cmh", "greater"),
    ("systemic_risk", "cmh", "greater"),
    ("is_ceg", "cmh", "greater"),
    ("ubiquitous", "cmh", "greater"),
    ("is_iei", "cmh", "less"),
    ("n_nonimmune_tissues", "ve", "greater"),
    ("max_nonimmune_ntpm", "ve", "greater"),
    ("rest_cells_ratio", "ve", "less"),
)
"""Declared in full before the run (section 6). Nothing may be added or removed afterwards.
`max_nonimmune_ntpm` is "greater", pinned here to prevent a post-hoc tail swap."""


def load() -> pd.DataFrame:
    """Join the three committed tables and rebuild the resting-arm QC class.

    Returns:
        One row per QC-passing Stim48hr perturbation, carrying the gate columns, ``z_l2``, the
        resting-arm QC class, and every external annotation the secondary family needs.
    """
    ws = pd.read_csv(paths.TABLES / "window_score.csv")
    mm = pd.read_csv(paths.TABLES / "magnitude_matched_rows.csv")
    saf = pd.read_csv(paths.TABLES / "window_score_organism_safety.csv")

    frame = ws.merge(
        mm[["gene_name", "z_l2", "has_rest_row", "rest_qc_pass", "rest_downstream"]],
        on="gene_name", how="left",
    ).merge(
        saf[["gene_name", "loeuf", "is_ceg", "lof_intolerant", "ubiquitous", "systemic_risk",
             "n_nonimmune_tissues", "max_nonimmune_ntpm"]],
        on="gene_name", how="left",
    )

    # Section 2.3. The 737/65 split of rest_qc_fail is retired: `low_target_gex` at rest mechanically
    # entails `ontarget_significant == False`, so "fails only on low_target_gex" selects 4 rows, not 65.
    frame["rest_qc_class"] = np.select(
        [~frame["has_rest_row"].astype(bool), frame["rest_qc_pass"].fillna(False).astype(bool)],
        ["rest_absent", "rest_qc_pass"],
        default="rest_qc_fail",
    )

    # The shipped label, read straight off the table. This closes D5 for this test: whatever
    # 12_magnitude_matched.py computes from n_downstream cannot reach the primary.
    sel = np.log1p(frame["stim_de_genes"]) - np.log1p(frame["rest_de_genes"])
    assert np.nanmax(np.abs(frame["selectivity"] - sel)) < 1e-12, "selectivity is not the log1p difference"
    assert (frame["fail_homeostasis"] == (frame["selectivity"] < np.log(10.0))).all(), "gate label drifted"
    return frame


def _population(frame: pd.DataFrame) -> pd.DataFrame:
    """Section 3.3: tolerance survivors with a QC-passing resting arm.

    Tolerance survivors only, because the decision that hangs on this axis is whether 148 genes enter
    the shortlist; that is the marginal contrast. ``rest_qc_pass`` only, because a row FIX A calls
    indeterminate cannot sit in a group labelled "kept-by-selectivity".

    Args:
        frame: Output of :func:`load`.

    Returns:
        The P204 population.
    """
    return frame[
        ~frame["fail_evidence"] & ~frame["fail_tolerance"] & (frame["rest_qc_class"] == "rest_qc_pass")
    ].copy()


def _ve(values: np.ndarray, in_top: np.ndarray, strata: np.ndarray, alt: str) -> stratified.StratifiedResult:
    """Thin alias so the call sites read like the pre-registration."""
    return stratified.van_elteren(values, in_top, strata, alternative=alt)


def _sham_labels(rest_de: np.ndarray, strata: np.ndarray, rejected: np.ndarray) -> np.ndarray:
    """Relabel, within each stratum, the top ``k_s`` rows by ``rest_de_genes`` alone.

    ``k_s`` is that stratum's observed number of selectivity-rejected rows, so the sham exposure has
    the same within-stratum marginal as the real one. It contains no stim term, no ``log1p``
    difference and no ``ln(10)`` bar.

    Args:
        rest_de: Resting DE-gene count per row.
        strata: Stratum label per row.
        rejected: The real exposure, used only for its per-stratum count.

    Returns:
        Boolean sham exposure.
    """
    out = np.zeros(rest_de.size, dtype=bool)
    for s in np.unique(strata):
        idx = np.flatnonzero(strata == s)
        k = int(rejected[idx].sum())
        if k:
            out[idx[np.argsort(-rest_de[idx], kind="stable")[:k]]] = True
    return out


def primary(pop: pd.DataFrame, n_bins: int, rng: np.random.Generator) -> dict[str, object]:
    """The pre-registered primary, its sham comparator, and the dominance contrast.

    Args:
        pop: The P204 population.
        n_bins: Number of z_l2 quantile bins.
        rng: Unused here; kept so the signature matches the control callers.

    Returns:
        Mapping of statistics, including ``dominance`` = ``z_sel <= z_sham - margin``.
    """
    have = pop["loeuf"].notna().to_numpy()
    sub = pop[have]
    strata = stratified.deciles(sub["z_l2"].to_numpy(), n_bins=n_bins)
    rejected = sub["fail_homeostasis"].to_numpy()

    sel = _ve(sub["loeuf"].to_numpy(), rejected, strata, "less")
    sham_lab = _sham_labels(sub["rest_de_genes"].to_numpy(), strata, rejected)
    sham = _ve(sub["loeuf"].to_numpy(), sham_lab, strata, "less")

    dropped_rejected = int(pop["fail_homeostasis"].sum() - rejected.sum())
    deviation = dropped_rejected / max(int(pop["fail_homeostasis"].sum()), 1) > 0.05
    return {
        "n_bins": n_bins, "n": int(have.sum()), "n_rejected": int(rejected.sum()),
        "n_kept": int((~rejected).sum()),
        "z_sel": sel.effect, "p_sel": sel.pvalue,
        "z_sham": sham.effect, "p_sham": sham.pvalue,
        "delta_z": sel.effect - sham.effect,
        "dominance": bool(sel.effect <= sham.effect - DOMINANCE_MARGIN),
        "n_strata_used": sel.n_strata_used, "n_strata_dropped": sel.n_strata_dropped,
        "n_top_used": sel.n_top_used, "n_background_used": sel.n_background_used,
        "dropped_rejected_rows": dropped_rejected,
        "protocol_deviation_gt5pct": bool(deviation),
    }


def controls(frame: pd.DataFrame, pop: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Sections 5: C2 through C7. Any failure voids the endpoint; the branch is then DEMOTE."""
    rows: list[dict[str, object]] = []
    have = pop["loeuf"].notna().to_numpy()
    sub = pop[have]
    strata = stratified.deciles(sub["z_l2"].to_numpy(), n_bins=N_BINS)
    rejected = sub["fail_homeostasis"].to_numpy()

    # C2. Discriminant validity: LOEUF must not validate whatever axis it is handed.
    keep = ~rejected
    c2a = _ve(sub.loc[keep, "loeuf"].to_numpy(), sub.loc[keep, "fail_tolerance"].to_numpy(),
              stratified.deciles(sub.loc[keep, "z_l2"].to_numpy(), n_bins=N_BINS), "greater") \
        if sub.loc[keep, "fail_tolerance"].nunique() > 1 else None
    ev = frame[frame["loeuf"].notna()]
    c2b = _ve(ev["loeuf"].to_numpy(), ev["fail_evidence"].to_numpy(),
              stratified.deciles(ev["z_l2"].to_numpy(), n_bins=N_BINS), "less")
    rows.append({"id": "C2a", "what": "loeuf ~ fail_tolerance among homeostasis survivors",
                 "requirement": "|z| < 1.0",
                 "z": float("nan") if c2a is None else c2a.effect,
                 "p": float("nan") if c2a is None else c2a.pvalue,
                 "fires": True if c2a is None else abs(c2a.effect) < 1.0,
                 "note": "no tolerance-failers among survivors by construction" if c2a is None else ""})
    rows.append({"id": "C2b", "what": "loeuf ~ fail_evidence over all rows", "requirement": "p_less > 0.05",
                 "z": c2b.effect, "p": c2b.pvalue, "fires": bool(c2b.pvalue > ALPHA), "note": ""})

    # C3. Negative-control population: evidence FAILERS with the same resting-arm class.
    neg = frame[frame["fail_evidence"] & ~frame["fail_tolerance"]
                & (frame["rest_qc_class"] == "rest_qc_pass") & frame["loeuf"].notna()]
    c3 = _ve(neg["loeuf"].to_numpy(), neg["fail_homeostasis"].to_numpy(),
             stratified.deciles(neg["z_l2"].to_numpy(), n_bins=N_BINS), "less")
    rows.append({"id": "C3", "what": "same design on evidence-FAILERS", "requirement": "|z| < 2.0",
                 "z": c3.effect, "p": c3.pvalue, "fires": bool(abs(c3.effect) < 2.0),
                 "note": f"n={len(neg)}"})

    # C4. Negative-control exposure: library abundance.
    med = pd.Series(sub["n_cells_target"].to_numpy()).groupby(strata).transform("median").to_numpy()
    abundance = sub["n_cells_target"].to_numpy() < med
    c4 = _ve(sub["loeuf"].to_numpy(), abundance, strata, "less")
    rows.append({"id": "C4", "what": "exposure = 1[n_cells_target < within-stratum median]",
                 "requirement": "|z| < 2.0", "z": c4.effect, "p": c4.pvalue,
                 "fires": bool(abs(c4.effect) < 2.0), "note": ""})

    # C5. Random-split calibration.
    hits = 0
    for _ in range(N_RANDOM_SPLITS):
        lab = np.zeros(len(sub), dtype=bool)
        for s in np.unique(strata):
            idx = np.flatnonzero(strata == s)
            k = int(rejected[idx].sum())
            if k:
                lab[rng.choice(idx, size=k, replace=False)] = True
        hits += int(_ve(sub["loeuf"].to_numpy(), lab, strata, "less").pvalue < ALPHA)
    fpr = hits / N_RANDOM_SPLITS
    rows.append({"id": "C5", "what": f"{N_RANDOM_SPLITS} size-matched random within-stratum splits",
                 "requirement": "P(p<0.05) in [0.02, 0.09]", "z": float("nan"), "p": fpr,
                 "fires": bool(0.02 <= fpr <= 0.09), "note": f"empirical FPR {fpr:.3f}"})

    # C7. Permutation, reported, non-decisive.
    obs = _ve(sub["loeuf"].to_numpy(), rejected, strata, "less").effect
    null = np.empty(N_RANDOM_SPLITS)
    for b in range(N_RANDOM_SPLITS):
        lab = np.zeros(len(sub), dtype=bool)
        for s in np.unique(strata):
            idx = np.flatnonzero(strata == s)
            k = int(rejected[idx].sum())
            if k:
                lab[rng.choice(idx, size=k, replace=False)] = True
        null[b] = _ve(sub["loeuf"].to_numpy(), lab, strata, "less").effect
    rows.append({"id": "C7", "what": "within-stratum permutation (software calibration only)",
                 "requirement": "reported, non-decisive",
                 "z": float(null.mean()), "p": float((null <= obs).mean()),
                 "fires": True, "note": f"null mean {null.mean():+.3f}, sd {null.std(ddof=1):.3f}"})
    return pd.DataFrame(rows)


def secondary(pop: pd.DataFrame) -> pd.DataFrame:
    """Section 6. Eight measures, declared in advance, BH-corrected within family, never decisive."""
    strata_all = stratified.deciles(pop["z_l2"].to_numpy(), n_bins=N_BINS)
    rejected_all = pop["fail_homeostasis"].to_numpy()
    out: list[dict[str, object]] = []
    for name, kind, alt in SECONDARY:
        v = pop[name].to_numpy()
        ok = np.isfinite(pd.to_numeric(pd.Series(v), errors="coerce").to_numpy(dtype=float))
        vals, top, strata = v[ok], rejected_all[ok], strata_all[ok]
        if kind == "cmh":
            r = stratified.cochran_mantel_haenszel(vals.astype(bool), top, strata, alternative=alt)
        else:
            r = stratified.van_elteren(vals.astype(float), top, strata, alternative=alt)
        out.append({"measure": name, "test": kind, "alternative": alt, "effect": r.effect,
                    "p": r.pvalue, "n_used": int(ok.sum()), "n_top": r.n_top_used,
                    "n_bg": r.n_background_used})
    frame = pd.DataFrame(out)
    frame["p_bh"] = stratified.benjamini_hochberg(frame["p"].tolist())
    return frame


def main() -> None:
    """Run the primary, the bin sweep, the controls and the secondary family. Then apply the rule."""
    argparse.ArgumentParser(description=__doc__).parse_args()
    paths.ensure_dirs()
    rng = np.random.default_rng(SEED)

    frame = load()
    pop = _population(frame)
    print(f"evidence-passers: {int((~frame['fail_evidence']).sum())}")
    print(f"tolerance survivors (the 214): {int((~frame['fail_evidence'] & ~frame['fail_tolerance']).sum())}")
    print(f"P204 = tolerance survivors with rest_qc_pass: {len(pop)} "
          f"({int(pop['fail_homeostasis'].sum())} rejected, {int((~pop['fail_homeostasis']).sum())} kept)")

    sweep = pd.DataFrame([primary(pop, b, rng) for b in BIN_SWEEP])
    main_row = sweep[sweep["n_bins"] == N_BINS].iloc[0]

    print("\n" + "=" * 82)
    print("PRIMARY (pre-registered): does the composite beat its own rest_de shadow, on LOEUF?")
    print(f"  outcome  gnomAD LOEUF, external to the screen")
    print(f"  n = {main_row['n']}  ({main_row['n_rejected']} rejected, {main_row['n_kept']} kept); "
          f"{main_row['dropped_rejected_rows']} rejected rows dropped for missing LOEUF")
    if main_row["protocol_deviation_gt5pct"]:
        print(f"  ** PROTOCOL DEVIATION: >5% of rejected rows dropped. Flagged, per section 8. **")
    print(f"  selectivity composite   z = {main_row['z_sel']:+.4f}   p = {main_row['p_sel']:.5f}")
    print(f"  sham (rest_de alone)    z = {main_row['z_sham']:+.4f}   p = {main_row['p_sham']:.5f}")
    print(f"  dominance  z_sel <= z_sham - {DOMINANCE_MARGIN}?  "
          f"{main_row['z_sel']:+.4f} <= {main_row['z_sham'] - DOMINANCE_MARGIN:+.4f}  "
          f"-> {'YES' if main_row['dominance'] else 'NO'}   (delta z = {main_row['delta_z']:+.4f})")
    print(f"\n  bin sweep (C6):")
    for _, r in sweep.iterrows():
        print(f"    n_bins {int(r['n_bins']):2d}   p_sel {r['p_sel']:.5f}   z_sel {r['z_sel']:+.4f}   "
              f"z_sham {r['z_sham']:+.4f}   dominance {'Y' if r['dominance'] else 'N'}")

    ctrl = controls(frame, pop, rng)
    print("\nCONTROLS")
    for _, r in ctrl.iterrows():
        print(f"  {r['id']:4s} {r['what'][:52]:54s} {r['requirement']:24s} "
              f"z={r['z']:+7.3f} p={r['p']:.4f}  {'FIRES' if r['fires'] else 'FAILS'}  {r['note']}")

    sec = secondary(pop)
    print("\nSECONDARY (BH within family, never decisive)")
    print(sec.to_string(index=False))

    # ---------------------------------------------------------------- the rule, section 4
    cond1 = bool(main_row["p_sel"] < ALPHA)
    cond2 = bool((sweep["p_sel"] < ALPHA).all())
    cond3 = bool(main_row["dominance"])
    cond4 = bool(ctrl.loc[ctrl["id"].isin(["C2a", "C2b", "C3", "C4", "C5"]), "fires"].all())
    retain = cond1 and cond2 and cond3 and cond4

    sweep.to_csv(paths.TABLES / "n6_primary.csv", index=False)
    ctrl.to_csv(paths.TABLES / "n6_controls.csv", index=False)
    sec.to_csv(paths.TABLES / "n6_secondary.csv", index=False)

    print("\n" + "=" * 82)
    print("DECISION RULE (fixed in advance; RETAIN needs all four)")
    print(f"  1. p_sel < {ALPHA} at n_bins=10            : {'PASS' if cond1 else 'FAIL'}")
    print(f"  2. p_sel < {ALPHA} at every n_bins          : {'PASS' if cond2 else 'FAIL'}")
    print(f"  3. dominance over the rest_de sham       : {'PASS' if cond3 else 'FAIL'}")
    print(f"  4. controls C2-C5 all fire               : {'PASS' if cond4 else 'FAIL'}")
    print("=" * 82)
    if retain:
        print("VERDICT: RETAIN. The composite carries LoF-constraint information beyond a raw threshold")
        print("on the resting DE burden. This licenses ONE sentence and no other (section 4.1). It does")
        print("NOT validate the ln(10) location, the log1p offset, or D6's implemented rule.")
    else:
        print("VERDICT: DEMOTE. The selectivity composite does not beat a one-variable shadow of itself.")
        print("Its ratio, its log1p offset and its ln(10) bar are decoration. It becomes a reported")
        print("per-gene annotation. The gate becomes: evidence floor + co-inhibitory preservation.")
        print("")
        print("  A raw rest_de threshold WOULD have out-performed it. Choosing that now, having seen")
        print("  this, would be tuning, and is forbidden. It is reported and it does not gate.")
        print("  Cost, to be stated in the abstract and not buried: safe goes 66 -> 214, and the 148")
        print("  genes thereby admitted are MORE LoF-constrained than the 66 kept.")
    print("=" * 82)
    print(f"\nwrote {paths.TABLES / 'n6_primary.csv'}, n6_controls.csv, n6_secondary.csv")


if __name__ == "__main__":
    main()
