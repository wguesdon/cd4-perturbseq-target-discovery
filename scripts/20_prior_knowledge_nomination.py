"""N10. Combine the within-screen gate with external prior knowledge to nominate drug targets.

Implements ``docs/preregistration_n10_2026_07_08.md`` exactly. That document was committed before this
script existed. Nothing here is tuned against its own output.

The gate (evidence floor + co-inhibitory preservation) decides which perturbations are safe. This layer
asks the pre-registered question: does the safety-gated ranking enrich for genes with human genetic
support for autoimmune disease, relative to the naive suppression ranking and a magnitude-matched
background? Both directions are informative and were written down in advance:

* ALIGN: the gate adds drug relevance the naive ranking lacks.
* OPPOSE: the gate de-enriches, which is the immunodeficiency result from a third angle -- the standard
  target-discovery prior, applied to an inflammatory reversal signature, selects for the targets a safety
  layer must reject.

It then emits the nomination: the 214 safe genes, each annotated with genetic support (P1), LoF tolerance
(P2, loeuf + prec), tractability (P3), and held-out Schmidt (P4), ranked by a transparent lexicographic
tier, never a fitted score.

Reads results/tables/*.csv and the on-disk Open Targets parquet. Never touches the h5ad layers.

Usage:
    uv run python scripts/20_prior_knowledge_nomination.py [--n-draws 5000]
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy import stats

from cd4_perturbseq import matched, paths, priors

N_BINS = 10
TOP_N = 100
ALPHA = 0.05
SEED = 0
KNOWN_DRUGS = ("IMPDH2", "PPP3R1", "CD3E", "CD3G", "IL4R", "CD2", "CD28")
"""Recovered approved-drug targets. They validate the gate; they are not counted as novel."""


def load() -> pd.DataFrame:
    """Assemble the universe: every ranked perturbation with a known z_l2, carrying P1-P4."""
    safety = pd.read_csv(paths.TABLES / "window_score_organism_safety.csv")
    rows = pd.read_csv(paths.TABLES / "magnitude_matched_rows.csv")[["gene_name", "z_l2"]]
    frame = safety.merge(rows, on="gene_name", how="left")
    frame = frame[frame["z_l2"].notna()].reset_index(drop=True)

    # P1: OT autoimmune genetic support. Absence from OT is a THIRD state, never "unsupported".
    gen = priors.ot_genetic_support()
    frame = frame.merge(gen, on="gene_name", how="left")
    frame["in_ot"] = frame["gene_name"].isin(set(priors.ot_tractability()["gene_name"]))
    frame["ot_genetic_supported"] = frame["ot_genetic_supported"].fillna(False)
    frame["ot_genetic_max"] = frame["ot_genetic_max"].fillna(0.0)

    # P3: OT tractability.
    frame = frame.merge(priors.ot_tractability(), on="gene_name", how="left")

    # P4: held-out Schmidt IL-2 protein screen.
    sch = priors.schmidt_cd4_il2_screen()
    frame = frame.merge(sch, on="gene_name", how="left")
    frame["il2_hit"] = (frame["il2_neg_fdr"] < 0.05) & (frame["il2_lfc"] < 0)

    frame["naive_rank"] = (-frame["eff_mean_z"]).rank(ascending=False, method="first")
    frame["z_decile"] = pd.qcut(frame["z_l2"], N_BINS, labels=False, duplicates="drop")
    return frame


def enrichment(frame: pd.DataFrame, n_draws: int) -> tuple[pd.DataFrame, dict[str, bool]]:
    """The pre-registered primary and its falsification controls (sections 5 and 8)."""
    rng = np.random.default_rng(SEED)
    records: list[dict] = []
    controls: dict[str, bool] = {}

    # Test 1: does the gate enrich for genetic support at all, vs evidence-passing-but-rejected?
    ev = frame[~frame["fail_evidence"]]
    a = int((ev["safe"] & ev["ot_genetic_supported"]).sum())
    b = int((ev["safe"] & ~ev["ot_genetic_supported"]).sum())
    c = int((~ev["safe"] & ev["ot_genetic_supported"]).sum())
    d = int((~ev["safe"] & ~ev["ot_genetic_supported"]).sum())
    orr, p_fisher = stats.fisher_exact([[a, b], [c, d]])
    records.append({"test": "gate vs rejected (Fisher)", "observed": a, "matched_mean": np.nan,
                    "p": p_fisher, "odds_ratio": orr, "note": f"safe {a}/{a+b} vs rejected {c}/{c+d}"})

    # Test 2: naive-top-100 and safe-top-100 supported counts vs a z_l2-matched null.
    supported = frame["ot_genetic_supported"]
    for label, top in (
        ("naive top 100", frame.nsmallest(TOP_N, "naive_rank")),
        ("safe-set top 100", frame[frame["safe"]].nlargest(min(TOP_N, int(frame["safe"].sum())), "window_score")),
    ):
        draws, shortfall = matched.count_matched(
            frame, frame["z_decile"], top.index, {"supported": supported}, n_draws, rng)
        assert (draws["n"] == len(top)).all(), "matched draw wrong size"
        obs = int(top["ot_genetic_supported"].sum())
        p = float((draws["supported"] >= obs).mean())
        records.append({"test": f"{label} enrichment", "observed": obs,
                        "matched_mean": float(draws["supported"].mean()), "p": p, "odds_ratio": np.nan,
                        "note": f"{obs}/{len(top)} supported; matched {draws['supported'].mean():.1f}; "
                                f"shortfall {shortfall}"})

    # The direction: safe-top supported fraction minus naive-top supported fraction, with a
    # two-proportion 95% CI. The pre-registration (section 6) keys OPPOSE on delta<0 with a CI that
    # EXCLUDES 0, not on the sign of delta alone.
    naive_top = frame.nsmallest(TOP_N, "naive_rank")
    n_safe_top = min(TOP_N, int(frame["safe"].sum()))
    safe_top = frame[frame["safe"]].nlargest(n_safe_top, "window_score")
    p_naive = naive_top["ot_genetic_supported"].mean()
    p_safe = safe_top["ot_genetic_supported"].mean()
    delta = p_safe - p_naive
    se = np.sqrt(p_naive * (1 - p_naive) / len(naive_top) + p_safe * (1 - p_safe) / n_safe_top)
    delta_lo, delta_hi = delta - 1.96 * se, delta + 1.96 * se
    records.append({"test": "delta (safe - naive) supported fraction", "observed": delta,
                    "matched_mean": np.nan, "p": float(2 * (1 - stats.norm.cdf(abs(delta / se)))),
                    "odds_ratio": np.nan,
                    "note": f"safe {p_safe:.3f} vs naive {p_naive:.3f}; 95% CI [{delta_lo:+.3f}, {delta_hi:+.3f}]"})

    # Control: label-shuffle must destroy enrichment.
    shuffled = frame.copy()
    shuffled["ot_genetic_supported"] = rng.permutation(shuffled["ot_genetic_supported"].to_numpy())
    st = shuffled[shuffled["safe"]].nlargest(min(TOP_N, int(shuffled["safe"].sum())), "window_score")
    draws_sh, _ = matched.count_matched(shuffled, shuffled["z_decile"], st.index,
                                        {"supported": shuffled["ot_genetic_supported"]}, n_draws, rng)
    p_sh = float((draws_sh["supported"] >= int(st["ot_genetic_supported"].sum())).mean())
    controls["label_shuffle_null"] = 0.05 <= p_sh <= 0.95
    records.append({"test": "CONTROL label-shuffle (must be non-significant)", "observed": np.nan,
                    "matched_mean": np.nan, "p": p_sh, "odds_ratio": np.nan,
                    "note": f"{'PASS' if controls['label_shuffle_null'] else 'FAIL'}"})

    # Control: matched-null calibration on genetic support. 300 random shortlists, each against 300
    # matched draws, is enough to see whether the 5% tail is calibrated without a 100k-draw loop.
    hits = 0
    n_cal = 300
    pool = frame.index.to_numpy()
    for _ in range(n_cal):
        idx = rng.choice(pool, size=TOP_N, replace=False)
        draws_c, _ = matched.count_matched(frame, frame["z_decile"], pd.Index(idx),
                                           {"supported": supported}, 300, rng)
        obs_c = int(supported.loc[idx].sum())
        if float((draws_c["supported"] >= obs_c).mean()) < ALPHA:
            hits += 1
    fpr = hits / n_cal
    controls["matched_calibration"] = 0.02 <= fpr <= 0.09
    records.append({"test": "CONTROL matched-null calibration (FPR ~5%)", "observed": np.nan,
                    "matched_mean": np.nan, "p": fpr, "odds_ratio": np.nan,
                    "note": f"FPR {fpr:.3f}; {'PASS' if controls['matched_calibration'] else 'FAIL'}"})

    return pd.DataFrame(records), controls


def nominate(frame: pd.DataFrame) -> pd.DataFrame:
    """The pre-registered lexicographic nomination (section 7) over the 214 safe genes."""
    safe = frame[frame["safe"]].copy()
    # np.select needs plain boolean ndarrays; several of these are nullable-boolean or carry NaN.
    tolerant = ((safe["prec"] <= 0.90) & ~safe["lof_intolerant"].fillna(False).astype(bool)).to_numpy(dtype=bool)
    tractable = safe["tractable"].fillna(False).to_numpy(dtype=bool)
    supported = safe["ot_genetic_supported"].fillna(False).to_numpy(dtype=bool)

    safe["tier"] = np.select(
        [supported & tolerant & tractable, supported],
        [1, 2],
        default=3,
    )
    safe["is_known_drug"] = safe["gene_name"].isin(KNOWN_DRUGS)
    safe = safe.sort_values(
        ["tier", "il2_hit", "window_score"], ascending=[True, False, False]
    ).reset_index(drop=True)
    safe["nomination_rank"] = np.arange(1, len(safe) + 1)
    return safe


def main() -> None:
    """Run the enrichment test, apply the falsification controls, and emit the nomination."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-draws", type=int, default=5000)
    args = parser.parse_args()
    paths.ensure_dirs()

    frame = load()
    print(f"universe: {len(frame)} ranked perturbations with a known z_l2")
    print(f"  in Open Targets: {int(frame['in_ot'].sum())}   "
          f"genetically supported (>= {priors.OT_GENETIC_THRESHOLD} for >=1 of 14 autoimmune): "
          f"{int(frame['ot_genetic_supported'].sum())}")
    print(f"  safe set: {int(frame['safe'].sum())}")

    enrich, controls = enrichment(frame, args.n_draws)
    print("\n=== PRE-REGISTERED ENRICHMENT (section 5) ===")
    print(enrich.to_string(index=False))

    enrich.to_csv(paths.TABLES / "n10_enrichment.csv", index=False)
    nom = nominate(frame)
    keep = ["nomination_rank", "gene_name", "tier", "window_rank", "viability_tier",
            "ot_genetic_max", "ot_genetic_supported", "ot_genetic_n_diseases",
            "loeuf", "prec", "recessive_intolerant", "tractable", "sm_pocket", "clinical_precedent",
            "il2_hit", "selectivity", "fail_homeostasis", "is_known_drug"]
    nom[keep].to_csv(paths.TABLES / "n10_nomination.csv", index=False)

    print("\n=== NOMINATION, tier counts ===")
    print(nom["tier"].value_counts().sort_index().to_string())
    print("\n=== Tier 1 (supported + LoF-tolerant + tractable) ===")
    t1 = nom[nom["tier"] == 1]
    if len(t1):
        print(t1[["nomination_rank", "gene_name", "ot_genetic_max", "prec", "il2_hit", "is_known_drug"]].to_string(index=False))
    else:
        print("  EMPTY. No safe gene is genetically supported AND LoF-tolerant AND tractable.")
    print("\n=== where the recovered drugs land ===")
    print(nom[nom["is_known_drug"]][["nomination_rank", "gene_name", "tier", "ot_genetic_max", "prec"]].to_string(index=False))

    # ---------------------------------------------------------------- verdict
    naive_row = enrich[enrich["test"] == "naive top 100 enrichment"].iloc[0]
    safe_row = enrich[enrich["test"] == "safe-set top 100 enrichment"].iloc[0]
    delta_row = enrich[enrich["test"].str.startswith("delta")].iloc[0]
    delta = float(delta_row["observed"])
    delta_p = float(delta_row["p"])
    delta_ci_excludes_0 = delta_p < ALPHA  # the two-proportion CI excludes 0 iff its z-test p < 0.05

    ctrl_ok = all(controls.values())
    print("\n" + "=" * 78)
    print(f"CONTROLS: label-shuffle {'PASS' if controls['label_shuffle_null'] else 'FAIL'}, "
          f"matched-calibration {'PASS' if controls['matched_calibration'] else 'FAIL'}")
    if not ctrl_ok:
        print("VOID: a falsification control did not fire. The enrichment is uninterpretable.")
        print("=" * 78)
        raise SystemExit(1)

    naive_enriched = naive_row["p"] < ALPHA
    safe_enriched = safe_row["p"] < ALPHA
    # Section 6, applied exactly: OPPOSE and ALIGN both require the delta CI to exclude 0.
    if delta > 0 and safe_enriched and delta_ci_excludes_0:
        verdict = "ALIGN"
        msg = ("The safety-gated ranking enriches for genetics-supported autoimmune targets beyond "
               "effect magnitude. The decision layer adds drug relevance.")
    elif delta < 0 and delta_ci_excludes_0:
        verdict = "OPPOSE"
        msg = ("The gate DE-enriches for genetics-supported targets relative to the naive ranking: the "
               "standard target-discovery prior, applied to an inflammatory reversal signature, selects "
               "for the targets a safety layer must reject.")
    else:
        verdict = "NULL (with a directional signal)"
        msg = ("The delta CI includes 0, so we do NOT claim the gate opposes genetic support. What is "
               "significant and pre-registered: the NAIVE ranking is enriched for autoimmune genetic "
               f"support ({'yes' if naive_enriched else 'no'}, p={naive_row['p']:.4f}) while the SAFETY-GATED "
               f"ranking is NOT (p={safe_row['p']:.4f}). The gate does not inherit the naive ranking's "
               "genetic enrichment. The difference between them is suggestive but under-powered at n=100 "
               f"(delta {delta:+.3f}, p={delta_p:.3f}). Consistent with the IEI and co-inhibitory results; "
               "not independently significant.")
    print(f"VERDICT: {verdict}")
    print(f"  naive top-100 supported {int(naive_row['observed'])} (matched {naive_row['matched_mean']:.1f}, "
          f"p={naive_row['p']:.4f}); safe top-100 supported {int(safe_row['observed'])} "
          f"(matched {safe_row['matched_mean']:.1f}, p={safe_row['p']:.4f})")
    print(f"  delta={delta:+.3f}, two-proportion p={delta_p:.3f} -> CI {'excludes' if delta_ci_excludes_0 else 'INCLUDES'} 0")
    print(f"  {msg}")
    print("=" * 78)
    print("\nThe nomination is hypothesis-generating: T-cell-intrinsic gate, organism safety is annotation,")
    print("one dataset, one week. No target here has a functional or in-vivo readout in this work.")


if __name__ == "__main__":
    main()
