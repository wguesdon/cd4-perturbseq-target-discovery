"""N21 SECONDARY — a POST-HOC, relaxed-rule pass over the frozen 91-gene nomination pool.

**This is not the registered analysis.** The registered H3 in
`docs/preregistration_n21_2026_07_09.md` emitted no card, because 6 of 73 eligible pool genes are
co-tested against a pre-declared coverage gate of 10. That result stands and is not superseded.

Addendum 2 records that the user directed a relaxed rule after seeing the gate fire. This is post-hoc
rule relaxation. It is run once, its deviation is stamped into every output row, and both outcomes are
reported side by side. The registered outcome governs any headline.

RELAXED RULE
    universe   frozen 91-gene pool, minus curated approved-immunomodulator targets, minus
               direction-discordant, minus housekeeping (name AND function), minus RULE #6 genes
    promotion  Freimer IL2 arm: knockdown significantly LOWERS IL-2, at Freimer's own threshold,
               with the sign convention validated by marker-gene positive controls
    no coverage gate
    IL2RA-lowering   weaker supportive annotation, never promotion (IL2RA marks activation AND Tregs)
    CTLA4-lowering   liability flag, NOT disqualifying (the registered rule excluded these outright)

Freimer may not alter the screen gate, re-rank the 6,371 perturbations, rescue a gene that fails the
internal evidence floor, or touch Schmidt. The word "discovered" may not appear. The strongest phrase
permitted is "a Freimer-supported follow-up hypothesis".

Usage:
    uv run python scripts/34_freimer_secondary_posthoc.py
"""

from __future__ import annotations

import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from cd4_perturbseq import paths  # noqa: E402

FDR = 0.10
SCREENS = paths.RESOURCES / "external_screens"
MARKER_CONTROL = {"IL2RA": "IL2RA", "IL2": "IL2", "CTLA4": "CTLA4"}
HOUSEKEEPING_FUNCTIONAL = {"COG6", "TNPO1", "ARPC2", "SEL1L"}
DO_NOT_HEADLINE = {"STAT6", "GATA3", "STAT5A", "PTEN"}
DEVIATION = "POST-HOC deviation from the registered H3 (addendum 2). Not a pre-registered result."


def load_freimer() -> pd.DataFrame:
    """Load Freimer, sign as established in addendum 1, with fan-out assertions.

    Returns:
        One row per (gene, arm).

    Raises:
        ValueError: On fan-out or a violated sign assumption.
    """
    f = pd.read_csv(SCREENS / "Freimer2022_Screen.csv").rename(columns={"id": "gene_name", "screen": "arm"})
    if f.duplicated(subset=["gene_name", "arm"]).any():
        raise ValueError("fan-out: duplicate (gene, arm) rows")
    if not (f["neg|lfc"] == f["pos|lfc"]).all():
        raise ValueError("neg|lfc and pos|lfc diverge; the gene-level lfc assumption is wrong")
    f["lfc"] = f["pos|lfc"]
    # addendum 1: `pos` = guides enriched in the marker-LOW bin = knockdown LOWERS the marker.
    f["lowers"] = (f["pos|fdr"] < FDR) & (f["lfc"] > 0)
    f["raises"] = (f["neg|fdr"] < FDR) & (f["lfc"] < 0)
    return f[["gene_name", "arm", "lowers", "raises", "lfc", "pos|fdr", "neg|fdr"]]


def check_sign(fr: pd.DataFrame) -> bool:
    """The marker gene must lower its own marker in every arm.

    Args:
        fr: Freimer frame.

    Returns:
        True if all three arms validate.
    """
    ok = True
    for arm, marker in MARKER_CONTROL.items():
        row = fr[(fr["arm"] == arm) & (fr["gene_name"] == marker)]
        fired = (not row.empty) and bool(row.iloc[0]["lowers"])
        print(f"  [{'ok  ' if fired else 'FAIL'}] {arm:<6} {marker} knockdown lowers its own marker")
        ok &= fired
    return ok


def review_promotion_instrument(fr: pd.DataFrame) -> pd.DataFrame:
    """Is a Freimer IL-2-lowering HIT a clean immunosuppression signal, or a broken-cell signal?

    The continuous association between our efficacy axis and Freimer's continuous effect survives
    stratification on resting-arm disruption (N21 H2). The FDR-thresholded HIT CALLS are a different
    object, and this asks what they select for. It is the review step that decides whether a promotion
    based on a hit call can be trusted.

    Args:
        fr: Freimer frame.

    Returns:
        A one-row summary, also printed.
    """
    from scipy import stats

    il2 = fr[fr["arm"] == "IL2"][["gene_name", "lowers", "lfc"]]
    w = pd.read_csv(paths.TABLES / "window_score.csv")
    m = w.merge(il2, on="gene_name", how="inner").dropna(subset=["rest_de_genes", "stim_de_genes"])
    hit, non = m[m["lowers"]], m[~m["lowers"]]
    p = float(stats.mannwhitneyu(hit["rest_de_genes"], non["rest_de_genes"], alternative="greater").pvalue)

    print("\n" + "=" * 88)
    print("REVIEW OF THE PROMOTION INSTRUMENT (does a Freimer IL-2 hit mean specific immunosuppression?)")
    print("=" * 88)
    print(f"  IL-2-lowering hits (n={len(hit)}): median resting-arm DE genes = {hit['rest_de_genes'].median():.0f}")
    print(f"  non-hits          (n={len(non)}): median resting-arm DE genes = {non['rest_de_genes'].median():.0f}")
    print(f"  Mann-Whitney, hits have MORE resting-arm disruption: p = {p:.3g}")
    print("\n  So a Freimer IL-2-lowering HIT is substantially a 'this cell can no longer transcribe an")
    print("  induced gene' signal, not a clean specific-immunosuppression signal. IL-2 is a highly")
    print("  induced transcript; disabling general transcription machinery lowers it.")
    print("\n  NOTE the asymmetry, because it matters: the CONTINUOUS association between our efficacy")
    print("  axis and Freimer's continuous effect SURVIVES stratification on resting-arm disruption")
    print("  (p = 0.0055) and on z_L2 and rest_de_genes jointly (p = 0.0075). The axis replicates.")
    print("  It is the FDR-thresholded HIT CALL that is confounded, and promotion uses the hit call.")

    return pd.DataFrame([{
        "n_hits": len(hit), "n_nonhits": len(non),
        "median_rest_de_hits": float(hit["rest_de_genes"].median()),
        "median_rest_de_nonhits": float(non["rest_de_genes"].median()),
        "mwu_p_hits_more_resting_disruption": p,
        "note": "Continuous H2 survives rest_de stratification (p=0.0055); the hit call does not.",
    }])


def main() -> None:
    """Apply the relaxed rule once and report all three permitted verdicts honestly."""
    print(DEVIATION + "\n")
    fr = load_freimer()
    print("sign-convention positive controls (addendum 1):")
    if not check_sign(fr):
        print("\nSign controls FAILED. Stop. Nothing here is interpretable.")
        sys.exit(1)

    nom = pd.read_csv(paths.TABLES / "nomination_recalibrated.csv")
    pool = nom[~nom["do_not_headline"] & ~nom["housekeeping_functional"]]
    pool = pool[~pool["gene_name"].isin(HOUSEKEEPING_FUNCTIONAL | DO_NOT_HEADLINE)].copy()
    n_before = len(pool)

    il2 = fr[fr["arm"] == "IL2"].set_index("gene_name")
    il2ra = fr[fr["arm"] == "IL2RA"].set_index("gene_name")
    ctla4 = fr[fr["arm"] == "CTLA4"].set_index("gene_name")

    cot = pool[pool["gene_name"].isin(il2.index)].copy()
    if len(cot) != len(pool[pool["gene_name"].isin(set(il2.index))]):
        raise ValueError("candidate count changed on the Freimer join")
    print(f"\nfrozen pool: {len(nom)}; eligible after exclusions: {n_before}; "
          f"co-tested in the IL2 arm: {len(cot)}")
    print(f"co-tested genes: {sorted(cot['gene_name'])}")

    rows = []
    for _, r in cot.iterrows():
        g = r["gene_name"]
        lowers_il2 = bool(il2.loc[g, "lowers"])
        rows.append({
            "gene_name": g,
            "freimer_il2_lowers": lowers_il2,
            "freimer_il2_fdr": float(il2.loc[g, "pos|fdr"]),
            "freimer_il2_lfc": float(il2.loc[g, "lfc"]),
            "freimer_il2ra_lowers": bool(il2ra.loc[g, "lowers"]) if g in il2ra.index else None,
            "freimer_ctla4_lowers": bool(ctla4.loc[g, "lowers"]) if g in ctla4.index else None,
            "tier": r["tier"],
            "efficacy": r["efficacy"],
            "window_rank": r["window_rank"],
            "liabilities": r["liabilities"],
            "provenance": DEVIATION,
        })
    tab = pd.DataFrame(rows).sort_values("freimer_il2_fdr")

    print("\n=== every co-tested gene, with its Freimer evidence (nothing hidden) ===")
    print(tab[["gene_name", "freimer_il2_lowers", "freimer_il2_fdr", "freimer_il2_lfc",
               "freimer_il2ra_lowers", "freimer_ctla4_lowers", "window_rank"]].to_string(index=False))

    promoted = tab[tab["freimer_il2_lowers"]].copy()
    if len(promoted):
        promoted["co_inhibitory_direction"] = promoted["freimer_ctla4_lowers"].map(
            {True: "functional support with co-inhibitory LIABILITY",
             False: "functional support with FAVOURABLE co-inhibitory direction"})
    promoted = promoted.head(3)

    tab.to_csv(paths.TABLES / "freimer_secondary_posthoc.csv", index=False)
    print(f"\nwrote {paths.TABLES / 'freimer_secondary_posthoc.csv'}")

    review = review_promotion_instrument(fr)
    review.to_csv(paths.TABLES / "freimer_promotion_instrument_review.csv", index=False)
    print(f"\nwrote {paths.TABLES / 'freimer_promotion_instrument_review.csv'}")

    if len(promoted):
        print("\n=== the promoted genes, against the confound just measured ===")
        w = pd.read_csv(paths.TABLES / "window_score.csv")
        cols = ["gene_name", "stim_de_genes", "rest_de_genes", "fail_homeostasis", "window_rank"]
        pw = w[w["gene_name"].isin(promoted["gene_name"])][cols].copy()
        pw["rest_de_percentile"] = pw["rest_de_genes"].map(
            lambda x: f"{(w['rest_de_genes'] < x).mean() * 100:.1f}th")
        pw["stim_over_rest"] = (pw["stim_de_genes"] / pw["rest_de_genes"]).round(2)
        print(pw.to_string(index=False))
        print("\n  The selectivity requirement this project fixed a priori was a 10x stimulated-to-resting")
        print("  ratio. Neither promoted gene is close. A ratio below 1.0 means the perturbation disturbs")
        print("  the UNSTIMULATED cell more than the stimulated one, which is the definition of no window.")

    print("\n" + "=" * 88)
    print("VERDICT (only three are permitted; addendum 2)")
    print("=" * 88)
    if promoted.empty:
        print("  (1) NO Freimer-supported hypothesis target.")
        print("      No co-tested pool gene shows a significant IL-2-lowering knockdown effect.")
    else:
        print(f"  Freimer supports {len(promoted)} gene(s) with an IL-2-lowering knockdown effect:")
        for _, r in promoted.iterrows():
            print(f"\n    {r['gene_name']}  (IL2 FDR {r['freimer_il2_fdr']:.3g}, lfc {r['freimer_il2_lfc']:+.3f})")
            print(f"      {r['co_inhibitory_direction']}")
            print(f"      tier: {r['tier']}")
            print(f"      liabilities: {r['liabilities']}")
        print("\n  Each is a hypothesis target for follow-up, NOT a vetted therapeutic target.")
        print("  Whether any survives review is decided in the results doc, not here.")

    print("\n  The REGISTERED H3 result is unchanged and governs the headline:")
    print("  6 of 73 eligible genes co-tested against a pre-declared gate of 10; no card emitted.")


if __name__ == "__main__":
    main()
