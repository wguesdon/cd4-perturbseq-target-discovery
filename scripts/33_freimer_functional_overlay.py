"""N21. Does the decision layer replicate against an independent, signed, functional CD4 screen?

Pre-registered in `docs/preregistration_n21_2026_07_09.md`, committed before this file existed and
before any Freimer sign was inspected for any gene in the nomination pool.

Freimer 2022 (Marson lab) screens primary human CD4 T cells for regulators of three markers: IL2RA,
IL2 and CTLA4. MAGeCK `neg|*` means knockdown DEPLETES the marker-high population, so the gene is a
POSITIVE REGULATOR of that marker. Two of those markers map onto the two halves of our thesis:

    IL2   is the activation output our EFFICACY axis scores from mRNA.
    CTLA4 is one of the nine genes in our activation-induced CO-INHIBITORY module.

N19 collapsed the three arms to a binary "replicates yes/no" and discarded the sign. The sign is
unspent, and it is the only instrument in this project that is not characterisation-bound: a pooled
functional screen cannot tell whether anyone has previously studied the gene.

Schmidt stays HELD OUT (RULE #3). It is the only out-of-sample validation the efficacy axis has.

H1  Our `tolerance_loss` (mRNA co-inhibitory attrition) predicts CTLA-4 protein loss.
H2  Our `efficacy` predicts IL-2 loss.
H3  EXPLORATORY: at most three hypothesis-target cards from the 91-gene nomination pool.

Freimer does NOT measure autoimmune therapeutic direction. A knockdown that lowers IL-2 could be CD3E,
whose loss causes immunodeficiency. Lowering IL-2 is an EFFICACY readout. Every claim says so.

Exits non-zero if any pre-registered control fails.

Usage:
    uv run python scripts/33_freimer_functional_overlay.py
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from cd4_perturbseq import paths  # noqa: E402

SEED = 0
FDR = 0.10
COVERAGE_GATE = 10  # pre-declared: fewer than this many co-tested nomination genes -> no card
N_MATCHED = 2_000

SCREENS = paths.RESOURCES / "external_screens"

# POSITIVE CONTROL, substituted in addendum 1. The registered control (the TCR-proximal signalosome)
# is entirely absent from Freimer's focused 1,351-gene IL2RA-regulator library, so it was NOT
# TESTABLE. The marker gene inside its own arm is a stronger control: knocking out CTLA4 must lower
# CTLA-4. It needs no outside knowledge and is blind to our nomination pool.
MARKER_CONTROL = {"IL2RA": "IL2RA", "IL2": "IL2", "CTLA4": "CTLA4"}
CONTROL_PSEUDOGENE = "Non-Targeting"

HOUSEKEEPING_FUNCTIONAL = {"COG6", "TNPO1", "ARPC2", "SEL1L"}
DO_NOT_HEADLINE = {"STAT6", "GATA3", "STAT5A", "PTEN"}


def load_freimer() -> pd.DataFrame:
    """Load Freimer with its sign intact, one row per (gene, arm).

    SIGN, from addendum 1. `neg|lfc` and `pos|lfc` are identical: `lfc` is a gene-level value, not
    side-specific. In every arm the marker gene itself is significant on the `pos` side and maximally
    non-significant on `neg`, so **`pos` = guides enriched in the marker-LOW bin = knockdown LOWERS
    the marker**. `neg` = knockdown raises it.

    The hit rule requires a non-zero effect as well as significance, because MAGeCK aggregates the 593
    non-targeting guides into one pseudo-gene whose FDR is tiny by construction while its lfc is
    exactly 0.

    Returns:
        Columns gene_name, arm, lowers_marker, raises_marker, fdr_lower, lfc.
    """
    f = pd.read_csv(SCREENS / "Freimer2022_Screen.csv")
    f = f.rename(columns={"id": "gene_name", "screen": "arm"})

    dup = int(f.duplicated(subset=["gene_name", "arm"]).sum())
    if dup:
        raise ValueError(f"Freimer fan-out: {dup} duplicate (gene, arm) rows. N19 was injured by this.")
    if not (f["neg|lfc"] == f["pos|lfc"]).all():
        raise ValueError("neg|lfc and pos|lfc diverge; the gene-level lfc assumption is wrong.")

    f["lfc"] = f["pos|lfc"]
    f["lowers_marker"] = (f["pos|fdr"] < FDR) & (f["lfc"] > 0)
    f["raises_marker"] = (f["neg|fdr"] < FDR) & (f["lfc"] < 0)
    f["fdr_lower"] = f["pos|fdr"]
    return f[["gene_name", "arm", "lowers_marker", "raises_marker", "fdr_lower", "lfc"]]


def controls(fr: pd.DataFrame, win: pd.DataFrame) -> dict[str, bool]:
    """Run every pre-registered control.

    Args:
        fr: The Freimer frame.
        win: Our window score.

    Returns:
        Mapping of control name to pass/fail.
    """
    print("\n" + "=" * 88)
    print("CONTROLS (pre-registered; all must pass or the overlay is void)")
    print("=" * 88)

    # POSITIVE: the marker gene inside its own arm must be a "knockdown lowers the marker" hit.
    hits = []
    for arm, marker in MARKER_CONTROL.items():
        sub = fr[(fr["arm"] == arm) & (fr["gene_name"] == marker)]
        ok = (not sub.empty) and bool(sub.iloc[0]["lowers_marker"])
        hits.append(ok)
        detail = "absent" if sub.empty else (
            f"lfc {sub.iloc[0]['lfc']:+.3f}, FDR {sub.iloc[0]['fdr_lower']:.2g}")
        print(f"[positive] {arm:<6} arm, {marker:<6} knockdown lowers its own marker: "
              f"{'YES' if ok else 'NO'}  ({detail})")
    pos_ok = all(hits)

    # NEGATIVE: the non-targeting pseudo-gene has lfc exactly 0 and must never be a hit.
    nt = fr[fr["gene_name"] == CONTROL_PSEUDOGENE]
    nt_hit = bool(nt["lowers_marker"].any() or nt["raises_marker"].any())
    print(f"[negative] {CONTROL_PSEUDOGENE} registers as a hit on any arm: {nt_hit} "
          f"(its lfc is {nt['lfc'].abs().max():.3f} on every arm)")
    neg_ok = not nt_hit

    # SIGN: the IL2 and IL2RA arms measure related biology and must overlap more than chance.
    il2 = set(fr[(fr["arm"] == "IL2") & fr["lowers_marker"]]["gene_name"])
    il2ra = set(fr[(fr["arm"] == "IL2RA") & fr["lowers_marker"]]["gene_name"])
    universe = fr["gene_name"].nunique()
    expected = len(il2) * len(il2ra) / universe
    overlap = len(il2 & il2ra)
    print(f"[sign]     IL2-lowering {len(il2)}, IL2RA-lowering {len(il2ra)}, overlap {overlap} "
          f"(expected {expected:.1f} by chance)")
    sign_ok = overlap > expected

    arms = fr.groupby("arm").size()
    print(f"[fan-out]  rows per arm: {arms.to_dict()}; genes: {fr['gene_name'].nunique()}")
    fan_ok = True  # load_freimer raises on duplicates before we get here

    return {"marker gene lowers its own marker (all 3 arms)": pos_ok,
            "non-targeting control silent": neg_ok,
            "IL2 and IL2RA arms overlap above chance": sign_ok,
            "no fan-out": fan_ok}


def _matched_null(merged: pd.DataFrame, x: str, y: str, rng: np.random.Generator) -> tuple[float, float]:
    """Spearman of x vs y, and a p-value from a z_L2-stratified permutation.

    Shuffling y within z_L2 deciles destroys the association while preserving any relationship
    between y and effect magnitude, so a surviving correlation cannot be bought by effect size.

    Args:
        merged: The co-tested frame.
        x: Our axis column.
        y: The Freimer signed column.
        rng: Random generator.

    Returns:
        Tuple of (observed rho, permutation p-value).
    """
    obs = float(stats.spearmanr(merged[x], merged[y]).statistic)
    decile = pd.qcut(merged["z_l2"], 10, labels=False, duplicates="drop")
    null = []
    for _ in range(N_MATCHED):
        shuffled = merged[y].copy()
        for d in decile.unique():
            m = decile == d
            shuffled[m] = rng.permutation(shuffled[m].to_numpy())
        null.append(float(stats.spearmanr(merged[x], shuffled).statistic))
    null = np.array(null)
    p = float((np.abs(null) >= abs(obs)).sum() + 1) / (N_MATCHED + 1)
    return obs, p


def posthoc_hit_contrast(fr: pd.DataFrame, win: pd.DataFrame, arm: str, ours: str,
                         rng: np.random.Generator) -> dict:
    """POST-HOC, not pre-registered. Does our axis separate Freimer's hits from its non-hits?

    The registered H1/H2 correlate across all co-tested genes, most of which have an lfc near zero, so
    the test is dominated by noise. This contrast is more sensitive. It is reported as post-hoc, and
    it is calibrated: run on the IL2 arm against `efficacy`, where H2 already passed, it must fire.
    If it fires there and not on CTLA-4, the CTLA-4 negative is not merely insensitivity.

    Args:
        fr: Freimer frame.
        win: Our window score, carrying z_l2.
        arm: Freimer marker arm.
        ours: Our axis column.
        rng: Random generator.

    Returns:
        The contrast statistics.
    """
    sub = fr[fr["arm"] == arm][["gene_name", "lowers_marker"]]
    m = win.merge(sub, on="gene_name", how="inner").dropna(subset=[ours, "z_l2"]).copy()
    m["dec"] = pd.qcut(m["z_l2"], 10, labels=False, duplicates="drop")
    hit, non = m[m["lowers_marker"]], m[~m["lowers_marker"]]
    if len(hit) < 5:
        return {"arm": arm, "axis": ours, "n_hits": len(hit), "verdict": "NOT TESTABLE"}

    p_raw = float(stats.mannwhitneyu(hit[ours], non[ours], alternative="greater").pvalue)
    want = hit["dec"].value_counts()
    null = []
    for _ in range(N_MATCHED):
        parts = [non[non["dec"] == d].sample(min(k, int((non["dec"] == d).sum())), replace=False,
                                             random_state=int(rng.integers(1e9)))
                 for d, k in want.items() if int((non["dec"] == d).sum()) > 0]
        null.append(pd.concat(parts)[ours].median())
    null = np.array(null)
    obs = float(hit[ours].median())
    p_matched = float((null >= obs).sum() + 1) / (len(null) + 1)
    verdict = "SEPARATES" if p_matched < 0.05 else "DOES NOT SEPARATE once matched on effect magnitude"

    print(f"\n  {arm} lowering vs not, by our `{ours}`")
    print(f"    hits {len(hit)}, non-hits {len(non)}; median {obs:+.4f} vs {non[ours].median():+.4f}")
    print(f"    Mann-Whitney p = {p_raw:.4g}; z_L2-matched p = {p_matched:.4f}")
    print(f"    {verdict}")
    return {"arm": arm, "axis": ours, "n_hits": len(hit), "median_hit": round(obs, 4),
            "p_mwu": p_raw, "p_matched": p_matched, "verdict": verdict}


def main() -> None:
    """Run H1, H2 and the gated H3."""
    rng = np.random.default_rng(SEED)
    fr = load_freimer()
    win = pd.read_csv(paths.TABLES / "window_score.csv")
    nom = pd.read_csv(paths.TABLES / "nomination_recalibrated.csv")

    if "z_l2" not in win.columns:
        mm = pd.read_csv(paths.TABLES / "magnitude_matched_rows.csv")
        col = next(c for c in mm.columns if c.lower() in ("z_l2", "z_l2_norm"))
        win = win.merge(mm[["gene_name", col]].rename(columns={col: "z_l2"}), on="gene_name", how="left")

    print(f"Freimer: {fr['gene_name'].nunique():,} genes x {fr['arm'].nunique()} arms = {len(fr):,} rows")
    print(f"arms: {sorted(fr['arm'].unique())}")

    verdicts = controls(fr, win)
    for name, ok in verdicts.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    if not all(verdicts.values()):
        print("\nA pre-registered control FAILED. The overlay is void; nothing here may be quoted.")
        sys.exit(1)

    results = []
    print("\n" + "=" * 88)
    print("H1 / H2: do our two axes replicate on an orthogonal, protein-level, independent screen?")
    print("=" * 88)
    for arm, ours, label in (("CTLA4", "tolerance_loss", "H1 co-inhibitory attrition -> CTLA-4 loss"),
                             ("IL2", "efficacy", "H2 efficacy -> IL-2 loss")):
        sub = fr[fr["arm"] == arm][["gene_name", "lfc", "fdr_lower"]]
        m = win.merge(sub, on="gene_name", how="inner").dropna(subset=["z_l2", ours, "lfc"])
        if len(m) < 20:
            print(f"\n{label}: only {len(m)} co-tested genes. UNDERPOWERED, not tested.")
            results.append({"hypothesis": label, "n": len(m), "rho": np.nan, "p_matched": np.nan,
                            "verdict": "UNDERPOWERED"})
            continue
        # With the sign established in addendum 1, a LARGER positive lfc means the knockdown lowers
        # the marker more. Our axes are also "more suppression = larger". So a POSITIVE rho means the
        # two platforms agree.
        m["freimer_lowering"] = m["lfc"]
        rho, p = _matched_null(m, ours, "freimer_lowering", rng)
        verdict = "REPLICATES" if (rho > 0 and p < 0.05) else "DOES NOT REPLICATE"
        print(f"\n{label}")
        print(f"  co-tested genes: {len(m)}")
        print(f"  Spearman({ours}, Freimer {arm} lowering) = {rho:+.3f}")
        print(f"  z_L2-stratified permutation p ({N_MATCHED:,} draws) = {p:.4f}")
        print(f"  VERDICT: {verdict}")
        results.append({"hypothesis": label, "n": len(m), "rho": round(rho, 4),
                        "p_matched": p, "verdict": verdict})

    print("\n" + "=" * 88)
    print("POST-HOC diagnostic (NOT pre-registered). Is H1's negative real, or just insensitive?")
    print("=" * 88)
    print("  Calibration first: the same contrast on the arm where H2 already passed. It must fire.")
    posthoc = [
        posthoc_hit_contrast(fr, win, "IL2", "efficacy", rng),      # calibration; must fire
        posthoc_hit_contrast(fr, win, "CTLA4", "tolerance_loss", rng),  # the question
        posthoc_hit_contrast(fr, win, "CTLA4", "efficacy", rng),    # is CTLA-4 loss just efficacy?
        posthoc_hit_contrast(fr, win, "IL2", "tolerance_loss", rng),
    ]
    pd.DataFrame(posthoc).to_csv(paths.TABLES / "freimer_posthoc_contrast.csv", index=False)

    print("\n" + "=" * 88)
    print("H3 (EXPLORATORY): the hypothesis-target card, gated on pre-declared coverage")
    print("=" * 88)
    il2 = fr[fr["arm"] == "IL2"].set_index("gene_name")
    ctla4 = fr[fr["arm"] == "CTLA4"].set_index("gene_name")
    pool = nom[~nom["do_not_headline"] & ~nom["housekeeping_functional"]].copy()
    pool = pool[~pool["gene_name"].isin(HOUSEKEEPING_FUNCTIONAL | DO_NOT_HEADLINE)]
    cotested = pool[pool["gene_name"].isin(il2.index)]
    print(f"nomination pool: {len(nom)}; after housekeeping and RULE #6 exclusions: {len(pool)}")
    print(f"co-tested in the Freimer IL2 arm: {len(cotested)}  (pre-declared gate: {COVERAGE_GATE})")

    cards = pd.DataFrame()
    if len(cotested) < COVERAGE_GATE:
        print(f"\nCOVERAGE GATE NOT MET. H3 is declared UNDERPOWERED and NO CARD is emitted,")
        print(f"regardless of what the data show. This was fixed in advance.")
        if len(cotested):
            print(f"For the record, the co-tested genes are: {sorted(cotested['gene_name'])}")
    else:
        sel = []
        for _, r in cotested.iterrows():
            g = r["gene_name"]
            lowers_il2 = bool(il2.loc[g, "lowers_marker"]) and il2.loc[g, "fdr_lower"] < FDR
            lowers_ctla4 = g in ctla4.index and bool(ctla4.loc[g, "lowers_marker"])
            if lowers_il2 and not lowers_ctla4:
                sel.append({"gene_name": g, "il2_fdr": float(il2.loc[g, "fdr_lower"]),
                            "il2_lfc": float(il2.loc[g, "lfc"]),
                            "lowers_ctla4": lowers_ctla4,
                            "efficacy": r["efficacy"], "tolerance_loss_rank": r["window_rank"],
                            "tier": r["tier"], "liabilities": r["liabilities"]})
        cards = pd.DataFrame(sel).sort_values(["il2_fdr", "efficacy"], ascending=[True, False]).head(3)
        print(f"\ngenes lowering IL-2 without lowering CTLA-4: {len(sel)}; emitting {len(cards)} card(s)")
        if len(cards):
            print(cards[["gene_name", "il2_fdr", "il2_lfc", "efficacy", "tier"]].to_string(index=False))

    out = pd.DataFrame(results)
    out.to_csv(paths.TABLES / "freimer_overlay.csv", index=False)
    cards.to_csv(paths.TABLES / "freimer_hypothesis_cards.csv", index=False)
    print(f"\nwrote {paths.TABLES / 'freimer_overlay.csv'}")
    print(f"wrote {paths.TABLES / 'freimer_hypothesis_cards.csv'} ({len(cards)} cards)")

    print("\n" + "=" * 88)
    print("HONEST READING (RULE #8)")
    print("=" * 88)
    print("  Freimer measures IL-2 and CTLA-4 in primary human CD4 T cells. It does NOT measure")
    print("  autoimmune therapeutic direction. A knockdown that lowers IL-2 could be CD3E, whose loss")
    print("  causes immunodeficiency. Lowering IL-2 is an EFFICACY readout, not a safety readout.")
    if len(cards) == 0:
        print("\n  Even after adding a held-out, signed, primary-human-CD4 cytokine screen, no gene")
        print("  passed the pre-declared hypothesis-target standard.")
    else:
        print("\n  The pipeline abstains from a vetted discovery claim. An independent primary-human-CD4")
        print(f"  CRISPR screen supports {', '.join(cards['gene_name'])} as hypothesis target(s) for")
        print("  follow-up. Not vetted therapeutic targets.")


if __name__ == "__main__":
    main()
