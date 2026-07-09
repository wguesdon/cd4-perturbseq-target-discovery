"""N19. Cross-screen concordance: do our candidates replicate in independent T-cell CRISPR screens?

RULE #9 build, the discovery push the user chose. N16 gave zero VETTED novel targets; N17 showed direction
of effect is intrinsically hard. This asks a different, replication question: are our within-screen
candidates (effector-module suppressors that pass the safety gate) reproduced as genuine functional
regulators in INDEPENDENT T-cell CRISPR screens? A gene that is a hit in our Perturb-seq screen AND in a
separate functional screen is far less likely to be an artifact of our one dataset, which is the strongest
lightweight way "combining datasets" can raise a candidate's credibility.

Independent screens (vendored, MIT, provenance in resources/external_screens/PROVENANCE.md):
  - Freimer 2022: regulators of IL2RA surface expression (Marson lab).
  - Arce 2025: fitness in Resting_Teff / Stimulated_Teff / Resting_Treg.
Schmidt & Steinhart 2022 is NOT used here (RULE #3, held out for validation only).

Design (RULE #1 controls first):
  - Universe = genes tested in BOTH our screen and the independent screen (no collider on our own filters).
  - Enrichment: are our SAFE genes enriched for independent-screen hits vs the co-tested background, with a
    magnitude-matched null so the signal is not just "big-effect genes hit everywhere"? Positive control:
    recovered approved-drug targets should be enriched. Negative control: a label shuffle must destroy it.
  - Replicated candidates: SAFE + genetically-supported + independent-screen hit, ranked. Honest caveat:
    replication means "real regulator", NOT "right therapeutic direction" (still UNRESOLVED, N17).

Reads results/tables/*.csv and resources/external_screens/*.csv. Never touches the h5ad layers.

Usage:
    uv run python scripts/27_cross_screen_concordance.py [--fdr 0.05] [--n-draws 20000]
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy import stats

from cd4_perturbseq import matched, paths

SCREENS = paths.REPO_ROOT / "resources" / "external_screens"
KNOWN_DRUGS = frozenset({"IMPDH2", "PPP3R1", "CD3E", "CD3G", "IL4R", "CD2", "CD28"})
SEED = 0


def _clean(df: pd.DataFrame) -> pd.DataFrame:
    """Strip a UTF-8 BOM from the id column and header if present."""
    df.columns = [c.replace("﻿", "") for c in df.columns]
    df.iloc[:, 0] = df.iloc[:, 0].astype(str).str.replace("﻿", "", regex=False)
    return df


def load_freimer(fdr: float) -> pd.DataFrame:
    """Per-gene hit/direction for the Freimer 2022 screen, aggregated over its 3 marker arms.

    The file has one row per (gene, marker) for markers IL2RA/IL2/CTLA4 (3 x 1351). We collapse to one row
    per gene by taking the most significant arm, so the merge does not fan out to duplicate rows.
    """
    f = _clean(pd.read_csv(SCREENS / "Freimer2022_Screen.csv")).rename(columns={"id": "gene_name"})
    f["row_min_fdr"] = f[["neg|fdr", "pos|fdr"]].min(axis=1)
    # depleted (neg) = knockdown lowers the marker -> positive regulator of the activation marker
    f["row_dir"] = np.where(f["neg|fdr"] <= f["pos|fdr"], "pos_regulator", "neg_regulator")
    f = f.sort_values("row_min_fdr").groupby("gene_name", as_index=False).first()
    f["freimer_min_fdr"] = f["row_min_fdr"]
    f["freimer_hit"] = f["freimer_min_fdr"] < fdr
    f["freimer_dir"] = f["row_dir"]
    f["freimer_marker"] = f["screen"]
    return f[["gene_name", "freimer_hit", "freimer_min_fdr", "freimer_dir", "freimer_marker"]]


def load_arce(fdr: float) -> pd.DataFrame:
    """Per-gene hit/direction for the Arce 2025 fitness screen (Stimulated_Teff primary)."""
    a = _clean(pd.read_csv(SCREENS / "Arce2025_Screen.csv")).rename(columns={"id": "gene_name"})
    cond = "Stimulated_Teff"
    a["arce_min_fdr"] = a[[f"neg|fdr.{cond}", f"pos|fdr.{cond}"]].min(axis=1)
    a["arce_hit"] = a["arce_min_fdr"] < fdr
    # depleted (neg) = knockdown reduces stimulated Teff fitness -> gene required for effector expansion
    a["arce_dir"] = np.where(a[f"neg|fdr.{cond}"] <= a[f"pos|fdr.{cond}"], "required_for_fitness",
                             "suppresses_fitness")
    return a[["gene_name", "arce_hit", "arce_min_fdr", "arce_dir"]]


def _matched_enrichment(frame: pd.DataFrame, hit_col: str, n_draws: int, rng) -> dict:
    """Are SAFE genes enriched for independent-screen hits vs a z_l2-decile-matched null of co-tested genes?"""
    co = frame[frame[hit_col].notna()].copy()  # co-tested universe
    co["z_decile"] = pd.qcut(co["z_l2"], 10, labels=False, duplicates="drop")
    safe = co[co["safe"]]
    obs = int(safe[hit_col].sum())
    draws, shortfall = matched.count_matched(co, co["z_decile"], safe.index,
                                             {"hit": co[hit_col].astype(bool)}, n_draws, rng)
    p = float((draws["hit"] >= obs).mean())
    return {"universe": len(co), "safe": len(safe), "safe_hits": obs,
            "matched_mean": float(draws["hit"].mean()), "p": p, "shortfall": shortfall}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--n-draws", type=int, default=20000)
    args = ap.parse_args()
    paths.ensure_dirs()
    rng = np.random.default_rng(SEED)

    saf = pd.read_csv(paths.TABLES / "window_score_organism_safety.csv")
    rows = pd.read_csv(paths.TABLES / "magnitude_matched_rows.csv")[["gene_name", "z_l2"]]
    nom = pd.read_csv(paths.TABLES / "n10_nomination.csv")[
        ["gene_name", "ot_genetic_supported", "ot_genetic_max", "ot_genetic_n_diseases", "tractable"]]
    frame = saf.merge(rows, on="gene_name", how="left").merge(nom, on="gene_name", how="left")
    frame = frame[frame["z_l2"].notna()].reset_index(drop=True)

    frame = frame.merge(load_freimer(args.fdr), on="gene_name", how="left")
    frame = frame.merge(load_arce(args.fdr), on="gene_name", how="left")
    frame["any_hit"] = frame[["freimer_hit", "arce_hit"]].any(axis=1)
    # co-tested in at least one independent screen
    frame["cotested"] = frame["freimer_hit"].notna() | frame["arce_hit"].notna()
    frame["any_hit"] = np.where(frame["cotested"], frame["any_hit"], np.nan)

    print(f"our tested genes: {len(frame)}")
    for name, col in (("Freimer", "freimer_hit"), ("Arce", "arce_hit")):
        n_co = int(frame[col].notna().sum())
        n_hit = int(frame[col].fillna(False).sum())
        print(f"  {name}: co-tested {n_co}, hits (FDR<{args.fdr}) {n_hit}")

    print("\n=== ENRICHMENT: are SAFE genes enriched for independent-screen hits? (magnitude-matched null) ===")
    controls = {}
    for name, col in (("Freimer", "freimer_hit"), ("Arce", "arce_hit"), ("either", "any_hit")):
        r = _matched_enrichment(frame, col, args.n_draws, rng)
        sig = r["p"] < 0.05
        print(f"  {name}: safe {r['safe_hits']}/{r['safe']} hits vs matched {r['matched_mean']:.1f}; "
              f"p={r['p']:.4f} {'*' if sig else ''} (universe {r['universe']})")

    # Positive control: recovered drug targets should be enriched among independent-screen hits.
    co = frame[frame["cotested"]]
    drugs = co[co["gene_name"].isin(KNOWN_DRUGS)]
    drug_hit = int(drugs["any_hit"].fillna(False).sum())
    base = float(co["any_hit"].fillna(False).mean())
    p_drug = stats.binomtest(drug_hit, len(drugs), base, alternative="greater").pvalue if len(drugs) else np.nan
    controls["pos_drugs"] = (p_drug < 0.10) if not np.isnan(p_drug) else False
    print(f"\n  CONTROL+ recovered drugs: {drug_hit}/{len(drugs)} are independent-screen hits "
          f"(base rate {base:.2f}, binom p={p_drug:.3f}) {'PASS' if controls['pos_drugs'] else 'weak'}")

    # Negative control: label shuffle must destroy the 'safe' enrichment.
    sh = co.copy()
    sh["safe"] = rng.permutation(sh["safe"].to_numpy())
    sh["z_decile"] = pd.qcut(sh["z_l2"], 10, labels=False, duplicates="drop")
    st = sh[sh["safe"]]
    d_sh, _ = matched.count_matched(sh, sh["z_decile"], st.index,
                                    {"hit": sh["any_hit"].fillna(False).astype(bool)}, args.n_draws, rng)
    p_sh = float((d_sh["hit"] >= int(st["any_hit"].fillna(False).sum())).mean())
    controls["neg_shuffle"] = 0.05 <= p_sh <= 0.95
    print(f"  CONTROL- label-shuffle: p={p_sh:.3f} {'PASS' if controls['neg_shuffle'] else 'FAIL'}")

    # ---- Replicated candidates: safe + genetically supported + independent-screen hit
    frame["ot_genetic_supported"] = frame["ot_genetic_supported"].fillna(False)
    rep = frame[(frame["safe"]) & (frame["ot_genetic_supported"]) & (frame["any_hit"].fillna(False))].copy()
    rep["n_screens_hit"] = frame[["freimer_hit", "arce_hit"]].fillna(False).sum(axis=1)
    rep = rep.sort_values(["ot_genetic_n_diseases", "efficacy"], ascending=[False, False])
    cols = ["gene_name", "efficacy", "window_rank", "ot_genetic_max", "ot_genetic_n_diseases",
            "tractable", "is_iei", "freimer_hit", "freimer_marker", "freimer_dir", "arce_hit", "arce_dir"]
    print("\n=== REPLICATED CANDIDATES: safe + genetically supported + independent-screen hit ===")
    print(f"  n = {len(rep)} (of {int(frame['safe'].sum())} safe genes)")
    if len(rep):
        print(rep[cols].to_string(index=False))

    keep = ["gene_name", "safe", "efficacy", "window_rank", "ot_genetic_supported", "ot_genetic_max",
            "ot_genetic_n_diseases", "tractable", "is_iei", "freimer_hit", "freimer_min_fdr", "freimer_dir",
            "freimer_marker", "arce_hit", "arce_min_fdr", "arce_dir", "any_hit"]
    frame[keep].to_csv(paths.TABLES / "cross_screen_concordance.csv", index=False)
    print(f"\nwrote {paths.TABLES / 'cross_screen_concordance.csv'}")

    ctrl_ok = controls.get("neg_shuffle", False)
    print("\n" + "=" * 80)
    print(f"CONTROLS: negative shuffle {'PASS' if controls['neg_shuffle'] else 'FAIL'}, "
          f"positive drugs {'PASS' if controls['pos_drugs'] else 'weak'}")
    if not ctrl_ok:
        print("VOID: the negative control failed; enrichment is uninterpretable.")
    print("Replication raises confidence a gene is a REAL functional regulator across datasets. It does NOT")
    print("resolve therapeutic direction (N17): a replicated candidate is a stronger hypothesis, not a")
    print("vetted nomination.")
    print("=" * 80)


if __name__ == "__main__":
    main()
