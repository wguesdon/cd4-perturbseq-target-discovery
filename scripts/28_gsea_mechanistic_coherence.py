"""N18. GSEA mechanistic coherence: does the efficacy axis suppress canonical activation programs?

RULE #9, confirmatory (told the user this will not find a new target). The N16 skeptic lens asked whether a
top-efficacy perturbation genuinely suppresses the T-cell activation program or just moves a lot of genes
(a magnitude effect). This tests it: preranked GSEA on the stimulated (Stim48hr) DE z-scores against MSigDB
Hallmark, comparing the TOP-EFFICACY safe perturbations with a MAGNITUDE-MATCHED low-efficacy control set.

If the efficacy axis measures real activation-program suppression, the top-efficacy set should show strongly
NEGATIVE NES for activation/inflammatory Hallmark sets (IL2-STAT5, TNFA-NFKB, inflammatory, IFN-gamma,
allograft), the magnitude-matched control should show this much more weakly (specificity, not magnitude),
and an unrelated housekeeping set (adipogenesis, myogenesis) should be null in both (RULE #1 controls).

The ranking uses only the 10,282 panel genes, so the panel is the background by construction (never the
genome). GSEA statistics come from gseapy prerank (gene-set permutation).

Reads the h5ad zscore layer (via h5py, the panel-wide per-perturbation z), results/tables/*.csv, and a
cached Hallmark GMT (fetched once from Enrichr). Writes results/tables/gsea_*.csv.

Usage:
    uv run --group gsea python scripts/28_gsea_mechanistic_coherence.py

`gseapy` lives in the `gsea` dependency group rather than the core dependencies: it is the only
script that needs it, and it pulls in `requests`. The group is declared in `pyproject.toml`, so a
clean checkout reproduces this script from `uv sync --group gsea` alone.
"""

from __future__ import annotations

import argparse
import urllib.request

import h5py
import numpy as np
import pandas as pd

from cd4_perturbseq import paths

GMT_URL = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020"
GMT_PATH = paths.EXTERNAL / "gene_sets" / "msigdb_hallmark_2020.gmt"
CONDITION = "Stim48hr"
N_TOP = 25
SEED = 0

# Hallmark sets by expected behaviour, for the RULE #1 controls.
ACTIVATION_SETS = ["IL-2/STAT5 Signaling", "TNF-alpha Signaling via NF-kB", "Inflammatory Response",
                   "Interferon Gamma Response", "Allograft Rejection"]
HOUSEKEEPING_SETS = ["Adipogenesis", "Myogenesis", "Spermatogenesis", "Pancreas Beta Cells"]


def fetch_gmt() -> str:
    """Download the Hallmark GMT once and cache it. Returns the local path."""
    GMT_PATH.parent.mkdir(parents=True, exist_ok=True)
    if not GMT_PATH.exists():
        req = urllib.request.Request(GMT_URL, headers={"User-Agent": "curl/8"})
        with urllib.request.urlopen(req, timeout=60) as r:
            GMT_PATH.write_bytes(r.read())
    return str(GMT_PATH)


def select_groups() -> tuple[list[str], list[str], pd.DataFrame]:
    """Top-efficacy safe genes, and a z_l2-decile-matched low-efficacy safe control set."""
    saf = pd.read_csv(paths.TABLES / "window_score_organism_safety.csv")
    z = pd.read_csv(paths.TABLES / "magnitude_matched_rows.csv")[["gene_name", "z_l2"]]
    df = saf[saf["safe"]].merge(z, on="gene_name", how="inner").dropna(subset=["efficacy", "z_l2"])
    df["z_decile"] = pd.qcut(df["z_l2"], 10, labels=False, duplicates="drop")
    top = df.nlargest(N_TOP, "efficacy")
    # magnitude-matched control: same z_l2 decile counts, drawn from the bottom half of efficacy
    low = df[df["efficacy"] <= df["efficacy"].median()]
    rng = np.random.default_rng(SEED)
    ctrl_idx: list[int] = []
    for dec, k in top["z_decile"].value_counts().items():
        pool = low[low["z_decile"] == dec]
        if len(pool):
            ctrl_idx += list(rng.choice(pool.index, size=min(k, len(pool)), replace=False))
    ctrl = df.loc[ctrl_idx]
    return top["gene_name"].tolist(), ctrl["gene_name"].tolist(), df


def read_zscores(genes: list[str]) -> pd.Series:
    """Mean panel-wide z-score per gene-symbol across the given perturbations, condition = Stim48hr."""
    import anndata as ad
    A = ad.read_h5ad(paths.DE_STATS_H5AD, backed="r")
    obs = A.obs
    sym = A.var["gene_name"].to_numpy()
    mask = (obs["culture_condition"] == CONDITION) & (obs["target_contrast_gene_name"].isin(genes))
    row_idx = np.where(mask.to_numpy())[0]
    if len(row_idx) == 0:
        raise SystemExit("no matching perturbation rows")
    order = np.argsort(row_idx)
    with h5py.File(paths.DE_STATS_H5AD, "r") as f:
        z = f["layers/zscore"][row_idx[order], :]  # (n_pert, 10282)
    mean_z = np.nanmean(z, axis=0)
    s = pd.Series(mean_z, index=sym)
    s = s.groupby(level=0).mean().dropna()  # collapse duplicate symbols
    return s


def run_gsea(ranked: pd.Series, gmt: str, label: str) -> pd.DataFrame:
    """Preranked GSEA of one ranked list against Hallmark. Returns the res2d table with a label."""
    import gseapy as gp
    rnk = ranked.sort_values(ascending=False).reset_index()
    rnk.columns = ["gene", "score"]
    pre = gp.prerank(rnk=rnk, gene_sets=gmt, min_size=15, max_size=500,
                     permutation_num=1000, seed=SEED, threads=4, outdir=None, no_plot=True)
    res = pre.res2d.copy()
    res["group"] = label
    for c in ("NES", "NOM p-val", "FDR q-val"):
        if c in res:
            res[c] = pd.to_numeric(res[c], errors="coerce")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args()
    paths.ensure_dirs()
    gmt = fetch_gmt()

    top_genes, ctrl_genes, df = select_groups()
    print(f"top-efficacy safe (n={len(top_genes)}): {top_genes[:8]} ...")
    print(f"magnitude-matched low-efficacy control (n={len(ctrl_genes)}): {ctrl_genes[:8]} ...")
    a = df.set_index("gene_name")
    print(f"  efficacy: top mean {a.loc[top_genes,'efficacy'].mean():.3f} vs ctrl {a.loc[ctrl_genes,'efficacy'].mean():.3f}")
    print(f"  |z_l2|:   top mean {a.loc[top_genes,'z_l2'].mean():.3f} vs ctrl {a.loc[ctrl_genes,'z_l2'].mean():.3f} "
          "(matched)")

    top_rank = read_zscores(top_genes)
    ctrl_rank = read_zscores(ctrl_genes)
    print(f"\nranked genes: top {len(top_rank)}, ctrl {len(ctrl_rank)} (panel background)")

    res_top = run_gsea(top_rank, gmt, "top_efficacy")
    res_ctrl = run_gsea(ctrl_rank, gmt, "magnitude_control")
    both = pd.concat([res_top, res_ctrl], ignore_index=True)
    both.to_csv(paths.TABLES / "gsea_hallmark.csv", index=False)

    def show(res: pd.DataFrame, label: str) -> None:
        r = res.set_index("Term")
        print(f"\n=== {label}: activation sets (expect NEGATIVE NES) ===")
        for s in ACTIVATION_SETS:
            hit = [t for t in r.index if s.split("/")[0][:8].lower() in t.lower() or s.lower() in t.lower()]
            for t in hit[:1]:
                print(f"  {t}: NES={r.loc[t,'NES']:+.2f}  FDR={r.loc[t,'FDR q-val']:.3f}")
        print(f"  housekeeping controls (expect ~0 / non-sig):")
        for s in HOUSEKEEPING_SETS:
            hit = [t for t in r.index if s.lower() in t.lower()]
            for t in hit[:1]:
                print(f"    {t}: NES={r.loc[t,'NES']:+.2f}  FDR={r.loc[t,'FDR q-val']:.3f}")

    show(res_top, "TOP-EFFICACY")
    show(res_ctrl, "MAGNITUDE-MATCHED CONTROL")

    # ---- RULE #1 verdicts
    rt = res_top.set_index("Term")
    def nes(res_idx, key):
        hit = [t for t in res_idx.index if key.lower() in t.lower()]
        return res_idx.loc[hit[0], "NES"] if hit else np.nan
    il2_top = nes(rt, "IL-2/STAT5")
    infl_top = nes(rt, "Inflammatory Response")
    house_top = nes(rt, "Adipogenesis")
    rc = res_ctrl.set_index("Term")
    il2_ctrl = nes(rc, "IL-2/STAT5")
    pos_ok = (il2_top < 0 and rt.loc[[t for t in rt.index if 'il-2/stat5' in t.lower()][0], "FDR q-val"] < 0.25)
    neg_ok = abs(house_top) < 1.5 or nes(rt, "Adipogenesis") != nes(rt, "Adipogenesis")  # null-ish
    spec_ok = (not np.isnan(il2_ctrl)) and (il2_top < il2_ctrl)  # more suppressed in top than control
    print("\n" + "=" * 78)
    print(f"POSITIVE CONTROL (IL-2/STAT5 down in top-efficacy): NES={il2_top:+.2f} -> {'PASS' if pos_ok else 'weak'}")
    print(f"NEGATIVE CONTROL (housekeeping ~null in top): Adipogenesis NES={house_top:+.2f} -> "
          f"{'PASS' if neg_ok else 'FAIL'}")
    print(f"SPECIFICITY (IL-2/STAT5 more suppressed in top than magnitude-matched): "
          f"top {il2_top:+.2f} vs ctrl {il2_ctrl:+.2f} -> {'PASS' if spec_ok else 'weak/FAIL'}")
    print(f"\nwrote {paths.TABLES / 'gsea_hallmark.csv'}")
    print("Confirmatory only: coherence of the efficacy axis with canonical activation programs. It does not")
    print("nominate a target and does not resolve direction (N17).")
    print("=" * 78)


if __name__ == "__main__":
    main()
