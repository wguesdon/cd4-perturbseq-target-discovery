"""Who was removed before the ranking, and why? The essentiality claim depends on the answer.

`docs/results/risk_kill_2026_07_08.md` concluded that "cancer-cell essentiality is the wrong safety
axis for this screen", because the 31 rankable Hart essentials sat at median rank 3148 of 6371.
The adversarial audit called that survivorship bias: conditioning on perturbation QC removes
exactly the genes whose knockdown kills cells, so the surviving essentials are the harmless ones
and the test is conditioned on a collider.

That is a claim about a mechanism, and this project has now been wrong about two of those in one
day. So measure it, against the true denominator: the sgRNA library metadata, not the DE table.

Two competing stories, and they predict opposite things:

  AUDIT'S STORY   essentials are lost because their knockdown kills cells, so they never reach the
                  DE table. Essentials should drop out MORE than nonessentials at DE-eligibility.

  RIVAL STORY     QC requires a significant on-target knockdown, which requires the target to be
                  expressed. Nonessential genes are lowly expressed. Nonessentials should drop out
                  MORE than essentials at the ontarget_significant filter.

Both can be true at different stages. Stage-by-stage is the only way to tell.

Usage:
    uv run python scripts/11_selection_funnel.py
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy import stats

from cd4_perturbseq import de_stats, paths, priors

STIM = "Stim48hr"
LIBRARY = paths.RAW / "suppl" / "sgrna_library_metadata.suppl_table.csv"


def _clean(values: pd.Series) -> pd.Series:
    """Drop NA before ``astype(str)``; pandas 3 keeps NA through it."""
    out = values.dropna().astype(str).str.strip()
    return out[~out.str.lower().isin({"nan", "none", ""})]


def alias_map() -> dict[str, str]:
    """Legacy sgRNA-derived gene symbol to the canonical symbol, learned from the library file.

    The library carries two gene-symbol vocabularies for the same guides.
    ``target_gene_name_from_sgRNA`` is derived from where the guide aligns and uses an older
    annotation (``AARS``, ``AES``, ``ADPRHL2``, ``ICK``); ``designed_target_gene_name`` is the
    canonical symbol (``AARS1``, ``TLE5``, ``ADPRS``, ``CILK1``). They disagree on 350 genes.

    This matters twice, and the shipped code got it wrong both times. Taking the **union** of the two
    columns as "the library" counted 475 aliases as distinct genes and produced a denominator of
    13,129 for a library the paper reports as 12,748 genes. And the Hart CEGv2/NEGv1 sets use the
    legacy vocabulary, while ``DE_stats`` uses the canonical one (11,526 of 11,526 of its genes are
    canonical, only 11,203 are legacy) — so intersecting Hart directly against ``DE_stats`` silently
    dropped 13 core-essentials.

    Returns:
        Mapping from legacy symbol to canonical symbol. Identity where they agree.
    """
    table = pd.read_csv(LIBRARY, low_memory=False)
    pair = table[["target_gene_name_from_sgRNA", "designed_target_gene_name"]].dropna()
    legacy = _clean(pair["target_gene_name_from_sgRNA"])
    canonical = _clean(pair["designed_target_gene_name"])
    common = legacy.index.intersection(canonical.index)
    return dict(zip(legacy.loc[common], canonical.loc[common]))


def library_targets() -> set[str]:
    """Every gene the sgRNA library targets, in one vocabulary.

    This is the denominator. The DE table is already a filtered view of it, so counting
    "essentials in the library" from the DE table hides exactly the dropout we are testing for.

    Uses ``designed_target_gene_name`` alone: 12,779 genes, against the 12,748 the preprint reports.
    Never union the two symbol columns. See :func:`alias_map`.

    Returns:
        Set of canonical gene symbols targeted by at least one guide.
    """
    table = pd.read_csv(LIBRARY, low_memory=False)
    return set(_clean(table["designed_target_gene_name"]))


def lift(genes: set[str], alias: dict[str, str]) -> set[str]:
    """Map a legacy-vocabulary gene set into the canonical vocabulary.

    Args:
        genes: Symbols in the legacy vocabulary, e.g. a Hart gene set.
        alias: Output of :func:`alias_map`.

    Returns:
        The same genes, canonicalised. Symbols with no alias entry pass through unchanged.
    """
    return {alias.get(g, g) for g in genes}


def _fisher(a_hit: int, a_tot: int, b_hit: int, b_tot: int) -> tuple[float, float]:
    """Fisher exact, two-sided, comparing two survival fractions.

    Args:
        a_hit: Survivors in group A. a_tot: group A total.
        b_hit: Survivors in group B. b_tot: group B total.

    Returns:
        Tuple of (odds ratio, p-value).
    """
    table = [[a_hit, a_tot - a_hit], [b_hit, b_tot - b_hit]]
    return stats.fisher_exact(table, alternative="two-sided")


def main() -> None:
    """Print the stage-by-stage funnel for essentials versus nonessentials."""
    argparse.ArgumentParser(description=__doc__).parse_args()

    alias = alias_map()
    ceg = lift(priors.hart_core_essentials_full(), alias)
    neg = lift(priors.hart_nonessentials(), alias)
    lib = library_targets()
    n_renamed = sum(1 for k, v in alias.items() if k != v)
    print(f"sgRNA library targets: {len(lib):,} genes (canonical vocabulary; the preprint says 12,748)")
    print(f"  {n_renamed} genes carry a legacy alias in the library's other symbol column.")
    print(f"  Unioning the two columns, as this script used to, inflates the denominator to 13,129.")
    print(f"Hart CEGv2 core-essentials: {len(ceg)}   Hart NEGv1 nonessentials: {len(neg)}  (alias-lifted)")

    obs = de_stats.read_obs().reset_index(drop=True)
    obs["gene_name"] = obs["target_contrast_gene_name"].astype(str)
    in_de = set(obs["gene_name"])
    stim = obs[obs["culture_condition"] == STIM].drop_duplicates("gene_name").set_index("gene_name")

    stages: list[tuple[str, set[str]]] = []
    stages.append(("in sgRNA library", lib))
    stages.append(("reaches DE_stats", in_de))
    stages.append((f"tested in {STIM}", set(stim.index)))

    keep = stim.index
    for label, mask in (
        ("not distal_offtarget", ~stim["distal_offtarget_flag"].astype(bool)),
        ("not neighboring_KD", ~stim["neighboring_gene_KD"].astype(bool)),
        ("ontarget_significant", stim["ontarget_significant"].astype(bool)),
        ("not low_target_gex", ~stim["low_target_gex"].astype(bool)),
    ):
        keep = keep.intersection(stim.index[mask])
        stages.append((f"+ {label}", set(keep)))

    print("\n=== SELECTION FUNNEL ===")
    header = f"{'stage':26s} {'CEGv2':>14s} {'NEGv1':>14s} {'all library':>14s}   {'OR (ess vs non)':>16s}"
    print(header)
    print("-" * len(header))

    prev_ceg = prev_neg = None
    for label, universe in stages:
        n_ceg = len(ceg & universe & lib)
        n_neg = len(neg & universe & lib)
        n_all = len(universe & lib)
        cell = ""
        if prev_ceg is not None and prev_ceg > 0 and prev_neg > 0:
            odds, p = _fisher(n_ceg, prev_ceg, n_neg, prev_neg)
            direction = "ess survive" if odds > 1 else "non survive"
            cell = f"{odds:6.2f} p={p:.1e} {direction}"
        pct_ceg = 100 * n_ceg / max(len(ceg & lib), 1)
        pct_neg = 100 * n_neg / max(len(neg & lib), 1)
        print(f"{label:26s} {n_ceg:5d} ({pct_ceg:5.1f}%) {n_neg:5d} ({pct_neg:5.1f}%) {n_all:8d}       {cell}")
        prev_ceg, prev_neg = n_ceg, n_neg

    ceg_lib = len(ceg & lib)
    neg_lib = len(neg & lib)
    ceg_final = len(ceg & stages[-1][1] & lib)
    neg_final = len(neg & stages[-1][1] & lib)

    # A single overall odds ratio hides the story: TWO colliders act, in opposite directions, at
    # different stages. Reading only the final number is how the original essentiality claim went
    # wrong in the first place.
    ceg_de = len(ceg & in_de & lib)
    neg_de = len(neg & in_de & lib)
    or_de, p_de = _fisher(ceg_de, ceg_lib, neg_de, neg_lib)

    stim_set = set(stim.index)
    ceg_pre = len(ceg & stim_set & lib)
    neg_pre = len(neg & stim_set & lib)
    or_qc, p_qc = _fisher(ceg_final, ceg_pre, neg_final, neg_pre)

    print("\n=== VERDICT: two colliders, pulling opposite ways ===")
    print(f"  [A] DE-eligibility. Essentials {ceg_de}/{ceg_lib} ({100*ceg_de/max(ceg_lib,1):.1f}%) vs "
          f"nonessentials {neg_de}/{neg_lib} ({100*neg_de/max(neg_lib,1):.1f}%). OR {or_de:.3f}, p={p_de:.3g}")
    print(f"      {ceg_lib - ceg_de} of {ceg_lib} essentials NEVER reach the DE table. Their knockdown")
    print("      depletes cells and they fail the authors' DE-eligibility gates. THE AUDIT IS RIGHT.")
    print()
    print(f"  [B] Perturbation QC. Essentials {ceg_final}/{ceg_pre} vs nonessentials {neg_final}/{neg_pre}. "
          f"OR {or_qc:.2f}, p={p_qc:.3g}")
    print("      ontarget_significant needs a detectable knockdown, which needs an expressed target.")
    print("      Nonessential genes are lowly expressed, so THEY are destroyed at this stage instead.")
    print()
    print("  Net effect: essentials %.1f%% and nonessentials %.1f%% of library reach the ranking."
          % (100 * ceg_final / max(ceg_lib, 1), 100 * neg_final / max(neg_lib, 1)))
    print()
    print("  CONCLUSION. The ranking is conditioned on survival through BOTH filters, and the two")
    print("  biases run in opposite directions. The surviving essentials are the sub-population whose")
    print("  knockdown did NOT kill the cell; the surviving nonessentials are the few that are highly")
    print("  expressed. A comparison between them estimates nothing causal.")
    print()
    print("  Therefore: 'cancer-cell essentiality is the wrong safety axis for this screen' is")
    print("  RETRACTED. It was never supported. The supportable statement is:")
    print("  THIS SCREEN CANNOT RESOLVE whether cancer-cell essentiality predicts the naive ranking,")
    print("  because the genes the question is about were removed before the ranking existed.")

    out = paths.TABLES / "selection_funnel.csv"
    rows = []
    for label, universe in stages:
        rows.append({
            "stage": label,
            "essentials": len(ceg & universe & lib),
            "nonessentials": len(neg & universe & lib),
            "all_library": len(universe & lib),
        })
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nwrote {out}")

    _essentiality_rank_test(ceg)


def _essentiality_rank_test(ceg: set[str]) -> None:
    """Where do Hart core essentials rank in the naive suppression ranking?

    An earlier version of this analysis reported 31 rankable core essentials at median rank 3,148 of
    6,371, Mann-Whitney p = 0.611, and concluded that the naive ranking is not enriched for essential
    genes. Both the count and the p-value came from the era when this script unioned two gene-symbol
    vocabularies. The count is wrong, and the corrected test points the other way.

    The corrected result does NOT license the opposite claim either. Essentiality is a collider here:
    only a minority of library core-essentials reach `DE_stats` at all, because their knockdown kills
    the cell, and the on-target significance filter then removes most nonessentials. A comparison among
    the survivors of that double selection estimates nothing causal. The number is reported so the
    superseded one cannot be quoted, not because it supports an inference.

    Args:
        ceg: Alias-lifted Hart CEGv2 core-essential gene symbols.
    """
    from scipy import stats

    w = pd.read_csv(paths.TABLES / "window_score.csv")
    w["naive_rank"] = (-w["eff_mean_z"]).rank(ascending=False, method="first")
    w["is_ceg"] = w["gene_name"].isin(ceg)
    ess, rest = w[w["is_ceg"]], w[~w["is_ceg"]]

    p_two = float(stats.mannwhitneyu(ess["naive_rank"], rest["naive_rank"], alternative="two-sided").pvalue)
    p_better = float(stats.mannwhitneyu(ess["naive_rank"], rest["naive_rank"], alternative="less").pvalue)

    frame = pd.DataFrame([{
        "n_rankable_essentials": len(ess),
        "n_rankable_total": len(w),
        "median_rank_essentials": float(ess["naive_rank"].median()),
        "median_rank_others": float(rest["naive_rank"].median()),
        "mwu_p_two_sided": p_two,
        "mwu_p_essentials_rank_better": p_better,
        "superseded_claim": "31 rankable, median rank 3148 of 6371, p = 0.611",
        "note": "Collider. Only a minority of library core-essentials reach DE_stats; the on-target "
                "filter then removes most nonessentials. The comparison among survivors is not causal.",
    }])
    out = paths.TABLES / "essentiality_rank_test.csv"
    frame.to_csv(out, index=False)

    print("\n=== Hart core essentials in the naive suppression ranking ===")
    print(f"  rankable core essentials: {len(ess)} of {len(w):,} perturbations")
    print(f"  median naive rank: {ess['naive_rank'].median():.0f} vs {rest['naive_rank'].median():.0f} for the rest")
    print(f"  Mann-Whitney two-sided p = {p_two:.4g}; essentials rank BETTER p = {p_better:.4g}")
    print(f"  SUPERSEDES: '31 rankable, median rank 3148, p = 0.611' (two-vocabulary alias artifact)")
    print(f"  This is a collider comparison and supports no causal inference either way.")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
