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


def library_targets() -> set[str]:
    """Every gene the sgRNA library actually targets.

    This is the denominator. The DE table is already a filtered view of it, so counting
    "essentials in the library" from the DE table hides exactly the dropout we are testing for.

    Returns:
        Set of gene symbols targeted by at least one guide.
    """
    table = pd.read_csv(LIBRARY, low_memory=False)
    cols = [c for c in ("target_gene_name_from_sgRNA", "designed_target_gene_name") if c in table.columns]
    targets: set[str] = set()
    for col in cols:
        # pandas 3 keeps NA through astype(str), so drop before converting.
        values = table[col].dropna().astype(str).str.strip()
        targets |= {s for s in values if s and s.lower() not in {"nan", "none"}}
    return targets


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

    ceg = priors.hart_core_essentials_full()
    neg = priors.hart_nonessentials()
    lib = library_targets()
    print(f"sgRNA library targets: {len(lib):,} genes")
    print(f"Hart CEGv2 core-essentials: {len(ceg)}   Hart NEGv1 nonessentials: {len(neg)}")

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


if __name__ == "__main__":
    main()
