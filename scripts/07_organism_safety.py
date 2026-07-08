"""Add the organism-level safety layer, and validate the axes before trusting them.

Our safety gate was T-cell-intrinsic. "Rest" meant resting CD4 T cells, so a mitochondrial or
one-carbon gene whose knockdown spares a quiescent T cell but destroys gut epithelium looked
perfectly safe. The Tier A shortlist is exactly that class, so this is not hypothetical.

Three additions, in the order a sceptic would demand them.

**First, validate our own viability measure.** We claimed that resting cell count is a
context-free viability signal, on the strength of a 283-gene subset. Hart CEGv2 (684 core
essentials) and Hart NEGv1 (927 nonessentials) are a positive and negative control pair. Core
essentials must deplete at rest; nonessentials must not. If both deplete, we are measuring
something other than essentiality and the tier is meaningless.

**Second, validate the new axes on the same controls.** Core essentials should be
loss-of-function intolerant in humans and ubiquitously expressed. Nonessentials should not be.
An axis that cannot separate CEGv2 from NEGv1 cannot be trusted to judge our shortlist.

**Third, only then apply them.** LOEUF asks whether living humans tolerate losing this gene.
Tissue breadth asks whether an inhibitor could ever spare the rest of the body.

Usage:
    uv run python scripts/07_organism_safety.py
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy import stats

from cd4_perturbseq import paths, priors

LOEUF_INTOLERANT = 0.60
UBIQUITOUS_TISSUES = 40


def _mwu(label: str, a: pd.Series, b: pd.Series, alternative: str, expect: str) -> bool | None:
    """Run a one-sided Mann-Whitney U test and print a verdict line.

    Returns None, never False, when a control arm is too small to test. An untested axis is not a
    failed axis, and conflating the two is how this project already reached one wrong conclusion:
    "zero essentials in the top 100" was read as evidence against essentiality when in truth the
    essentials had been removed by QC before ranking.

    Args:
        label: What is being compared.
        a: Values for the positive control (core essentials).
        b: Values for the negative control (nonessentials).
        alternative: ``"less"`` or ``"greater"``, describing ``a`` relative to ``b``.
        expect: Human-readable statement of the expected direction.

    Returns:
        True if the controls separate as expected, False if they do not, None if untestable.
    """
    a, b = a.dropna(), b.dropna()
    if len(a) < 5 or len(b) < 5:
        print(f"  {label:34s} NOT TESTABLE: n={len(a)} essential, {len(b)} nonessential")
        return None
    _, p = stats.mannwhitneyu(a, b, alternative=alternative)
    ok = p < 0.05
    print(
        f"  {label:34s} essential median {a.median():8.3f}  nonessential {b.median():8.3f}  "
        f"p={p:.2e}  {'PASS' if ok else 'FAIL'}   ({expect})"
    )
    return ok


def main() -> None:
    """Validate the safety axes on controls, then annotate the window score."""
    argparse.ArgumentParser(description=__doc__).parse_args()
    paths.ensure_dirs()

    window = pd.read_csv(paths.TABLES / "window_score.csv")
    constraint = priors.gnomad_constraint()
    breadth = priors.hpa_tissue_breadth()
    ceg = priors.hart_core_essentials_full()
    neg = priors.hart_nonessentials()

    print(f"gnomAD constraint: {len(constraint):,} genes")
    print(f"HPA tissue breadth: {len(breadth):,} genes across 51 consensus tissues")
    print(f"Hart CEGv2 core-essentials: {len(ceg)}   Hart NEGv1 nonessentials: {len(neg)}")

    frame = window.merge(constraint, on="gene_name", how="left").merge(breadth, on="gene_name", how="left")
    frame["is_ceg"] = frame["gene_name"].isin(ceg)
    frame["is_neg"] = frame["gene_name"].isin(neg)
    print(f"\nof {len(frame):,} analysable perturbations: "
          f"{int(frame['is_ceg'].sum())} are CEGv2 core-essential, {int(frame['is_neg'].sum())} are NEGv1 nonessential")
    print(f"  (the bundled 283-gene subset only ever had 31 rankable, which is why it could not adjudicate)")

    # ---------------------------------------------------------------- control validation
    ess = frame[frame["is_ceg"]]
    non = frame[frame["is_neg"]]

    print("\n=== CONTROL VALIDATION: can each axis separate CEGv2 from NEGv1? ===")
    verdicts = {
        "our viability axis": _mwu(
            "resting cells retained", ess["rest_cells_ratio"], non["rest_cells_ratio"],
            "less", "essentials should be DEPLETED at rest",
        ),
        "gnomAD LOEUF": _mwu(
            "LOEUF (LoF intolerance)", ess["loeuf"], non["loeuf"],
            "less", "essentials should be LoF-INTOLERANT, low LOEUF",
        ),
        "HPA breadth": _mwu(
            "non-immune tissues expressed", ess["n_nonimmune_tissues"], non["n_nonimmune_tissues"],
            "greater", "essentials should be UBIQUITOUS",
        ),
    }
    print()
    for name, ok in verdicts.items():
        state = {True: "validated on controls", False: "DID NOT SEPARATE THE CONTROLS",
                 None: "UNTESTED, control arm too small"}[ok]
        print(f"  {name:22s} {state}")

    if any(v is None for v in verdicts.values()):
        print(f"\n  Only {int(frame['is_neg'].sum())} of {len(neg)} Hart nonessential controls survive")
        print("  perturbation QC. Nonessential genes are lowly expressed, so they fail the")
        print("  ontarget_significant filter. This is a SELECTION problem, not a failed test:")
        print("  conditioning on QC conditions on a collider. Run the controls on the full obs")
        print("  (all genes with a resting row) before concluding anything about essentiality.")
        print("  See docs/results/adversarial_audit_2026_07_08.md, finding 2.")

    # ---------------------------------------------------------------- the organism-level layer
    frame["lof_intolerant"] = frame["loeuf"] < LOEUF_INTOLERANT
    frame["ubiquitous"] = frame["n_nonimmune_tissues"] >= UBIQUITOUS_TISSUES
    frame["systemic_risk"] = frame["lof_intolerant"] | frame["ubiquitous"]

    known = frame["loeuf"].notna() & frame["n_nonimmune_tissues"].notna()
    print(f"\ngenes with both annotations: {int(known.sum()):,} of {len(frame):,}")
    print(f"  LoF-intolerant (LOEUF < {LOEUF_INTOLERANT}): {int(frame['lof_intolerant'].sum()):,}")
    print(f"  ubiquitous (expressed in >= {UBIQUITOUS_TISSUES} of 45 non-immune tissues): {int(frame['ubiquitous'].sum()):,}")

    # ---------------------------------------------------------------- apply to the shortlist
    tier_a = frame[frame["safe"] & (frame["viability_tier"] == "non-depleting")].nlargest(20, "window_score")
    cols = ["gene_name", "window_score", "loeuf", "pli", "n_nonimmune_tissues",
            "max_nonimmune_ntpm", "lof_intolerant", "ubiquitous", "systemic_risk"]
    print("\n=== TIER A SHORTLIST, now judged at the level of a whole human ===")
    print(tier_a[cols].round(3).to_string(index=False))

    n_flagged = int(tier_a["systemic_risk"].sum())
    print(f"\n  {n_flagged} of {len(tier_a)} Tier A hits carry a systemic-risk flag.")
    survivors = tier_a[~tier_a["systemic_risk"]]["gene_name"].tolist()
    print(f"  Survive organism-level safety: {', '.join(survivors) if survivors else 'NONE'}")

    print("\n=== approved-drug targets, for calibration ===")
    truth = pd.read_csv(paths.GROUND_TRUTH / "immunomodulator_targets.csv")
    positives = set(truth.loc[truth["include_as_positive"], "gene_symbol"])
    drugs = frame[frame["gene_name"].isin(positives) & ~frame["fail_evidence"]]
    print(drugs[["gene_name", "loeuf", "pli", "n_nonimmune_tissues", "lof_intolerant", "ubiquitous", "safe"]].round(3).to_string(index=False))

    out = paths.TABLES / "window_score_organism_safety.csv"
    frame.to_csv(out, index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
