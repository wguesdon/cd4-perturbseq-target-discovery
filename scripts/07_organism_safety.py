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
"""Dominant LoF-intolerance cut. gnomAD's own convention is ~0.35 for strong constraint; 0.60 is a
more permissive "flag for review" bar. Documented, fixed a priori, never tuned against the shortlist.
This annotates; it does not gate."""

RECESSIVE_PREC = 0.90
"""gnomAD posterior of recessive constraint above which a potent inhibitor, which phenocopies
biallelic loss, is flagged. 0.90 is a strong-recessive bar. Documented, fixed a priori."""

BROAD_NTPM = 25.0
"""A meaningful transcript level in at least one non-immune tissue (nTPM). "Expressed" is usually
nTPM >= 1; 25 is a "substantially expressed" bar. Documented, fixed a priori, replaces the degenerate
tissue-count flag."""

UBIQUITOUS_TISSUES = 40
"""Legacy. The n_nonimmune_tissues >= 40 flag had no information (5,716 of 6,371; IQR all at 45).
Kept only so the legacy `ubiquitous` column still exists; it is out of systemic_risk."""


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
    # Three organism-level warnings, and NaN means UNKNOWN, never SAFE. The previous version wrote
    # `loeuf < 0.60`, and `NaN < 0.60` is False, so 365 genes with no LOEUF were silently filed as
    # LoF-tolerant. `pd.NA`-typed comparisons keep the unknowns unknown.
    frame["lof_intolerant"] = (frame["loeuf"] < LOEUF_INTOLERANT).where(frame["loeuf"].notna(), other=pd.NA).astype("boolean")

    # The recessive axis. LOEUF and pLI measure HAPLOinsufficiency. A potent inhibitor phenocopies
    # BIALLELIC loss, which is what gnomAD's prec captures. This is the axis that catches the
    # recessive mitochondrial-disease genes in Tier A that LOEUF waves through (POLG prec 0.9996,
    # VARS2 0.9999, ACAD9 0.9933, IMPDH2 0.9988, all with LOEUF that reads as tolerant).
    frame["recessive_intolerant"] = (frame["prec"] > RECESSIVE_PREC).where(frame["prec"].notna(), other=pd.NA).astype("boolean")

    # Breadth, on the continuous measure. `n_nonimmune_tissues >= 40` flagged 5,716 of 6,371 and
    # carried no information: its 25th, 50th and 75th percentiles are all 45.0. max_nonimmune_ntpm
    # spans 0.1 to 34,000 and actually separates. "Broadly expressed" here means a meaningful
    # transcript level (>= 25 nTPM) in at least one non-immune tissue, a documented a-priori cut.
    frame["broadly_expressed"] = (frame["max_nonimmune_ntpm"] >= BROAD_NTPM).where(
        frame["max_nonimmune_ntpm"].notna(), other=pd.NA).astype("boolean")
    # Legacy column, kept so nothing downstream breaks, but no longer part of systemic_risk.
    frame["ubiquitous"] = frame["n_nonimmune_tissues"] >= UBIQUITOUS_TISSUES

    # systemic_risk fires if ANY axis fires; it is UNKNOWN only when every axis is unknown.
    flags = frame[["lof_intolerant", "recessive_intolerant", "broadly_expressed"]]
    frame["systemic_risk"] = np.where(
        flags.eq(True).any(axis=1), True,
        np.where(flags.isna().all(axis=1), pd.NA, False),
    )
    frame["systemic_risk"] = frame["systemic_risk"].astype("boolean")

    print(f"\norganism-level annotation (REPORTED, never a gate; the gate is in 04):")
    print(f"  LoF-intolerant, dominant  (LOEUF < {LOEUF_INTOLERANT}): {int(frame['lof_intolerant'].eq(True).sum()):,}")
    print(f"  LoF-intolerant, recessive (prec > {RECESSIVE_PREC}):   {int(frame['recessive_intolerant'].eq(True).sum()):,}")
    print(f"  broadly expressed         (max non-immune nTPM >= {BROAD_NTPM}): {int(frame['broadly_expressed'].eq(True).sum()):,}")
    print(f"  any systemic-risk flag: {int(frame['systemic_risk'].eq(True).sum()):,}   "
          f"unknown (no annotation at all): {int(frame['systemic_risk'].isna().sum()):,}")

    # ---------------------------------------------------------------- apply to the shortlist
    tier_a = frame[frame["safe"] & (frame["viability_tier"] == "non-depleting")].nlargest(20, "window_score")
    cols = ["gene_name", "window_score", "loeuf", "prec", "max_nonimmune_ntpm",
            "lof_intolerant", "recessive_intolerant", "broadly_expressed", "systemic_risk"]
    print("\n=== TIER A SHORTLIST, now judged at the level of a whole human ===")
    print(tier_a[cols].round(3).to_string(index=False))

    n_flagged = int(tier_a["systemic_risk"].eq(True).sum())
    n_recessive = int(tier_a["recessive_intolerant"].eq(True).sum())
    print(f"\n  {n_flagged} of {len(tier_a)} Tier A hits carry a systemic-risk flag; "
          f"{n_recessive} are recessive-intolerant (a potent inhibitor phenocopies biallelic loss).")
    survivors = tier_a[tier_a["systemic_risk"].eq(False)]["gene_name"].tolist()
    print(f"  Clear ALL organism-level axes: {', '.join(survivors) if survivors else 'NONE'}")
    print("  systemic_risk saturates because a T-cell-intrinsic screen surfaces core cellular")
    print("  machinery; the informative content is the per-axis resolution, above all the recessive")
    print("  axis, which LOEUF and pLI are blind to. Nobody may act on Tier A on this screen alone.")

    print("\n=== approved-drug targets, for calibration ===")
    truth = pd.read_csv(paths.GROUND_TRUTH / "immunomodulator_targets.csv")
    positives = set(truth.loc[truth["include_as_positive"], "gene_symbol"])
    drugs = frame[frame["gene_name"].isin(positives) & ~frame["fail_evidence"]]
    print(drugs[["gene_name", "loeuf", "prec", "max_nonimmune_ntpm",
                 "lof_intolerant", "recessive_intolerant", "safe"]].round(3).to_string(index=False))
    print("  IMPDH2 (mycophenolate) is recessive-intolerant (prec 0.999): the gate already files it")
    print("  as antiproliferative, and mycophenolate does need monitoring. PPP3R1 (ciclosporin) is")
    print("  not (prec 0.28). The recessive axis tracks the clinical distinction, on n=few.")

    out = paths.TABLES / "window_score_organism_safety.csv"
    frame.to_csv(out, index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
