"""Hard stop: the nomination path may not read any field derived from approved-drug status.

`clinical_precedent` leaked into the nomination gate once already. It looked like an obvious bug fix,
it was committed, and it would have made the drug-recovery validation circular, because that field is
true for 32 of the 36 curated positives. It was reverted. This script exists so it cannot happen again
without the pipeline failing loudly.

Two families of check, both of which must pass or the script exits non-zero.

**Leakage.** The rebuilt nomination gate may use ONLY: screen-gate passage, structural tractability,
and the direction-discordant veto. It may NOT use approved-drug labels, curated-positive labels,
`clinical_precedent`, known-immunomodulator status, drug-recovery benchmark labels, or
`maxClinicalStage`. The check is behavioural rather than lexical: for each forbidden field, we verify
that permuting it across genes leaves the nominated set bit-for-bit identical. A field the gate cannot
see cannot change its output.

**Fan-out.** Every external join in this project has one row per gene. This repo has already been
injured twice by fan-out: a gene-symbol union inflated the library denominator, and Freimer's three
marker arms per gene inflated a drug control from p = 0.17 to p = 0.001. Assert it, do not hope.

Usage:
    uv run python scripts/32_leakage_guard.py
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from cd4_perturbseq import paths, priors  # noqa: E402

SEED = 0

# Fields that encode "somebody already made a drug against this gene". The gate may not read these at
# all: if it did, a gene could be PROMOTED because a drug exists, and the drug-recovery validation
# would be circular.
FORBIDDEN = ["clinical_precedent", "tractable_with_precedent", "il2_hit"]

# The only inputs the nomination gate is permitted to read for promotion.
PERMITTED = ["safe", "tractable", "direction_verdict"]

# Read, but only to REMOVE a gene from a novel-candidate list, never to add one. An earlier version of
# this script listed `is_known_drug` as forbidden and reported it clean, because `nominate()` read a
# `known` column derived BEFORE the permutation, so shuffling the source field could not move anything.
# The check was vacuous. It is now tested for the property that actually matters: monotone exclusion.
EXCLUSION_ONLY = ["is_known_drug"]


def nominate(frame: pd.DataFrame) -> set[str]:
    """The nomination gate, reading its inputs directly so a permutation propagates.

    Args:
        frame: One row per screen-passing gene.

    Returns:
        The nominated gene symbols.
    """
    keep = frame["tractable"].fillna(False).astype(bool) & ~frame["direction_verdict"].eq("DISCORDANT")
    known = frame["is_known_drug"].fillna(False).astype(bool)
    return set(frame.loc[keep & ~known, "gene_name"])


def check_exclusion_only(f: pd.DataFrame) -> bool:
    """`is_known_drug` may only shrink the nominated set, never grow it.

    Blanking the label must yield a SUPERSET, and the genes it adds back must be exactly the known
    drugs that otherwise pass the gate. If blanking the label ever ADDS a gene that is not a known
    drug, the label is doing promotion work and that is leakage.

    Args:
        f: The joined frame.

    Returns:
        True if the label is monotone-exclusionary.
    """
    real = nominate(f)
    blanked = nominate(f.assign(is_known_drug=False))
    added = blanked - real
    known = set(f.loc[f["is_known_drug"].fillna(False).astype(bool), "gene_name"])
    superset = real <= blanked
    only_known = added <= known
    print(f"\n[exclusion-only] blanking `is_known_drug` adds {len(added)} genes: {sorted(added)}")
    print(f"                 superset holds: {superset}; every added gene is a known drug: {only_known}")
    return superset and only_known


def build() -> pd.DataFrame:
    """Join everything the nomination path could conceivably see.

    Returns:
        One row per screen-passing gene, carrying permitted and forbidden fields alike.
    """
    nom = pd.read_csv(paths.TABLES / "n10_nomination.csv")
    doe = pd.read_csv(paths.TABLES / "direction_of_effect.csv")[["gene_name", "direction_verdict"]]
    tr = priors.ot_tractability()[["gene_name", "tractable", "tractable_with_precedent",
                                   "clinical_precedent"]]
    f = nom.drop(columns=["tractable", "clinical_precedent"], errors="ignore")
    f = f.merge(doe, on="gene_name", how="left").merge(tr, on="gene_name", how="left")
    f["safe"] = True  # every row of n10_nomination is already screen-passing
    return f


def check_leakage(f: pd.DataFrame, rng: np.random.Generator) -> list[str]:
    """Permute each forbidden field; the nominated set must not move.

    Args:
        f: The joined frame.
        rng: Random generator.

    Returns:
        Names of fields whose permutation changed the output. Empty means no leakage.
    """
    baseline = nominate(f)
    leaks = []
    print(f"baseline nomination: {len(baseline)} genes")
    for col in FORBIDDEN:
        if col not in f.columns:
            print(f"  [skip] {col:<26} not present")
            continue
        shuffled = f.assign(**{col: rng.permutation(f[col].to_numpy())})
        after = nominate(shuffled)
        ok = after == baseline
        print(f"  [{'ok  ' if ok else 'LEAK'}] {col:<26} permuting it changes "
              f"{len(after ^ baseline)} genes")
        if not ok:
            leaks.append(col)
    return leaks


def check_positive_control(f: pd.DataFrame, rng: np.random.Generator) -> bool:
    """The guard must be able to FAIL. Permute a PERMITTED field and require the output to move.

    A leakage test that passes on everything is worthless. `tractable` is permitted and load-bearing,
    so shuffling it must change the nominated set.

    Args:
        f: The joined frame.
        rng: Random generator.

    Returns:
        True if the guard demonstrably detects a real dependency.
    """
    baseline = nominate(f)
    shuffled = f.assign(tractable=rng.permutation(f["tractable"].to_numpy()))
    moved = len(nominate(shuffled) ^ baseline)
    print(f"\n[positive control] permuting the PERMITTED field `tractable` moves {moved} genes")
    print("                   (a guard that cannot fire would report 0 here)")
    return moved > 0


def check_fanout() -> list[str]:
    """Every external annotation must be one row per gene.

    Returns:
        Names of sources with duplicate gene rows.
    """
    sources = {
        "ot_tractability": priors.ot_tractability(),
        "gnomad_constraint": priors.gnomad_constraint(),
        "hpa_tissue_breadth": priors.hpa_tissue_breadth(),
        "schmidt_cd4_il2_screen": priors.schmidt_cd4_il2_screen(),
        "ot_genetic_support": priors.ot_genetic_support(),
    }
    bad = []
    print("\nfan-out check: one row per gene in every external join")
    for name, frame in sources.items():
        key = "gene_name" if "gene_name" in frame.columns else frame.columns[0]
        dupes = int(frame[key].duplicated().sum())
        print(f"  [{'ok  ' if dupes == 0 else 'FAN-OUT'}] {name:<24} {len(frame):>7,} rows, "
              f"{dupes} duplicate {key}")
        if dupes:
            bad.append(name)
    return bad


def main() -> None:
    """Run both families of check and gate the exit code on them."""
    rng = np.random.default_rng(SEED)
    f = build()
    print(f"screen-passing genes: {len(f)}")
    print(f"permitted for promotion: {PERMITTED}")
    print(f"forbidden entirely:      {FORBIDDEN}")
    print(f"exclusion-only:          {EXCLUSION_ONLY}\n")

    leaks = check_leakage(f, rng)
    fires = check_positive_control(f, rng)
    monotone = check_exclusion_only(f)
    fanout = check_fanout()

    print("\n" + "=" * 76)
    ok = not leaks and fires and monotone and not fanout
    if leaks:
        print(f"  [FAIL] the nomination gate READS {leaks}. That is drug-label leakage.")
    else:
        print("  [PASS] no forbidden field can change the nominated set")
    print(f"  [{'PASS' if fires else 'FAIL'}] the guard is falsifiable: a permitted field does move the output")
    print(f"  [{'PASS' if monotone else 'FAIL'}] `is_known_drug` only excludes; it never promotes")
    if fanout:
        print(f"  [FAIL] fan-out in {fanout}")
    else:
        print("  [PASS] every external join is one row per gene")
    print("=" * 76)

    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
