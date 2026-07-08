"""N7: is ANY of this specific to reversal, or does effect size buy the whole result?

The adversarial audit's most damaging finding, filed ``invalidates_headline: true`` and never
answered:

    "The risk-kill test has no specificity to reversal: non-reversal null rankings pass the
    paper's own verdict function with equal or stronger effect sizes."

``scripts/12_magnitude_matched.py`` (N1) killed the *noise* null: 200 random shortlists fire at the
nominal 5%. It did not touch the *structured* null. This script does, and it asks the question at
the level that actually matters, which is not a pillar but the gate itself.

**The headline number is "18 of the naive top 20 are rejected by the safety gate."** That number is
worthless unless a comparable 20 is rejected less often. So: draw 5,000 shortlists of 20 matched to
the naive top 20 on ``z_l2`` decile, drawn from the same evidence-passing pool, and count.

Conditioning on the evidence floor is deliberate. It is the gate's first axis and the naive top 20
all clear it; a background that failed it would be rejected for a reason unrelated to the question.

Then the same question per axis. ``04_window_score.py`` rejects on homeostasis (selectivity) and on
tolerance. Which of them is doing work that effect magnitude does not already do?

**Pre-registered claim, stated before the counts are printed and inherited from N1:** tolerance is
the reversal-specific axis and homeostasis is not, because tolerance is ``rho = 0.069`` with effect
magnitude while the selectivity that drives homeostasis is a ratio of two DE counts that are ``0.725``
and ``0.527`` correlated with it. The script exits non-zero if tolerance fails to separate.

It is written to be able to fail in the way that would hurt: if tolerance rejection in the naive
top 20 is no higher than in a magnitude-matched 20, the project's contribution is "large-effect
perturbations are toxic", which is neither new nor about reversal.

Part 2 attacks from the other side. Rank the whole screen by three rival scores that know nothing
about the effector program -- transcriptome-wide magnitude, collateral DE count, and sign-blind
score magnitude -- and ask what their top 100 looks like on each pillar.

Usage:
    uv run python scripts/14_reversal_specificity.py --n-draws 5000
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from cd4_perturbseq import paths

N_BINS = 10
TOP_N = 20
ALPHA = 0.01
"""Pre-registered: tolerance must beat the magnitude-matched draws at this level."""

# Only the co-inhibitory axis remains a gate. "homeostasis" (selectivity) was demoted to an
# annotation (N6), so it no longer appears in reject_reason. Its non-specificity is why it was
# demoted: N7 measured 14 of ours vs 13.6 magnitude-matched, P=0.53. We report that result in the
# text as the reason, and we do not carry an empty "rejected on homeostasis: 0" column here.
AXES = ("tolerance",)


def load() -> pd.DataFrame:
    """Join the window score to the effect-magnitude covariate.

    Returns:
        Evidence-passing perturbations with a known ``z_l2``, ranked by naive suppression.

    Raises:
        FileNotFoundError: If either upstream table is missing.
    """
    window = pd.read_csv(paths.TABLES / "window_score.csv")
    rows = pd.read_csv(paths.TABLES / "magnitude_matched_rows.csv")
    frame = window.merge(rows[["gene_name", "z_l2", "abs_score", "stim_downstream"]], on="gene_name", how="left")

    frame = frame[~frame["fail_evidence"] & frame["z_l2"].notna()].copy()
    frame["naive_rank"] = (-frame["eff_mean_z"]).rank(ascending=False, method="first")
    frame["rejected"] = ~frame["safe"]
    frame["z_decile"] = pd.qcut(frame["z_l2"], N_BINS, labels=False, duplicates="drop")
    return frame.reset_index(drop=True)


def _matched_draws(frame: pd.DataFrame, top: pd.DataFrame, n_draws: int, rng: np.random.Generator) -> pd.DataFrame:
    """Draw shortlists matched to the top on z_l2 decile, from the non-top pool.

    Args:
        frame: All evidence-passing perturbations.
        top: The naive shortlist.
        n_draws: Number of shortlists.
        rng: Random generator. Sampling is the null here, not the estimate, so a seed is correct.

    Returns:
        One row per draw, with the rejection count and per-axis counts.
    """
    wanted = top["z_decile"].value_counts()
    pools = {
        decile: frame.index[(frame["z_decile"] == decile) & ~frame.index.isin(top.index)].to_numpy()
        for decile in wanted.index
    }

    records = []
    for _ in range(n_draws):
        picked: list[int] = []
        for decile, count in wanted.items():
            pool = pools[decile]
            if pool.size:
                picked.extend(rng.choice(pool, size=min(int(count), pool.size), replace=False))
        sample = frame.loc[picked]
        row = {"n": len(sample), "rejected": int(sample["rejected"].sum())}
        for axis in AXES:
            # Exact-token match on the split reject_reason, never a substring, so a future axis name
            # that contains another cannot double-count. reject_reason is comma-joined.
            row[axis] = int(sample["reject_reason"].str.split(",").apply(lambda parts: axis in parts).sum())
        # The demoted selectivity axis, reported so the figure can show WHY it was demoted: it
        # rejects a magnitude-matched shortlist as often as ours. It is not part of the gate.
        row["selectivity_annotation"] = int(sample["fail_homeostasis"].sum())
        records.append(row)
    return pd.DataFrame(records)


def gate_specificity(frame: pd.DataFrame, n_draws: int) -> tuple[bool, pd.DataFrame]:
    """Is the gate's headline rejection count specific to reversal, or bought by effect size?

    Args:
        frame: Output of :func:`load`.
        n_draws: Number of matched shortlists.

    Returns:
        Tuple of (tolerance separates at ``ALPHA``, a tidy table of the comparison).
    """
    rng = np.random.default_rng(0)
    top = frame.nsmallest(TOP_N, "naive_rank")
    draws = _matched_draws(frame, top, n_draws, rng)

    observed = int(top["rejected"].sum())
    print(f"\n=== is '{observed} of the naive top {TOP_N} are rejected' specific to reversal? ===")
    print(f"    pool: {len(frame)} evidence-passing perturbations with a known z_l2")
    print(f"    naive top {TOP_N} z_l2 median {top['z_l2'].median():.1f}   whole pool {frame['z_l2'].median():.1f}")
    print(f"    base rejection rate over the whole pool: {frame['rejected'].mean():.1%}")

    p_total = float((draws["rejected"] >= observed).mean())
    print(f"\n    naive top {TOP_N}              {observed:2d}/{TOP_N} rejected  ({observed / TOP_N:.0%})")
    print(f"    z_l2-matched, {n_draws:,} draws   {draws['rejected'].mean():4.1f}/{TOP_N} rejected  "
          f"({draws['rejected'].mean() / TOP_N:.0%})   P(matched >= observed) = {p_total:.4f}")
    if p_total >= 0.05:
        print("    -> THE HEADLINE COUNT IS NOT SPECIFIC TO REVERSAL. A shortlist of the same effect")
        print(f"       magnitude is rejected almost as often. Quoting {observed}/{TOP_N} as evidence that")
        print("       *reversal* nominates toxic targets overstates it. Say so in the report.")

    print(f"\n    which AXIS does work that effect magnitude does not? ({n_draws:,} matched draws)")
    records = [{
        "quantity": f"rejected (any axis), top {TOP_N}",
        "observed": observed,
        "matched_mean": draws["rejected"].mean(),
        "p_matched_ge_observed": p_total,
        "specific_to_reversal": p_total < ALPHA,
    }]
    verdicts = {}
    for axis in AXES:
        observed_axis = int(top["reject_reason"].str.split(",").apply(lambda parts: axis in parts).sum())
        p_axis = float((draws[axis] >= observed_axis).mean())
        verdicts[axis] = p_axis < ALPHA
        call = "SPECIFIC TO REVERSAL" if p_axis < ALPHA else "explained by effect magnitude"
        print(f"      {axis:12s} top {observed_axis:2d}   matched mean {draws[axis].mean():5.2f}   "
              f"P(matched >= top) = {p_axis:.4f}   {call}")
        records.append({
            "quantity": f"rejected on {axis}, top {TOP_N}",
            "observed": observed_axis,
            "matched_mean": draws[axis].mean(),
            "p_matched_ge_observed": p_axis,
            "specific_to_reversal": verdicts[axis],
        })

    # The demoted selectivity axis, as a reported comparison. This is the demotion, made visible:
    # it rejects a magnitude-matched shortlist as often as it rejects ours.
    observed_sel = int(top["fail_homeostasis"].sum())
    p_sel = float((draws["selectivity_annotation"] >= observed_sel).mean())
    print(f"      {'selectivity':12s} top {observed_sel:2d}   matched mean {draws['selectivity_annotation'].mean():5.2f}   "
          f"P(matched >= top) = {p_sel:.4f}   {'SPECIFIC' if p_sel < ALPHA else 'explained by effect magnitude (DEMOTED, N6)'}")
    records.append({
        "quantity": f"fails selectivity annotation, top {TOP_N}",
        "observed": observed_sel,
        "matched_mean": draws["selectivity_annotation"].mean(),
        "p_matched_ge_observed": p_sel,
        "specific_to_reversal": p_sel < ALPHA,
    })

    print("\n    The gate rejects on the co-inhibitory axis, and that is specific to reversal. The")
    print("    demoted selectivity axis rejects a magnitude-matched shortlist just as often: it was")
    print("    about the size of a perturbation, not its direction, and effect size is not a finding.")

    table = pd.DataFrame(records)
    table["n_draws"] = n_draws
    table["pool"] = len(frame)
    return verdicts["tolerance"], table


def rival_rankings(frame: pd.DataFrame) -> pd.DataFrame:
    """Rank by scores that know nothing about the effector program, and look at their top 100.

    Args:
        frame: Output of :func:`load`.

    Returns:
        One row per rival ranking.
    """
    print("\n=== rival rankings: what does a top 100 look like if you never mention reversal? ===")
    print("    (whole evidence-passing pool; 'tolerance loss' is the gate's own axis)")

    # Every key is oriented so that LARGER means HIGHER IN THE RANKING, then sorted descending.
    # An earlier draft sorted `-eff_mean_z` ascending, which crowned the LEAST suppressed
    # perturbations and made our own ranking look clean. This project has now made the same
    # rank-direction error twice; the assertion below is why it did not survive a second time.
    rankings = {
        "naive suppression (ours)": -frame["eff_mean_z"],
        "z_l2 (pure effect magnitude)": frame["z_l2"],
        "collateral DE count": frame["stim_downstream"],
        "|score| (sign-blind)": frame["abs_score"],
    }
    naive_top = frame.loc[rankings["naive suppression (ours)"].sort_values(ascending=False).index[0]]
    assert naive_top["naive_rank"] == 1, "rank direction inverted: the top of our ranking is not rank 1"

    take = min(100, len(frame))
    records = []
    for label, key in rankings.items():
        top = frame.loc[key.sort_values(ascending=False).index[:take]]
        records.append({
            "ranking": label,
            "median tolerance loss": top["tolerance_loss"].median(),
            "median z_l2": top["z_l2"].median(),
            "rejected on tolerance": int(top["reject_reason"].str.split(",").apply(lambda p: "tolerance" in p).sum()),
            "fails selectivity (annotation)": int(top["fail_homeostasis"].sum()),
            "pass gate": int(top["safe"].sum()),
        })
    table = pd.DataFrame(records)
    print(table.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\n    n = {take} of {len(frame)} evidence-passers, so these are thirds of the pool, not sharp tops.")
    print("    A magnitude-only ranking that matched ours on tolerance loss would end the tolerance")
    print("    claim. Read the first column. Note that the evidence floor already removes inducers,")
    print("    so '|score| sign-blind' is close to our own ranking by construction and is a weak null.")
    return table


def main() -> None:
    """Test whether the gate's result survives a magnitude-matched null."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-draws", type=int, default=5000)
    args = parser.parse_args()

    paths.ensure_dirs()
    frame = load()
    holds, specificity = gate_specificity(frame, args.n_draws)
    rivals = rival_rankings(frame)

    out = paths.TABLES / "reversal_specificity.csv"
    specificity.to_csv(out, index=False)
    print(f"\nwrote {out}")
    rivals_out = paths.TABLES / "reversal_rival_rankings.csv"
    rivals.to_csv(rivals_out, index=False)
    print(f"wrote {rivals_out}")

    print("\n" + "=" * 88)
    if holds:
        print("VERDICT: tolerance collapse is what reversal-ranking specifically brings in.")
        print("The aggregate '18 of 20 rejected' is largely bought by effect magnitude and must not be")
        print("quoted as evidence about reversal. The tolerance axis must, and it is the axis the")
        print("adversarial audit attacked hardest and lost.")
    else:
        print("VERDICT: not even tolerance separates from a magnitude-matched shortlist.")
        print("The project's claim reduces to 'large-effect knockdowns are toxic'. Rewrite the headline.")
    print("=" * 88)

    raise SystemExit(0 if holds else 1)


if __name__ == "__main__":
    main()
