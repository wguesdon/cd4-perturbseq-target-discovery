"""Validate the stratified tests before any conclusion is built on them.

``src/cd4_perturbseq/stratified.py`` is hand-rolled, because pulling in statsmodels for two
closed forms would be the only heavyweight dependency in the project. Hand-rolled statistics are
exactly the kind of thing that is quietly wrong, so this script tries to break them.

Three tiers, in increasing severity.

**Against scipy.** With a single stratum, van Elteren must equal the tie-corrected
Mann-Whitney normal approximation, and the CMH statistic must equal the Yates-corrected Pearson
chi-square scaled by ``(n - 1) / n``. That scaling is not a fudge: Mantel-Haenszel uses the
hypergeometric variance, whose denominator is ``n^2 (n - 1)``, while Yates uses ``n^3``. Writing
this check the naive way makes it fail by 0.1%, which is how it was caught.

**Against brute force.** Shuffle the group label *within* each stratum, twenty thousand times.
That is the exact null both tests claim to test, and a closed form that disagrees with it is
wrong however elegant it looks.

**Against the null.** Simulate four hundred datasets with a stratum effect and no group effect.
A test whose type-I error is not near its nominal 0.05 cannot be used to reject anything.

Usage:
    uv run python scripts/12a_validate_stratified.py
"""

from __future__ import annotations

import argparse
from collections.abc import Callable

import numpy as np
from scipy import stats

from cd4_perturbseq.stratified import benjamini_hochberg, cochran_mantel_haenszel, van_elteren

RNG = np.random.default_rng(12345)
PERM_DRAWS = 20_000
NULL_TRIALS = 400


class Ledger:
    """Accumulate pass/fail outcomes and print them uniformly."""

    def __init__(self) -> None:
        self.outcomes: list[bool] = []

    def close(self, name: str, got: float, want: float, tol: float = 1e-6) -> None:
        """Assert two numbers agree.

        Args:
            name: What is being checked.
            got: Observed value.
            want: Reference value.
            tol: Absolute tolerance.
        """
        ok = bool(abs(got - want) < tol)
        self.outcomes.append(ok)
        print(f"  {name:54s} got {got:.8g}  want {want:.8g}  {'PASS' if ok else 'FAIL'}")

    def claim(self, name: str, ok: bool) -> None:
        """Assert a boolean claim.

        Args:
            name: What is being checked.
            ok: Whether it holds.
        """
        self.outcomes.append(bool(ok))
        print(f"  {name:54s} {'PASS' if ok else 'FAIL'}")


def _permutation_p(
    values: np.ndarray,
    in_top: np.ndarray,
    strata: np.ndarray,
    statistic: Callable[[np.ndarray, np.ndarray, np.ndarray], float],
    n_perm: int,
) -> float:
    """Monte-Carlo p-value under within-stratum permutation of the group label.

    Args:
        values: Outcome per row.
        in_top: Boolean group indicator.
        strata: Stratum label per row.
        statistic: Maps ``(values, in_top, strata)`` to a scalar, larger meaning "top is bigger".
        n_perm: Number of permutations.

    Returns:
        The one-sided Monte-Carlo p-value, with the usual plus-one correction.
    """
    observed = statistic(values, in_top, strata)
    rng = np.random.default_rng(7)
    shuffled = in_top.copy()
    blocks = [np.flatnonzero(strata == s) for s in np.unique(strata)]
    at_least = 0
    for _ in range(n_perm):
        for block in blocks:
            shuffled[block] = rng.permutation(in_top[block])
        if statistic(values, shuffled, strata) >= observed:
            at_least += 1
    return (at_least + 1) / (n_perm + 1)


def _van_elteren_z(values: np.ndarray, in_top: np.ndarray, strata: np.ndarray) -> float:
    """Van Elteren standard normal deviate, for the permutation harness."""
    return van_elteren(values, in_top, strata, "greater").statistic


def _cmh_z(values: np.ndarray, in_top: np.ndarray, strata: np.ndarray) -> float:
    """CMH standard normal deviate, for the permutation harness."""
    return cochran_mantel_haenszel(values, in_top, strata, "greater").statistic


def against_scipy(ledger: Ledger) -> None:
    """Single-stratum reductions must reproduce scipy exactly.

    Args:
        ledger: Outcome accumulator.
    """
    print("=== van Elteren, one stratum, no ties, vs the Mann-Whitney normal approximation ===")
    x = RNG.normal(size=300)
    top = np.zeros(300, dtype=bool)
    top[:60] = True
    x[:60] += 0.5
    result = van_elteren(x, top, np.zeros(300), "greater")
    _, reference = stats.mannwhitneyu(
        x[top], x[~top], alternative="greater", use_continuity=False, method="asymptotic"
    )
    ledger.close("van Elteren p == mannwhitneyu p", result.pvalue, reference)

    print("\n=== van Elteren, one stratum, heavy ties, vs tie-corrected Mann-Whitney ===")
    xt = RNG.integers(0, 4, size=400).astype(float)
    topt = np.zeros(400, dtype=bool)
    topt[:100] = True
    xt[:100] += 1
    tied = van_elteren(xt, topt, np.zeros(400), "greater")
    _, reference_tied = stats.mannwhitneyu(
        xt[topt], xt[~topt], alternative="greater", use_continuity=False, method="asymptotic"
    )
    ledger.close("van Elteren p == tie-corrected mannwhitneyu p", tied.pvalue, reference_tied)

    print("\n=== uninformative strata are dropped, not silently counted ===")
    xs = RNG.normal(size=200)
    tops = np.zeros(200, dtype=bool)
    tops[:50] = True
    strata = np.concatenate([np.zeros(100), np.ones(100)])
    dropped = van_elteren(xs, tops, strata, "greater")
    ledger.close("a stratum with zero top rows is dropped", float(dropped.n_strata_dropped), 1.0)
    ledger.close("only the informative stratum is used", float(dropped.n_strata_used), 1.0)

    print("\n=== CMH, one stratum: the MH statistic is Yates chi-square scaled by (n-1)/n ===")
    print("    Mantel-Haenszel uses the hypergeometric variance, denominator n^2 (n-1).")
    print("    Yates-corrected Pearson uses n^3. They are not equal, and that is correct.")
    flag = np.array([1] * 25 + [0] * 75 + [1] * 100 + [0] * 800, dtype=float)
    group = np.array([True] * 100 + [False] * 900)
    n = flag.size
    cmh = cochran_mantel_haenszel(flag, group, np.zeros(n), "two-sided")
    chi2_yates, _, _, _ = stats.chi2_contingency(np.array([[25, 75], [100, 800]]), correction=True)
    ledger.close("CMH z^2 == Yates chi2 * (n-1)/n", cmh.statistic**2, chi2_yates * (n - 1) / n)
    ledger.close("MH odds ratio == crude OR when there is one stratum", cmh.effect, (25 * 800) / (75 * 100), 1e-9)

    print("\n=== CMH recovers a common odds ratio where the crude one is confounded ===")
    flags: list[float] = []
    groups: list[bool] = []
    labels: list[float] = []
    for stratum, (a, b, c, d) in enumerate([(20, 20, 10, 30), (10, 90, 2, 98)]):
        flags += [1.0] * a + [0.0] * b + [1.0] * c + [0.0] * d
        groups += [True] * (a + b) + [False] * (c + d)
        labels += [float(stratum)] * (a + b + c + d)
    f, g, s = np.array(flags), np.array(groups), np.array(labels)
    crude = (f[g].sum() * (1 - f[~g]).sum()) / ((1 - f[g]).sum() * f[~g].sum())
    pooled = cochran_mantel_haenszel(f, g, s, "greater")
    stratum_ors = [(20 * 30) / (20 * 10), (10 * 98) / (90 * 2)]
    print(
        f"    crude OR {crude:.3f}   stratum ORs {stratum_ors[0]:.3f}, {stratum_ors[1]:.3f}   "
        f"MH OR {pooled.effect:.3f}"
    )
    ledger.claim("MH OR lies between the stratum ORs", min(stratum_ors) <= pooled.effect <= max(stratum_ors))


def against_brute_force(ledger: Ledger, n_perm: int) -> None:
    """Both closed forms must agree with a within-stratum permutation null.

    Args:
        ledger: Outcome accumulator.
        n_perm: Number of permutations.
    """
    print(f"\n=== permutation validation: shuffle the group label WITHIN each stratum, {n_perm:,}x ===")
    print("    This is the exact null. A wrong closed form cannot survive it.")

    size = 240
    strata = np.repeat([0.0, 1.0, 2.0], size // 3)
    group = np.zeros(size, dtype=bool)
    for index, fraction in enumerate([0.1, 0.3, 0.6]):
        rows = np.flatnonzero(strata == index)
        group[rows[: int(len(rows) * fraction)]] = True

    # The stratum shifts the outcome AND the group prevalence: exactly the confounding the
    # stratified tests exist to absorb.
    shift = np.array([0.0, 0.7, 1.4])[strata.astype(int)]
    continuous = RNG.normal(loc=shift + 0.45 * group, scale=1.0)
    binary = (RNG.random(size) < 0.15 + 0.35 * group).astype(float)

    closed = van_elteren(continuous, group, strata, "greater").pvalue
    brute = _permutation_p(continuous, group, strata, _van_elteren_z, n_perm)
    print(f"    van Elteren  closed-form p = {closed:.5f}   permutation p = {brute:.5f}")
    ledger.claim("van Elteren agrees with permutation (within 0.02)", abs(closed - brute) < 0.02)

    closed_cmh = cochran_mantel_haenszel(binary, group, strata, "greater").pvalue
    brute_cmh = _permutation_p(binary, group, strata, _cmh_z, n_perm)
    print(f"    CMH          closed-form p = {closed_cmh:.5f}   permutation p = {brute_cmh:.5f}")
    ledger.claim("CMH agrees with permutation (within 0.02)", abs(closed_cmh - brute_cmh) < 0.02)


def against_the_null(ledger: Ledger, trials: int) -> None:
    """Type-I error must sit near the nominal level when there is no group effect.

    Args:
        ledger: Outcome accumulator.
        trials: Number of simulated null datasets.
    """
    print(f"\n=== null calibration: {trials} datasets with a stratum effect and NO group effect ===")
    size = 240
    strata = np.repeat([0.0, 1.0, 2.0], size // 3)
    shift = np.array([0.0, 0.7, 1.4])[strata.astype(int)]
    rng = np.random.default_rng(99)

    false_positives_ve = 0
    false_positives_cmh = 0
    for _ in range(trials):
        group = np.zeros(size, dtype=bool)
        for index in range(3):
            rows = np.flatnonzero(strata == index)
            group[rng.choice(rows, size=len(rows) // 4, replace=False)] = True
        false_positives_ve += van_elteren(rng.normal(loc=shift), group, strata, "greater").pvalue < 0.05
        binary = (rng.random(size) < 0.2).astype(float)
        false_positives_cmh += cochran_mantel_haenszel(binary, group, strata, "greater").pvalue < 0.05

    rate_ve = false_positives_ve / trials
    rate_cmh = false_positives_cmh / trials
    print(f"    van Elteren type-I at alpha = 0.05: {rate_ve:.3f}")
    print(f"    CMH         type-I at alpha = 0.05: {rate_cmh:.3f}   (continuity correction makes it conservative)")
    ledger.claim("van Elteren type-I within [0.02, 0.09]", 0.02 <= rate_ve <= 0.09)
    ledger.claim("CMH type-I at or below 0.09, never anti-conservative", rate_cmh <= 0.09)


def multiplicity(ledger: Ledger) -> None:
    """The Benjamini-Hochberg wrapper must not deviate from scipy.

    Args:
        ledger: Outcome accumulator.
    """
    print("\n=== Benjamini-Hochberg ===")
    raw = [0.001, 0.02, 0.03, 0.5]
    adjusted = benjamini_hochberg(raw)
    reference = stats.false_discovery_control(np.asarray(raw), method="bh")
    ledger.close("BH matches scipy", float(np.abs(adjusted - reference).max()), 0.0, 1e-12)
    ledger.claim("BH passes NaN through", bool(np.isnan(benjamini_hochberg([0.01, float("nan"), 0.2])[1])))


def main() -> None:
    """Run every validation tier and fail the process if any check fails."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-perm", type=int, default=PERM_DRAWS)
    parser.add_argument("--n-null", type=int, default=NULL_TRIALS)
    args = parser.parse_args()

    ledger = Ledger()
    against_scipy(ledger)
    against_brute_force(ledger, args.n_perm)
    against_the_null(ledger, args.n_null)
    multiplicity(ledger)

    passed = sum(ledger.outcomes)
    total = len(ledger.outcomes)
    print("\n" + "=" * 76)
    if passed == total:
        print(f"VERDICT: the stratified tests are sound. {passed}/{total} checks passed.")
        print("They may be used to draw conclusions.")
    else:
        print(f"VERDICT: {total - passed} of {total} checks FAILED. Do not trust any result built on this module.")
    print("=" * 76)

    raise SystemExit(0 if passed == total else 1)


if __name__ == "__main__":
    main()
