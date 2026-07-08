"""Stratified tests that use every control row, instead of one seeded draw.

``scripts/02_risk_kill_reversal.py`` matches its background by sampling, once, with
``np.random.default_rng(0)``, one background row per top-K row per decile. That throws away
6,171 of 6,271 usable controls and makes the reported p-value a random variable: redrawing the
same design 2,000 times moves the IEI odds ratio over 1.87-16.12 and the p-value over
0.00032-0.129, with 11.5% of seeds landing above 0.05.

The fix is not a better seed. It is to stop sampling. Both tests below condition on the
matching variable as a stratum and pool the evidence across strata, so every control row is
used and no random number is drawn.

- :func:`cochran_mantel_haenszel` for a binary liability. The Mantel-Haenszel odds ratio, with
  the continuity-corrected CMH statistic.
- :func:`van_elteren` for a continuous liability. A stratified Wilcoxon rank-sum, weighting each
  stratum by ``1 / (n + 1)``, which is the optimal weight when the effect is a location shift on
  the probability scale rather than on the data scale.

Both drop strata that carry no information (fewer than two rows, or all rows in one group). A
stratum in which every row is a top-K row contributes nothing to a within-stratum comparison,
and silently pretending otherwise is how a magnitude-matched test can look powered when it is
not. :attr:`StratifiedResult.n_strata_used` and :attr:`StratifiedResult.n_top_used` are reported
so the caller can see how much of the top-K actually got compared.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


@dataclass(frozen=True)
class StratifiedResult:
    """Outcome of a stratified two-group test.

    Attributes:
        effect: Mantel-Haenszel odds ratio, or the rank-biserial-like z for van Elteren.
        statistic: The standard normal deviate.
        pvalue: One- or two-sided p-value, per the ``alternative`` requested.
        n_strata_used: Strata that carried information.
        n_strata_dropped: Strata dropped as uninformative.
        n_top_used: Top-group rows inside the informative strata.
        n_background_used: Background rows inside the informative strata.
    """

    effect: float
    statistic: float
    pvalue: float
    n_strata_used: int
    n_strata_dropped: int
    n_top_used: int
    n_background_used: int


def _clean(values: np.ndarray, in_top: np.ndarray, strata: np.ndarray) -> tuple[np.ndarray, ...]:
    """Drop rows with a missing value, a missing stratum, or a missing group label.

    Args:
        values: The measured quantity, or a boolean flag.
        in_top: Boolean group indicator.
        strata: Stratum label per row; NaN means unassigned.

    Returns:
        Tuple of the three arrays, filtered to complete rows.
    """
    values = np.asarray(values, dtype=np.float64)
    in_top = np.asarray(in_top, dtype=bool)
    strata = np.asarray(strata, dtype=np.float64)
    keep = np.isfinite(values) & np.isfinite(strata)
    return values[keep], in_top[keep], strata[keep]


def _one_sided(z: float, alternative: str) -> float:
    """Convert a standard normal deviate into a p-value.

    Args:
        z: Standard normal deviate, positive when the top group is larger.
        alternative: ``"greater"``, ``"less"``, or ``"two-sided"``.

    Returns:
        The p-value.

    Raises:
        ValueError: If ``alternative`` is not recognised.
    """
    if alternative == "greater":
        return float(stats.norm.sf(z))
    if alternative == "less":
        return float(stats.norm.cdf(z))
    if alternative == "two-sided":
        return float(2.0 * stats.norm.sf(abs(z)))
    raise ValueError(f"unknown alternative {alternative!r}")


def van_elteren(
    values: np.ndarray,
    in_top: np.ndarray,
    strata: np.ndarray,
    alternative: str = "greater",
) -> StratifiedResult:
    """Stratified Wilcoxon rank-sum test with van Elteren's ``1 / (n + 1)`` weights.

    Ranks are computed within each stratum, so the test asks only whether a top-group row tends
    to exceed a background row drawn from the same stratum. Ties get midranks and the variance is
    the exact permutation variance of the within-stratum rank sum, which is tie-corrected.

    This is the continuous analogue of :func:`cochran_mantel_haenszel`, and it replaces
    ``scipy.stats.mannwhitneyu`` on a subsampled background.

    Args:
        values: Continuous measurement, one per row. NaNs are dropped.
        in_top: Boolean, True for the top-K group.
        strata: Stratum label per row, for example a decile index. NaNs are dropped.
        alternative: Direction of the test for the top group.

    Returns:
        A :class:`StratifiedResult` whose ``effect`` is the standardised rank-sum excess, and
        whose ``statistic`` is the standard normal deviate.

    Raises:
        ValueError: If no stratum carries information.
    """
    values, in_top, strata = _clean(values, in_top, strata)

    numerator = 0.0
    denominator = 0.0
    used = dropped = n_top_used = n_bg_used = 0

    for label in np.unique(strata):
        rows = strata == label
        x = values[rows]
        g = in_top[rows]
        n = x.size
        m = int(g.sum())
        if n < 2 or m == 0 or m == n:
            dropped += 1
            continue

        ranks = stats.rankdata(x, method="average")
        rank_sum = float(ranks[g].sum())
        expected = m * (n + 1.0) / 2.0
        # Exact permutation variance of the rank sum, valid under ties because it is written
        # in terms of the realised midranks rather than assuming 1..n.
        spread = float(((ranks - (n + 1.0) / 2.0) ** 2).sum())
        variance = m * (n - m) / (n * (n - 1.0)) * spread

        weight = 1.0 / (n + 1.0)
        numerator += weight * (rank_sum - expected)
        denominator += weight**2 * variance

        used += 1
        n_top_used += m
        n_bg_used += n - m

    if used == 0 or denominator <= 0:
        raise ValueError("no stratum carried information; every stratum was empty or one-sided")

    z = numerator / np.sqrt(denominator)
    return StratifiedResult(
        effect=float(z),
        statistic=float(z),
        pvalue=_one_sided(z, alternative),
        n_strata_used=used,
        n_strata_dropped=dropped,
        n_top_used=n_top_used,
        n_background_used=n_bg_used,
    )


def cochran_mantel_haenszel(
    flag: np.ndarray,
    in_top: np.ndarray,
    strata: np.ndarray,
    alternative: str = "greater",
) -> StratifiedResult:
    """Cochran-Mantel-Haenszel test of a binary liability, stratified on a matching variable.

    Pools the 2x2 tables across strata. The Mantel-Haenszel odds ratio is
    ``sum(a_i d_i / n_i) / sum(b_i c_i / n_i)``, and the test statistic is the continuity-corrected
    standardised excess of ``sum(a_i)`` over its hypergeometric expectation.

    Args:
        flag: Boolean liability, one per row.
        in_top: Boolean, True for the top-K group.
        strata: Stratum label per row, for example a decile index.
        alternative: Direction of the test for the top group.

    Returns:
        A :class:`StratifiedResult` whose ``effect`` is the Mantel-Haenszel odds ratio. The odds
        ratio is ``inf`` when no discordant background pair exists, and ``nan`` when the pooled
        variance is zero.

    Raises:
        ValueError: If no stratum carries information.
    """
    flag_f, in_top, strata = _clean(np.asarray(flag, dtype=float), in_top, strata)
    flag_b = flag_f.astype(bool)

    r_sum = s_sum = 0.0
    a_sum = e_sum = v_sum = 0.0
    used = dropped = n_top_used = n_bg_used = 0

    for label in np.unique(strata):
        rows = strata == label
        f = flag_b[rows]
        g = in_top[rows]
        n = int(f.size)
        if n < 2 or g.all() or not g.any():
            dropped += 1
            continue

        a = float((f & g).sum())
        b = float((~f & g).sum())
        c = float((f & ~g).sum())
        d = float((~f & ~g).sum())

        r_sum += a * d / n
        s_sum += b * c / n
        a_sum += a
        e_sum += (a + b) * (a + c) / n
        v_sum += (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1.0))

        used += 1
        n_top_used += int(g.sum())
        n_bg_used += n - int(g.sum())

    if used == 0:
        raise ValueError("no stratum carried information; every stratum was empty or one-sided")

    odds = float("inf") if s_sum == 0 else r_sum / s_sum
    if v_sum <= 0:
        return StratifiedResult(odds, float("nan"), float("nan"), used, dropped, n_top_used, n_bg_used)

    # Continuity correction, applied toward the null so it can never inflate significance.
    excess = a_sum - e_sum
    corrected = np.sign(excess) * max(abs(excess) - 0.5, 0.0)
    z = corrected / np.sqrt(v_sum)
    return StratifiedResult(
        effect=odds,
        statistic=float(z),
        pvalue=_one_sided(z, alternative),
        n_strata_used=used,
        n_strata_dropped=dropped,
        n_top_used=n_top_used,
        n_background_used=n_bg_used,
    )


def benjamini_hochberg(pvalues: list[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values.

    A thin wrapper over :func:`scipy.stats.false_discovery_control`, present so that callers
    declare the correction rather than reaching past it. NaN p-values pass through as NaN.

    Args:
        pvalues: Raw p-values.

    Returns:
        Array of adjusted p-values, in the caller's order.
    """
    raw = np.asarray(pvalues, dtype=np.float64)
    finite = np.isfinite(raw)
    adjusted = np.full(raw.shape, np.nan)
    if finite.any():
        adjusted[finite] = stats.false_discovery_control(raw[finite], method="bh")
    return adjusted


def deciles(values: np.ndarray, n_bins: int = 10) -> np.ndarray:
    """Assign quantile-bin labels, leaving non-finite values unassigned.

    Args:
        values: The matching variable.
        n_bins: Number of quantile bins requested. Collapsed if ties prevent it.

    Returns:
        Float array of bin indices, NaN where ``values`` is not finite.
    """
    values = np.asarray(values, dtype=np.float64)
    labels = np.full(values.shape, np.nan)
    finite = np.isfinite(values)
    if finite.sum() < n_bins:
        labels[finite] = 0.0
        return labels
    import pandas as pd

    labels[finite] = pd.qcut(values[finite], n_bins, labels=False, duplicates="drop")
    return labels
