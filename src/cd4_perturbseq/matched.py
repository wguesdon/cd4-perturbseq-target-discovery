"""Effect-magnitude-matched shortlist draws, in one place.

Both the reversal-specificity test (``scripts/14``) and the prior-knowledge nomination
(``scripts/20``) ask the same question: is a property of an observed shortlist more common than in
shortlists of the same effect-magnitude profile? The null is a set of shortlists drawn to match the
observed one on a stratum (a ``z_l2`` decile), from the rest of the universe.

The earlier inline version in ``scripts/14`` truncated silently: it drew
``size=min(count, pool.size)`` per stratum, so when a stratum's pool was smaller than the observed
count the draw came back short and the comparison was against a smaller shortlist than intended. At
a top-20 shortlist against a 286-row pool that never bit; at a top-100 shortlist it does. This
version keeps the requested per-stratum count, sampling with replacement only where the pool forces
it, and reports the shortfall so a caller can assert the realised size.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


def draw_matched_indices(
    universe_index: pd.Index,
    strata: pd.Series,
    top_index: pd.Index,
    n_draws: int,
    rng: np.random.Generator,
) -> tuple[list[np.ndarray], int]:
    """Draw ``n_draws`` stratum-matched index sets from the universe, excluding the observed top.

    Each draw contains, for every stratum the observed shortlist occupies, the same number of rows
    the shortlist has in that stratum. Rows are drawn without replacement where the eligible pool is
    large enough, and with replacement only for a stratum whose pool is smaller than the requested
    count (recorded as a shortfall).

    Args:
        universe_index: Index of every eligible row (the full ranked universe).
        strata: Stratum label per row, indexed by ``universe_index``.
        top_index: Index of the observed shortlist (excluded from every pool).
        n_draws: Number of matched shortlists to draw.
        rng: Random generator. Sampling is the null, not the estimate, so a fixed seed is correct.

    Returns:
        Tuple of (list of ``n_draws`` index arrays each summing to the observed shortlist size, total
        number of stratum-slots that had to be drawn with replacement across all draws).
    """
    wanted = strata.loc[top_index].value_counts()
    pools = {
        stratum: universe_index[(strata == stratum) & ~universe_index.isin(top_index)].to_numpy()
        for stratum in wanted.index
    }

    draws: list[np.ndarray] = []
    replacement_slots = 0
    for _ in range(n_draws):
        picked: list = []
        for stratum, count in wanted.items():
            pool = pools[stratum]
            count = int(count)
            if pool.size == 0:
                # No eligible rows in this stratum at all; the count cannot be honoured. This is a
                # genuine coverage gap and the shortfall is recorded so the caller can assert on it.
                replacement_slots += count
                continue
            replace = pool.size < count
            if replace:
                replacement_slots += count - pool.size
            picked.extend(rng.choice(pool, size=count, replace=replace))
        draws.append(np.array(picked))
    return draws, replacement_slots


def count_matched(
    frame: pd.DataFrame,
    strata: pd.Series,
    top_index: pd.Index,
    predicates: dict[str, Sequence[bool] | pd.Series],
    n_draws: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, int]:
    """Count, over matched draws, how many rows satisfy each named boolean predicate.

    Args:
        frame: The universe, indexed to match ``strata``.
        strata: Stratum label per row.
        top_index: The observed shortlist.
        predicates: Mapping of name to a boolean mask over ``frame`` (aligned to its index).
        n_draws: Number of matched draws.
        rng: Random generator.

    Returns:
        Tuple of (one row per draw with the size ``n`` and a count column per predicate, the total
        replacement-slot shortfall from :func:`draw_matched_indices`).
    """
    draws, shortfall = draw_matched_indices(frame.index, strata, top_index, n_draws, rng)
    masks = {name: pd.Series(mask, index=frame.index).astype(bool) for name, mask in predicates.items()}
    records = []
    for idx in draws:
        row = {"n": len(idx)}
        for name, mask in masks.items():
            row[name] = int(mask.loc[idx].sum())
        records.append(row)
    return pd.DataFrame(records), shortfall
