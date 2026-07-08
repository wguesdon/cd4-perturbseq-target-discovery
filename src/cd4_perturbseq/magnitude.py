"""Transcriptome-wide effect magnitude, the covariate the risk-kill test should have matched on.

``scripts/02_risk_kill_reversal.py`` originally matched its background on ``n_cells_target``,
reasoning that perturbations assayed in more cells have more statistical power and therefore more
significant DE genes. The data say the opposite: ``Spearman(n_cells_target, stim_de_genes) =
-0.243``. Cell count is a viability readout, not a power proxy, and matching on it controls the
confound backwards.

The confound that actually exists is effect magnitude. A perturbation that moves the whole
transcriptome gets a large ``|mean z|`` over *any* gene set, the effector module included, and it
crosses the DE threshold on many genes for the same reason.

:func:`effect_magnitude` measures that directly, and excludes three things that would otherwise put
the score inside the variable meant to control for it:

- the perturbed gene's own column, whose knockdown z is large by construction because
  ``ontarget_significant`` is a QC requirement,
- the effector-module genes the suppression score is built from,
- the tolerance-module genes the tolerance pillar is built from.

The result is ``rho = +0.198`` with ``|naive_suppression|`` and ``rho = +0.069`` with tolerance
suppression, so it is a genuine covariate rather than a restatement of either.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd

from . import de_stats
from .paths import INTERIM, ensure_dirs

DEFAULT_CACHE = INTERIM / "z_l2_stim48.csv"
"""Where :func:`effect_magnitude` memoises itself. Gitignored; regenerated from the h5ad."""


def effect_magnitude(
    obs: pd.DataFrame,
    var: pd.DataFrame,
    rows: np.ndarray,
    exclude: Iterable[str] = (),
    cache: Path | None = DEFAULT_CACHE,
    verbose: bool = True,
) -> pd.DataFrame:
    """Per-perturbation transcriptome-wide effect magnitude, off-module and off-target.

    ``z_l2 = sqrt(sum of z^2)`` over the measured genes. This is the only step in the project that
    reads full-width rows of a 16.8 GB layer, so it is cached.
    :func:`~cd4_perturbseq.de_stats.read_layer_rows` downcasts to float32, which costs 262 MB for
    the 6,371 QC-passing stimulated perturbations.

    Args:
        obs: Full ``.obs`` frame, positionally indexed.
        var: Full ``.var`` frame.
        rows: Positional indices into ``.obs`` of the perturbations to score.
        exclude: Gene symbols dropped from the norm, beyond each row's own target gene.
        cache: Where to memoise. Reused only when its row count matches ``rows``. None disables.
        verbose: Print progress.

    Returns:
        DataFrame indexed by gene symbol with ``z_l2`` (off-module, off-target), ``z_l2_raw``
        (nothing excluded), and ``n_genes_used``.
    """
    if cache is not None and cache.exists():
        cached = pd.read_csv(cache)
        if len(cached) == len(rows):
            if verbose:
                print(f"  z_l2: reusing cache {cache} ({len(cached):,} rows)")
            return cached.set_index("gene_name")
        if verbose:
            print(f"  z_l2: cache has {len(cached):,} rows, need {len(rows):,}; recomputing")

    names = var["gene_name"].astype(str).to_numpy()
    column_of = {gene: index for index, gene in enumerate(names)}

    if verbose:
        print(f"  z_l2: reading {len(rows):,} full-width rows of the zscore layer (float32, ~262 MB) ...")
    block = de_stats.read_layer_rows(rows, layer="zscore")
    finite = np.isfinite(block)
    squared = np.where(finite, block, 0.0).astype(np.float64) ** 2
    z_l2_raw = np.sqrt(squared.sum(axis=1))

    module_cols = np.array(sorted({column_of[g] for g in exclude if g in column_of}), dtype=np.int64)
    if module_cols.size:
        squared[:, module_cols] = 0.0
        finite[:, module_cols] = False

    targets = obs.loc[rows, "target_contrast_gene_name"].astype(str).to_numpy()
    self_rows = np.array([i for i, gene in enumerate(targets) if gene in column_of], dtype=np.int64)
    self_cols = np.array([column_of[gene] for gene in targets if gene in column_of], dtype=np.int64)
    if self_rows.size:
        squared[self_rows, self_cols] = 0.0
        finite[self_rows, self_cols] = False
    if verbose:
        print(f"  z_l2: excluded {module_cols.size} module columns and {self_rows.size:,} self-target columns")

    frame = pd.DataFrame(
        {
            "gene_name": targets,
            "z_l2": np.sqrt(squared.sum(axis=1)),
            "z_l2_raw": z_l2_raw,
            "n_genes_used": finite.sum(axis=1),
        }
    )
    if cache is not None:
        ensure_dirs()
        frame.to_csv(cache, index=False)
        if verbose:
            print(f"  z_l2: wrote cache {cache}")
    return frame.set_index("gene_name")
