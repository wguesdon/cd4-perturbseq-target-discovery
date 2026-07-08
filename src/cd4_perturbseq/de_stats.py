"""Lazy accessors for the 16.8 GB `GWCD4i.DE_stats.h5ad` effect matrix.

The file holds six float64 layers of shape (33983 perturbation-condition pairs,
10282 measured genes). Materialising even one layer costs ~2.8 GB, so we never load
a whole layer. AnnData's ``backed`` mode only lazily backs ``.X``, not ``.layers``,
so we read the layers through ``h5py`` directly and use ``anndata.io.read_elem``
only for the small ``.obs`` and ``.var`` frames.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

from .paths import DE_STATS_H5AD

LAYERS = ("log_fc", "p_value", "adj_p_value", "baseMean", "lfcSE", "zscore")


def read_obs(path: Path = DE_STATS_H5AD) -> pd.DataFrame:
    """Read the perturbation-condition annotation frame.

    Args:
        path: Path to the DE stats h5ad.

    Returns:
        DataFrame with one row per (perturbed gene, culture condition) pair.
    """
    with h5py.File(path, "r") as handle:
        return read_elem(handle["obs"])


def read_var(path: Path = DE_STATS_H5AD) -> pd.DataFrame:
    """Read the measured-gene annotation frame.

    Args:
        path: Path to the DE stats h5ad.

    Returns:
        DataFrame with one row per measured gene.
    """
    with h5py.File(path, "r") as handle:
        return read_elem(handle["var"])


def qc_mask(obs: pd.DataFrame) -> pd.Series:
    """Routine perturbation QC, the filter every ranking in this project is built on.

    Drops perturbations with an off-target flag or a knocked-down neighbouring gene, and keeps
    only those with a significant on-target knockdown against an expressed target. Ranking a
    perturbation that never knocked its target down would be ranking noise.

    This was copy-pasted into ``scripts/04_window_score.py`` and ``scripts/10_tolerance_is_real.py``
    and lives here so the three agree by construction. Note that it conditions on a collider:
    ``ontarget_significant`` requires the target to be measured and expressed, which removes
    nonessential genes preferentially. See ``scripts/11_selection_funnel.py``.

    Args:
        obs: The perturbation-condition annotation frame, or any subset of its rows.

    Returns:
        Boolean series, True where the perturbation passes QC.
    """
    keep = ~obs["distal_offtarget_flag"].astype(bool)
    keep &= ~obs["neighboring_gene_KD"].astype(bool)
    keep &= obs["ontarget_significant"].astype(bool)
    keep &= ~obs["low_target_gex"].astype(bool)
    return keep


def describe(path: Path = DE_STATS_H5AD) -> dict[str, object]:
    """Summarise the on-disk structure without loading any layer.

    Args:
        path: Path to the DE stats h5ad.

    Returns:
        Mapping with the layer names, matrix shape, dtype, and HDF5 chunk shape.
    """
    with h5py.File(path, "r") as handle:
        layers = sorted(handle["layers"].keys())
        probe = handle["layers"][layers[0]]
        return {
            "layers": layers,
            "shape": probe.shape,
            "dtype": str(probe.dtype),
            "chunks": probe.chunks,
            "compression": probe.compression,
            "varm": sorted(handle["varm"].keys()) if "varm" in handle else [],
        }


def read_layer_columns(
    gene_indices: Sequence[int],
    layer: str = "zscore",
    path: Path = DE_STATS_H5AD,
) -> np.ndarray:
    """Read selected gene columns of a layer for every perturbation.

    HDF5 fancy indexing requires a strictly increasing selection, so the indices are
    sorted internally and the result is permuted back to the caller's order.

    Args:
        gene_indices: Positional indices into ``.var`` of the genes to read.
        layer: Layer name, one of :data:`LAYERS`.
        path: Path to the DE stats h5ad.

    Returns:
        Array of shape ``(n_obs, len(gene_indices))`` in the caller's column order.

    Raises:
        ValueError: If ``layer`` is not a known layer or ``gene_indices`` is empty.
    """
    if layer not in LAYERS:
        raise ValueError(f"unknown layer {layer!r}; expected one of {LAYERS}")
    indices = np.asarray(gene_indices, dtype=np.int64)
    if indices.size == 0:
        raise ValueError("gene_indices must be non-empty")

    order = np.argsort(indices)
    ascending = indices[order]
    if np.any(np.diff(ascending) == 0):
        raise ValueError("gene_indices must be unique")

    with h5py.File(path, "r") as handle:
        block = handle["layers"][layer][:, ascending]

    inverse = np.empty_like(order)
    inverse[order] = np.arange(order.size)
    return block[:, inverse]


def read_layer_rows(
    obs_indices: Sequence[int],
    layer: str = "zscore",
    path: Path = DE_STATS_H5AD,
    block_size: int = 512,
    dtype: str = "float32",
) -> np.ndarray:
    """Read whole rows of a layer for selected perturbations, downcast to save memory.

    The dataset is contiguous and float64, so a full row is 10,282 * 8 bytes of adjacent
    storage and row selection is efficient. Rows are read in blocks and downcast, because
    holding 6,371 rows as float64 costs 524 MB against 262 MB as float32, and the z-scores
    carry nowhere near float64 precision.

    Args:
        obs_indices: Positional indices into ``.obs`` of the perturbations to read.
        layer: Layer name, one of :data:`LAYERS`.
        path: Path to the DE stats h5ad.
        block_size: Number of rows to pull from HDF5 per read.
        dtype: Output dtype.

    Returns:
        Array of shape ``(len(obs_indices), n_vars)`` in the caller's row order.

    Raises:
        ValueError: If ``layer`` is unknown or ``obs_indices`` is empty.
    """
    if layer not in LAYERS:
        raise ValueError(f"unknown layer {layer!r}; expected one of {LAYERS}")
    rows = np.asarray(obs_indices, dtype=np.int64)
    if rows.size == 0:
        raise ValueError("obs_indices must be non-empty")

    order = np.argsort(rows)
    ascending = rows[order]

    with h5py.File(path, "r") as handle:
        dataset = handle["layers"][layer]
        out = np.empty((rows.size, dataset.shape[1]), dtype=dtype)
        for start in range(0, ascending.size, block_size):
            chunk = ascending[start : start + block_size]
            out[start : start + chunk.size] = dataset[chunk, :].astype(dtype, copy=False)

    inverse = np.empty_like(order)
    inverse[order] = np.arange(order.size)
    return out[inverse]
