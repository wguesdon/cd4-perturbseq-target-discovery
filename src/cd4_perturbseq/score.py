"""Calibrated gene-module scoring for perturbation effect matrices.

The naive way to ask "does knocking down gene G suppress the effector program?" is to average
the z-scores of the program's genes in G's perturbation. That statistic is wrong twice.

**It is uncalibrated.** The effector genes are co-regulated, so the effective number of
independent observations is far below the module size. Averaging 32 correlated z-scores has a
null standard deviation much larger than the 1/sqrt(32) an independence assumption implies. A
weakly-powered perturbation whose z-scores drift mildly negative therefore scores as highly as a
real hit. Empirically, `CAST` scored 1.19 on two significant DE genes.

**It is not competitive.** A perturbation that collapses the whole transcriptome drags the module
down with everything else. A self-contained test calls that a hit. It is not a hit; it is
cytotoxicity.

Both are the classic failures of self-contained gene-set testing, and CAMERA (Wu and Smyth,
Nucleic Acids Research 2012) fixes both. It compares the module against the rest of the same
perturbation's transcriptome, and it inflates the null variance by

    VIF = 1 + (m - 1) * rho_bar

where ``rho_bar`` is the mean pairwise inter-gene correlation inside the module. We estimate
``rho_bar`` across perturbations from the z-score matrix itself.

A third correction, specific to Perturb-seq: the perturbed gene is dropped from both the module
and the background. Its own on-target knockdown z-score is large, negative, and required to be
significant by QC. Leaving it in lets a gene inflate its own module score. Twenty-one effector
genes are themselves perturbed in this library, so this is not hypothetical.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def inter_gene_vif(z_module: np.ndarray) -> tuple[float, float]:
    """Variance inflation factor from the mean pairwise correlation inside a module.

    Args:
        z_module: Array of shape ``(n_perturbations, m)`` of z-scores for the module's genes.

    Returns:
        Tuple of (vif, rho_bar). ``rho_bar`` is the mean off-diagonal Pearson correlation.

    Raises:
        ValueError: If the module has fewer than two genes.
    """
    m = z_module.shape[1]
    if m < 2:
        raise ValueError("VIF requires at least two genes in the module")

    corr = np.corrcoef(z_module, rowvar=False)
    off_diagonal = corr[~np.eye(m, dtype=bool)]
    rho_bar = float(np.nanmean(off_diagonal))
    vif = 1.0 + (m - 1) * rho_bar
    # A strongly anti-correlated module could drive VIF below 1, which would make the test
    # anti-conservative. Clamp at 1: never claim more power than independence would give.
    return max(vif, 1.0), rho_bar


@dataclass(frozen=True)
class CameraResult:
    """Per-perturbation competitive gene-set statistics.

    Attributes:
        suppression_t: Positive when the module is suppressed relative to the rest of the
            transcriptome in that perturbation. This is the calibrated efficacy statistic.
        mean_module: Mean z-score across the module's genes, perturbed gene excluded.
        mean_rest: Mean z-score across all other measured genes.
        vif: Variance inflation factor used.
        rho_bar: Mean inter-gene correlation used to derive the VIF.
    """

    suppression_t: np.ndarray
    mean_module: np.ndarray
    mean_rest: np.ndarray
    vif: float
    rho_bar: float


def camera_suppression(
    z: np.ndarray,
    module_indices: np.ndarray,
    self_indices: np.ndarray,
    vif: float | None = None,
    rho_bar: float | None = None,
) -> CameraResult:
    """Competitive, correlation-aware test of module suppression, per perturbation.

    For each perturbation the perturbed gene is removed from both the module and the
    background, so an on-target knockdown cannot inflate its own score.

    Args:
        z: Array of shape ``(n_perturbations, n_genes)`` of z-scores.
        module_indices: Positional gene indices belonging to the module.
        self_indices: For each perturbation, the positional index of its own perturbed gene in
            the gene axis, or -1 when that gene is not measured.
        vif: Precomputed variance inflation factor. Computed from ``z`` if omitted.
        rho_bar: Precomputed mean inter-gene correlation, reported alongside ``vif``.

    Returns:
        A :class:`CameraResult`. ``suppression_t`` is positive when the module is suppressed
        more than the rest of the transcriptome.

    Raises:
        ValueError: If shapes are inconsistent.
    """
    n_pert, n_genes = z.shape
    if self_indices.shape[0] != n_pert:
        raise ValueError("self_indices must have one entry per perturbation")

    module_mask = np.zeros(n_genes, dtype=bool)
    module_mask[module_indices] = True
    m_full = int(module_mask.sum())

    if vif is None:
        vif, rho_bar = inter_gene_vif(z[:, module_indices])

    z64 = z.astype(np.float64, copy=False)

    sum_all = z64.sum(axis=1)
    sumsq_all = np.einsum("ij,ij->i", z64, z64)
    sum_module = z64[:, module_mask].sum(axis=1)

    has_self = self_indices >= 0
    rows = np.arange(n_pert)
    self_z = np.zeros(n_pert, dtype=np.float64)
    self_z[has_self] = z64[rows[has_self], self_indices[has_self]]
    self_in_module = np.zeros(n_pert, dtype=bool)
    self_in_module[has_self] = module_mask[self_indices[has_self]]

    # Drop the perturbed gene from every total, module and background alike.
    n_used = n_genes - has_self.astype(np.int64)
    sum_all_used = sum_all - self_z
    sumsq_all_used = sumsq_all - self_z**2

    m_used = m_full - self_in_module.astype(np.int64)
    sum_module_used = sum_module - np.where(self_in_module, self_z, 0.0)

    n_rest = n_used - m_used
    sum_rest = sum_all_used - sum_module_used

    mean_module = sum_module_used / m_used
    mean_rest = sum_rest / n_rest

    variance = (sumsq_all_used - sum_all_used**2 / n_used) / (n_used - 1)
    sd = np.sqrt(np.maximum(variance, 1e-12))

    standard_error = sd * np.sqrt(vif / m_used + 1.0 / n_rest)
    # Positive when the module sits BELOW the rest of the transcriptome.
    suppression_t = (mean_rest - mean_module) / standard_error

    return CameraResult(
        suppression_t=suppression_t,
        mean_module=mean_module,
        mean_rest=mean_rest,
        vif=float(vif),
        rho_bar=float(rho_bar if rho_bar is not None else np.nan),
    )


def zscore(values: np.ndarray) -> np.ndarray:
    """Standardise a vector, ignoring NaNs.

    Args:
        values: Input array.

    Returns:
        Array with mean 0 and unit standard deviation over the finite entries.
    """
    finite = np.isfinite(values)
    mu = np.nanmean(values[finite])
    sigma = np.nanstd(values[finite])
    if sigma == 0:
        return np.zeros_like(values)
    return (values - mu) / sigma
