"""Inspect the structure of `GWCD4i.DE_stats.h5ad` without loading any layer.

Confirms the shape, dtype, chunking, layer names, and the `.obs` columns actually
present, which differ between the December 2025 supplementary CSV and the May 2026
h5ad (the latter adds the guide and donor concordance columns).

Usage:
    uv run python scripts/00_inspect_de_stats.py
"""

from __future__ import annotations

from cd4_perturbseq import de_stats, paths


def main() -> None:
    """Print the on-disk structure and the obs/var schema."""
    path = paths.DE_STATS_H5AD
    print(f"file: {path}")
    print(f"size: {path.stat().st_size / 1e9:.2f} GB\n")

    info = de_stats.describe(path)
    for key, value in info.items():
        print(f"{key:14s} {value}")

    obs = de_stats.read_obs(path)
    var = de_stats.read_var(path)

    print(f"\nobs: {obs.shape}")
    for col in obs.columns:
        print(f"  {col:32s} {str(obs[col].dtype):>10s}  n_missing={int(obs[col].isna().sum())}")

    print(f"\nvar: {var.shape}  columns={list(var.columns)}")
    print(f"var head: {list(var.index[:5])}")

    print("\nculture_condition counts:")
    print(obs["culture_condition"].value_counts().to_string())

    print(f"\nunique perturbed genes: {obs['target_contrast_gene_name'].nunique()}")


if __name__ == "__main__":
    main()
