"""Canonical filesystem paths for the project.

All paths are resolved relative to the repository root so that scripts can be run
from any working directory.
"""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

DATA = REPO_ROOT / "data"
RAW = DATA / "raw"
INTERIM = DATA / "interim"
PROCESSED = DATA / "processed"
EXTERNAL = DATA / "external"

PRIORS = EXTERNAL / "gwt_priors"
"""Gene lists and reference screens bundled by the source paper's analysis repo."""

DE_STATS_H5AD = RAW / "GWCD4i.DE_stats.h5ad"
"""Per-perturbation x per-condition differential expression effect matrix."""

DE_STATS_OBS_CSV = RAW / "suppl" / "DE_stats.suppl_table.csv"
"""Tabular form of the h5ad `.obs` (December 2025 snapshot; lacks confidence columns)."""

RESULTS = REPO_ROOT / "results"
FIGURES = RESULTS / "figures"
TABLES = RESULTS / "tables"

RESOURCES = REPO_ROOT / "resources"
GROUND_TRUTH = RESOURCES / "ground_truth"

HANDOFF_INPUTS = RESOURCES / "handoff_inputs"
"""Small, COMMITTED artifacts for Claude Science handoffs.

Claude Science runs in a browser UI, potentially on a different machine, and will never have
the 16.8 GB h5ad. Handoff inputs therefore live here rather than under the gitignored
`data/interim`, so that a `git pull` is enough to run any handoff. They double as the
reproducibility bundle for a reviewer who does not want to download the effect matrix.
"""


def ensure_dirs() -> None:
    """Create the derived-output directories if they do not already exist."""
    for path in (INTERIM, PROCESSED, FIGURES, TABLES, HANDOFF_INPUTS):
        path.mkdir(parents=True, exist_ok=True)
