"""Gene modules that define the therapeutic-window objective.

The window we are after is: suppress the pro-inflammatory effector program in
stimulated cells, while sparing the tolerance program, the resting transcriptome,
and viability. That means the effector module is the objective and the tolerance
module is a penalty, even though both are induced by stimulation.
"""

from __future__ import annotations

import pandas as pd

from .paths import INTERIM

TOLERANCE_GENES: frozenset[str] = frozenset(
    {
        # Regulatory T cell identity and suppressive effectors.
        "FOXP3",
        "IKZF2",
        "IL10",
        "TGFB1",
        "LRRC32",
        "ENTPD1",
        "NT5E",
        # Co-inhibitory checkpoints. Knocking these down releases the brakes.
        "CTLA4",
        "PDCD1",
        "LAG3",
        "HAVCR2",
        "TIGIT",
    }
)
"""Stimulation-induced genes whose suppression is a liability, not a benefit."""


def load_activation_program() -> pd.DataFrame:
    """Load the activation program written by ``scripts/01_build_activation_program.py``.

    Returns:
        DataFrame with columns ``gene_name``, ``curated``, ``external``,
        ``arce_log2fc``, ``source``.

    Raises:
        FileNotFoundError: If the program has not been built yet.
    """
    path = INTERIM / "activation_program.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found; run scripts/01_build_activation_program.py first"
        )
    return pd.read_csv(path)


def effector_core(program: pd.DataFrame | None = None) -> list[str]:
    """High-confidence pro-inflammatory effector genes.

    Genes supported by both the curated effector list and the external Arce
    stimulation contrast, with the tolerance module removed.

    Args:
        program: Optional pre-loaded program frame; loaded from disk if omitted.

    Returns:
        Sorted gene symbols.
    """
    frame = load_activation_program() if program is None else program
    core = frame.loc[frame["source"] == "both", "gene_name"]
    return sorted(set(core) - TOLERANCE_GENES)


def effector_broad(program: pd.DataFrame | None = None) -> list[str]:
    """The broad stimulation-induced program, tolerance genes removed.

    Args:
        program: Optional pre-loaded program frame; loaded from disk if omitted.

    Returns:
        Sorted gene symbols.
    """
    frame = load_activation_program() if program is None else program
    return sorted(set(frame["gene_name"]) - TOLERANCE_GENES)


def tolerance_module(program: pd.DataFrame | None = None) -> list[str]:
    """Tolerance genes that are present in the activation program.

    Args:
        program: Optional pre-loaded program frame; loaded from disk if omitted.

    Returns:
        Sorted gene symbols.
    """
    frame = load_activation_program() if program is None else program
    return sorted(set(frame["gene_name"]) & TOLERANCE_GENES)
