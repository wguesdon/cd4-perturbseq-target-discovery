"""Loaders for the prior gene lists and reference screens.

Every loader returns gene symbols. Provenance for all of these files is the source
paper's analysis repository, `emdann/GWT_perturbseq_analysis_2025@master` (MIT),
fetched by ``scripts/fetch_priors.sh``.
"""

from __future__ import annotations

import pandas as pd

from .paths import PRIORS


def _read_symbol_list(filename: str) -> set[str]:
    """Read a headerless one-symbol-per-line gene list.

    Args:
        filename: Basename of the file inside the priors directory.

    Returns:
        Set of gene symbols, whitespace-stripped and empty entries dropped.
    """
    series = pd.read_csv(PRIORS / filename, header=None).iloc[:, 0]
    return {s for s in series.astype(str).str.strip() if s and s.lower() != "nan"}


def core_essential_genes() -> set[str]:
    """Hart core-essential genes as bundled by the authors.

    Note:
        This file contains 283 symbols, whereas the published Hart CEG2 set has
        roughly 684. It is therefore a pre-filtered subset, presumably intersected
        with genes measured or perturbed in this screen. Treat it as a
        high-precision, low-recall essentiality prior, not as complete CEG2.

    Returns:
        Set of core-essential gene symbols.
    """
    return _read_symbol_list("core_essentials_hart.tsv")


def druggable_classes() -> dict[str, set[str]]:
    """Druggable-genome protein classes.

    Returns:
        Mapping from class name to the set of gene symbols in that class.
    """
    files = {
        "kinase": "kinases.tsv",
        "gpcr": "gpcr_union.tsv",
        "ion_channel": "ion_channels.tsv",
        "enzyme": "enzymes.tsv",
        "transporter": "transporters.tsv",
        "nuclear_receptor": "nuclear_receptors.tsv",
        "catalytic_receptor": "catalytic_receptors.tsv",
    }
    return {name: _read_symbol_list(f) for name, f in files.items()}


def iei_genes() -> set[str]:
    """Genes whose loss of function causes an inborn error of immunity (IUIS 2024).

    The `Genetic defect` column mixes gene symbols with cytogenetic lesions (for
    example `11q23del`). We keep only entries that look like HGNC symbols.

    Returns:
        Set of gene symbols implicated in human immunodeficiency.
    """
    table = pd.read_csv(PRIORS / "IUIS-IEI-list-July-2024V2.csv")
    defects = table["Genetic defect"].astype(str).str.strip()
    symbol_like = defects.str.fullmatch(r"[A-Z][A-Z0-9\-]{1,14}")
    return set(defects[symbol_like.fillna(False)])


def immune_effector_genes() -> pd.DataFrame:
    """Curated immune effector genes with their category.

    Returns:
        DataFrame with columns ``gene_name`` and ``category`` (Cytokine, Receptor,
        TF, Others). Category whitespace is stripped, which merges the duplicate
        ``TF`` level present in the source file.
    """
    table = pd.read_csv(PRIORS / "immune_effector_genes.csv")
    table.columns = [c.strip().lower() for c in table.columns]
    table["gene_name"] = table["gene_name"].astype(str).str.strip()
    table["category"] = table["category"].astype(str).str.strip()
    return table


def arce_stim_vs_rest() -> pd.DataFrame:
    """Arce 2024 bulk RNA-seq DESeq2 results, Teff stimulation versus resting.

    An external, non-perturbational definition of the T cell activation program.

    Returns:
        DataFrame with columns ``gene_name``, ``log2FoldChange``, ``padj``.
    """
    table = pd.read_csv(
        PRIORS / "Arce2024_20230130_DESeq2_output_AAVS1_Teff_Stimulation_vs_Resting.csv"
    )
    table["gene_name"] = table["gene_name"].astype(str).str.strip()
    return table[["gene_name", "log2FoldChange", "padj"]].dropna(subset=["padj"])


def schmidt_cd4_il2_screen() -> pd.DataFrame:
    """Schmidt & Steinhart 2022 genome-wide CRISPRi screen, CD4+ IL-2 production.

    Provides an orthogonal, CD4-native readout of which knockdowns suppress the
    activation program, measured by protein-level cytokine output rather than
    transcriptome.

    Returns:
        DataFrame with columns ``gene_name`` and ``il2_lfc``. Negative ``il2_lfc``
        means the knockdown reduces IL-2 production.
    """
    table = pd.read_csv(PRIORS / "SchmidtSteinhart2022_CRISPRi_screen_gene_phenotypes.csv")
    cd4 = table[table["phenotype"] == "CD4+ IL2"].copy()
    cd4["gene_name"] = cd4["id"].astype(str).str.strip()
    # MAGeCK reports the same lfc in the neg| and pos| blocks; take the neg| column.
    cd4 = cd4.rename(columns={"neg|lfc": "il2_lfc", "neg|fdr": "il2_neg_fdr"})
    return cd4[["gene_name", "il2_lfc", "il2_neg_fdr"]]
