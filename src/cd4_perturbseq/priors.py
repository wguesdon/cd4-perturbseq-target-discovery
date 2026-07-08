"""Loaders for the prior gene lists and reference screens.

Every loader returns gene symbols. Provenance for all of these files is the source
paper's analysis repository, `emdann/GWT_perturbseq_analysis_2025@master` (MIT),
fetched by ``scripts/fetch_priors.sh``.
"""

from __future__ import annotations

import pandas as pd

from .paths import EXTERNAL, PRIORS

SAFETY = EXTERNAL / "safety_priors"
"""Organism-level safety priors, fetched by ``scripts/06_fetch_safety_priors.sh``."""

IMMUNE_TISSUES = frozenset(
    {"appendix", "bone marrow", "lymph node", "spleen", "thymus", "tonsil"}
)
"""HPA consensus tissues that are lymphoid. Expression here is on-target, not a liability."""


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


def gnomad_constraint() -> pd.DataFrame:
    """Human loss-of-function intolerance from gnomAD v4.0.

    LOEUF is the upper bound of the 90% confidence interval on the observed-over-expected ratio
    of loss-of-function variants. A low value means living humans carry far fewer LoF variants in
    this gene than chance would give, so selection is removing them. It is the clearest
    organism-level warning that exists, and it is orthogonal to anything measurable in a CD4 T
    cell screen.

    The file carries one row per transcript. We keep MANE Select transcripts and take the first
    row per gene, since duplicates within a gene carry identical constraint values.

    Returns:
        DataFrame with columns ``gene_name``, ``loeuf``, ``pli``.
    """
    table = pd.read_csv(
        SAFETY / "gnomad.v4.0.constraint_metrics.tsv",
        sep="\t",
        usecols=["gene", "mane_select", "lof.oe_ci.upper", "lof.pLI"],
        low_memory=False,
    )
    table = table[table["mane_select"].astype(str).str.lower() == "true"]
    table = table.rename(
        columns={"gene": "gene_name", "lof.oe_ci.upper": "loeuf", "lof.pLI": "pli"}
    )
    table = table.dropna(subset=["gene_name"]).drop_duplicates("gene_name")
    return table[["gene_name", "loeuf", "pli"]]


def hpa_tissue_breadth() -> pd.DataFrame:
    """Expression breadth across the 51 Human Protein Atlas consensus tissues.

    A gene expressed in every tissue has no therapeutic window outside its target tissue: an
    inhibitor cannot be told to spare the gut. Expression in lymphoid tissue is on-target and is
    therefore excluded from the liability measure.

    Returns:
        DataFrame with ``gene_name``, ``n_tissues_expressed`` (nTPM >= 1, all 51 tissues),
        ``n_nonimmune_tissues`` (of 45), and ``max_nonimmune_ntpm``.
    """
    table = pd.read_csv(
        SAFETY / "rna_tissue_consensus.tsv",
        sep="\t",
        usecols=["Gene name", "Tissue", "nTPM"],
    )
    table = table.rename(columns={"Gene name": "gene_name", "Tissue": "tissue", "nTPM": "ntpm"})
    expressed = table[table["ntpm"] >= 1.0]

    total = expressed.groupby("gene_name")["tissue"].nunique().rename("n_tissues_expressed")

    nonimmune = table[~table["tissue"].isin(IMMUNE_TISSUES)]
    nonimmune_expressed = nonimmune[nonimmune["ntpm"] >= 1.0]
    n_nonimmune = nonimmune_expressed.groupby("gene_name")["tissue"].nunique().rename("n_nonimmune_tissues")
    max_nonimmune = nonimmune.groupby("gene_name")["ntpm"].max().rename("max_nonimmune_ntpm")

    return pd.concat([total, n_nonimmune, max_nonimmune], axis=1).reset_index()


def hart_core_essentials_full() -> set[str]:
    """The complete Hart CEGv2 core-essential set, 684 genes.

    The copy bundled by the source paper's analysis repo is a 283-symbol subset with only 68
    genes perturbed in this library, so it could never adjudicate essentiality. This is the
    original.

    Returns:
        Set of gene symbols.
    """
    table = pd.read_csv(SAFETY / "hart_CEGv2.txt", sep="\t")
    return set(table["GENE"].astype(str).str.strip())


def hart_nonessentials() -> set[str]:
    """Hart NEGv1, 927 genes believed nonessential in every cell line tested.

    The negative control our data-native viability axis has never had. Core essentials should be
    depleted at rest; these should not be. If both deplete, the axis is measuring something other
    than essentiality.

    Returns:
        Set of gene symbols.
    """
    table = pd.read_csv(SAFETY / "hart_NEGv1.txt", sep="\t")
    return set(table["GENE"].astype(str).str.strip())


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
