"""Loaders for the prior gene lists and reference screens.

Every loader returns gene symbols. Provenance for all of these files is the source
paper's analysis repository, `emdann/GWT_perturbseq_analysis_2025@master` (MIT),
fetched by ``scripts/fetch_priors.sh``.
"""

from __future__ import annotations

import pandas as pd

import glob

from .paths import EXTERNAL, PRIORS

SAFETY = EXTERNAL / "safety_priors"
"""Organism-level safety priors, fetched by ``scripts/06_fetch_safety_priors.sh``."""

OPEN_TARGETS = EXTERNAL / "open_targets"
"""Open Targets 26.06 bulk parquet, fetched by ``scripts/15`` and ``scripts/19`` (gitignored)."""

OT_GENETIC_THRESHOLD = 0.1
"""The source paper's autoimmune genetic-evidence cut (Zhu, Dann et al., Methods p.50). Adopted,
not tuned. N10 reports a sensitivity at 0.05 and 0.2."""

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


def iei_genes(lof_only: bool = True) -> set[str]:
    """Genes whose LOSS of function causes an inborn error of immunity (IUIS 2024).

    Three corrections over the naive read of this table.

    **Mechanism.** The table has a ``GOF/DN`` column. Sixteen entries are gain-of-function
    diseases, where loss of function is not the pathogenic mechanism and therefore says nothing
    about what CRISPRi does: ``CXCR4``, ``HCK``, ``LYN``, ``NFKBIA``, ``NLRC4``, ``NLRP12``,
    ``NLRP3``, ``PLCG1``, ``PLCG2``, ``SAMD9``, ``SAMD9L``, ``STAT4``, ``STAT6``, ``STING1``,
    ``SYK``, ``TLR8``. Entries annotated ``LOF GOF`` or ``LOF, DN`` do involve loss of function
    and are kept.

    **Composite entries.** ``C4A+C4B`` is two genes and was silently dropped by the symbol regex.

    **Cytogenetic lesions.** The ``Genetic defect`` column also holds things like ``11q23del``
    and ``14q32 deletion or mutation``, which are not genes. Those are dropped.

    Note:
        Membership here is NOT a safety gate. Verified 2026-07-08: the flag is more enriched among
        approved immunomodulator targets (25.0%, OR 8.31) than among the perturbations a naive
        reversal ranking calls toxic (14.0%, OR 4.16). A gene whose loss causes immunodeficiency is
        a gene that matters for immunity, which is what a good immune drug target is. Use this as a
        reported liability annotation. See ``docs/results/adversarial_audit_2026_07_08.md``.

    Args:
        lof_only: Drop entries whose stated mechanism is purely gain of function.

    Returns:
        Set of gene symbols implicated in human immunodeficiency.
    """
    table = pd.read_csv(PRIORS / "IUIS-IEI-list-July-2024V2.csv")
    table.columns = [c.strip() for c in table.columns]

    defect = table["Genetic defect"].astype(str).str.strip()
    mechanism = table["GOF/DN"].astype(str).str.strip().str.upper()

    # Some rows carry the mechanism inside the defect string, e.g. "PIK3CD GOF", so the GOF/DN
    # column alone does not catch them.
    inline_gof = defect.str.contains(r"\bGOF\b", regex=True)
    if lof_only:
        keep = (mechanism != "GOF") & ~inline_gof
        defect = defect[keep]

    # Drop parenthetical aliases: "CD40 (TNFRSF5)" -> "CD40", "KMT2D (MLL2)" -> "KMT2D".
    cleaned = defect.str.replace(r"\s*\([^)]*\)", "", regex=True)
    # Composite and multi-gene entries: "C4A+C4B", "CFHR1 CFHR2. CFHR3 CFHR4 CFHR5".
    tokens = cleaned.str.split(r"[+,.\s]+", regex=True).explode().str.strip()

    # A symbol is all-caps alphanumeric, or a CNorfNNN open reading frame. Cytogenetic lesions
    # ("11q23del", "Del10p13-p14", "Large deletion of 22q11.2") and "Unknown" fail both.
    is_symbol = tokens.str.fullmatch(r"[A-Z][A-Z0-9\-]{1,13}") | tokens.str.fullmatch(r"C\d+orf\d+")
    mechanism_words = {"GOF", "LOF", "DN", "AR", "AD", "XL"}
    symbols = {s for s in tokens[is_symbol.fillna(False)] if s not in mechanism_words}
    return symbols


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

    ``prec`` is gnomAD's posterior probability that the gene is recessive: heterozygous LoF is
    tolerated, homozygous LoF is not. It is the statistic LOEUF and pLI are blind to, and it is the
    one that matters here, because a potent inhibitor phenocopies BIALLELIC loss, not the
    heterozygous loss LOEUF measures. Several Tier A mitochondrial genes (POLG, VARS2, ACAD9) have
    LOEUF near 1 (so LOEUF calls them safe to inhibit) yet ``prec`` above 0.99.

    Returns:
        DataFrame with columns ``gene_name``, ``loeuf``, ``pli``, ``prec``.
    """
    table = pd.read_csv(
        SAFETY / "gnomad.v4.0.constraint_metrics.tsv",
        sep="\t",
        usecols=["gene", "mane_select", "lof.oe_ci.upper", "lof.pLI", "lof.pRec"],
        low_memory=False,
    )
    table = table[table["mane_select"].astype(str).str.lower() == "true"]
    table = table.rename(
        columns={"gene": "gene_name", "lof.oe_ci.upper": "loeuf", "lof.pLI": "pli", "lof.pRec": "prec"}
    )
    table = table.dropna(subset=["gene_name"]).drop_duplicates("gene_name")
    return table[["gene_name", "loeuf", "pli", "prec"]]


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


def _ot_ensembl_to_symbol() -> pd.Series:
    """Map Open Targets Ensembl gene ids to approved symbols, from the on-disk target parquet.

    Returns:
        Series indexed by Ensembl id, values are approved symbols.
    """
    parts = [pd.read_parquet(f, columns=["id", "approvedSymbol"])
             for f in sorted(glob.glob(str(OPEN_TARGETS / "target_*.parquet")))]
    table = pd.concat(parts, ignore_index=True).dropna(subset=["approvedSymbol"]).drop_duplicates("id")
    return table.set_index("id")["approvedSymbol"]


def ot_genetic_support(threshold: float = OT_GENETIC_THRESHOLD) -> pd.DataFrame:
    """Open Targets autoimmune genetic-association support, per gene symbol.

    Reads the filtered ``autoimmune_genetic_assoc.parquet`` produced by
    ``scripts/19_fetch_ot_genetic_evidence.sh`` (genetic_association datatype scores for the 14
    autoimmune diseases the source paper used), maps Ensembl ids to symbols, and summarises per gene.

    Args:
        threshold: A disease counts as supported if its association score is at least this. The
            source paper used 0.1.

    Returns:
        DataFrame with ``gene_name``, ``ot_genetic_max`` (max score across the 14 diseases),
        ``ot_genetic_n_diseases`` (how many clear ``threshold``), and ``ot_genetic_supported``
        (bool, at least one disease clears it). One row per gene present in the association table.
    """
    assoc = pd.read_parquet(OPEN_TARGETS / "autoimmune_genetic_assoc.parquet")
    assoc["gene_name"] = assoc["targetId"].map(_ot_ensembl_to_symbol())
    assoc = assoc.dropna(subset=["gene_name"])
    grouped = assoc.groupby("gene_name")
    out = pd.DataFrame({
        "ot_genetic_max": grouped["associationScore"].max(),
        "ot_genetic_n_diseases": grouped["associationScore"].apply(lambda s: int((s >= threshold).sum())),
    }).reset_index()
    out["ot_genetic_supported"] = out["ot_genetic_n_diseases"] > 0
    return out


def ot_tractability() -> pd.DataFrame:
    """Open Targets small-molecule and antibody tractability, summarised per gene symbol.

    The raw field is a list of ``{modality, id, value}`` buckets. We reduce it to a few booleans a
    triage layer can read: a small-molecule pocket exists, the gene is in a druggable protein family,
    it is antibody-accessible at the cell surface, and there is any clinical precedent for a drug
    against it in any modality.

    ``tractable`` is the disjunction of ALL FOUR flags, ``clinical_precedent`` included. Until
    2026-07-09 it was ``sm_pocket or sm_druggable_family or ab_accessible``, silently discarding
    ``clinical_precedent`` — the strongest tractability evidence that exists, since it means a drug
    against the gene has already reached the clinic. The omission inverted the flag on the genes that
    matter most: ``MALT1``, which has clinical-stage allosteric protease inhibitors, read as
    untractable, while ``ICAM2`` read as tractable on cell-surface localisation alone, with neither a
    pocket nor a precedent. See ``docs/preregistration_n20_2026_07_09.md``.

    Returns:
        DataFrame with ``gene_name``, ``sm_pocket``, ``sm_druggable_family``, ``ab_accessible``,
        ``clinical_precedent``, and ``tractable`` (any of the four).
    """
    parts = [pd.read_parquet(f, columns=["approvedSymbol", "tractability"])
             for f in sorted(glob.glob(str(OPEN_TARGETS / "target_*.parquet")))]
    table = pd.concat(parts, ignore_index=True).dropna(subset=["approvedSymbol"]).drop_duplicates("approvedSymbol")

    rows = []
    for symbol, buckets in zip(table["approvedSymbol"], table["tractability"], strict=True):
        flags = {(b["modality"], b["id"]): bool(b["value"]) for b in (buckets if buckets is not None else [])}
        sm_pocket = flags.get(("SM", "High-Quality Pocket"), False) or flags.get(("SM", "Med-Quality Pocket"), False)
        sm_family = flags.get(("SM", "Druggable Family"), False)
        ab = any(flags.get(("AB", k), False) for k in
                 ("UniProt loc high conf", "GO CC high conf", "UniProt loc med conf", "GO CC med conf",
                  "Human Protein Atlas loc", "UniProt SigP or TMHMM"))
        clinical = any(flags.get((m, k), False) for m in ("SM", "AB", "PR", "OC")
                       for k in ("Approved Drug", "Advanced Clinical", "Phase 1 Clinical"))
        rows.append({
            "gene_name": symbol, "sm_pocket": sm_pocket, "sm_druggable_family": sm_family,
            "ab_accessible": ab, "clinical_precedent": clinical,
            "tractable": sm_pocket or sm_family or ab or clinical,
        })
    return pd.DataFrame(rows)
