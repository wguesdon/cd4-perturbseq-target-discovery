"""Build the CD4+ T cell activation program used as the therapeutic-window objective.

The program is the set of genes induced when a CD4+ T cell is stimulated. A candidate
target is interesting if knocking it down suppresses this program in stimulated cells.

Two independent definitions are combined so the score does not rest on one source:

1. `curated`  -- the authors' 117 immune effector genes (cytokines, receptors, TFs).
2. `external` -- genes significantly up in Arce 2024 Teff stimulation versus rest,
   a bulk RNA-seq experiment outside this perturb-seq screen.

Writing both, plus their intersection, lets the downstream score be reported against
a narrow high-confidence core and a broad program.

Usage:
    uv run python scripts/01_build_activation_program.py
"""

from __future__ import annotations

import argparse

import pandas as pd

from cd4_perturbseq import paths, priors

ARCE_PADJ = 0.05
ARCE_LFC = 1.0


def build_program(arce_padj: float, arce_lfc: float) -> pd.DataFrame:
    """Assemble the activation program from the curated and external sources.

    Args:
        arce_padj: Adjusted p-value ceiling for the Arce stim-versus-rest contrast.
        arce_lfc: Log2 fold change floor for calling a gene stimulation-induced.

    Returns:
        DataFrame with columns ``gene_name``, ``curated``, ``external``,
        ``arce_log2fc``, and ``source``.
    """
    effectors = priors.immune_effector_genes()
    curated = set(effectors["gene_name"])

    arce = priors.arce_stim_vs_rest()
    up = arce[(arce["padj"] < arce_padj) & (arce["log2FoldChange"] > arce_lfc)]
    # A symbol can appear more than once; keep its strongest induction.
    up = up.sort_values("log2FoldChange", ascending=False).drop_duplicates("gene_name")
    external = set(up["gene_name"])

    all_genes = sorted(curated | external)
    lfc_by_gene = dict(zip(up["gene_name"], up["log2FoldChange"], strict=True))

    frame = pd.DataFrame({"gene_name": all_genes})
    frame["curated"] = frame["gene_name"].isin(curated)
    frame["external"] = frame["gene_name"].isin(external)
    frame["arce_log2fc"] = frame["gene_name"].map(lfc_by_gene)
    frame["source"] = frame.apply(
        lambda r: "both" if r["curated"] and r["external"] else ("curated" if r["curated"] else "external"),
        axis=1,
    )
    return frame.sort_values(["source", "gene_name"]).reset_index(drop=True)


def main() -> None:
    """Build the program and write it to ``data/interim/activation_program.csv``."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arce-padj", type=float, default=ARCE_PADJ)
    parser.add_argument("--arce-lfc", type=float, default=ARCE_LFC)
    args = parser.parse_args()

    paths.ensure_dirs()
    program = build_program(args.arce_padj, args.arce_lfc)

    out = paths.INTERIM / "activation_program.csv"
    program.to_csv(out, index=False)

    counts = program["source"].value_counts()
    print(f"activation program: {len(program)} genes -> {out}")
    for source in ("both", "curated", "external"):
        print(f"  {source:9s} {counts.get(source, 0):5d}")
    core = program[program["source"] == "both"]
    print(f"\nhigh-confidence core (curated AND external), {len(core)} genes:")
    print("  " + ", ".join(core["gene_name"].head(30)))


if __name__ == "__main__":
    main()
