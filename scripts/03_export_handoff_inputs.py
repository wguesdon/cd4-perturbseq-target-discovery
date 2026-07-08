"""Export the small, machine-portable inputs a Claude Science handoff needs.

Claude Science runs in a browser UI, possibly on a different machine, and will not have the
16.8 GB h5ad. Every handoff must therefore operate on small artifacts. This script writes
them from the already-computed ranking, so a handoff never needs the effect matrix.

Outputs, all small and all COMMITTED under resources/handoff_inputs/ so that a plain
`git pull` on any machine is enough to run a handoff:

    rankable_genes.txt         the 6,371 genes that survived perturbation QC
    measured_genes.txt         the 10,282 measured transcriptome genes
    perturbed_genes.txt        the 11,526 genes perturbed anywhere in the library
    ground_truth_coverage.csv  per-positive coverage: perturbed, measured, rankable, naive rank

Usage:
    uv run python scripts/03_export_handoff_inputs.py
"""

from __future__ import annotations

import shutil

import pandas as pd

from cd4_perturbseq import paths

RANKING = paths.TABLES / "risk_kill_naive_reversal.csv"
GROUND_TRUTH = paths.GROUND_TRUTH / "immunomodulator_targets.csv"


def main() -> None:
    """Write the rankable gene list and the per-positive coverage table."""
    paths.ensure_dirs()

    if not RANKING.exists():
        raise FileNotFoundError(f"{RANKING} missing; run scripts/02_risk_kill_reversal.py first")

    ranked = pd.read_csv(RANKING)
    rankable = sorted(ranked["target_contrast_gene_name"].astype(str).unique())

    out_rankable = paths.HANDOFF_INPUTS / "rankable_genes.txt"
    out_rankable.write_text("\n".join(rankable) + "\n")
    print(f"wrote {out_rankable}  ({len(rankable)} genes surviving QC)")

    # Promote the two gene lists out of gitignored data/interim into the committed bundle.
    for name in ("measured_genes.txt", "perturbed_genes.txt"):
        src = paths.INTERIM / name
        if not src.exists():
            raise FileNotFoundError(f"{src} missing; run scripts/00_inspect_de_stats.py first")
        shutil.copy2(src, paths.HANDOFF_INPUTS / name)
        print(f"wrote {paths.HANDOFF_INPUTS / name}")

    measured = set((paths.INTERIM / "measured_genes.txt").read_text().split())
    perturbed = set((paths.INTERIM / "perturbed_genes.txt").read_text().split())
    rankable_set = set(rankable)

    truth = pd.read_csv(GROUND_TRUTH)
    positives = truth[truth["include_as_positive"]].copy()
    positives["perturbed"] = positives["gene_symbol"].isin(perturbed)
    positives["measured"] = positives["gene_symbol"].isin(measured)
    positives["rankable"] = positives["gene_symbol"].isin(rankable_set)

    rank_of = dict(zip(ranked["target_contrast_gene_name"].astype(str), ranked["rank"], strict=True))
    positives["naive_rank"] = positives["gene_symbol"].map(rank_of)

    out_cov = paths.HANDOFF_INPUTS / "ground_truth_coverage.csv"
    cols = ["gene_symbol", "drug_examples", "mechanism", "perturbed", "measured", "rankable", "naive_rank"]
    positives[cols].to_csv(out_cov, index=False)
    print(f"wrote {out_cov}")

    n = len(positives)
    print(f"\ncoverage of {n} curated positives:")
    print(f"  perturbed in library : {int(positives['perturbed'].sum())}")
    print(f"  also measured        : {int(positives['measured'].sum())}")
    print(f"  also rankable        : {int(positives['rankable'].sum())}")
    print(f"\n  never perturbed: {', '.join(sorted(positives.loc[~positives['perturbed'], 'gene_symbol']))}")
    lost = positives[positives["perturbed"] & ~positives["rankable"]]["gene_symbol"]
    print(f"  perturbed but failed QC: {', '.join(sorted(lost))}")


if __name__ == "__main__":
    main()
