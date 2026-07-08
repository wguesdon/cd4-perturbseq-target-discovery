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
    cs2_blind_ranking.csv      the ranked table with OUR ANALYTICAL CHOICES STRIPPED OUT

The last one needs explaining. ``risk_kill_naive_reversal.csv`` carries ``z_l2_decile``,
``n_cells_decile`` and ``matched_background``. Those are not measurements; they are the record of
which covariate we decided to stratify on and which background we decided to draw. Handing that
table to an independent reviewer and asking it to choose a confounder is not a blind test, because
the answer is written in the column names.

``cs2_blind_ranking.csv`` keeps every measurement, including both candidate covariates
(``n_cells_target`` and ``z_l2``), and drops every column that encodes a conclusion. The reviewer
must pick its own confounder and justify the pick. See ``claude_life_science/RUN_CS1_AND_CS2.md``.

Usage:
    uv run python scripts/03_export_handoff_inputs.py
"""

from __future__ import annotations

import shutil

import pandas as pd

from cd4_perturbseq import paths

RANKING = paths.TABLES / "risk_kill_naive_reversal.csv"
GROUND_TRUTH = paths.GROUND_TRUTH / "immunomodulator_targets.csv"

BLIND_DROP = ("z_l2_decile", "n_cells_decile", "matched_background", "in_top")
"""Columns that record an analytical decision rather than a measurement. Stripped for CS2."""


def write_blind_ranking(ranked: pd.DataFrame) -> None:
    """Write the ranked table with our analytical decisions removed.

    A reviewer asked to find the confounder must not be handed a column called ``z_l2_decile``.
    Both candidate covariates survive, because without ``z_l2`` the reviewer cannot compute effect
    magnitude at all: it needs full rows of a 16.8 GB layer. Neither is labelled as the answer.

    Args:
        ranked: The ranked frame from ``results/tables/risk_kill_naive_reversal.csv``.
    """
    blind = ranked.drop(columns=[c for c in BLIND_DROP if c in ranked.columns])
    out = paths.HANDOFF_INPUTS / "cs2_blind_ranking.csv"
    blind.to_csv(out, index=False)

    leaked = [c for c in BLIND_DROP if c in blind.columns]
    if leaked:
        raise AssertionError(f"blind table still encodes our choices: {leaked}")
    print(f"wrote {out}  ({len(blind):,} rows, {len(blind.columns)} columns, "
          f"{len(BLIND_DROP)} decision columns stripped)")
    print(f"  kept both candidate covariates: n_cells_target, z_l2")


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

    write_blind_ranking(ranked)

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
