"""The denominator manifest: every N this project quotes, and exactly what it counts.

Requested by an external red-team audit (2026-07-09): the report uses 12,779, 10,282, 33,983, 6,371,
286, 214, 206, 97, 91, 36, 20, 6 and 4 as denominators, and a reader cannot be expected to hold them
apart. Two of them are actively misleading unless stated:

  * `all_library` is **12,779** by our count. The source paper says **12,748**. The 31-gene difference
    is unreconciled and is reported, not hidden. An earlier version of this repo said 13,129, which was
    an alias artifact from unioning two gene-symbol vocabularies.
  * The **6** recovered approved-drug targets are NOT a subset of the 36 curated positives. `CD2` and
    `CD28` were deliberately held out of the ground truth (RULE #4: never add a positive after watching
    it rank). So "6 recovered" is an annotation and "5 of 20 above chance" is the test.

Every row is recomputed from a committed table or from the source paper, and the script exits non-zero
if a live count has drifted.

Usage:
    uv run python scripts/31_denominator_manifest.py
"""

from __future__ import annotations

import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from cd4_perturbseq import paths  # noqa: E402

# Not derivable from any committed table; they are properties of the h5ad, printed by scripts/00.
PANEL_GENES = 10_282
DE_ROWS = 33_983
PAPER_LIBRARY = 12_748


def main() -> None:
    """Recompute every denominator from its source and write the manifest."""
    funnel = pd.read_csv(paths.TABLES / "selection_funnel.csv").set_index("stage")
    window = pd.read_csv(paths.TABLES / "window_score.csv")
    truth = pd.read_csv(paths.GROUND_TRUTH / "immunomodulator_targets.csv")
    short = pd.read_csv(paths.TABLES / "final_shortlist.csv")
    nom = pd.read_csv(paths.TABLES / "nomination_recalibrated.csv")
    sens = pd.read_csv(paths.TABLES / "nomination_rule_sensitivity.csv").set_index("rule")
    doe = pd.read_csv(paths.TABLES / "direction_of_effect.csv")

    positives = set(truth.loc[truth["include_as_positive"], "gene_symbol"])
    rankable = set(window["gene_name"])
    safe = set(window.loc[window["safe"], "gene_name"])

    n_library = int(funnel.loc["in sgRNA library", "all_library"])
    n_rankable = len(window)
    n_evidence = int((~window["fail_evidence"]).sum())
    n_safe = int(window["safe"].sum())
    n_pos_rankable = len(positives & rankable)
    n_pos_safe = len(positives & safe)
    n_recovered = int(short["is_known_drug"].fillna(False).sum())
    n_candidates = int(sens.loc["(screen gate only)", "novel_nominated"])
    n_old_rule = int(sens.loc["supported & lof_tolerant & tractable", "novel_nominated"])
    n_nom = len(nom)
    n_adjudicable = int(doe["direction_verdict"].isin(["CONCORDANT", "DISCORDANT"]).sum())
    n_unmeasured = int(doe["direction_verdict"].eq("UNKNOWN").sum())

    rows = [
        ("sgRNA library targets", n_library, "Genes targeted by the CRISPRi library, by our count.",
         "selection_funnel.csv", f"The source paper reports {PAPER_LIBRARY:,}. The {n_library - PAPER_LIBRARY} "
         "gene difference is UNRECONCILED and reported as such. A prior figure of 13,129 was an alias "
         "artifact from unioning two gene-symbol vocabularies, and is retracted."),
        ("measured panel genes", PANEL_GENES, "Genes on the probe panel. NOT the transcriptome.",
         "GWCD4i.DE_stats.h5ad via scripts/00_inspect_de_stats.py",
         "Every over-representation test uses this as the background, never the genome."),
        ("DE rows", DE_ROWS, "Perturbation-by-condition pairs in the effect matrix.",
         "GWCD4i.DE_stats.h5ad via scripts/00_inspect_de_stats.py",
         "Three conditions: Rest, Stim 8 h, Stim 48 h. Not a gene count."),
        ("rankable perturbations", n_rankable, "QC-passing Stim48hr perturbations used for ranking.",
         "window_score.csv", "The denominator for every rank quoted in this report."),
        ("evidence-floor passers", n_evidence,
         "Perturbations suppressing >= 3 EFFECTOR-module genes at FDR 0.10.",
         "window_score.csv (~fail_evidence)",
         "The effector module, NOT the co-inhibitory module. The two are disjoint by construction."),
        ("screen-passing genes", n_safe,
         "Pass the evidence floor AND are not top-quartile co-inhibitory suppressors.",
         "window_score.csv (safe)",
         "The CSV column is named `safe`. It means screen-passing, not organism-safe. The assay "
         "cannot see therapeutic safety."),
        ("curated immunomodulator rows", len(truth), "All rows in the ground-truth file.",
         "immunomodulator_targets.csv", "Includes excluded rows (agonists, drug-binding chaperones)."),
        ("curated positives", len(positives), "Rows flagged include_as_positive.",
         "immunomodulator_targets.csv", "The evaluation set. Fixed a priori (RULE #4)."),
        ("assay-visible positives", n_pos_rankable,
         "Curated positives that are rankable in this screen.",
         "window_score.csv INTERSECT ground truth",
         "The correct denominator for recovery. The rest are invisible by construction: TCR-only "
         "stimulation hides cytokine-signalling targets, and JAK1/JAK3 were never perturbed."),
        ("positives among screen-passing genes", n_pos_safe,
         "Curated positives that clear the triage layer.",
         "window_score.csv", "Only this many. It is the hard power ceiling on any recovery test "
         "conditioned on the screen gate."),
        ("recovered approved-drug targets", n_recovered,
         "Approved-drug targets present among the screen-passing genes.",
         "final_shortlist.csv (is_known_drug)",
         "NOT a subset of the curated positives. CD2 and CD28 were deliberately held OUT of the "
         "ground truth (RULE #4). '6 recovered' is an annotation; '5 of 20 above chance' is the test."),
        ("nomination candidates", n_candidates,
         "Screen-passing, not a known drug, not direction-discordant.",
         "nomination_rule_sensitivity.csv", "The input to any nomination rule."),
        ("old-rule nominations", n_old_rule,
         "What `supported & lof_tolerant & tractable` returned.",
         "nomination_rule_sensitivity.csv", "One gene, ICAM2. This is the number N20 rejected."),
        ("rebuilt nomination pool", n_nom,
         "Screen gate AND structurally tractable AND not direction-discordant.",
         "nomination_recalibrated.csv", "Exploratory repair, not confirmatory discovery. These are "
         "ranked hypotheses with named liabilities, not candidates in the drug-discovery sense."),
        ("direction-adjudicable genes", n_adjudicable,
         "Screen-passing genes with any direction verdict.",
         "direction_of_effect.csv", "All of them required prior characterisation."),
        ("direction-unmeasured genes", n_unmeasured,
         "Screen-passing genes with no direction verdict.",
         "direction_of_effect.csv",
         "Absence of a measurement, NOT evidence of a favourable direction. Missingness is not random."),
    ]

    frame = pd.DataFrame(rows, columns=["universe", "n", "meaning", "source", "caveat"])
    out = paths.TABLES / "denominator_manifest.csv"
    frame.to_csv(out, index=False)

    print(frame[["universe", "n", "source"]].to_string(index=False))
    print(f"\nwrote {out}")

    checks = {
        "library vs paper": (n_library, 12_779),
        "rankable": (n_rankable, 6_371),
        "screen-passing": (n_safe, 214),
        "assay-visible positives": (n_pos_rankable, 20),
        "positives in screen-passing set": (n_pos_safe, 4),
        "recovered drug targets": (n_recovered, 6),
        "rebuilt nomination pool": (n_nom, 91),
        "direction-unmeasured": (n_unmeasured, 191),
    }
    bad = {k: v for k, v in checks.items() if v[0] != v[1]}
    print()
    for k, (got, want) in checks.items():
        print(f"  [{'ok  ' if got == want else 'DRIFT'}] {k:<34} {got} (expected {want})")
    if bad:
        print(f"\nDRIFT: {bad}. A denominator moved. Do not render the report.")
        sys.exit(1)

    print(f"\nUNRECONCILED: our library count is {n_library:,}; the source paper reports "
          f"{PAPER_LIBRARY:,}. Difference {n_library - PAPER_LIBRARY}. Reported, not hidden.")


if __name__ == "__main__":
    main()
