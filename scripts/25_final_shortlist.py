"""N16. Synthesise the honest final novel-candidate shortlist across every committed axis.

Runs the RULE #9 loop as a synthesis, not a new measurement. The question is the competition's
literal one: given everything the pipeline has computed (the within-screen safety gate, direction of
effect, autoimmune genetics, LoF tolerance, tractability, held-out Schmidt, and the recovery audit),
**what new drug targets did we actually find, and how confident are we?**

Step 1 (literature). The field's credibility bar for a target nomination is not a single screen hit. It
is the convergence the target-discovery literature rewards: human genetic support roughly doubles the
odds a target reaches approval (Nelson 2015; King 2019; Minikel 2024), the therapeutic direction must be
loss-of-function for a knockdown-nominated gene (a knockdown mimics an inhibitor; a negative regulator is
the wrong way round -- Manguso 2017), the gene must tolerate loss of function in humans (LOEUF / pRec;
Minikel 2020), and it must be tractable. A candidate that clears the screen but fails any of these is not
a credible novel nomination. This script encodes that bar as an explicit, auditable bucketing -- it does
NOT fit a score (RULE #2: transparent lexicographic tiers, never a tuned combination).

Step 2 (implement). Join the committed tables and bucket the 214 safe genes into:
  - VALIDATION: recovered approved-drug targets (the pipeline works; not novel).
  - DEMOTED:    direction-DISCORDANT genes (negative regulators; an inhibitor is anti-therapeutic).
  - SHORTLIST:  genetically supported + LoF-tolerant + tractable + not a known drug + not discordant.
                This is the strict "new targets" answer the user asked for.
  - WIDER:      genetically supported + not a known drug + not discordant, that miss strict LoF-tolerant
                AND tractable together (reported as secondary, with the reason each misses).

Step 3 (critical review) and step 4 (report) happen in docs/results/final_shortlist_*.md; this script
prints every number that document cites and writes results/tables/final_shortlist.csv.

Honesty controls (must all pass, or the synthesis is void):
  - every recovered approved drug lands in VALIDATION, never in SHORTLIST;
  - PTPN2 and RC3H1 land in DEMOTED, never in SHORTLIST;
  - no SHORTLIST gene is DISCORDANT or a known drug (definitional, re-checked as a guard).

Reads results/tables/*.csv only. Never touches the h5ad layers and never calls an LLM.

Usage:
    uv run python scripts/25_final_shortlist.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from cd4_perturbseq import paths

# Recovered approved-drug targets (N10 KNOWN_DRUGS). They validate the gate; they are not novel.
KNOWN_DRUGS = frozenset({"IMPDH2", "PPP3R1", "CD3E", "CD3G", "IL4R", "CD2", "CD28"})


def load() -> pd.DataFrame:
    """Join the safety table, the N10 nomination, and the N11 direction verdict on gene_name.

    Returns:
        One row per safe gene, carrying every axis the shortlist is built from.
    """
    saf = pd.read_csv(paths.TABLES / "window_score_organism_safety.csv")
    nom = pd.read_csv(paths.TABLES / "n10_nomination.csv")
    doe = pd.read_csv(paths.TABLES / "direction_of_effect.csv")

    saf_cols = ["gene_name", "efficacy", "tolerance_loss", "n_module_down", "is_iei",
                "max_nonimmune_ntpm", "lof_intolerant", "loeuf", "prec", "recessive_intolerant", "safe"]
    frame = nom.merge(saf[saf_cols], on="gene_name", how="left", suffixes=("", "_saf"))
    frame = frame.merge(doe[["gene_name", "direction_verdict", "direction_basis",
                             "is_negative_regulator"]], on="gene_name", how="left")
    # loeuf/prec exist in both; keep the nomination copy, fill from safety if missing.
    for col in ("loeuf", "prec", "recessive_intolerant"):
        if f"{col}_saf" in frame:
            frame[col] = frame[col].fillna(frame[f"{col}_saf"])
    return frame


def bucket(frame: pd.DataFrame) -> pd.DataFrame:
    """Assign each safe gene to a credibility bucket by the pre-stated rules (no fitted score)."""
    f = frame.copy()
    f["is_known_drug"] = f["gene_name"].isin(KNOWN_DRUGS)
    # LoF-tolerant, identical rule to N10 nominate(): not pRec-recessive-intolerant and not LoF-constrained.
    f["lof_tolerant"] = (f["prec"].fillna(1.0) <= 0.90) & (~f["lof_intolerant"].fillna(True).astype(bool))
    f["tractable"] = f["tractable"].fillna(False).astype(bool)
    f["supported"] = f["ot_genetic_supported"].fillna(False).astype(bool)
    f["discordant"] = f["direction_verdict"].eq("DISCORDANT")
    f["concordant_or_unknown"] = f["direction_verdict"].isin(["CONCORDANT", "UNKNOWN"])

    strict = f["supported"] & f["lof_tolerant"] & f["tractable"] & ~f["is_known_drug"] & ~f["discordant"]
    wider = f["supported"] & ~f["is_known_drug"] & ~f["discordant"] & ~strict

    f["bucket"] = np.select(
        [f["is_known_drug"], f["discordant"], strict, wider],
        ["VALIDATION", "DEMOTED", "SHORTLIST", "WIDER"],
        default="GATE_ONLY",  # safe + passes the gate, but no genetic support (hypothesis pool)
    )
    return f


def rank_shortlist(short: pd.DataFrame) -> pd.DataFrame:
    """Transparent lexicographic order: direction, held-out support, genetic breadth, screen rank."""
    s = short.copy()
    s["direction_priority"] = s["direction_verdict"].map({"CONCORDANT": 0, "UNKNOWN": 1}).fillna(2)
    s = s.sort_values(
        ["direction_priority", "il2_hit", "ot_genetic_n_diseases", "ot_genetic_max", "window_rank"],
        ascending=[True, False, False, False, True],
    ).reset_index(drop=True)
    s["shortlist_rank"] = np.arange(1, len(s) + 1)
    return s


def _controls(f: pd.DataFrame) -> dict[str, bool]:
    """Falsifiable honesty controls. All must pass or the synthesis is void."""
    known_in_shortlist = f[(f["bucket"] == "SHORTLIST") & f["is_known_drug"]]
    disc_in_shortlist = f[(f["bucket"] == "SHORTLIST") & f["discordant"]]
    known_all_validation = f[f["is_known_drug"]]["bucket"].eq("VALIDATION").all()
    ptpn2_rc3h1 = f[f["gene_name"].isin(["PTPN2", "RC3H1"])]
    demoted_ok = set(ptpn2_rc3h1["bucket"]) <= {"DEMOTED"} and len(ptpn2_rc3h1) == 2
    return {
        "no_known_drug_in_shortlist": known_in_shortlist.empty,
        "no_discordant_in_shortlist": disc_in_shortlist.empty,
        "all_known_drugs_are_validation": bool(known_all_validation),
        "ptpn2_rc3h1_demoted": bool(demoted_ok),
    }


def main() -> None:
    """Bucket the safe genes, run the honesty controls, print the shortlist, write the table."""
    paths.ensure_dirs()
    frame = load()
    f = bucket(frame)

    print(f"universe: {len(f)} safe genes (gate-passing)")
    print("\n=== direction verdicts among safe genes ===")
    print(f["direction_verdict"].value_counts(dropna=False).to_string())
    print("\n=== buckets ===")
    print(f["bucket"].value_counts().to_string())

    controls = _controls(f)
    print("\n=== HONESTY CONTROLS ===")
    for k, v in controls.items():
        print(f"  {'PASS' if v else 'FAIL'}  {k}")
    if not all(controls.values()):
        print("\nVOID: an honesty control failed. The bucketing is wrong; fix before trusting the shortlist.")
        raise SystemExit(1)

    show = ["gene_name", "direction_verdict", "ot_genetic_max", "ot_genetic_n_diseases",
            "loeuf", "prec", "tractable", "il2_hit", "efficacy", "window_rank", "n_module_down", "is_iei"]

    print("\n" + "=" * 90)
    print("VALIDATION (recovered approved-drug targets -- the pipeline recovers known biology, not novel)")
    print("=" * 90)
    val = f[f["bucket"] == "VALIDATION"].sort_values("window_rank")
    print(val[show].to_string(index=False))

    print("\n" + "=" * 90)
    print("DEMOTED (direction-DISCORDANT: negative regulators; a knockdown/inhibitor is anti-therapeutic)")
    print("=" * 90)
    dem = f[f["bucket"] == "DEMOTED"]
    print(dem[show + ["direction_basis"]].to_string(index=False))

    print("\n" + "=" * 90)
    print("SHORTLIST (novel: genetically supported + LoF-tolerant + tractable + not a known drug + not discordant)")
    print("=" * 90)
    short = rank_shortlist(f[f["bucket"] == "SHORTLIST"])
    if short.empty:
        print("  EMPTY.")
    else:
        print(short[["shortlist_rank"] + show].to_string(index=False))

    print("\n" + "=" * 90)
    print("WIDER (genetically supported, novel, concordant-or-unknown, but missing strict LoF-tolerant AND tractable)")
    print("=" * 90)
    wider = f[f["bucket"] == "WIDER"].copy()
    wider["miss_reason"] = np.where(~wider["lof_tolerant"], "not LoF-tolerant",
                             np.where(~wider["tractable"], "not tractable", "?"))
    wider = wider.sort_values(["direction_verdict", "ot_genetic_n_diseases"], ascending=[True, False])
    if wider.empty:
        print("  EMPTY.")
    else:
        print(wider[["gene_name", "direction_verdict", "miss_reason", "ot_genetic_max",
                     "ot_genetic_n_diseases", "loeuf", "prec", "tractable", "il2_hit",
                     "efficacy", "window_rank"]].to_string(index=False))

    # ---- write the committed table (all buckets, so the report can cite any of them)
    out = f.copy()
    out["shortlist_rank"] = np.nan
    if not short.empty:
        out = out.merge(short[["gene_name", "shortlist_rank"]], on="gene_name", how="left",
                        suffixes=("", "_s"))
        out["shortlist_rank"] = out["shortlist_rank_s"].combine_first(out["shortlist_rank"])
        out = out.drop(columns=[c for c in out.columns if c.endswith("_s") or c.endswith("_saf")])
    keep = ["gene_name", "bucket", "shortlist_rank", "direction_verdict", "direction_basis",
            "is_known_drug", "supported", "lof_tolerant", "tractable",
            "ot_genetic_max", "ot_genetic_n_diseases", "loeuf", "prec", "recessive_intolerant",
            "il2_hit", "efficacy", "tolerance_loss", "window_rank", "n_module_down", "is_iei",
            "max_nonimmune_ntpm", "viability_tier"]
    keep = [c for c in keep if c in out.columns]
    out[keep].sort_values(["bucket", "window_rank"]).to_csv(paths.TABLES / "final_shortlist.csv", index=False)
    print(f"\nwrote {paths.TABLES / 'final_shortlist.csv'}")

    print("\n" + "=" * 90)
    print("HONEST HEADLINE (for the report, RULE #8):")
    n_short = len(short)
    n_wider = len(wider)
    print(f"  The pipeline recovers {len(val)} approved-drug targets (validation) and correctly demotes")
    print(f"  {len(dem)} direction-discordant genes (PTPN2, RC3H1). The strict novel shortlist is {n_short} gene(s);")
    print(f"  a wider genetically-supported novel set adds {n_wider}. Direction is UNKNOWN for the novel")
    print("  candidates (N11/N12): they are hypothesis-generating, not vetted nominations.")
    print("=" * 90)


if __name__ == "__main__":
    main()
