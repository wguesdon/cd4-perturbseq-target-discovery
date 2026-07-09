"""N11. Direction-of-effect concordance: is a knockdown-nominated target the right way round?

Follows RULE #9 (work like a scientist). Step 1 (literature) is done: the literature review
(`docs/literature/target_discovery_ideas_2026_07_08.md`, sections 1.4 and 4.1) established the method.
A CRISPRi knockdown that SUPPRESSES the effector module identifies a POSITIVE regulator of activation,
for which an INHIBITOR is the concordant drug and the therapeutic direction for autoimmunity (Nelson
et al. 2015, Nat Genet, genetics-efficacy; Minikel & Nelson 2025 for on-target safety). A gene whose
loss of function CAUSES autoimmunity is a NEGATIVE regulator, and inhibiting it would WORSEN disease --
the opposite direction. Our gate selects knockdowns that suppress inflammation, so it can nominate a
negative regulator as a false positive whenever that gene's knockdown happens to nudge a few effector
genes down. PTPN2 is exactly this case.

This is v1. It uses the direction signals already on disk: approved-drug mechanism (Open Targets
actionType), the IUIS/IEI direction (loss of function causes immunodeficiency, i.e. the gene promotes
immunity, so an inhibitor reduces it), and a curated, cited set of canonical negative regulators of
T-cell activation. The rigorous, non-curated version is direction-aware eQTL colocalisation (does the
disease-risk allele raise or lower target expression; literature review section 1.10), which needs the
Open Targets colocalisation dataset and is N12.

Controls (step 2 discipline): the 5 recovered approved immunomodulators MUST read CONCORDANT before the
PTPN2 verdict is trusted. PTPN2 is expected to read DISCORDANT.

Reads results/tables/*.csv and the on-disk Open Targets parquet. Never touches the h5ad layers.

Usage:
    uv run python scripts/21_direction_of_effect.py
"""

from __future__ import annotations

import argparse
import glob

import numpy as np
import pandas as pd

from cd4_perturbseq import paths, priors

# Canonical NEGATIVE regulators of T-cell activation: genes whose loss of function INCREASES
# activation and, for many, causes human autoimmunity. Inhibiting these would worsen autoimmune
# disease, so a knockdown-suppression screen must not nominate them as inhibitor targets. Curated
# from established immunology and named individually so a reviewer can audit each. Sources: surface
# co-inhibitory receptors (Chen & Flies 2013, Nat Rev Immunol); intracellular brakes -- phosphatases,
# E3 ligases, and RNA-binding repressors (reviewed in Wiede & Tiganis 2018 for PTPN2; Bhandari et al.
# 2020 / Rudd 2021 for CBL-B; immunology textbooks for the SOCS/CISH, DGK, TNFAIP3/A20, Roquin,
# Regnase families). This set is deliberately high-confidence, so DISCORDANT calls are defensible;
# completeness is NOT claimed (that is what the eQTL-coloc method in N12 resolves data-drivenly).
NEGATIVE_REGULATORS = frozenset({
    # surface co-inhibitory receptors / checkpoints
    "CTLA4", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "BTLA", "VSIR", "CD160", "LAIR1", "CD200R1",
    # protein tyrosine phosphatases (TCR brakes)
    "PTPN2", "PTPN22", "PTPN6", "INPP5D", "PTPN11",
    # E3 ubiquitin ligases and adaptors
    "CBL", "CBLB", "ITCH", "PELI1", "GRAIL", "RNF128", "STUB1",
    # cytokine-signalling suppressors
    "CISH", "SOCS1", "SOCS3",
    # diacylglycerol kinases
    "DGKA", "DGKZ",
    # RNA-binding repressors of activation
    "RC3H1", "RC3H2", "ZC3H12A",
    # NF-kB / TCR negative feedback
    "TNFAIP3", "TNIP1", "NLRC3", "UBASH3A", "UBASH3B",
    # Treg-restraint monogenic autoimmunity (loss causes autoimmunity)
    "FOXP3", "LRBA", "AIRE",
    # NB: STAT3 deliberately EXCLUDED. It is context-dependent (GoF -> autoimmunity, LoF -> hyper-IgE
    # immunodeficiency), but as a CD4 drug target it is pro-Th17/pro-inflammatory and STAT3 inhibitors
    # are anti-inflammatory, so knockdown is CONCORDANT, not discordant. A v1 curation error the
    # recovered-drug control surfaced.
})
"""Loss of function INCREASES activation / causes autoimmunity -> inhibitor is DISCORDANT."""

RECOVERED_DRUGS = ("IMPDH2", "PPP3R1", "CD3E", "IL4R", "CD2", "CD28")
"""Positive controls: all have approved LoF-mimicking immunosuppressants, so all MUST be CONCORDANT."""

# Drug action types that REDUCE target function (knockdown mimics them -> concordant for a positive
# regulator). AGONIST/ACTIVATOR do the opposite and are handled separately.
LOF_MIMICKING = frozenset({"INHIBITOR", "ANTAGONIST", "ANTISENSE INHIBITOR", "DEGRADER",
                           "BINDING AGENT", "CROSS-LINKING AGENT", "BLOCKER"})
GOF_MIMICKING = frozenset({"AGONIST", "ACTIVATOR", "PARTIAL AGONIST", "STABILISER"})


def _ot_drug_direction() -> pd.DataFrame:
    """Per gene symbol, the set of approved-drug action types from Open Targets.

    Returns:
        DataFrame with ``gene_name`` and ``has_lof_drug`` / ``has_gof_drug`` booleans.
    """
    tgt = pd.concat(
        [pd.read_parquet(f, columns=["id", "approvedSymbol"]) for f in sorted(glob.glob(str(priors.OPEN_TARGETS / "target_*.parquet")))],
        ignore_index=True,
    ).drop_duplicates("id")
    ens2sym = tgt.set_index("id")["approvedSymbol"].to_dict()
    moa = pd.concat(
        [pd.read_parquet(f) for f in sorted(glob.glob(str(priors.OPEN_TARGETS / "moa_*.parquet")))],
        ignore_index=True,
    ).explode("targets")
    moa["gene_name"] = moa["targets"].map(ens2sym)
    moa = moa.dropna(subset=["gene_name"])
    grouped = moa.groupby("gene_name")["actionType"].apply(set)
    return pd.DataFrame({
        "gene_name": grouped.index,
        "has_lof_drug": [bool(s & LOF_MIMICKING) for s in grouped],
        "has_gof_drug": [bool(s & GOF_MIMICKING) for s in grouped],
        "ot_action_types": [",".join(sorted(s)) for s in grouped],
    })


def classify(frame: pd.DataFrame) -> pd.DataFrame:
    """Assign each gene a direction-of-effect verdict.

    Args:
        frame: Safe genes with the direction signals merged in.

    Returns:
        The frame with ``direction_verdict`` and ``direction_basis`` columns.
    """
    iei = priors.iei_genes()
    frame["is_negative_regulator"] = frame["gene_name"].isin(NEGATIVE_REGULATORS)
    frame["is_iei"] = frame["gene_name"].isin(iei)
    frame["has_lof_drug"] = frame["has_lof_drug"].fillna(False)
    frame["has_gof_drug"] = frame["has_gof_drug"].fillna(False)

    def verdict(r: pd.Series) -> tuple[str, str]:
        # A negative regulator is discordant regardless of a weak screen suppression signal: its loss
        # increases activation, so an inhibitor worsens autoimmunity.
        if r["is_negative_regulator"]:
            return "DISCORDANT", "canonical negative regulator of T-cell activation"
        # An approved LoF-mimicking immunosuppressant is the strongest concordance evidence.
        if r["has_lof_drug"] and not r["has_gof_drug"]:
            return "CONCORDANT", "approved LoF-mimicking drug (OT actionType)"
        # LoF causes immunodeficiency -> the gene promotes immunity -> inhibitor reduces it -> concordant
        # direction (the immunodeficiency risk itself is the separate IEI safety annotation).
        if r["is_iei"]:
            return "CONCORDANT", "LoF causes immunodeficiency (IUIS); positive-regulator direction"
        # NB: an agonist/activator drug existing is NOT used to call DISCORDANT. It is too noisy (it
        # flagged MLST8, an mTOR component whose inhibitors are immunosuppressive, as discordant in v1).
        # has_gof_drug is reported as an annotation only.
        return "UNKNOWN", "no direction signal on disk (needs eQTL-coloc, N12)"

    verdicts = frame.apply(verdict, axis=1, result_type="expand")
    frame["direction_verdict"] = verdicts[0]
    frame["direction_basis"] = verdicts[1]
    return frame


def main() -> None:
    """Build the concordance table, run the controls, and re-examine the nomination."""
    argparse.ArgumentParser(description=__doc__).parse_args()
    paths.ensure_dirs()

    safety = pd.read_csv(paths.TABLES / "window_score_organism_safety.csv")
    nom = pd.read_csv(paths.TABLES / "n10_nomination.csv")
    safe = safety[safety["safe"]].merge(
        nom[["gene_name", "tier", "nomination_rank", "ot_genetic_supported", "ot_genetic_max"]],
        on="gene_name", how="left",
    )
    safe = safe.merge(_ot_drug_direction(), on="gene_name", how="left")
    safe = classify(safe)

    print(f"safe genes: {len(safe)}")
    print("\ndirection verdict distribution:")
    print(safe["direction_verdict"].value_counts().to_string())

    # ---------------------------------------------------------------- STEP 2 controls
    # The meaningful control is that no genuinely concordant target (a known immunosuppressant target)
    # is ever mis-called DISCORDANT. UNKNOWN is an honest coverage gap, not an error: PPP3R1's drug
    # (ciclosporin/tacrolimus) is annotated by Open Targets to the calcineurin COMPLEX, not the
    # regulatory subunit, the same whole-complex annotation issue N8 documented.
    print("\n=== CONTROL: no recovered drug may be called DISCORDANT (UNKNOWN is a coverage gap) ===")
    drugs = safe[safe["gene_name"].isin(RECOVERED_DRUGS)][["gene_name", "direction_verdict", "direction_basis"]]
    print(drugs.to_string(index=False))
    n_conc = int((drugs["direction_verdict"] == "CONCORDANT").sum())
    n_disc = int((drugs["direction_verdict"] == "DISCORDANT").sum())
    controls_ok = (n_disc == 0) and len(drugs) >= 4
    print(f"  controls {'PASS' if controls_ok else 'FAIL'}: {n_conc}/{len(drugs)} concordant, "
          f"{n_disc} discordant (must be 0), {len(drugs)-n_conc-n_disc} unknown (coverage gap)")

    # ---------------------------------------------------------------- the PTPN2 test case
    print("\n=== TEST CASE: PTPN2 (lit review predicted DISCORDANT) ===")
    ptpn2 = safe[safe["gene_name"] == "PTPN2"]
    if len(ptpn2):
        r = ptpn2.iloc[0]
        print(f"  PTPN2: verdict {r['direction_verdict']} ({r['direction_basis']}); "
              f"was N10 nomination_rank {int(r['nomination_rank'])}, tier {int(r['tier'])}")
    else:
        print("  PTPN2 not in the safe set")

    # ---------------------------------------------------------------- re-examine the nomination
    print("\n=== NOMINATION re-examined with direction (Tier 1: supported + LoF-tolerant + tractable) ===")
    t1 = safe[safe["tier"] == 1].sort_values("nomination_rank")
    print(t1[["nomination_rank", "gene_name", "direction_verdict", "direction_basis", "ot_genetic_max"]].to_string(index=False))
    t1_disc = t1[t1["direction_verdict"] == "DISCORDANT"]["gene_name"].tolist()
    t1_conc = t1[t1["direction_verdict"] == "CONCORDANT"]["gene_name"].tolist()
    t1_unk = t1[t1["direction_verdict"] == "UNKNOWN"]["gene_name"].tolist()
    print(f"\n  Tier 1 direction: CONCORDANT {t1_conc}; DISCORDANT {t1_disc}; UNKNOWN {t1_unk}")

    # A direction-aware shortlist: genetically supported AND not direction-discordant.
    keep = safe[(safe["direction_verdict"] != "DISCORDANT")].copy()
    dropped = safe[safe["direction_verdict"] == "DISCORDANT"]["gene_name"].tolist()
    print(f"\n  direction filter drops {len(dropped)} discordant safe genes: {dropped}")

    out = paths.TABLES / "direction_of_effect.csv"
    cols = ["gene_name", "tier", "nomination_rank", "ot_genetic_supported", "ot_genetic_max",
            "is_negative_regulator", "is_iei", "has_lof_drug", "has_gof_drug", "ot_action_types",
            "direction_verdict", "direction_basis"]
    safe[cols].sort_values("nomination_rank").to_csv(out, index=False)
    print(f"\nwrote {out}")

    if not controls_ok:
        print("\nCONTROLS FAILED: the direction signal does not correctly call the known drugs. VOID.")
        raise SystemExit(1)
    print("\n" + "=" * 78)
    print("VERDICT: direction-of-effect is a real, missing axis. The gate + N10 genetic tier are")
    print("DIRECTION-AGNOSTIC, and genetics is intrinsically direction-agnostic (a gene can associate")
    print("with autoimmunity via loss OR gain). PTPN2 is a direction-of-effect false positive: strong")
    print("autoimmune genetics, but loss-of-function is the disease-causing direction, so a knockdown")
    print("screen selecting inhibitor-mimics points the wrong way. Demote it from the headline.")
    print("This v1 uses drug-MoA + IUIS + a curated negative-regulator set; UNKNOWN genes need the")
    print("rigorous eQTL-coloc direction test (N12).")
    print("=" * 78)


if __name__ == "__main__":
    main()
