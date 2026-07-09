"""N17. Direction of effect via Open Targets GWAS-vs-eQTL colocalisation (v2 proxy).

RULE #9 build. The question from N16: our knockdown screen mimics an inhibitor, so a novel hit is only a
valid autoimmune target if lowering the gene lowers disease (a positive regulator; CONCORDANT). Genetics is
direction-agnostic, so 191/214 safe genes are direction-UNKNOWN. The rigorous fix is colocalisation of the
autoimmune GWAS signal with a cis-eQTL: does the disease-risk direction track higher or lower expression?

Open Targets precomputes this. For each target it exposes QTL credible sets and their colocalisation with
GWAS credible sets, including ``betaRatioSignAverage`` (the harmonised relative sign of the two effects) and
``h4`` (posterior probability the two share one causal variant). We read direction from OT rather than
harmonise variants ourselves (which a feasibility probe showed is fragile: the exact GWAS risk variant is
often absent from served eQTL subsets, and cross-resource allele alignment is error-prone).

CRITICAL discipline (RULE #1): the sign convention is NOT hardcoded. A first read of PTPN2 suggested one
mapping; a calibration on TYK2 (known concordant, deucravacitinib) and TNFAIP3 (known discordant, A20
haploinsufficiency) showed the OPPOSITE. So we LEARN the convention from a panel of genes whose therapeutic
direction is established (approved inhibitor => concordant; monogenic loss-of-function autoimmunity =>
discordant), MEASURE accuracy on that panel, and only emit watchlist verdicts if the panel is classified
well. PTPN2 is allelically heterogeneous (its autoimmune colocalisations split ~evenly), which is itself a
finding: a single per-gene direction label is sometimes not well defined.

This is a v2 proxy, not a bespoke fine-mapping: it inherits OT's coloc calls and cell-type coverage, and a
gene with conflicting colocalisations is reported HETEROGENEOUS, not forced.

Reads results/tables/final_shortlist.csv. Calls the Open Targets GraphQL API (no bulk download). Never
touches the h5ad layers.

Usage:
    uv run python scripts/26_direction_eqtl.py [--h4 0.8] [--genes GENE ...]
"""

from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.request

import pandas as pd

from cd4_perturbseq import paths

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# Calibration panel: genes whose therapeutic direction is established, used to LEARN the sign convention
# and measure accuracy before trusting any novel call. CONCORDANT = an approved inhibitor treats the
# disease (so reducing the gene is therapeutic; higher expression => more disease). DISCORDANT = monogenic
# loss of function causes autoimmunity (so the gene is protective; higher expression => less disease).
PANEL = {
    # concordant, approved inhibitor (the drug IS the therapeutic direction)
    "TYK2": "CONCORDANT",     # deucravacitinib (psoriasis); protective LoF P1104A
    "S1PR1": "CONCORDANT",    # ozanimod/fingolimod (UC, MS)
    "IL6R": "CONCORDANT",     # tocilizumab (RA)
    "TNF": "CONCORDANT",      # infliximab/adalimumab (IBD, RA)
    "ITGA4": "CONCORDANT",    # natalizumab (MS, Crohn)
    # discordant, monogenic loss-of-function autoimmunity (the gene is a brake)
    "CTLA4": "DISCORDANT",    # CTLA4 haploinsufficiency (CHAI)
    "TNFAIP3": "DISCORDANT",  # A20 haploinsufficiency (HA20)
    "SH2B3": "DISCORDANT",    # LNK loss -> autoimmunity
    "LRBA": "DISCORDANT",     # LRBA deficiency -> immune dysregulation
}

AUTOIMMUNE_KEYS = (
    "crohn", "coliti", "bowel", "diabetes", "celiac", "arthrit", "lupus", "sclerosis", "immune",
    "psoria", "autoimmun", "vitiligo", "graves", "thyroid", "asthma", "allerg", "ankylos", "spondyl",
    "rheumat", "sjogren", "juvenile idiopathic",
)

RECOVERED_CONTROLS = ("PTPN2", "RC3H1", "IL4R", "CD28")
"""PTPN2/RC3H1 curated DISCORDANT (N11); IL4R/CD28 recovered inhibitor targets (expect concordant-ish)."""


def _post(query: str, retries: int = 3) -> dict:
    """POST a GraphQL query to Open Targets, with simple retry. Returns the parsed ``data`` dict."""
    body = json.dumps({"query": query}).encode()
    for attempt in range(retries):
        try:
            req = urllib.request.Request(OT_URL, data=body,
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                out = json.loads(resp.read())
            if "errors" in out:
                raise RuntimeError(out["errors"][0].get("message", "graphql error"))
            return out["data"]
        except (urllib.error.URLError, RuntimeError, TimeoutError) as exc:
            if attempt == retries - 1:
                raise
            time.sleep(2 * (attempt + 1))
    return {}


def symbol_to_ensembl(symbol: str) -> str | None:
    """Resolve a gene symbol to its Ensembl ID via the OT search endpoint."""
    q = f'{{search(queryString:"{symbol}",entityNames:["target"]){{hits{{id entity name}}}}}}'
    try:
        hits = _post(q)["search"]["hits"]
    except Exception:
        return None
    for h in hits:
        if h.get("entity") == "target" and (h.get("name") == symbol or h.get("id", "").startswith("ENSG")):
            return h["id"]
    return hits[0]["id"] if hits and hits[0].get("id", "").startswith("ENSG") else None


def gene_direction(ensembl: str, h4_min: float) -> dict:
    """Tally autoimmune GWAS-vs-eQTL colocalisation sign for a gene.

    Returns a dict with the count of positive and negative ``betaRatioSignAverage`` colocalisations
    (h4 >= ``h4_min``, eQTL only, autoimmune traits), the set of traits, and a raw dominant sign.
    """
    q = (f'{{target(ensemblId:"{ensembl}"){{credibleSets{{rows{{studyType '
         f'colocalisation(studyTypes:[gwas]){{rows{{h4 betaRatioSignAverage '
         f'otherStudyLocus{{study{{traitFromSource}}}}}}}}}}}}}}}}')
    try:
        rows = _post(q)["target"]["credibleSets"]["rows"]
    except Exception:
        return {"n_pos": 0, "n_neg": 0, "traits": [], "raw": "NO_DATA"}
    n_pos = n_neg = 0
    traits: set[str] = set()
    for cs in rows:
        if cs.get("studyType") != "eqtl":
            continue
        for co in (cs.get("colocalisation") or {}).get("rows", []):
            trait = (((co.get("otherStudyLocus") or {}).get("study") or {}).get("traitFromSource") or "")
            if not any(k in trait.lower() for k in AUTOIMMUNE_KEYS):
                continue
            if (co.get("h4") or 0) < h4_min:
                continue
            b = co.get("betaRatioSignAverage")
            if b is None:
                continue
            if b > 0:
                n_pos += 1
            elif b < 0:
                n_neg += 1
            traits.add(trait)
    total = n_pos + n_neg
    if total == 0:
        raw = "NO_COLOC"
    elif max(n_pos, n_neg) / total >= 0.7 and total >= 2:
        raw = "POS" if n_pos > n_neg else "NEG"
    else:
        raw = "HETEROGENEOUS"
    return {"n_pos": n_pos, "n_neg": n_neg, "traits": sorted(traits), "raw": raw}


def calibrate(panel_results: dict[str, dict]) -> tuple[str, float, list[str]]:
    """Learn which raw sign (POS/NEG) maps to CONCORDANT from the known-direction panel.

    Returns (concordant_sign, accuracy, notes). Only panel genes with a resolved raw sign (POS/NEG,
    not HETEROGENEOUS/NO_COLOC) vote.
    """
    votes = {"POS": {"CONCORDANT": 0, "DISCORDANT": 0}, "NEG": {"CONCORDANT": 0, "DISCORDANT": 0}}
    notes = []
    for g, expected in PANEL.items():
        raw = panel_results.get(g, {}).get("raw", "NO_DATA")
        if raw in ("POS", "NEG"):
            votes[raw][expected] += 1
        notes.append(f"{g}: expected {expected}, raw {raw} (pos {panel_results.get(g,{}).get('n_pos',0)}/"
                     f"neg {panel_results.get(g,{}).get('n_neg',0)})")
    # If NEG predicts CONCORDANT (as the TYK2/TNFAIP3 probe suggested), score that mapping.
    neg_is_concordant = votes["NEG"]["CONCORDANT"] + votes["POS"]["DISCORDANT"]
    pos_is_concordant = votes["POS"]["CONCORDANT"] + votes["NEG"]["DISCORDANT"]
    resolved = sum(1 for g in PANEL if panel_results.get(g, {}).get("raw") in ("POS", "NEG"))
    if resolved == 0:
        return "UNKNOWN", 0.0, notes
    if neg_is_concordant >= pos_is_concordant:
        return "NEG", neg_is_concordant / resolved, notes
    return "POS", pos_is_concordant / resolved, notes


def verdict(raw: str, concordant_sign: str) -> str:
    """Map a raw dominant sign to a direction verdict using the calibrated convention."""
    if raw in ("NO_DATA", "NO_COLOC"):
        return "UNRESOLVED"
    if raw == "HETEROGENEOUS":
        return "HETEROGENEOUS"
    return "CONCORDANT" if raw == concordant_sign else "DISCORDANT"


def main() -> None:
    """Calibrate the sign convention on the known panel, then classify the watchlist if calibration holds."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h4", type=float, default=0.8)
    ap.add_argument("--genes", nargs="*", help="override gene set (default: panel + watchlist + controls)")
    a = ap.parse_args()
    paths.ensure_dirs()

    short = pd.read_csv(paths.TABLES / "final_shortlist.csv")
    watchlist = short[short["bucket"].isin(["SHORTLIST", "WIDER"])]["gene_name"].tolist()
    genes = a.genes or sorted(set(list(PANEL) + list(RECOVERED_CONTROLS) + watchlist))

    print(f"resolving direction for {len(genes)} genes (h4>={a.h4}) via OT coloc ...\n")
    results: dict[str, dict] = {}
    for i, g in enumerate(genes, 1):
        ens = symbol_to_ensembl(g)
        if not ens:
            results[g] = {"ensembl": None, "raw": "NO_DATA", "n_pos": 0, "n_neg": 0, "traits": []}
            print(f"  [{i}/{len(genes)}] {g}: no Ensembl id")
            continue
        d = gene_direction(ens, a.h4)
        d["ensembl"] = ens
        results[g] = d
        print(f"  [{i}/{len(genes)}] {g} ({ens}): raw {d['raw']} (pos {d['n_pos']}/neg {d['n_neg']})")

    panel_results = {g: results[g] for g in PANEL}
    concordant_sign, accuracy, notes = calibrate(panel_results)
    print("\n=== CALIBRATION on known-direction panel ===")
    for n in notes:
        print("  " + n)
    print(f"\n  learned convention: raw '{concordant_sign}' => CONCORDANT; panel accuracy {accuracy:.2f}")
    trusted = accuracy >= 0.8 and concordant_sign in ("POS", "NEG")
    print(f"  method {'TRUSTED' if trusted else 'NOT TRUSTED'} (need >=0.80 panel accuracy)")

    # Assemble the full table.
    recs = []
    for g in genes:
        d = results[g]
        v = verdict(d["raw"], concordant_sign) if trusted else "UNRESOLVED(method-uncalibrated)"
        role = ("panel:" + PANEL[g]) if g in PANEL else ("control" if g in RECOVERED_CONTROLS else "watchlist")
        recs.append({"gene_name": g, "role": role, "ensembl": d.get("ensembl"),
                     "n_coloc_pos": d["n_pos"], "n_coloc_neg": d["n_neg"], "raw": d["raw"],
                     "direction_eqtl": v, "n_traits": len(d.get("traits", [])),
                     "traits": "; ".join(d.get("traits", [])[:6])})
    tab = pd.DataFrame(recs)
    tab.to_csv(paths.TABLES / "direction_eqtl.csv", index=False)

    print("\n=== CONTROLS ===")
    for g in ("TYK2", "TNFAIP3", "PTPN2"):
        row = tab[tab["gene_name"] == g]
        if len(row):
            r = row.iloc[0]
            print(f"  {g}: {r['direction_eqtl']} (pos {r['n_coloc_pos']}/neg {r['n_coloc_neg']}, raw {r['raw']})")

    if trusted:
        print("\n=== WATCHLIST direction verdicts ===")
        w = tab[tab["role"] == "watchlist"].sort_values("direction_eqtl")
        print(w[["gene_name", "direction_eqtl", "n_coloc_pos", "n_coloc_neg", "n_traits"]].to_string(index=False))
        from collections import Counter
        print("\n  verdict counts:", dict(Counter(w["direction_eqtl"])))

    print(f"\nwrote {paths.TABLES / 'direction_eqtl.csv'}")
    print("\nv2 proxy: inherits OT coloc calls + cell-type coverage; a gene with conflicting colocalisations")
    print("is HETEROGENEOUS, not forced; PTPN2's split is a real finding, not a failure.")


if __name__ == "__main__":
    main()
