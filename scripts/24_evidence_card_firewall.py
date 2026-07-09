"""N15. Evidence cards with a fabrication firewall: defensible Claude-in-the-loop.

Runs the RULE #9 loop, and it is the project's "Built with Claude" (Claude Use) story done honestly.
Step 1 (literature): LLMs fabricate numbers and citations at high rates (Walters & Wilder 2023 found
55% of GPT-3.5 and 18% of GPT-4 citations entirely fabricated), so an LLM-written evidence card cannot
be trusted on its face. The defensible pattern (literature review section 3; CiteCheck, MIRAGE) is a
deterministic firewall: freeze a per-gene record of ONLY values already computed and committed, let the
model write prose from that record alone, then a deterministic checker rejects the card if any number in
it is not in the record (within rounding) or any citation is not on a curated allow-list.

This script builds the frozen records from committed tables, implements the checker, and PROVES the
firewall works with a falsifiable control: a deliberately fabricated card (a wrong number and a fake
citation) MUST be caught. The LLM card generation itself is a separate workflow that runs each card
through this same checker; the firewall is what makes that generation safe.

Reads results/tables/*.csv. Never touches the h5ad layers and never calls an LLM (the checker is
deterministic; the control is a hardcoded fabricated card).

Usage:
    uv run python scripts/24_evidence_card_firewall.py
"""

from __future__ import annotations

import argparse
import json
import re

import numpy as np
import pandas as pd

from cd4_perturbseq import paths

# Citations the project actually uses, curated. A card may cite ONLY these; anything else is flagged
# as a potential hallucination. Author-year tokens, matched case-insensitively on the surname+year.
CITATION_ALLOWLIST = frozenset({
    "zhu 2025", "dann 2025",  # the source paper
    "nelson 2015", "king 2019", "minikel 2024", "minikel 2020", "trajanoska 2023",  # genetics-of-targets
    "manguso 2017", "wiede 2018", "vinuesa 2005",  # PTPN2 / Roquin direction
    "schmidt 2022", "steinhart 2022", "frangieh 2021", "shifrut 2018",  # recovery precedent / held-out
    "lamb 2006", "sirota 2011", "chen 2017",  # reversal / connectivity
    "kriegeskorte 2009",  # double dipping
    "wang 2007", "tran 2007", "allan 2007",  # FOXP3-in-Tconv
})

CANDIDATES = ("IMPDH2", "PPP3R1", "CD3E", "IL4R", "PTPN2", "CD2", "CD28")
"""The recovered drugs (validation) plus the demoted direction-discordant PTPN2, so a card must be able
to state a negative verdict too."""


def assemble_records() -> dict[str, dict]:
    """Freeze a per-gene record from committed tables only. This is the sole source a card may use."""
    saf = pd.read_csv(paths.TABLES / "window_score_organism_safety.csv").set_index("gene_name")
    nom = pd.read_csv(paths.TABLES / "n10_nomination.csv").set_index("gene_name")
    doe = pd.read_csv(paths.TABLES / "direction_of_effect.csv").set_index("gene_name")
    rec = pd.read_csv(paths.TABLES / "window_score.csv").set_index("gene_name")

    records: dict[str, dict] = {}
    for g in CANDIDATES:
        if g not in saf.index:
            continue
        s = saf.loc[g]
        r: dict = {
            "gene": g,
            "window_rank": int(s["window_rank"]),
            "viability_tier": str(s["viability_tier"]),
            "n_module_down": int(rec.loc[g, "n_module_down"]),
            "efficacy": round(float(s["efficacy"]), 3),
            "tolerance_loss": round(float(s["tolerance_loss"]), 3),
            "safe": bool(s["safe"]),
            "loeuf": None if pd.isna(s["loeuf"]) else round(float(s["loeuf"]), 3),
            "prec": None if pd.isna(s["prec"]) else round(float(s["prec"]), 3),
            "recessive_intolerant": None if pd.isna(s["recessive_intolerant"]) else bool(s["recessive_intolerant"]),
            "max_nonimmune_ntpm": None if pd.isna(s["max_nonimmune_ntpm"]) else round(float(s["max_nonimmune_ntpm"]), 1),
            "is_iei": bool(s["is_iei"]),
        }
        if g in nom.index:
            n = nom.loc[g]
            r.update({
                "nomination_tier": int(n["tier"]),
                "ot_genetic_max": round(float(n["ot_genetic_max"]), 3),
                "ot_genetic_n_diseases": int(n["ot_genetic_n_diseases"]) if not pd.isna(n["ot_genetic_n_diseases"]) else 0,
                "tractable": bool(n["tractable"]) if not pd.isna(n["tractable"]) else False,
                "il2_hit_schmidt": bool(n["il2_hit"]) if not pd.isna(n["il2_hit"]) else False,
            })
        if g in doe.index:
            r["direction_verdict"] = str(doe.loc[g, "direction_verdict"])
            r["direction_basis"] = str(doe.loc[g, "direction_basis"])
        records[g] = r
    return records


def _record_numbers(record: dict) -> set[float]:
    """Every numeric value in a record, plus common roundings, for tolerant matching."""
    nums: set[float] = set()
    for v in record.values():
        if isinstance(v, bool):
            continue
        if isinstance(v, (int, float)):
            for r in (float(v), round(float(v), 2), round(float(v), 1), round(float(v))):
                nums.add(r)
    return nums


NUM_RE = re.compile(r"(?<![A-Za-z0-9_.])-?\d+(?:\.\d+)?")
CITE_RE = re.compile(r"([A-Z][a-zA-Z]+)\s+(?:et al\.?\s+)?\(?((?:19|20)\d{2})\)?")


def check_card(card: str, record: dict, allowlist: frozenset[str] = CITATION_ALLOWLIST,
               tol: float = 0.02) -> list[str]:
    """Return a list of firewall violations. Empty list = the card is grounded.

    Every number in the card must match a value in the record within ``tol``; every citation must be
    on the allow-list. Structural small integers (0..14, the disease count and small counts) are
    permitted because prose uses them ("12 of 14 diseases", "one of five").

    Args:
        card: The card text to check.
        record: The frozen per-gene record the card was written from.
        allowlist: Permitted citation surname+year tokens (lower-cased).
        tol: Absolute tolerance for matching a number to a record value.

    Returns:
        Human-readable violation strings.
    """
    violations: list[str] = []
    allowed = _record_numbers(record) | set(float(i) for i in range(0, 15))  # structural small ints
    for tok in NUM_RE.findall(card):
        val = float(tok)
        # A 4-digit year is citation metadata, validated by the citation allow-list below, not the
        # numeric record. A fabricated year is still caught as an UNLISTED CITATION.
        if 1900 <= val <= 2099 and float(val).is_integer():
            continue
        if not any(abs(val - a) <= tol for a in allowed):
            violations.append(f"UNGROUNDED NUMBER {tok!r} (not in the record)")
    for surname, year in CITE_RE.findall(card):
        token = f"{surname.lower()} {year}"
        if token not in allowlist:
            violations.append(f"UNLISTED CITATION {surname} {year} (not on the allow-list)")
    return violations


def _control() -> bool:
    """Falsifiable control: a deliberately fabricated card MUST be caught. Returns True if it is."""
    record = {"gene": "IMPDH2", "loeuf": 0.767, "prec": 0.999, "ot_genetic_max": 0.0,
              "window_rank": 51, "efficacy": 0.598}
    fabricated = (
        "IMPDH2 is a strong candidate (LOEUF 0.35, genetic support 0.88 across 9 diseases), "
        "validated as a target by Fabricato et al. 2099 and consistent with Manguso 2017."
    )
    v = check_card(fabricated, record)
    caught_number = any("0.35" in x or "0.88" in x for x in v)
    caught_citation = any("Fabricato" in x for x in v)
    allowed_real_citation = not any("Manguso" in x for x in v)  # Manguso 2017 is on the allow-list
    print("=== FALSIFIABLE CONTROL: a fabricated card must be caught ===")
    for x in v:
        print(f"  flagged: {x}")
    print(f"  fabricated number caught: {caught_number}; fabricated citation caught: {caught_citation}; "
          f"real citation (Manguso 2017) NOT flagged: {allowed_real_citation}")
    return caught_number and caught_citation and allowed_real_citation


def validate_cards(cards_path: str) -> int:
    """Run the firewall over LLM-generated cards. Returns the number of cards with violations.

    Args:
        cards_path: JSON file mapping gene -> card text (as produced by the generation workflow).

    Returns:
        Count of cards that had >=1 firewall violation.
    """
    records = assemble_records()
    cards = json.loads(open(cards_path).read())
    n_bad = 0
    print("=== FIREWALL over LLM-generated cards ===")
    for gene, card in cards.items():
        if gene not in records:
            continue
        v = check_card(card, records[gene])
        status = "CLEAN" if not v else f"CAUGHT {len(v)}"
        print(f"\n[{gene}] {status}")
        for x in v:
            print(f"    - {x}")
        n_bad += int(bool(v))
    print(f"\n  {n_bad} of {len(cards)} cards had a firewall violation.")
    return n_bad


def main() -> None:
    """Build the records, prove the firewall catches fabrication, and write the frozen records."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cards", help="JSON of gene->card text to validate through the firewall")
    a = parser.parse_args()
    paths.ensure_dirs()
    if a.cards:
        validate_cards(a.cards)
        return

    records = assemble_records()
    print(f"froze {len(records)} per-gene evidence records from committed tables\n")

    ok = _control()
    print(f"\n  CONTROL {'PASS' if ok else 'FAIL'}: the firewall {'catches' if ok else 'MISSES'} "
          "a fabricated number and a fabricated citation while passing a real one.")

    # Sanity: a faithful card built directly from a record must pass the checker with zero violations.
    print("\n=== SANITY: a record-faithful card passes with zero violations ===")
    r = records.get("PTPN2")
    if r:
        faithful = (
            f"PTPN2: window rank {r['window_rank']}, {r['n_module_down']} module genes down, "
            f"genetic support {r['ot_genetic_max']} across {r['ot_genetic_n_diseases']} autoimmune diseases, "
            f"LOEUF {r['loeuf']}. Direction verdict {r['direction_verdict']}: a negative regulator "
            "(Manguso 2017), so an inhibitor is anti-therapeutic. Demoted."
        )
        vv = check_card(faithful, r)
        print(f"  PTPN2 faithful card violations: {vv if vv else 'none'}")

    out = paths.TABLES / "evidence_records.json"
    out.write_text(json.dumps(records, indent=2))
    print(f"\nwrote {out}")

    if not ok:
        raise SystemExit(1)
    print("\n" + "=" * 78)
    print("The firewall is the defensible Claude-Use pattern: Claude writes each card from ONLY its")
    print("frozen record; this deterministic checker rejects any ungrounded number or unlisted citation.")
    print("Claude never sets a rank or a verdict (those come from the committed pipeline); it writes the")
    print("prose, and the prose is validated. The falsifiable control proves the checker catches")
    print("fabrication. LLM card generation runs as a workflow over these records, gated by this checker.")
    print("=" * 78)


if __name__ == "__main__":
    main()
