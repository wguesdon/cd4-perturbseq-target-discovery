"""Scratch verification: check the peer-critic's factual claims about ICAM2 against committed tables.

Read-only. Confirms every number before it enters the shortlist report (RULE #8 / N15 firewall lesson).
"""
from __future__ import annotations
import pandas as pd
from cd4_perturbseq import paths

saf = pd.read_csv(paths.TABLES / "window_score_organism_safety.csv").set_index("gene_name")
nom = pd.read_csv(paths.TABLES / "n10_nomination.csv").set_index("gene_name")

cols = ["viability_tier", "homeostasis_verdict", "fail_homeostasis", "fail_evidence", "fail_tolerance",
        "safe", "max_nonimmune_ntpm", "n_nonimmune_tissues", "broadly_expressed", "efficacy",
        "tolerance_loss", "n_module_down", "window_rank", "loeuf", "prec", "recessive_intolerant"]
print("=== ICAM2 record (window_score_organism_safety.csv) ===")
print(saf.loc["ICAM2", cols].to_string())
print("\n=== ICAM2 (n10_nomination.csv) ===")
print(nom.loc["ICAM2", ["tier", "fail_homeostasis", "ot_genetic_max", "ot_genetic_n_diseases",
                        "tractable", "sm_pocket", "clinical_precedent", "il2_hit", "is_known_drug"]].to_string())

print("\n=== GATE_ONLY top window ranks (safe, not genetically supported) — critic named PRKAR1A/KRR1/ENO1 ===")
safe = saf[saf["safe"]].copy()
supported = set(nom[nom["ot_genetic_supported"]].index) if "ot_genetic_supported" in nom else set()
gate_only = safe[~safe.index.isin(supported)].sort_values("window_rank").head(6)
print(gate_only[["window_rank", "viability_tier", "fail_homeostasis", "efficacy", "n_module_down"]].to_string())

print("\n=== depleting-at-rest among the 6 VALIDATION genes? ===")
for g in ["IMPDH2", "CD2", "CD3E", "CD28", "PPP3R1", "IL4R"]:
    if g in saf.index:
        print(f"  {g}: {saf.loc[g,'viability_tier']}, fail_homeostasis={saf.loc[g,'fail_homeostasis']}")

print("\n=== recovery statistics (recovery_pvalue.csv) ===")
try:
    rec = pd.read_csv(paths.TABLES / "recovery_pvalue.csv")
    print(rec.to_string(index=False))
except Exception as e:
    print("  could not read:", e)
