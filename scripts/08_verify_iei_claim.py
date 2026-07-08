"""Verify the audit's headline-invalidating claim, myself, from the data.

CLAIM: the IUIS inborn-errors-of-immunity flag is MORE enriched among our own approved-drug
positives than among the naive top-100. If so, "is an IEI gene" cannot separate a toxic
knockdown from a good drug target, and our safety gate is partly tautological.

Also test whether the IUIS parsing is the problem. The IUIS table has a `GOF/DN` column. Genes
whose disease mechanism is GAIN of function should never have counted: our gate assumes LOSS of
function causes immunodeficiency, which is the direction CRISPRi produces.
"""

from __future__ import annotations

import sys

import pandas as pd
from scipy import stats

REPO = "/home/will/Documents/Github/cd4-perturbseq-target-discovery"
sys.path.insert(0, f"{REPO}/src")

from cd4_perturbseq import priors  # noqa: E402

pd.set_option("display.width", 200)

# ---------------------------------------------------------------- the IUIS table, properly read
raw = pd.read_csv(f"{REPO}/data/external/gwt_priors/IUIS-IEI-list-July-2024V2.csv")
raw.columns = [c.strip() for c in raw.columns]
print("IUIS columns:", [c for c in raw.columns][:6])
print(f"rows: {len(raw)}")

gof = raw["GOF/DN"].astype(str).str.strip().str.upper()
print("\nGOF/DN values:")
print(gof.value_counts().head(8).to_string())

defects = raw["Genetic defect"].astype(str).str.strip()
symbol_like = defects.str.fullmatch(r"[A-Z][A-Z0-9\-]{1,14}")

current = set(defects[symbol_like.fillna(False)])
print(f"\ncurrent priors.iei_genes() equivalent: {len(current)} symbols")

# Exclude entries whose stated mechanism is gain of function or dominant negative only.
is_gof_only = gof.isin({"GOF"})
lof_like = symbol_like.fillna(False) & ~is_gof_only
lof_only = set(defects[lof_like])
print(f"after dropping GOF-only entries:        {len(lof_only)} symbols  (-{len(current)-len(lof_only)})")

dropped = sorted(current - lof_only)
print(f"  dropped as gain-of-function: {', '.join(dropped[:18])}{' ...' if len(dropped) > 18 else ''}")

# Which real symbols does the regex silently drop?
not_matched = sorted(set(defects[~symbol_like.fillna(False)]))
plausible = [d for d in not_matched if d.isupper() and 2 <= len(d) <= 15 and " " not in d]
print(f"\nregex drops {len(not_matched)} entries; plausible gene symbols among them:")
print("  " + ", ".join(plausible[:20]))

# ---------------------------------------------------------------- the actual claim
ranked = pd.read_csv(f"{REPO}/results/tables/risk_kill_naive_reversal.csv").rename(
    columns={"target_contrast_gene_name": "gene_name"}
)
truth = pd.read_csv(f"{REPO}/resources/ground_truth/immunomodulator_targets.csv")
positives = set(truth.loc[truth["include_as_positive"], "gene_symbol"])


def fisher(flag: pd.Series, group: pd.Series) -> tuple[float, float, int, int]:
    """Fisher exact, one-sided greater, of `flag` enrichment inside `group`."""
    a = int((flag & group).sum())
    b = int((~flag & group).sum())
    c = int((flag & ~group).sum())
    d = int((~flag & ~group).sum())
    odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    return odds, p, a, int(group.sum())


for label, iei_set in (("CURRENT (regex, GOF included)", current), ("LoF-only (GOF dropped)", lof_only)):
    ranked["is_iei"] = ranked["gene_name"].isin(iei_set)
    ranked["is_positive"] = ranked["gene_name"].isin(positives)
    ranked["in_top100"] = ranked["rank"] <= 100

    print(f"\n=== {label} ===")
    bg_rate = ranked.loc[~ranked["in_top100"] & ~ranked["is_positive"], "is_iei"].mean()
    print(f"  background IEI rate: {bg_rate:.2%}")

    odds_t, p_t, n_t, tot_t = fisher(ranked["is_iei"], ranked["in_top100"])
    print(f"  naive top-100    : {n_t}/{tot_t} are IEI ({n_t/tot_t:.1%})   OR={odds_t:.2f}  p={p_t:.3g}")

    odds_p, p_p, n_p, tot_p = fisher(ranked["is_iei"], ranked["is_positive"])
    print(f"  approved positives: {n_p}/{tot_p} are IEI ({n_p/max(tot_p,1):.1%})   OR={odds_p:.2f}  p={p_p:.3g}")

    verdict = "IEI CANNOT separate toxic from on-pathway" if odds_p >= odds_t else "IEI separates them"
    print(f"  -> {verdict}")

    ieipos = sorted(set(ranked.loc[ranked["is_positive"] & ranked["is_iei"], "gene_name"]))
    print(f"  approved drug targets our gate would REJECT as IEI: {', '.join(ieipos) if ieipos else 'none'}")

# ---------------------------------------------------------------- which positives are IUIS at all
print("\n=== every curated positive, IEI status (regardless of rankability) ===")
for name, iei_set in (("current", current), ("LoF-only", lof_only)):
    hits = sorted(positives & iei_set)
    print(f"  {name:9s}: {len(hits)}/{len(positives)} -> {', '.join(hits)}")
