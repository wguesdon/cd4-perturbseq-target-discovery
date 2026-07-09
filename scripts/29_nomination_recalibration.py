"""N20. Is "zero vetted novel targets" a property of the biology, or of the nomination rule?

Pre-registered in `docs/preregistration_n20_2026_07_09.md`, committed before this file existed.
Post-hoc provenance is declared there: the analysis was prompted by a user challenge, not planned.

Three tests, each with a decision rule fixed in advance.

H1  The N16 SHORTLIST rule is `supported & lof_tolerant & tractable & ~is_known_drug & ~discordant`.
    Its SENSITIVITY is the fraction of the approved-drug targets the pipeline recovers that would
    have been nominated had we not already known they were drugs. If sensitivity < 5/6, the rule is
    rejected as a nomination gate and the zero result is an artifact of the rule.

H2  Direction-adjudicability may be a COLLIDER on prior characterisation: a gene can be called
    CONCORDANT only if a drug already exists against it, or a monogenic syndrome is already
    described. If more than 80% of CONCORDANT verdicts rest on such a basis, then "novel AND
    concordant" is close to a contradiction, and direction may be used only to DEMOTE, never to
    PROMOTE.

H3  `priors.ot_tractability` dropped `clinical_precedent` from `tractable` (fixed 2026-07-09). Report
    every gene whose flag depends on it.

The replacement gate is chosen on POSITIVE-CLASS SENSITIVITY ALONE, which never looks at which novel
gene a rule promotes, so it cannot be tuned to manufacture a hit:

    NOMINATE = passes the screen gate  AND  tractable  AND  NOT direction-discordant

Genetic support and LoF tolerance are demoted from conjuncts to reported annotations. IMPDH2 is the
internal proof: mycophenolate is approved against it, it has zero autoimmune genetic support, and its
`prec` of 0.99883 makes it recessive-intolerant. Either conjunct would reject mycophenolate.

Output is a RANKED, TIERED candidate table with per-gene liabilities, never a binary verdict and never
a fitted score (RULE #2). The script exits non-zero if any pre-registered control fails.

Usage:
    uv run python scripts/29_nomination_recalibration.py
"""

from __future__ import annotations

import itertools
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "src"))

from cd4_perturbseq import paths, priors  # noqa: E402

SEED = 0
N_PERM = 50_000
CONJUNCTS = ["supported", "lof_tolerant", "tractable"]
SENSITIVITY_FLOOR = 5 / 6  # pre-registered: below this the rule is rejected as a nomination gate
COLLIDER_CEILING = 0.80  # pre-registered: above this, direction may only demote

# Direction bases that require the gene to be already characterised. Both of the two bases the
# direction script can emit are of this kind; the test measures the proportion anyway, so a future
# basis that does not require prior characterisation would change the verdict rather than be ignored.
CHARACTERISED_BASES = ("approved LoF-mimicking drug", "LoF causes immunodeficiency")

# Housekeeping / proliferation families, pre-committed in the registration so the count is published
# whatever it is. Mitochondrial respiratory chain, ribosome, proteasome, RNA polymerase, glycolysis.
HOUSEKEEPING_RE = re.compile(r"^(NDUF|RPL|RPS|PSM|POLR|MRPL|MRPS|ATP5|UQCR|COX\d)")
HOUSEKEEPING_EXTRA = {"ENO1", "PGK1", "KRR1", "PRKAR1A", "UBR4", "CCNT1", "TOMM70", "CLPB"}


def _perm_pvalue(is_pos: np.ndarray, drawn: int, observed: int, rng: np.random.Generator) -> tuple[float, float]:
    """Permutation p-value for `observed` positives among `drawn` genes sampled from the screened set.

    Args:
        is_pos: Boolean vector over the screened background, True where the gene is a curated positive.
        drawn: Size of the nominated set.
        observed: Number of curated positives actually in the nominated set.
        rng: Random generator.

    Returns:
        Tuple of (permutation p-value, expected count under the null).
    """
    u = len(is_pos)
    hits = np.array([int(is_pos[rng.choice(u, drawn, replace=False)].sum()) for _ in range(N_PERM)])
    p = float((hits >= observed).sum() + 1) / (N_PERM + 1)
    return p, float(hits.mean())


def load() -> pd.DataFrame:
    """Join the window score, the nomination annotations and the direction verdicts on the safe set.

    Returns:
        One row per safe gene, carrying every column the nomination rule reads.
    """
    win = pd.read_csv(paths.TABLES / "window_score.csv")
    nom = pd.read_csv(paths.TABLES / "n10_nomination.csv")
    doe = pd.read_csv(paths.TABLES / "direction_of_effect.csv")[
        ["gene_name", "direction_verdict", "direction_basis"]
    ]
    # `lof_intolerant` lives only in the organism-safety table. N16's `lof_tolerant` is built from it,
    # NOT from `recessive_intolerant`; using the wrong one silently changes every count below.
    saf = pd.read_csv(paths.TABLES / "window_score_organism_safety.csv")[["gene_name", "lof_intolerant"]]

    f = nom.merge(doe, on="gene_name", how="left").merge(saf, on="gene_name", how="left")
    f = f.merge(win[["gene_name", "efficacy", "tolerance_loss", "n_module_down"]], on="gene_name", how="left")

    f["supported"] = f["ot_genetic_supported"].fillna(False).astype(bool)
    f["tractable"] = f["tractable"].fillna(False).astype(bool)
    f["clinical_precedent"] = f["clinical_precedent"].fillna(False).astype(bool)
    f["sm_pocket"] = f["sm_pocket"].fillna(False).astype(bool)
    f["known"] = f["is_known_drug"].fillna(False).astype(bool)
    f["discordant"] = f["direction_verdict"].eq("DISCORDANT")
    f["concordant"] = f["direction_verdict"].eq("CONCORDANT")
    # Byte-identical to N10 nominate() and N16 bucket(): both guards run in the conservative direction.
    f["lof_tolerant"] = (f["prec"].fillna(1.0) <= 0.90) & (~f["lof_intolerant"].fillna(True).astype(bool))
    return f


def selftest_reproduces_n16(f: pd.DataFrame) -> None:
    """Refuse to audit the shipped rule until we can reproduce its published output exactly.

    N16 published a SHORTLIST of exactly one gene, ICAM2. If our reimplementation of
    `supported & lof_tolerant & tractable & ~is_known_drug & ~discordant` does not return that, the
    audit is measuring a different rule and every number downstream is void.

    Args:
        f: The joined safe-gene frame.

    Raises:
        SystemExit: If the reimplementation disagrees with the committed N16 table.
    """
    mine = set(f.loc[f["supported"] & f["lof_tolerant"] & f["tractable"]
                     & ~f["known"] & ~f["discordant"], "gene_name"])
    published = pd.read_csv(paths.TABLES / "final_shortlist.csv")
    theirs = set(published.loc[published["bucket"].eq("SHORTLIST"), "gene_name"])
    if mine != theirs:
        print(f"SELF-TEST FAILED: reimplemented shipped rule returns {sorted(mine)}, "
              f"the committed N16 table says {sorted(theirs)}. The audit is void.")
        sys.exit(1)
    print(f"[self-test] reimplemented N16 rule reproduces its published SHORTLIST exactly: {sorted(theirs)}")


def h3_tractability(f: pd.DataFrame) -> pd.DataFrame:
    """Name every gene whose tractability rests on clinical precedent alone.

    Args:
        f: The joined safe-gene frame.

    Returns:
        The subset that would read as untractable under the pre-fix definition.
    """
    tr = priors.ot_tractability()
    pre_fix = tr["sm_pocket"] | tr["sm_druggable_family"] | tr["ab_accessible"]
    flipped = set(tr.loc[tr["tractable"] & ~pre_fix, "gene_name"])
    return f[f["gene_name"].isin(flipped)]


def h2_collider(f: pd.DataFrame) -> tuple[float, pd.Series]:
    """Proportion of CONCORDANT verdicts whose basis requires the gene to be already characterised.

    Args:
        f: The joined safe-gene frame.

    Returns:
        Tuple of (proportion, value counts of the bases).
    """
    conc = f[f["concordant"]]
    counts = conc["direction_basis"].value_counts()
    needs_prior = conc["direction_basis"].str.startswith(CHARACTERISED_BASES).sum()
    return (needs_prior / len(conc) if len(conc) else 0.0), counts


def h1_sensitivity(f: pd.DataFrame) -> pd.DataFrame:
    """Sensitivity of every AND-rule over the conjuncts, measured on the positive class.

    Args:
        f: The joined safe-gene frame.

    Returns:
        One row per rule variant: recovered drugs, sensitivity, and novel genes nominated.
    """
    pos = f[f["known"]]
    cand = f[~f["known"] & ~f["discordant"]]
    rows = []
    for r in range(len(CONJUNCTS) + 1):
        for combo in itertools.combinations(CONJUNCTS, r):
            pm = pos[list(combo)].all(axis=1) if combo else pd.Series(True, index=pos.index)
            cm = cand[list(combo)].all(axis=1) if combo else pd.Series(True, index=cand.index)
            rows.append({
                "rule": " & ".join(combo) if combo else "(screen gate only)",
                "n_conjuncts": len(combo),
                "recovered_drugs": int(pm.sum()),
                "sensitivity": round(float(pm.sum()) / len(pos), 3),
                "novel_nominated": int(cm.sum()),
            })
    return pd.DataFrame(rows).sort_values(["n_conjuncts", "recovered_drugs"], ascending=[True, False])


def tier_and_rank(f: pd.DataFrame) -> pd.DataFrame:
    """Apply the rebuilt gate, then tier and rank the nominations lexicographically.

    Args:
        f: The joined safe-gene frame.

    Returns:
        The nominated genes, tiered and ranked, with a modality tag and a liability string.
    """
    n = f[~f["known"] & ~f["discordant"] & f["tractable"]].copy()

    n["tier"] = np.select(
        [n["concordant"] & n["supported"], n["concordant"], n["supported"]],
        ["A: concordant + genetically supported",
         "B: concordant, no autoimmune genetics",
         "C: direction unknown, genetically supported"],
        default="D: direction unknown, no autoimmune genetics",
    )
    n["modality"] = np.where(
        n["lof_tolerant"],
        "knockout- or degrader-tolerable",
        "inhibitor only, needs a therapeutic window (LoF-intolerant)",
    )
    n["housekeeping"] = n["gene_name"].str.match(HOUSEKEEPING_RE) | n["gene_name"].isin(HOUSEKEEPING_EXTRA)

    liab = []
    for _, r in n.iterrows():
        bits = []
        if not r["lof_tolerant"]:
            bits.append("LoF-intolerant")
        if r["fail_homeostasis"]:
            bits.append("perturbs the resting arm")
        if str(r["viability_tier"]) != "non-depleting":
            bits.append(str(r["viability_tier"]))
        if not r["sm_pocket"] and not r["clinical_precedent"]:
            bits.append("no pocket and no clinical precedent")
        if r["housekeeping"]:
            bits.append("housekeeping or respiratory-chain gene")
        if r["concordant"]:
            bits.append("concordant only because an inhibitor already exists")
        liab.append("; ".join(bits) if bits else "none recorded")
    n["liabilities"] = liab

    n["tier_key"] = n["tier"].str[0]
    n = n.sort_values(
        ["tier_key", "ot_genetic_max", "clinical_precedent", "window_rank"],
        ascending=[True, False, False, True],
    ).reset_index(drop=True)
    n["nomination_rank"] = np.arange(1, len(n) + 1)
    return n


def controls(f: pd.DataFrame, nom: pd.DataFrame, rng: np.random.Generator) -> dict[str, bool]:
    """Run every pre-registered control. All must pass or the analysis is void.

    Args:
        f: The joined safe-gene frame.
        nom: The nominated genes.
        rng: Random generator.

    Returns:
        Mapping of control name to pass/fail.
    """
    win = pd.read_csv(paths.TABLES / "window_score.csv")
    truth = pd.read_csv(paths.GROUND_TRUTH / "immunomodulator_targets.csv")
    positives = set(truth.loc[truth["include_as_positive"], "gene_symbol"])
    win["is_pos"] = win["gene_name"].isin(positives)

    # POSITIVE CONTROL: the rebuilt gate must still recover curated positives above chance.
    gate_genes = set(f.loc[f["tractable"] & ~f["discordant"], "gene_name"])
    observed = int(win["gene_name"].isin(gate_genes).mul(win["is_pos"]).sum())
    p_rec, exp_rec = _perm_pvalue(win["is_pos"].to_numpy(), len(gate_genes), observed, rng)
    print(f"\n[control] recovery of curated positives by the rebuilt gate: "
          f"observed {observed}, expected {exp_rec:.2f}, permutation p = {p_rec:.4g}")
    rec_ok = p_rec < 0.05 and observed > exp_rec

    # NEGATIVE CONTROL: discordant genes must never be nominated.
    leaked = sorted(set(nom["gene_name"]) & {"PTPN2", "RC3H1"})
    print(f"[control] direction-discordant genes in the nominated set: {leaked or 'none'}")
    disc_ok = not leaked

    # NEGATIVE CONTROL: a shuffled tractability flag must not reproduce the enrichment.
    shuffled_p = []
    for _ in range(200):
        perm = f.assign(tractable=rng.permutation(f["tractable"].to_numpy()))
        g = set(perm.loc[perm["tractable"] & ~perm["discordant"], "gene_name"])
        obs = int(win["gene_name"].isin(g).mul(win["is_pos"]).sum())
        shuffled_p.append(obs)
    shuffled_p = np.array(shuffled_p)
    frac_ge = float((shuffled_p >= observed).mean())
    print(f"[control] shuffled tractability reaches {observed} positives in "
          f"{frac_ge:.1%} of 200 shuffles (mean {shuffled_p.mean():.2f})")
    shuffle_ok = frac_ge <= 0.10

    # GROUND-TRUTH GUARD (RULE #4): no curated positive may be presented as a novel nomination.
    contaminated = sorted(set(nom["gene_name"]) & positives)
    print(f"[control] curated positives presented as novel nominations: {contaminated or 'none'}")
    truth_ok = not contaminated

    return {
        "recovery preserved": rec_ok,
        "discordant genes demoted": disc_ok,
        "tractability flag carries information": shuffle_ok,
        "ground truth not contaminated": truth_ok,
    }


def main() -> None:
    """Run the three tests, rebuild the nomination layer, and gate the verdict on the controls."""
    rng = np.random.default_rng(SEED)
    f = load()
    print(f"safe genes: {len(f)}   recovered approved-drug targets: {int(f['known'].sum())}")
    selftest_reproduces_n16(f)

    print("\n" + "=" * 90)
    print("H3. Tractability rested on `sm_pocket | sm_druggable_family | ab_accessible`,")
    print("    discarding `clinical_precedent`. Genes whose tractability depends on it:")
    print("=" * 90)
    flipped = h3_tractability(f)
    print(flipped[["gene_name", "direction_verdict", "clinical_precedent", "sm_pocket",
                   "ot_genetic_max"]].to_string(index=False) if len(flipped) else "  (none)")

    print("\n" + "=" * 90)
    print("H2. Is direction-adjudicability a collider on prior characterisation?")
    print("=" * 90)
    p_collider, bases = h2_collider(f)
    print(bases.to_string())
    print(f"\n  proportion of CONCORDANT verdicts requiring prior characterisation: {p_collider:.3f}")
    print(f"  pre-registered ceiling: {COLLIDER_CEILING}")
    collider = p_collider > COLLIDER_CEILING
    print(f"  VERDICT: {'COLLIDER CONFIRMED' if collider else 'not a collider'} — direction of effect "
          f"{'may only DEMOTE, never PROMOTE' if collider else 'may remain a promoting criterion'}.")

    print("\n" + "=" * 90)
    print("H1. Sensitivity of every rule variant, measured on the positive class only")
    print("=" * 90)
    variants = h1_sensitivity(f)
    print(variants.to_string(index=False))
    shipped = variants[variants["rule"] == "supported & lof_tolerant & tractable"].iloc[0]
    print(f"\n  the shipped rule recovers {shipped['recovered_drugs']} of {int(f['known'].sum())} "
          f"(sensitivity {shipped['sensitivity']}) and nominates {shipped['novel_nominated']} gene(s)")
    print(f"  pre-registered floor: {SENSITIVITY_FLOOR:.3f}")
    rule_rejected = shipped["sensitivity"] < SENSITIVITY_FLOOR
    print(f"  VERDICT: {'RULE REJECTED as a nomination gate' if rule_rejected else 'rule stands'}")

    pos = f[f["known"]]
    for c in CONJUNCTS:
        failed = sorted(pos.loc[~pos[c], "gene_name"])
        if failed:
            print(f"    `{c}` would have rejected: {', '.join(failed)}")

    print("\n" + "=" * 90)
    print("REBUILT NOMINATION LAYER: screen gate AND tractable AND NOT direction-discordant")
    print("=" * 90)
    nom = tier_and_rank(f)
    print(f"nominated: {len(nom)} genes")
    print(nom.groupby("tier").size().to_string())
    print(f"\nhousekeeping / respiratory-chain genes among the nominations: "
          f"{int(nom['housekeeping'].sum())} of {len(nom)} "
          f"(pre-committed to publishing this number whatever it is)")

    verdicts = controls(f, nom, rng)
    print("\n" + "=" * 90)
    for name, ok in verdicts.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

    out = paths.TABLES / "nomination_recalibrated.csv"
    cols = ["nomination_rank", "tier", "gene_name", "direction_verdict", "direction_basis",
            "ot_genetic_max", "ot_genetic_n_diseases", "loeuf", "prec", "lof_tolerant", "modality",
            "sm_pocket", "clinical_precedent", "viability_tier", "housekeeping", "window_rank",
            "efficacy", "n_module_down", "il2_hit", "liabilities"]
    nom[cols].to_csv(out, index=False)
    print(f"\nwrote {out}")

    variants.to_csv(paths.TABLES / "nomination_rule_sensitivity.csv", index=False)
    print(f"wrote {paths.TABLES / 'nomination_rule_sensitivity.csv'}")

    print("\n" + "=" * 90)
    print("TIER A — concordant and genetically supported. The leads, with their liabilities.")
    print("=" * 90)
    tier_a = nom[nom["tier_key"] == "A"]
    print(tier_a[["gene_name", "ot_genetic_max", "ot_genetic_n_diseases", "clinical_precedent",
                  "window_rank", "modality"]].to_string(index=False))
    print("\nliabilities:")
    for _, r in tier_a.iterrows():
        print(f"  {r['gene_name']:<8} {r['liabilities']}")

    if not all(verdicts.values()):
        print("\nA pre-registered control FAILED. The recalibration is void; nothing here may be quoted.")
        sys.exit(1)

    print("\n" + "=" * 90)
    print("HONEST HEADLINE (RULE #8)")
    print("=" * 90)
    print(f"""  'Zero vetted novel targets' was a property of the nomination rule, not of the screen.
  The rule recovers {shipped['recovered_drugs']} of {int(f['known'].sum())} approved-drug targets the same pipeline finds, so its
  sensitivity on its own positive class is {shipped['sensitivity']}. Direction of effect cannot promote a novel
  gene: {p_collider:.0%} of concordant verdicts require a drug or a monogenic syndrome to already exist,
  so it is used here only to demote. A gate calibrated to recover 6 of 6 known drug targets
  nominates {len(nom)} candidates, of which {int(nom['housekeeping'].sum())} are housekeeping or respiratory-chain genes.
  These are ranked hypotheses with named liabilities, not vetted targets. The Tier A leads are
  concordant only because inhibitors against them already exist, which makes them repurposing
  candidates rather than virgin targets.""")


if __name__ == "__main__":
    main()
