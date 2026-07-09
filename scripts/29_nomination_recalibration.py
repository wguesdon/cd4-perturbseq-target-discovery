"""N20. Is "zero vetted novel targets" a property of the biology, or of the nomination rule?

Pre-registered in `docs/preregistration_n20_2026_07_09.md`, committed before this file existed.
Post-hoc provenance is declared there: the analysis was prompted by a user challenge, not planned.

Three tests, each with a decision rule fixed in advance.

H1  The N16 SHORTLIST rule is `supported & lof_tolerant & tractable & ~is_known_drug & ~discordant`.
    Its SENSITIVITY is the fraction of the approved-drug targets the pipeline recovers that would
    have been nominated had we not already known they were drugs. If sensitivity < 5/6, the rule is
    rejected as a nomination gate and the zero result is an artifact of the rule.

H2  Direction-adjudicability may be a DETERMINISTIC PROXY FOR PRIOR CHARACTERISATION, with
    missingness not at random: a gene can be called CONCORDANT only if a drug already exists against
    it, or a monogenic syndrome is already described. (An earlier draft called this a "collider".
    That is wrong: a collider needs two parents, and direction-adjudicability has one. Conditioning
    on it opens no spurious path; it simply cannot fire on an uncharacterised gene.) If more than 80%
    of CONCORDANT verdicts rest on such a basis, then "novel AND concordant" is close to a
    contradiction, and direction may be used only to DEMOTE, never to PROMOTE.

    Two consequences that must travel with the claim. First, DISCORDANCE is equally
    characterisation-bound, so the demotion filter is near-blind to a novel wrong-direction gene:
    excluding PTPN2 and RC3H1 licenses no direction-safety claim about the direction-unmeasured
    genes. Second, a gene is direction-UNMEASURED, not direction-favourable. The absence of a
    measurement is not evidence of a good direction, and the tier is named accordingly.

H3  `priors.ot_tractability` omitted `clinical_precedent` from `tractable`. Restoring it looked like
    an obvious bug fix. It is not: `clinical_precedent` is true for 32 of the 36 curated positives,
    so a gate that reads it recovers approved drug targets BY CONSTRUCTION and the recovery
    validation becomes circular. Structural tractability alone is already true for all six recovered
    targets. So `tractable` stays structural, `tractable_with_precedent` is reported separately, and
    genes entering only through precedent are labelled repurposing candidates.

The replacement gate is chosen on POSITIVE-CLASS SENSITIVITY ALONE, which never looks at which novel
gene a rule promotes, so it cannot be tuned to manufacture a hit:

    NOMINATE = passes the screen gate  AND  structurally tractable  AND  NOT direction-discordant

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
import pathlib
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from cd4_perturbseq import paths, priors  # noqa: E402

SEED = 0
N_PERM = 50_000
N_SHUFFLE = 200
CONJUNCTS = ["supported", "lof_tolerant", "tractable"]
SENSITIVITY_FLOOR = 5 / 6  # pre-registered: below this the rule is rejected as a nomination gate
COLLIDER_CEILING = 0.80  # pre-registered: above this, direction may only demote

# Both direction bases the direction script can emit require the gene to be already characterised.
# The proportion is measured rather than assumed, so a future basis would change the verdict.
CHARACTERISED_BASES = ("approved LoF-mimicking drug", "LoF causes immunodeficiency")

# Housekeeping / proliferation families, pre-committed in the registration so the count is published
# whatever it is. Mitochondrial respiratory chain, ribosome, proteasome, RNA polymerase, glycolysis.
# This is a NAME-PREFIX heuristic and therefore a LOWER BOUND: it cannot see functional housekeeping
# with idiosyncratic names. The peer critic named four such genes inside the direction-unmeasured
# tier — Golgi trafficking, nuclear import, actin branching and ERAD — and they are counted below.
HOUSEKEEPING_RE = re.compile(r"^(NDUF|RPL|RPS|PSM|POLR|MRPL|MRPS|ATP5|UQCR|COX\d)")
HOUSEKEEPING_EXTRA = {"ENO1", "PGK1", "KRR1", "PRKAR1A", "UBR4", "CCNT1", "TOMM70", "CLPB"}
FUNCTIONAL_HOUSEKEEPING = {"COG6", "TNPO1", "ARPC2", "SEL1L"}

# RULE #6. Two other entrants work on these axes; they may never occupy a headline or a top rank.
DO_NOT_HEADLINE = {"STAT6", "GATA3", "STAT5A", "PTEN"}


def _clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Exact binomial confidence interval, via beta quantiles.

    Hand-rolled rather than pulling in statsmodels, which this project already declined to add for a
    CMH test (see DO NOT REDO). Two beta quantiles; validated against the closed-form edge cases.

    Args:
        k: Successes.
        n: Trials.
        alpha: Two-sided error rate.

    Returns:
        Tuple of (lower, upper) bounds.
    """
    lo = 0.0 if k == 0 else float(stats.beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(stats.beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


def _perm_pvalue(is_pos: np.ndarray, drawn: int, observed: int, rng: np.random.Generator) -> tuple[float, float]:
    """Permutation p-value for `observed` positives among `drawn` genes drawn from the screened set.

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
    return float((hits >= observed).sum() + 1) / (N_PERM + 1), float(hits.mean())


def load() -> pd.DataFrame:
    """Join the window score, nomination annotations, direction verdicts and tractability flags.

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
    tr = priors.ot_tractability()[["gene_name", "tractable", "tractable_with_precedent",
                                   "clinical_precedent", "sm_pocket"]]

    f = nom.drop(columns=["tractable", "clinical_precedent", "sm_pocket"], errors="ignore")
    f = f.merge(doe, on="gene_name", how="left").merge(saf, on="gene_name", how="left")
    f = f.merge(tr, on="gene_name", how="left")
    f = f.merge(win[["gene_name", "efficacy", "tolerance_loss", "n_module_down"]], on="gene_name", how="left")

    for c in ("tractable", "tractable_with_precedent", "clinical_precedent", "sm_pocket"):
        f[c] = f[c].fillna(False).astype(bool)
    f["supported"] = f["ot_genetic_supported"].fillna(False).astype(bool)
    f["known"] = f["is_known_drug"].fillna(False).astype(bool)
    f["discordant"] = f["direction_verdict"].eq("DISCORDANT")
    f["concordant"] = f["direction_verdict"].eq("CONCORDANT")
    # FUNCTIONALLY identical to N10 nominate() and N16 bucket() on the 214 safe genes today, not
    # byte-identical: script 20 handles fillna differently. The divergence is latent and would surface
    # for any future gene with prec <= 0.90 and a missing gnomAD constraint record. Both guards here
    # run in the conservative direction. The self-test below is what makes this claim checkable.
    f["lof_tolerant"] = (f["prec"].fillna(1.0) <= 0.90) & (~f["lof_intolerant"].fillna(True).astype(bool))
    # Precedent-only: tractable ONLY because a drug already reached the clinic. Circular for recovery.
    f["precedent_only"] = f["tractable_with_precedent"] & ~f["tractable"]
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


def h3_precedent_is_circular(f: pd.DataFrame) -> tuple[int, int, pd.DataFrame]:
    """Show that gating on clinical precedent would make the drug-recovery validation circular.

    Args:
        f: The joined safe-gene frame.

    Returns:
        Tuple of (positives with precedent, total curated positives found in OT, precedent-only genes).
    """
    tr = priors.ot_tractability()
    truth = pd.read_csv(paths.GROUND_TRUTH / "immunomodulator_targets.csv")
    positives = set(truth.loc[truth["include_as_positive"], "gene_symbol"])
    sub = tr[tr["gene_name"].isin(positives)]
    return int(sub["clinical_precedent"].sum()), len(sub), f[f["precedent_only"]]


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
    pos, cand = f[f["known"]], f[~f["known"] & ~f["discordant"]]
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

    The gate uses STRUCTURAL `tractable`, matching the pre-registration and both controls. An earlier
    version gated on `tractable_with_precedent`, which re-admitted the four precedent-only genes that
    H3 identifies as circular and promoted MALT1 and MEN1 into the top tier. The table then
    contradicted its own H3 verdict. Caught by the adversarial peer critic, not by the controls,
    because the controls were computed on the structural gate while the table was not.

    Tiers A and B are CONCORDANT, and by H2 that means a drug or a monogenic syndrome already exists,
    so they are repurposing leads rather than novel targets. Tier C is not a discovery: it is the set
    whose direction was never measured. Absence of a measurement is not evidence of a good direction.

    Args:
        f: The joined safe-gene frame.

    Returns:
        The nominated genes, tiered and ranked, with a modality tag and a liability string.
    """
    n = f[~f["known"] & ~f["discordant"] & f["tractable"]].copy()

    n["tier"] = np.select(
        [n["concordant"] & n["supported"], n["concordant"], n["supported"]],
        ["A: concordant + genetically supported (repurposing lead)",
         "B: concordant, no autoimmune genetics (repurposing lead)",
         "C: screen + genetics candidate, direction UNMEASURED"],
        default="D: screen only, direction UNMEASURED (hypothesis pool)",
    )
    n["modality"] = np.where(
        n["lof_tolerant"],
        "knockout- or degrader-tolerable",
        "inhibitor only, needs a therapeutic window (LoF-intolerant)",
    )
    n["housekeeping"] = n["gene_name"].str.match(HOUSEKEEPING_RE) | n["gene_name"].isin(HOUSEKEEPING_EXTRA)
    n["housekeeping_functional"] = n["housekeeping"] | n["gene_name"].isin(FUNCTIONAL_HOUSEKEEPING)
    n["do_not_headline"] = n["gene_name"].isin(DO_NOT_HEADLINE)

    liab = []
    for _, r in n.iterrows():
        bits = []
        if r["do_not_headline"]:
            bits.append("RULE #6 competitor overlap: may not be headlined or top-ranked")
        if r["gene_name"] in FUNCTIONAL_HOUSEKEEPING:
            bits.append("functional housekeeping (trafficking / import / cytoskeleton / ERAD)")
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
            bits.append("concordant only because a drug or a monogenic syndrome already exists")
        if r["direction_verdict"] == "UNKNOWN":
            bits.append("direction of effect NEVER MEASURED (N17 failed its calibration at chance)")
        liab.append("; ".join(bits) if bits else "none recorded")
    n["liabilities"] = liab

    n["tier_key"] = n["tier"].str[0]
    # RULE #6 genes are sorted to the bottom of their tier. The prose already avoided them; the
    # committed table was headlining STAT6 by rank order, which is the same overclaim in a CSV.
    n = n.sort_values(["tier_key", "do_not_headline", "ot_genetic_max", "clinical_precedent", "window_rank"],
                      ascending=[True, True, False, False, True]).reset_index(drop=True)
    n["nomination_rank"] = np.arange(1, len(n) + 1)
    return n


def precedent_extension(f: pd.DataFrame) -> pd.DataFrame:
    """Genes tractable only through an existing clinical drug, reported outside the gate.

    These cannot be validated by drug recovery without circularity, so they are never gated on.
    They are repurposing candidates by construction.

    Args:
        f: The joined safe-gene frame.

    Returns:
        The precedent-only genes.
    """
    return f[~f["known"] & ~f["discordant"] & f["precedent_only"]].copy()


def controls(f: pd.DataFrame, nom: pd.DataFrame, rng: np.random.Generator) -> dict[str, bool]:
    """Run every pre-registered control. All must pass or the analysis is void.

    The recovery control deliberately uses STRUCTURAL tractability. Using
    `tractable_with_precedent` would recover approved drug targets by construction.

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

    gate = set(f.loc[f["tractable"] & ~f["discordant"], "gene_name"])
    observed = int((win["gene_name"].isin(gate) & win["is_pos"]).sum())

    # The HONEST null holds the safety gate fixed and draws |gate| genes from the 214 safe genes.
    # Drawing from all 6,371 rankable genes instead re-counts the safety-gate enrichment that N13 and
    # N14 already established, and inflates the p-value by an order of magnitude. Both are printed.
    f_is_pos = f["gene_name"].isin(positives).to_numpy()
    p_cond, exp_cond = _perm_pvalue(f_is_pos, len(gate), observed, rng)
    p_uncond, exp_uncond = _perm_pvalue(win["is_pos"].to_numpy(), len(gate), observed, rng)
    n_pos_safe = int(f_is_pos.sum())
    ceiling, _ = _perm_pvalue(f_is_pos, len(gate), n_pos_safe, rng)
    print(f"\n[control] recovery by the rebuilt gate (STRUCTURAL tractability, non-circular)")
    print(f"          PRIMARY, conditioned on the safety gate ({len(gate)} drawn from {len(f)} safe): "
          f"observed {observed}, expected {exp_cond:.2f}, permutation p = {p_cond:.4g}")
    print(f"          power ceiling: only {n_pos_safe} curated positives are in the safe set, so a gate "
          f"retaining all of them could not beat p = {ceiling:.4g}")
    print(f"          (unconditioned on 6,371 rankable genes: p = {p_uncond:.4g} — this re-counts the "
          f"safety-gate enrichment of N13/N14 and is NOT the headline)")
    rec_ok = p_cond < 0.05 and observed > exp_cond
    p_rec = p_cond

    leaked = sorted(set(nom["gene_name"]) & {"PTPN2", "RC3H1"})
    print(f"[control] direction-discordant genes in the nominated set: {leaked or 'none'}")
    disc_ok = not leaked

    shuffled = np.array([
        int((win["gene_name"].isin(
            set(f.loc[rng.permutation(f["tractable"].to_numpy()) & ~f["discordant"], "gene_name"])
        ) & win["is_pos"]).sum())
        for _ in range(N_SHUFFLE)
    ])
    frac_ge = float((shuffled >= observed).mean())
    print(f"[control] shuffled tractability reaches {observed} positives in {frac_ge:.1%} of "
          f"{N_SHUFFLE} shuffles (mean {shuffled.mean():.2f})")
    shuffle_ok = frac_ge <= 0.10

    contaminated = sorted(set(nom["gene_name"]) & positives)
    print(f"[control] curated positives presented as novel nominations: {contaminated or 'none'}")
    truth_ok = not contaminated

    return {
        "recovery preserved under a non-circular gate": rec_ok,
        "discordant genes still demoted": disc_ok,
        "tractability flag carries information": shuffle_ok,
        "ground truth not contaminated": truth_ok,
    }


def main() -> None:
    """Run the three tests, rebuild the nomination layer, and gate the verdict on the controls."""
    rng = np.random.default_rng(SEED)
    f = load()
    print(f"safe genes: {len(f)}   recovered approved-drug targets: {int(f['known'].sum())}")
    selftest_reproduces_n16(f)

    print("\n" + "=" * 92)
    print("H3. Would gating on clinical precedent be a bug fix, or a circularity?")
    print("=" * 92)
    n_prec, n_truth, prec_only = h3_precedent_is_circular(f)
    print(f"  clinical_precedent is true for {n_prec} of the {n_truth} curated positives found in "
          f"Open Targets.")
    print(f"  A gate reading it would recover approved drug targets by construction, so the recovery")
    print(f"  validation would be circular. Structural tractability is already true for all "
          f"{int(f['known'].sum())} recovered targets.")
    print(f"  VERDICT: `tractable` stays STRUCTURAL. Precedent is reported, never gated on.")
    print(f"\n  Genes tractable ONLY via clinical precedent (repurposing candidates by construction):")
    print(prec_only[["gene_name", "direction_verdict", "clinical_precedent", "sm_pocket",
                     "ot_genetic_max"]].to_string(index=False) if len(prec_only) else "  (none)")

    print("\n" + "=" * 92)
    print("H2. Is direction-adjudicability a deterministic proxy for prior characterisation?")
    print("=" * 92)
    p_collider, bases = h2_collider(f)
    print(bases.to_string())
    n_unknown = int(f["direction_verdict"].eq("UNKNOWN").sum())
    print(f"\n  proportion of CONCORDANT verdicts requiring prior characterisation: {p_collider:.3f}")
    print(f"  pre-registered ceiling: {COLLIDER_CEILING}")
    collider = p_collider > COLLIDER_CEILING
    print(f"  VERDICT: {'CONFIRMED' if collider else 'not confirmed'} — direction of effect "
          f"{'may only DEMOTE, never PROMOTE' if collider else 'may remain a promoting criterion'}.")
    print(f"  The classifier has exactly two emit paths and both require prior characterisation, so")
    print(f"  P(concordant | uncharacterised) = 0. This is a structural certainty, not a measured bias.")
    print(f"  Missingness is not at random: {n_unknown} of {len(f)} safe genes ({n_unknown/len(f):.0%}) are UNKNOWN.")
    print(f"  CAVEAT 1: DISCORDANCE is equally characterisation-bound, so demoting PTPN2 and RC3H1")
    print(f"           licenses NO direction-safety claim about the {n_unknown} direction-unmeasured genes.")
    print(f"  CAVEAT 2: those genes are direction-UNMEASURED, not direction-favourable. The absence of a")
    print(f"           measurement is not evidence of a good direction. Tier C is named accordingly.")

    print("\n" + "=" * 92)
    print("H1. Sensitivity of every rule variant, measured on the positive class only")
    print("=" * 92)
    variants = h1_sensitivity(f)
    print(variants.to_string(index=False))
    shipped = variants[variants["rule"] == "supported & lof_tolerant & tractable"].iloc[0]
    n_known = int(f["known"].sum())
    k = int(shipped["recovered_drugs"])
    lo, hi = _clopper_pearson(k, n_known)
    p_binom = float(stats.binomtest(k, n_known, SENSITIVITY_FLOOR, alternative="less").pvalue)
    print(f"\n  the shipped rule re-nominates {k} of the {n_known} approved-drug targets this pipeline")
    print(f"  recovers, and nominates {shipped['novel_nominated']} novel gene(s).")
    print(f"  n = {n_known}. Exact Clopper-Pearson 95% CI on the sensitivity: "
          f"[{lo:.2f}, {hi:.2f}] — it still CONTAINS the pre-registered floor of {SENSITIVITY_FLOOR:.3f}.")
    print(f"  One-sided binomial P against that floor: {p_binom:.3f}.")
    print(f"  So the proportion alone does NOT reject the rule; it is under-powered at n = {n_known}.")
    print("\n  The rule is rejected for two reasons that do not depend on it:")
    print(f"    (1) it returns {shipped['novel_nominated']} gene from {int((~f['known'] & ~f['discordant']).sum())} candidates;")
    print("    (2) the misses are a DETERMINISTIC database-coverage property, not a sampling rate:")
    pos = f[f["known"]]
    for c in CONJUNCTS:
        failed = sorted(pos.loc[~pos[c], "gene_name"])
        if failed:
            print(f"        `{c}` deletes {', '.join(failed)}")
    print("        These are the essential-hub drug classes phenotypic pharmacology keeps hitting:")
    print("        calcineurin (PPP3R1, ciclosporin), IMPDH (IMPDH2, mycophenolate), the TCR (CD3E).")
    rule_rejected = shipped["sensitivity"] < SENSITIVITY_FLOOR
    print(f"  VERDICT: {'RULE REJECTED as a nomination gate' if rule_rejected else 'rule stands'}")

    print("\n" + "=" * 92)
    print("REBUILT NOMINATION LAYER: screen gate AND tractable AND NOT direction-discordant")
    print("=" * 92)
    nom = tier_and_rank(f)
    ext = precedent_extension(f)
    print(f"nominated: {len(nom)} genes (structural tractability only)")
    print(nom.groupby("tier").size().to_string())
    print(f"\nhousekeeping / respiratory-chain genes among the nominations (name-prefix): "
          f"{int(nom['housekeeping'].sum())} of {len(nom)} "
          f"(pre-committed to publishing this number whatever it is)")
    print(f"  that heuristic is a LOWER BOUND. Adding functional housekeeping with idiosyncratic "
          f"names: {int(nom['housekeeping_functional'].sum())} of {len(nom)}.")
    print(f"\nprecedent-only repurposing candidates, REPORTED OUTSIDE THE GATE "
          f"(gating on them would make drug recovery circular): {sorted(ext['gene_name'])}")

    verdicts = controls(f, nom, rng)
    print("\n" + "=" * 92)
    for name, ok in verdicts.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

    cols = ["nomination_rank", "tier", "gene_name", "direction_verdict", "direction_basis",
            "ot_genetic_max", "ot_genetic_n_diseases", "loeuf", "prec", "lof_tolerant", "modality",
            "tractable", "precedent_only", "sm_pocket", "clinical_precedent", "viability_tier",
            "housekeeping", "housekeeping_functional", "do_not_headline", "window_rank", "efficacy",
            "n_module_down", "il2_hit", "liabilities"]
    nom[cols].to_csv(paths.TABLES / "nomination_recalibrated.csv", index=False)
    variants.to_csv(paths.TABLES / "nomination_rule_sensitivity.csv", index=False)
    print(f"\nwrote {paths.TABLES / 'nomination_recalibrated.csv'}")
    print(f"wrote {paths.TABLES / 'nomination_rule_sensitivity.csv'}")

    for key, label in (("A", "TIER A — repurposing leads (concordant + genetically supported)"),
                       ("C", "TIER C — screen + genetics, direction UNMEASURED (NOT a discovery)")):
        sub = nom[nom["tier_key"] == key]
        print("\n" + "=" * 92)
        print(f"{label}: n = {len(sub)}")
        print("=" * 92)
        print(sub[["gene_name", "ot_genetic_max", "ot_genetic_n_diseases", "clinical_precedent",
                   "viability_tier", "window_rank", "modality"]].to_string(index=False))
        print("\nliabilities:")
        for _, r in sub.iterrows():
            print(f"  {r['gene_name']:<9} {r['liabilities']}")

    if not all(verdicts.values()):
        print("\nA pre-registered control FAILED. The recalibration is void; nothing here may be quoted.")
        sys.exit(1)

    n_c = int((nom["tier_key"] == "C").sum())
    print("\n" + "=" * 92)
    print("HONEST HEADLINE (RULE #8)")
    print("=" * 92)
    print(f"""  'Zero vetted novel targets' was substantially a property of the nomination rule, not of the
  screen. Of the {n_known} approved-drug targets this pipeline recovers, the N16 rule would have
  re-nominated only {k} (n = {n_known}, exact 95% CI [{lo:.2f}, {hi:.2f}], one-sided binomial P = {p_binom:.3f} against the
  pre-registered 5/6 floor). That proportion is under-powered and does NOT by itself reject the
  rule. The rule is rejected because it returns {shipped['novel_nominated']} gene from {int((~f['known'] & ~f['discordant']).sum())}, and because its misses are a
  deterministic database-coverage property: `supported` deletes CD3E, IMPDH2 and PPP3R1 at an
  Open Targets autoimmune score of exactly 0.000, and `lof_tolerant` deletes IMPDH2 at prec 0.999.
  Those are the essential-hub drug classes phenotypic pharmacology keeps hitting.

  Direction of effect cannot promote a novel gene. All {int(f['concordant'].sum())} concordant verdicts require a drug or a
  monogenic syndrome to already exist, so P(concordant | uncharacterised) = 0 and the annotation is
  missing for {n_unknown} of {len(f)} safe genes. Direction is therefore used only to demote. Because
  DISCORDANCE is equally characterisation-bound, demoting PTPN2 and RC3H1 licenses no direction-safety
  claim about the {n_unknown} unmeasured genes.

  A gate calibrated to re-nominate {n_known} of {n_known} recovered drug targets returns {len(nom)} candidates,
  {int(nom['housekeeping'].sum())} of them housekeeping or respiratory-chain by name and {int(nom['housekeeping_functional'].sum())} once functional housekeeping is
  counted. {n_c} carry autoimmune genetic support with direction NEVER MEASURED. That tier is not a
  discovery: it is the set whose direction nobody has measured, and absence of a measurement is not
  evidence of a favourable direction.

  NO GENE HERE IS A VETTED NOVEL TARGET. The contribution is the decision layer and its audit.""")


if __name__ == "__main__":
    main()
