"""N8: how many approved immunomodulator targets could this benchmark EVER recover?

The drug-target-recovery benchmark is failing, and not because the ranking is bad. Of 36
hand-curated positives, 32 are perturbed, 30 are measured, and only 20 survive perturbation QC.
``JAK1`` and ``JAK3`` were never perturbed. AUROC is 0.542 with a bootstrap 95% CI of
[0.373, 0.707], which covers chance.

Before we torture that number, we should find out whether the positive set is simply too small.
That is a question about the world, not about our ranking, and Open Targets answers it.

    **Count 3 is the number this script exists to produce**: how many approved, immune-indication,
    loss-of-function-mimicking drug targets are among the 6,371 genes our ranking can even see.
    That is the maximum achievable positive-set size. Above roughly 60 the benchmark is powered and
    we rebuild it. Near 20 it cannot be rebuilt, and we retire drug-target recovery as the headline
    and validate the efficacy axis on the held-out Schmidt IL-2 screen instead.

Either answer is publishable. The failure mode is a number that looks fine because nobody checked
its denominator.

**The direction-of-effect constraint.** CRISPRi knockdown removes protein function, so a
perturbation can phenocopy an inhibitor, antagonist, blocker or degrader, and can never phenocopy
an agonist. ``NR3C1`` is the sharpest trap: glucocorticoids are the most prescribed
immunosuppressants in the world and they are *agonists*. ``FKBP1A`` is the second: tacrolimus must
*bind* it to inhibit calcineurin, so knocking it down abolishes the drug rather than mimicking it.
Any positive set that ignores this is silently wrong. Open Targets gives us ``actionType`` and we
gate on it.

**Two judgement calls are exposed rather than buried**, because count 3 depends on both and a single
number would hide that.

1. *Which diseases are immune?* Open Targets has an "immune system disorder" therapeutic area, and
   **asthma is not in it**. Neither is atopic eczema. Dupilumab targets ``IL4R``, is indicated for
   both, and is one of our own curated positives, so a strict therapeutic-area filter would silently
   drop it. We therefore report a STRICT set (the therapeutic area alone) and a BROAD set (adding
   autoimmune, allergic, hypersensitivity, asthma, atopic eczema and inflammatory bowel disease).
2. *Which mechanisms mimic a knockdown?* ``INHIBITOR`` and ``ANTAGONIST`` obviously do.
   ``MODULATOR`` and ``BINDING AGENT`` might. We report a CORE set and a PERMISSIVE set.

3. *Is the annotated target the real target?* Open Targets lists metformin's mechanism against all
   forty subunits of mitochondrial complex I. A benchmark positive must be a gene a drug was
   designed against, not a subunit of a complex it happens to bind.

The result is a 2x2x2 of ceilings rather than one number, and the decision is read off the whole box.
A fourth trap is not a judgement call but an outright error, and it is described on
:func:`approved_immune_targets`: clinical stage must be read per INDICATION, not per drug-target pair.

Version-pinned to the Open Targets release in ``data/external/open_targets/VERSION``. Target-drug
links move between quarterly releases, so a ceiling quoted without a release version is not a number.

RULE #4 is enforced here, not just stated. ``CD2`` (alefacept) ranks second on our shortlist and is
NOT in our curated ground truth, deliberately. If Open Targets returns it independently, that is
evidence. If we had typed it in, it would be nothing. This script reports whether it comes back.

Usage:
    bash scripts/15_fetch_open_targets.sh
    uv run python scripts/16_open_targets_benchmark.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from cd4_perturbseq import de_stats, paths

OT_DIR = paths.EXTERNAL / "open_targets"

APPROVED_STAGE = "APPROVAL"
"""Open Targets encodes clinical stage as a string. There is no ``PHASE_4``; approval is ``APPROVAL``."""

IMMUNE_STRICT = {"MONDO_0005046"}
"""The "immune system disorder" therapeutic area. 1,189 diseases. Excludes asthma."""

IMMUNE_BROAD = IMMUNE_STRICT | {
    "MONDO_0007179",  # autoimmune disease
    "MONDO_0005271",  # allergic disease
    "MONDO_0000605",  # hypersensitivity reaction disease
    "MONDO_0004979",  # asthma
    "MONDO_0004980",  # atopic eczema
    "MONDO_0005265",  # inflammatory bowel disease
}
"""Adds the inflammatory and allergic indications that sit outside the immune therapeutic area.
Without asthma and atopic eczema, dupilumab's target ``IL4R`` is not an immune drug target."""

CANCER_ROOT = {"MONDO_0045024"}
"""The "cancer or benign tumor" therapeutic area, 3,661 diseases, SUBTRACTED from every immune set.

This is not a nicety. MONDO classifies multiple myeloma, acute myeloid leukaemia, diffuse large
B-cell lymphoma and chronic lymphocytic leukaemia as descendants of "immune system disorder",
because they are malignancies of immune cells. Without this subtraction, bortezomib counts as an
approved immune drug and the proteasome subunit ``PSMB5`` becomes an immunomodulator target, as do
the histone deacetylases (vorinostat) and beta-tubulin (vincristine). The first run of this script
returned 173 "rankable approved immune targets" and most of them were cytotoxic chemotherapy."""

MAX_TARGETS_PER_MECHANISM = 3
"""Reject mechanism rows that name more than this many targets.

Open Targets annotates metformin's mechanism against every subunit of mitochondrial complex I, and
vincristine's against the tubulin family. Those rows make ``NDUFA4`` an approved drug target for
type 1 diabetes. A benchmark positive must be a gene a drug was designed against, not a subunit of
a complex the drug happens to bind."""

LOF_CORE = frozenset({
    "INHIBITOR", "ANTAGONIST", "BLOCKER", "DEGRADER", "ANTISENSE INHIBITOR",
    "ALLOSTERIC ANTAGONIST", "NEGATIVE MODULATOR", "NEGATIVE ALLOSTERIC MODULATOR",
})
"""Mechanisms a knockdown unambiguously phenocopies."""

LOF_PERMISSIVE = LOF_CORE | {"INVERSE AGONIST", "DISRUPTING AGENT", "MODULATOR", "BINDING AGENT"}
"""An upper bound. ``MODULATOR`` and ``BINDING AGENT`` are ambiguous in Open Targets and may be
either direction, so including them can only inflate the ceiling, never deflate it."""

POWERED_AT = 60
"""Count 3 at or above this and the benchmark is worth rebuilding. Fixed before the count was run.

The decision is taken on the UPPER BOUND: broadest disease set, most permissive mechanism set, no
specificity filter. If even that cannot reach 60, no defensible definition can."""


def _aslist(value) -> list:
    """Coerce a parquet list column to a Python list.

    Arrow list columns arrive as numpy arrays, and ``array or []`` raises rather than falling back,
    because an array's truth value is ambiguous.

    Args:
        value: A cell from a list-typed column, possibly None or NaN.

    Returns:
        A list, empty when the cell holds nothing.
    """
    if value is None:
        return []
    try:
        return list(value)
    except TypeError:
        return []


def _read_parts(pattern: str) -> pd.DataFrame:
    """Read every parquet part matching a glob under the Open Targets directory.

    Args:
        pattern: Glob relative to ``OT_DIR``.

    Returns:
        The concatenated frame.

    Raises:
        FileNotFoundError: If nothing matches. Run ``scripts/15_fetch_open_targets.sh``.
    """
    files = sorted(OT_DIR.glob(pattern))
    if not files:
        raise FileNotFoundError(f"no {pattern} under {OT_DIR}; run bash scripts/15_fetch_open_targets.sh")
    return pd.concat([pq.read_table(f).to_pandas() for f in files], ignore_index=True)


def _under(disease: pd.DataFrame, roots: set[str]) -> pd.Series:
    """Boolean mask of diseases that are, or descend from, any of the roots.

    Args:
        disease: The Open Targets disease frame.
        roots: Ontology ids.

    Returns:
        Boolean series aligned to ``disease``.
    """
    return disease.apply(
        lambda row: bool(roots & (set(_aslist(row["ancestors"])) | {row["id"]} | set(_aslist(row["therapeuticAreas"])))),
        axis=1,
    )


def immune_disease_ids(disease: pd.DataFrame, roots: set[str], drop_cancer: bool = True) -> set[str]:
    """Immune diseases, with malignancies of immune cells removed.

    MONDO places multiple myeloma, AML, DLBCL and CLL under "immune system disorder". They are
    immune-cell cancers, and the drugs approved for them are cytotoxic chemotherapy, not
    immunomodulators. Leaving them in turns the proteasome and beta-tubulin into immunomodulator
    targets.

    Args:
        disease: The Open Targets disease frame.
        roots: Ontology ids treated as immune.
        drop_cancer: Subtract everything under the cancer therapeutic area.

    Returns:
        The set of matching disease ids.
    """
    keep = _under(disease, roots)
    if drop_cancer:
        keep &= ~_under(disease, CANCER_ROOT)
    return set(disease.loc[keep, "id"])


def lof_drug_target_pairs(
    moa: pd.DataFrame, action_types: frozenset[str], max_targets: int | None
) -> set[tuple[str, str]]:
    """Drug-target pairs whose mechanism a loss-of-function perturbation can phenocopy.

    Args:
        moa: The drug mechanism-of-action frame.
        action_types: Accepted ``actionType`` values.
        max_targets: Reject mechanism rows naming more targets than this. None disables the filter,
            which readmits complex-I subunits and the tubulin family.

    Returns:
        Set of ``(drugId, targetId)`` tuples.
    """
    subset = moa[moa["actionType"].isin(action_types)]
    pairs: set[tuple[str, str]] = set()
    for _, row in subset.iterrows():
        targets = _aslist(row["targets"])
        if max_targets is not None and len(targets) > max_targets:
            continue
        for drug in _aslist(row["chemblIds"]):
            for target in targets:
                pairs.add((str(drug), str(target)))
    return pairs


def approved_immune_targets(
    indication: pd.DataFrame, immune: set[str], lof_pairs: set[tuple[str, str]]
) -> pd.DataFrame:
    """Targets with at least one drug APPROVED FOR an immune indication, LoF-mimicking.

    The stage must be linked to the indication. ``clinical_target.maxClinicalStage`` is the maximum
    stage a drug-target pair reached across *all* indications, so joining it to an immune indication
    conjoins an approval earned somewhere else with a trial that never finished. Vorinostat is
    approved for cutaneous T-cell lymphoma and merely *trialled* in graft-versus-host disease and
    Crohn's; read from ``clinical_target`` it looks like an approved immunomodulator, and ``HDAC1``
    becomes a benchmark positive. ``clinical_indication`` carries the stage per indication, so the
    approval and the disease are the same fact.

    Args:
        indication: The ``clinical_indication`` frame: one row per drug-disease pair.
        immune: Immune disease ids, cancer already subtracted.
        lof_pairs: Accepted ``(drugId, targetId)`` tuples.

    Returns:
        One row per surviving ``(targetId, drugId, diseaseId)``.
    """
    approved = indication[
        (indication["maxClinicalStage"] == APPROVED_STAGE) & (indication["diseaseId"].isin(immune))
    ]
    targets_of: dict[str, set[str]] = {}
    for drug, target in lof_pairs:
        targets_of.setdefault(drug, set()).add(target)

    rows = []
    for _, row in approved.iterrows():
        for target in targets_of.get(str(row["drugId"]), ()):
            rows.append({"targetId": target, "drugId": row["drugId"], "diseaseId": row["diseaseId"]})
    return pd.DataFrame(rows, columns=["targetId", "drugId", "diseaseId"])


def ceiling(frame: pd.DataFrame, symbol_of: dict[str, str], measured: set[str], rankable: set[str]) -> dict:
    """The three counts the benchmark decision turns on.

    Args:
        frame: Output of :func:`approved_immune_targets`.
        symbol_of: Ensembl id to approved gene symbol, over ALL human targets.
        measured: The 10,282 measured transcriptome genes.
        rankable: The 6,371 genes surviving perturbation QC.

    Returns:
        Mapping with counts 1, 2 and 3, and the rankable symbols.
    """
    targets = set(frame["targetId"])
    symbols = {symbol_of[t] for t in targets if t in symbol_of}
    return {
        "count1_targets": len(targets),
        "count2_measured": len(symbols & measured),
        "count3_rankable": len(symbols & rankable),
        "symbols": symbols,
        "rankable_symbols": sorted(symbols & rankable),
    }


VALIDATION = {
    "IMPDH2": "mycophenolate, an approved immunosuppressant, inhibits it",
    "IL17A": "secukinumab, approved for psoriasis, neutralises it",
    "IL12B": "ustekinumab, approved for psoriasis, neutralises it",
    "DHODH": "teriflunomide, approved for multiple sclerosis, inhibits it",
    "TNF": "adalimumab, approved for rheumatoid arthritis, neutralises it",
}
"""Textbook immunomodulator targets. If the filters cannot recover these, the filters are wrong."""

ANTI_VALIDATION = {
    "PSMB5": "bortezomib is chemotherapy for myeloma, not an immunomodulator",
    "TUBB4B": "vincristine is chemotherapy; its target is the tubulin family",
    "HDAC1": "vorinostat is chemotherapy for lymphoma",
    "NDUFA4": "metformin is annotated against all 40 subunits of complex I",
}
"""Things the first, broken run of this script called approved immunomodulator targets."""


def validate(symbols: set[str], symbol_of: dict[str, str]) -> bool:
    """Check the filters against things we already know, before trusting the count.

    Args:
        symbols: Gene symbols the filters returned as approved immune LoF targets.
        symbol_of: The full Ensembl to symbol map, used to confirm the gene exists at all.

    Returns:
        True if every positive control is recovered and every negative control is excluded.
    """
    known = set(symbol_of.values())
    print("\n=== VALIDATION: can these filters recover drugs we already know about? ===")
    ok = True
    for gene, why in VALIDATION.items():
        found = gene in symbols
        exists = gene in known
        ok &= found or not exists
        state = "recovered" if found else ("MISSING" if exists else "absent from Open Targets entirely")
        print(f"  +  {gene:8s} {state:34s} {why}")
    print("  (a positive control that Open Targets has never heard of cannot fail this test)")
    for gene, why in ANTI_VALIDATION.items():
        found = gene in symbols
        ok &= not found
        print(f"  -  {gene:8s} {'STILL PRESENT, filter failed' if found else 'correctly excluded':34s} {why}")
    print(f"\n  filters {'PASS' if ok else 'FAIL'} their own controls")
    return ok


def main() -> None:
    """Bound the benchmark, validate the filters, audit our curation, enforce RULE #4."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--powered-at", type=int, default=POWERED_AT)
    args = parser.parse_args()
    paths.ensure_dirs()

    version_file = OT_DIR / "VERSION"
    if not version_file.exists():
        raise FileNotFoundError(f"{version_file} missing; run bash scripts/15_fetch_open_targets.sh")
    version = version_file.read_text().strip()
    print(f"Open Targets release {version}")

    indication = _read_parts("clinical_indication.parquet")
    disease = _read_parts("disease.parquet")
    moa = _read_parts("moa_*.parquet")
    targets_tbl = pd.concat(
        [pq.read_table(f, columns=["id", "approvedSymbol"]).to_pandas() for f in sorted(OT_DIR.glob("target_*.parquet"))],
        ignore_index=True,
    )
    symbol_of = dict(zip(targets_tbl["id"].astype(str), targets_tbl["approvedSymbol"].astype(str), strict=True))
    print(f"  clinical_indication {len(indication):,}   disease {len(disease):,}   mechanisms {len(moa):,}   "
          f"targets {len(symbol_of):,}")

    measured = set((paths.HANDOFF_INPUTS / "measured_genes.txt").read_text().split())
    rankable = set((paths.HANDOFF_INPUTS / "rankable_genes.txt").read_text().split())

    # ---------------------------------------------------------------- the grid of ceilings
    print("\n=== THE CEILING. How many approved immune LoF targets can this ranking even see? ===")
    print(f"    'Powered' was fixed at count 3 >= {args.powered_at} before any count was computed, and is")
    print("    judged on the UPPER BOUND: broadest diseases, most permissive mechanisms, no specificity filter.")
    print("    Cancer is subtracted from every immune set. MONDO calls myeloma an immune system disorder.")

    grid = []
    primary = upper = None
    for d_label, roots in (("strict (immune TA)", IMMUNE_STRICT), ("broad (+asthma, eczema, IBD)", IMMUNE_BROAD)):
        immune = immune_disease_ids(disease, roots, drop_cancer=True)
        for m_label, actions in (("core", LOF_CORE), ("permissive", LOF_PERMISSIVE)):
            for s_label, max_t in ((f"<={MAX_TARGETS_PER_MECHANISM} targets", MAX_TARGETS_PER_MECHANISM),
                                   ("unfiltered", None)):
                pairs = lof_drug_target_pairs(moa, actions, max_t)
                frame = approved_immune_targets(indication, immune, pairs)
                counts = ceiling(frame, symbol_of, measured, rankable)
                grid.append({"disease_set": d_label, "n_immune_diseases": len(immune), "moa_set": m_label,
                             "specificity": s_label,
                             **{k: v for k, v in counts.items() if k not in ("symbols", "rankable_symbols")}})
                if d_label.startswith("broad") and m_label == "core" and max_t is not None:
                    primary = (frame, counts)
                if d_label.startswith("broad") and m_label == "permissive" and max_t is None:
                    upper = (frame, counts)
    table = pd.DataFrame(grid)
    print()
    print(table.to_string(index=False))

    frame_p, counts_p = primary
    frame_u, counts_u = upper
    print(f"\n    PRIMARY  (broad diseases, core mechanisms, specificity filter): count 3 = {counts_p['count3_rankable']}")
    print(f"    UPPER    (every judgement call made in the benchmark's favour):  count 3 = {counts_u['count3_rankable']}")

    filters_ok = validate(counts_p["symbols"], symbol_of)

    # ---------------------------------------------------------------- the decision
    count3 = counts_u["count3_rankable"]
    print("\n" + "=" * 92)
    if not filters_ok:
        print("VERDICT: WITHHELD. The filters fail their own positive and negative controls, so the")
        print("ceiling below is not trustworthy. Fix the filters before quoting any number from this run.")
    elif count3 >= args.powered_at:
        print(f"VERDICT: count 3 = {count3} >= {args.powered_at} at the upper bound. The benchmark is POWERED.")
        print(f"The defensible PRIMARY count is {counts_p['count3_rankable']}. Rebuild on the union set and")
        print("report AUROC against both the hand-curated and the Open Targets positive sets, so the reader")
        print("can see how much the answer depends on the gold standard.")
    else:
        print(f"VERDICT: count 3 = {count3} < {args.powered_at}, even with every call made in its favour.")
        print("Drug-target recovery CANNOT be rebuilt into a powered benchmark on this dataset.")
        print("Retire it as the headline. Report it honestly as an underpowered check, and validate the")
        print("efficacy axis on the held-out Schmidt IL-2 screen (AUROC 0.702, CI [0.591, 0.814], n=33).")
        print("The ceiling is a property of the perturbation library, not of our ranking.")
    print("=" * 92)

    # ---------------------------------------------------------------- audit our curation
    truth = pd.read_csv(paths.GROUND_TRUTH / "immunomodulator_targets.csv")
    ot_symbols = counts_p["symbols"]
    all_ot_genes = set(symbol_of.values())

    def state(gene: str) -> str:
        if gene in ot_symbols:
            return "OT_CONFIRMS"
        if gene not in all_ot_genes:
            return "NOT_IN_OPEN_TARGETS"
        if gene not in measured:
            return "OT_SILENT_not_measured_here"
        return "OT_SILENT"

    truth["ot_state"] = truth["gene_symbol"].map(state)
    truth["ot_rankable"] = truth["gene_symbol"].isin(counts_p["rankable_symbols"])
    truth["disagrees_with_curation"] = truth["include_as_positive"] & (truth["ot_state"] != "OT_CONFIRMS")

    print("\n=== Task A: Open Targets versus our 46 hand-curated rows ===")
    print("    OT_SILENT_not_measured_here is NOT a disagreement. It means our transcriptome cannot")
    print("    see the gene, so no ranking of ours could ever have recovered it.")
    print(truth.groupby(["include_as_positive", "ot_state"]).size().to_string())
    flagged = truth[truth["disagrees_with_curation"]]
    print(f"\n  positives Open Targets does not confirm: {len(flagged)}")
    print(flagged[["gene_symbol", "ot_state", "mechanism"]].to_string(index=False))

    print("\n  the direction-of-effect traps we excluded by hand:")
    for gene in ("NR3C1", "IL2", "FKBP1A", "PPIA"):
        row = truth[truth["gene_symbol"] == gene]
        if len(row):
            r = row.iloc[0]
            confirms = r["ot_state"] == "OT_CONFIRMS"
            print(f"    {gene:8s} we excluded it; OT {'RETURNS it -> disagreement' if confirms else 'does not return it -> agree'}")

    # ---------------------------------------------------------------- Task B
    curated = set(truth.loc[truth["include_as_positive"], "gene_symbol"])
    missed = sorted(set(counts_p["rankable_symbols"]) - curated)
    print(f"\n=== Task B: rankable approved immune LoF targets we did NOT curate: {len(missed)} ===")
    print(f"  {', '.join(missed) if missed else 'none'}")

    # ---------------------------------------------------------------- RULE #4
    print("\n=== RULE #4: CD2 (alefacept) ranks Tier A #2 and was deliberately NOT added by hand ===")
    if "CD2" in counts_p["rankable_symbols"]:
        print("  Open Targets returns CD2 INDEPENDENTLY, under the primary definition. That is evidence.")
        print("  It may now enter the positive set, sourced to Open Targets rather than to us watching it rank.")
    elif "CD2" in ot_symbols:
        print("  Open Targets knows CD2 as an approved immune LoF target, but it is not rankable here.")
    else:
        print("  Open Targets does NOT return CD2 under the primary definition. Do not add it. RULE #4 holds.")

    # ---------------------------------------------------------------- Task C
    print("\n=== Task C: the JAK family, and what the library lost ===")
    for gene in ("JAK1", "JAK2", "JAK3", "TYK2"):
        print(f"  {gene:6s} state={state(gene):28s} rankable_here={gene in rankable}")

    # ---------------------------------------------------------------- persist
    out = paths.TABLES / "open_targets_benchmark_ceiling.csv"
    table.insert(0, "ot_release", version)
    table["filters_pass_controls"] = filters_ok
    table.to_csv(out, index=False)
    print(f"\nwrote {out}")

    audit_out = paths.GROUND_TRUTH / "immunomodulator_targets.open_targets.csv"
    truth.to_csv(audit_out, index=False)
    print(f"wrote {audit_out}")

    extra = frame_p.copy()
    extra["gene_symbol"] = extra["targetId"].map(symbol_of)
    extra["in_measured"] = extra["gene_symbol"].isin(measured)
    extra["in_rankable"] = extra["gene_symbol"].isin(rankable)
    extra["already_curated"] = extra["gene_symbol"].isin(curated)
    extra_out = paths.GROUND_TRUTH / "open_targets_additional_positives.csv"
    extra.sort_values(["in_rankable", "gene_symbol"], ascending=[False, True]).to_csv(extra_out, index=False)
    print(f"wrote {extra_out}")

    raise SystemExit(0 if filters_ok else 1)


if __name__ == "__main__":
    main()
