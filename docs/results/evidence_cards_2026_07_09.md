# N15 — Evidence cards with a fabrication firewall (the Claude-in-the-loop story)

Run under RULE #9. Script `scripts/24_evidence_card_firewall.py`; frozen records
`results/tables/evidence_records.json`.

## 1. Question and method

The project's deliverable format is a sourced evidence card per candidate. An LLM is the natural writer,
but LLMs fabricate numbers and citations at high rates (Walters & Wilder 2023: 55% of GPT-3.5 and 18% of
GPT-4 citations entirely fabricated), so an LLM-written card cannot be trusted on its face. The literature
review (§3) gives the defensible pattern: freeze a per-gene record of only values already computed and
committed; let Claude write prose from that record alone; then a **deterministic checker** rejects the card
if any number in it is not in the record (within rounding) or any citation is not on a curated allow-list.
Claude writes the prose; it never sets a rank or a verdict, which come from the committed pipeline.

## 2. The firewall, and its falsifiable control

`scripts/24_evidence_card_firewall.py` builds the frozen records and implements `check_card`, which flags:
- any number in the card not matching a record value within tolerance (4-digit years are exempt and handled
  by the citation check), and
- any `Surname Year` citation not on the allow-list of the ~20 papers the project actually uses.

**Falsifiable control (required by RULE #9):** a deliberately fabricated card — "LOEUF 0.35, genetic support
0.88 across 9 diseases … Fabricato et al. 2099" against a record with LOEUF 0.767 — MUST be caught. It is:
the checker flags 0.35, 0.88 and the citation Fabricato 2099, while leaving a real citation (Manguso 2017)
untouched. A record-faithful card passes with zero violations. Control **PASS**.

## 3. Demonstrated on Claude-generated cards

Claude wrote a card for each of 7 genes from its frozen record only. All 7 pass the firewall with **zero
violations**. An adversarial variant of the CD2 card, embellished with an ungrounded genetic score (0.94,
the record says 0.68) and a fake citation (Smith 2088), is **caught** on both. This shows the firewall works
on real model output, not just the hardcoded control.

## 4. The cards (validated, every number from the committed tables)

**IMPDH2** — recovered positive control, target of mycophenolate. Knockdown suppresses the module (6 genes
down; efficacy 0.598), passes the gate in the depleting-at-rest antiproliferative tier. Direction concordant
(approved LoF-mimicking drug). No autoimmune genetic association (0.0). Safety flags recessive constraint
(prec 0.999) and broad expression (max nTPM 296.5) — mycophenolate needs monitoring. Recovery validates the
gate (Frangieh 2021; Schmidt 2022).

**PPP3R1** — calcineurin regulatory subunit, target of ciclosporin/tacrolimus. Passes the gate (3 genes
down; efficacy 0.255; non-depleting), LoF-tolerant (LOEUF 0.678; prec 0.283), tractable. Direction UNKNOWN
in v1 only because the drug is recorded against the calcineurin complex, not this subunit; the mechanism is
concordant. A recovery, not a novel nomination.

**CD3E** — anti-CD3 axis, recovered. Strong suppression (7 genes down; efficacy 0.981), confirmed by the
held-out Schmidt screen (Schmidt 2022). LoF causes immunodeficiency (IEI), the positive-regulator direction,
concordant. Validation.

**IL4R** — target of dupilumab, strongest autoimmune genetic support in this set (0.923 across 2 diseases),
direction concordant (approved antagonist), safe and LoF-tolerant (LOEUF 0.604). Its knockdown does not
suppress the module here (efficacy −0.086): a cytokine-signalling mechanism largely invisible to TCR-only
stimulation, so it ranks low by suppression while passing the annotations.

**PTPN2** — **demoted**. Strong autoimmune genetics (0.762 across 12 diseases), passes the gate at the
evidence-floor minimum (3 genes down; efficacy 0.405), LoF-tolerant (LOEUF 0.619), tractable — but it is a
canonical negative regulator whose loss of function causes autoimmunity, so an inhibitor is anti-therapeutic
(Manguso 2017). A knockdown screen models the wrong direction; genetics is direction-agnostic and cannot
rescue it. A worked example of the direction-of-effect check, not a nomination.

**CD2** — target of alefacept, recovered, autoimmune genetic support (0.68 across 3 diseases), Schmidt-
confirmed (Schmidt 2022). Passes the gate (14 genes down; efficacy 0.68), LoF-tolerant (LOEUF 1.144),
concordant. Validation.

**CD28** — costimulation axis (abatacept), recovered, genetic support (0.72 across 12 diseases), Schmidt-
confirmed. Passes the gate (9 genes down; efficacy 0.654). LoF causes immunodeficiency (IEI), concordant.
Validation.

## 5. Critical review

- **This is a defensible Claude-Use pattern, not autonomous discovery.** The published ranking, the tiers
  and the verdicts are 100% the deterministic pipeline. Claude writes prose over frozen numbers, and the
  prose is checked. It cannot move a rank, invent a statistic, or cite a paper we do not use.
- **The firewall prevents fabrication, not misinterpretation** (stated, per the literature review). A number
  can be grounded and still be interpreted wrongly in prose; the checker guarantees the first, not the
  second. Human review of the prose remains necessary.
- **The card set is honest about what it contains:** six of seven are recovered approved drugs (validation),
  and the one genetically-supported novel-looking gene, PTPN2, is demoted on direction. That is the true
  state of the shortlist after N11-N14, and the cards show it rather than manufacturing a hit.
- **Note on the workflow:** a first attempt to generate the cards via a fan-out of subagents with a strict
  structured-output schema ran away on retries and was stopped; the cards here were written directly and
  validated by the same firewall. The firewall, not the generation method, is the deliverable.

## 6. Verdict

The evidence-card firewall is a working, falsifiable, defensible Claude-in-the-loop component: Claude writes,
a deterministic checker grounds every number and citation, and a seeded fabrication is provably caught. It is
the honest form of the "Built with Claude" story — the model does the writing, not the deciding.
