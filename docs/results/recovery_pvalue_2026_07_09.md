# N14 — "The gate recovers approved drug targets", as a p-value

Run under RULE #9. Script `scripts/23_recovery_pvalue.py`; table `results/tables/recovery_pvalue.csv`.

## 1. Question and method

The project's most intuitive validation is that the gate, given only the knockdown transcriptome, recovers
targets of approved immunomodulators. Until now that was an anecdote ("it recovers mycophenolate,
calcineurin, …"). The literature review (§1.2) makes it rigorous only if two conditions hold: the recovery
is quantified against the **screened-gene background** (not the genome), and no threshold was **fit to the
known targets** (Kriegeskorte et al. 2009, double dipping). Recovering known targets to validate a screen
is itself standard (Shifrut et al. 2018; Frangieh et al. 2021; Schmidt & Steinhart 2022).

## 2. The recovery, decomposed (hypergeometric + 50,000-draw permutation, screened background)

Of 36 curated approved-immunomodulator targets, **20 are in the screened panel**.

| test | observed | expected by chance | enrichment | permutation p |
|---|---|---|---|---|
| positives passing the **evidence floor** (efficacy axis) | 5 of 20 | 0.90 | **5.6×** | **0.0017** |
| positives passing the **full gate** (safe) | 4 of 20 | 0.67 | 6.0× | 0.0037 |
| positives kept by **safety, given evidence** | 4 of 5 | 3.74 | 1.1× | p=0.63 (n.s.) |

**The recovery is driven by the efficacy axis, not the safety axes.** The evidence floor — a perturbation
must significantly suppress ≥3 effector-module genes — is what pulls in known inhibitor targets, at 5.6×
chance (p=0.0017). Conditional on passing it, the safety axes do **not** preferentially retain drugs
(p=0.63), which is exactly correct: the gate does not know which genes are drugs. The safety axes simply
apply the therapeutic-window filter, keeping the 4 with a window and correctly rejecting **CD3G**, a
narrow-therapeutic-index anti-CD3, on the co-inhibitory axis.

The recovered set is `CD3E`, `IL4R`, `IMPDH2`, `PPP3R1` (all pass the gate) plus `CD3G` (efficacy-visible,
safety-rejected). `IMPDH2` is mycophenolate; `PPP3R1` is calcineurin (the non-redundant regulatory
subunit); `IL4R` is dupilumab; `CD3E`/`CD3G` are the anti-CD3 axis.

## 3. Recall is bounded by assay blindness, not by the gate

**15 of 20 screened positives never pass the evidence floor** — their knockdown does not suppress the
module in this assay. These are not gate failures; they are invisible to the screen by design:
cytokine-signalling targets with no polarising cytokine applied (`JAK2`, `TYK2`, `S1PR1`, `MTOR`, `IL2RA`,
`PDE4B`, `CCR4`), and the calcineurin catalytic subunits (`PPP3CA`, `PPP3CB`) silenced by paralogue
compensation so only the regulatory subunit `PPP3R1` is visible. Recall is a property of the perturbation
library and the assay condition (N8's benchmark ceiling), not of the ranking. So the precise claim is: the
efficacy axis significantly recovers the **subset of known targets the assay can see**, not all of them.

## 4. Double-dipping / leakage audit (Kriegeskorte 2009)

| quantity | provenance |
|---|---|
| evidence floor = 3 module genes down | a priori from the DE structure (`04` docstring), not fit to the positives |
| tolerance gate = p75 of evidence-passers | a data-defined quantile, not fit to the positives |
| effector & co-inhibitory module gene sets | from activation biology (`01`) + curation, not from the drugs |
| ground-truth positives | used only to **report** where drugs land (`04` `validate()`), never to set a threshold |
| Schmidt & Steinhart held-out screen | RULE #3: never enters the score, gate, or any threshold |
| CD2 (alefacept) at Tier A | RULE #4: deliberately not added to the positives after watching it rank |

**The recovery is out of sample.** No gate threshold or module was fit to maximise it; the positives touch
only the reporting tables, not the decision rule.

## 5. Critical review

- **Against the literature:** recovery-of-known-targets is the accepted validation, and the rigorous form
  is exactly this — a permutation p against the screened background with a leakage audit. The result agrees
  with the field's expectation that a functional readout recovers direction-concordant inhibitor targets.
- **The decomposition is the honest part.** A weaker analysis would report "4 of 20 known targets are safe,
  6× chance" and imply the safety gate is validated by drug recovery. It is not: the efficacy axis carries
  the signal; the safety axes are not drug-enriched given evidence (p=0.63). Reporting the full gate's 6×
  without the conditional would overstate what the safety layer recovers.
- **n is small** (20 screened positives, 5 visible). The p=0.0017 is solid, but the estimate rests on 5
  genes; it is a direction-and-significance result, not a precision claim.
- **This is validation, not the headline.** The drug-recovery benchmark as a scored deliverable is retired
  (N8): the library ceiling is 38–53 of a 60-target set, a property of the library. Recovery here is the
  sanity check that the efficacy axis is finding real inhibitor targets, which licenses trusting the gate's
  novel nominations — subject to the direction-of-effect caveat (N11) that genetics and the screen cannot
  resolve for a novel gene.

## 6. Verdict

The efficacy axis recovers approved-immunomodulator targets far above chance (5.6×, permutation p=0.0017),
out of sample, with recall capped by assay blindness. The safety axes apply the therapeutic-window filter
without themselves being drug-enriched. This is a defensible, quantified validation that the gate is
sensible — and it is explicitly not the headline.
