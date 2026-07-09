# N17 — Direction of effect via eQTL colocalisation: a self-invalidating proxy (honest negative result)

Run under RULE #9. Script `scripts/26_direction_eqtl.py`; table `results/tables/direction_eqtl.csv`. This is
the build the user asked for after N16: resolve direction of effect for the novel watchlist by combining the
screen with immune eQTL and GWAS data. The feasibility step (N17 step 1) was green on data availability. This
step built the method, calibrated it on known-direction genes, and it **did not survive its own controls**.
That is the result, and it is reported as such.

## 1. The method, and why it is self-calibrating

Open Targets precomputes colocalisation between GWAS credible sets and cis-eQTL credible sets, exposing
`betaRatioSignAverage` (the harmonised relative sign of the disease and expression effects) and `h4` (the
posterior that they share one causal variant). In principle a positive relative sign means the
expression-raising allele is risk-raising, so higher expression means more disease, so a knockdown is
CONCORDANT; a negative sign means the gene is protective, so a knockdown is DISCORDANT.

The sign convention is **not** hardcoded (RULE #1). A first read of PTPN2 suggested one mapping; a two-gene
check (TYK2, TNFAIP3) suggested the opposite. So the script learns the convention from a panel of nine genes
whose therapeutic direction is established, and applies it to the watchlist **only if the panel is classified
with accuracy at least 0.80**:
- CONCORDANT anchors, an approved inhibitor treats the disease: TYK2 (deucravacitinib), S1PR1 (ozanimod),
  IL6R (tocilizumab), TNF (anti-TNF), ITGA4 (natalizumab).
- DISCORDANT anchors, monogenic loss of function causes autoimmunity: CTLA4 (CHAI), TNFAIP3 (HA20), SH2B3,
  LRBA.

## 2. Result: the method classifies known genes no better than chance

| panel gene | expected | eQTL-coloc raw sign (pos/neg) | evaluable? |
|---|---|---|---|
| TYK2 | CONCORDANT | NEG (3/15) | yes |
| CTLA4 | DISCORDANT | NEG (4/72) | yes |
| S1PR1 | CONCORDANT | NO_COLOC | no |
| TNF | CONCORDANT | NO_COLOC | no |
| SH2B3 | DISCORDANT | NO_COLOC | no |
| LRBA | DISCORDANT | NO_COLOC | no |
| IL6R | CONCORDANT | HETEROGENEOUS (2/2) | no |
| ITGA4 | CONCORDANT | HETEROGENEOUS (1/0) | no |
| TNFAIP3 | DISCORDANT | HETEROGENEOUS (1/0) | no |

**Panel accuracy 0.50.** Only two of nine anchors were even evaluable, and they are the decisive failure:
TYK2 (concordant) and CTLA4 (discordant) are **both dominant-negative**, so no sign convention can separate
them. The method reported **NOT TRUSTED** and left the entire watchlist `UNRESOLVED`. It did not manufacture
verdicts (RULE #1: a method that fails its controls is void).

## 3. Three failure modes, all real and mostly principled

1. **Sparse colocalisation coverage.** Most anchor and watchlist genes have no autoimmune eQTL coloc at
   h4 ≥ 0.8: S1PR1, TNF, SH2B3, LRBA, and among the watchlist ICAM2, HDAC7, MALT1, TBX21, S1PR1 all return
   `NO_COLOC`. The genes we would most want to vet are the ones with the least coloc data.
2. **Allelic and tissue heterogeneity.** PTPN2 splits almost evenly (21 positive, 17 negative) across its
   autoimmune colocalisations; SMAD4 similarly (21/13). A single per-gene direction label is not well defined
   when independent credible sets at the locus disagree. This is not noise; it is real biology the proxy
   cannot collapse to one call.
3. **eQTL direction captures only expression-mediated effects.** TYK2 is the clean example. Its causal
   autoimmune variant P1104A changes the *activity* of the protein, not its mRNA level, so the mRNA eQTL
   points the wrong way relative to the therapeutic direction. Any expression-only proxy will misclassify
   genes whose disease mechanism is through coding or activity change rather than transcript abundance. This
   is a limitation of the data type, not of the code.

## 4. What this means

- **For the watchlist and ICAM2:** direction stays UNRESOLVED. ICAM2, the single filter-intersection
  candidate, has no autoimmune eQTL coloc at all, so even this route says nothing about it. The N16 verdict
  is unchanged.
- **For the user's question (method flaw vs data gap):** it sharpens the answer. Direction of effect is
  indeed the missing axis, and the *data* to address it exists (N17 step 1). But a *reliable* resolution is
  genuinely hard: a lightweight coloc-sign proxy fails its controls, because the anchors lack coverage, loci
  are heterogeneous, and expression is not the whole mechanism. So "zero vetted novel targets" is not a
  weakness of our screen or a quick-fixable gap. It reflects the intrinsic difficulty of assigning direction,
  which is why the decision layer, not a single dataset, is the honest contribution.
- **A couple of genes gave a clean, biology-consistent signal but are not promoted:** CTLA4 read
  cleanly negative (consistent with a protective, discordant gene) and STAT3 cleanly positive (17/0,
  consistent with the concordant call in the N16 wider triage). These are encouraging but anecdotal; with the
  panel at 0.50 the method is uncalibrated and none of its calls are applied.

## 5. Critical review — does it survive literature comparison?

The failure modes match the known limitations of eQTL-based causal-direction inference. eQTL and protein/pQTL
or activity effects are frequently discordant, so expression direction can mislead about therapeutic
direction (the TYK2 P1104A case is textbook). Colocalisation is tissue- and context-specific and often
absent for a given gene-trait pair, capping coverage. Allelic heterogeneity at autoimmune loci is well
documented, and averaging a sign across independent credible sets is not meaningful. None of this contradicts
the field; it reproduces, on our exact candidate set, why direction of effect is treated as a hard,
dedicated problem (Mendelian randomisation with pQTLs, activity-aware fine-mapping) rather than a lookup.

## 6. Verdict and what a rigorous version needs

The lightweight eQTL-coloc direction proxy is **not trustworthy for these genes** and is not shipped as a
vetting axis. Resolving direction properly would need: tissue-matched (activated CD4) colocalisation, GWAS
betas explicitly oriented to disease risk, per-credible-set (not per-gene-averaged) handling of allelic
heterogeneity, and pQTL or coding-variant direction to cover activity-mediated genes. That is a dedicated
study, not an API proxy, and it is the honest next step for anyone who wants to vet this watchlist. Until
then the novel candidates remain direction-UNRESOLVED, hypothesis-generating, and not nominations — exactly
the N16 conclusion, now reinforced by a method that tried to overturn it and could not.

The value of this iteration is the near-miss it prevented: a first two-gene calibration looked like it
worked, and a per-gene direction label would have been shipped with an inverted, chance-level convention.
The panel and its 0.80 gate caught it before it reached the report.
