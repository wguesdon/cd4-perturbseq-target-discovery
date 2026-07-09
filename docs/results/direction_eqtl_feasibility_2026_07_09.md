# N17 (step 1) — Can combining this screen with eQTL data resolve direction of effect? Feasibility: YES

Run under RULE #9. This is the literature/method + feasibility step, prompted by a direct question: is our
inability to name a vetted novel target a flaw in the method, or a missing data axis that another dataset
supplies? The answer is the latter, and the missing axis is now obtainable.

## 1. Why this is the binding constraint (recap of N11/N16)

A CRISPRi knockdown mimics an inhibitor, so a hit is only a valid autoimmune target if lowering the gene is
therapeutic. Genetics says a gene matters; it does not say whether risk comes from too much or too little of
it. So for a novel gene with no drug and no curated mechanism, neither the screen nor the genetics fixes the
direction. 191 of 214 safe genes are therefore direction-UNKNOWN (N11), and the one filter-intersection
candidate, ICAM2, is unresolved with active counter-evidence (N16). This is a property of the assay, not a
bug. The fix is an external dataset that carries direction.

## 2. The rigorous method, and why N12 stalled

Direction of effect is read from the sign of the risk-allele-to-expression relationship: if an autoimmune
risk allele is associated with **higher** target expression, the gene is elevated in disease and a knockdown
is **concordant**; if the risk allele associates with **lower** expression, loss drives disease and an
inhibitor is **discordant** (the PTPN2 pattern). N12 planned this via the Open Targets colocalisation bulk
parquet and stalled: 21 GB over a slow link, out of budget for a ~44-gene question.

## 3. Feasibility probe: the data is API-accessible, in matched cell states (verdict GREEN)

Both halves of the signal are freely queryable per gene, bypassing the 21 GB download.

**eQTL direction — eQTL Catalogue v2 REST API** (`https://www.ebi.ac.uk/eqtl/api/v2`). Reachable from the
sandbox. A region query returns every cis association with a signed `beta`, `se`, `pvalue`, `rsid`,
`ref`/`alt`. Verified on the textbook ERAP2 eQTL (rs2248374, beta −1.33, p = 6e−30 in DICE naive CD4). The
catalogue holds **CD4 T cell datasets in the exact stimulation state of our screen**, which is the decisive
match:

| study | condition | n | why it matters |
|---|---|---|---|
| Schmiedel_2018 (DICE) | CD4+ T, anti-CD3/CD28 4h | 89 | activated CD4, matches our restimulation |
| Cytoimmgen | CD4+ T, anti-CD3/CD28 16h/40h/5D | ~90–100 each | activated CD4 time-course |
| Schmiedel_2018 (DICE) | CD4+ T, naive | 88 | resting comparator |
| OneK1K | CD4+ T (and TCM/TEM/CTL), naive | ~981 | large naive power |
| BLUEPRINT, Randolph, Perez, Nathan | CD4+ T, naive / flu | 89–249 | replication across cohorts |

Retrieval note: the endpoint paginates by genomic position; `molecular_trait_id`/`gene_id` filters return
"No results", so query the gene's cis window by `pos=chr:start-end` and filter client-side on
`molecular_trait_id`. Confirmed working.

**GWAS risk-allele direction — Open Targets Platform GraphQL API**
(`https://api.platform.opentargets.org/api/v4/graphql`). Reachable (returns PTPN2 for its Ensembl ID). This
supplies the autoimmune credible set / lead variant and its risk direction per target-disease, so no bulk
download is needed. (The GWAS Catalog by-gene REST endpoint returned no associations and is not the route.)

## 4. The N17 build (step 2, next iteration) and its controls

For each watchlist gene (44) plus recovered-drug positive controls:
1. Pull the lead cis-eQTL in activated CD4 (DICE/Cytoimmgen) and naive CD4 (OneK1K), record the
   expression-raising allele from the `beta` sign.
2. Pull the autoimmune GWAS lead variant and risk allele from Open Targets.
3. Harmonise alleles at the shared (or LD-proxy) variant and assign CONCORDANT / DISCORDANT / UNRESOLVED.

**Controls that must fire before any verdict is trusted (RULE #1):**
- **PTPN2 must return DISCORDANT.** Autoimmune risk at PTPN2 is loss of function (Manguso 2017; Wiede &
  Tiganis), so the risk allele should track lower expression/activity. If our proxy reproduces this, it is
  calibrated on the one gene whose direction we already know. This is the key validating control.
- A recovered inhibitor target with genetic support (e.g. IL4R or CD28) should return CONCORDANT.
- A gene with no colocalised eQTL must return UNRESOLVED, never a forced call.

## 5. Honest limits, stated before building

This is a **v2 direction proxy, not formal colocalisation.** It aligns a GWAS lead variant with an eQTL
lead variant; it does not prove they share one causal variant, so LD can produce a spurious sign match.
Allele harmonisation across resources is error-prone and will be checked. The eQTL cell state (4h/40h
anti-CD3/CD28) approximates but is not identical to the disease-relevant state. The proxy can **promote** a
watchlist gene toward a candidate (concordant) or **kill** it (discordant, like PTPN2), pending formal
coloc/MR; it cannot by itself make a nomination. It is exactly the axis the knockdown screen cannot supply,
sourced by combining datasets — which is the answer to the question that prompted it.

## 6. Verdict

The failure to name a vetted novel target is **predominantly a missing-data-axis problem, not a method
flaw.** Direction of effect is unobservable in a single-condition knockdown screen. It becomes observable by
combining the screen with immune-cell eQTL and GWAS direction, and that combination is feasible per gene via
public APIs in matched activated-CD4 cell states. N17 builds it, calibrated on PTPN2. The N12 "infeasible"
verdict stands for the full 21 GB colocalisation dump but is **superseded** for the watchlist by this lighter
per-gene route.
