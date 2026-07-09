# N11 — Direction of effect: is a knockdown-nominated target the right way round?

Run under RULE #9 (work like a scientist): literature → implement → critical review → report.
Script `scripts/21_direction_of_effect.py`; table `results/tables/direction_of_effect.csv`.

## 1. The question and why it matters

Our gate selects perturbations whose **knockdown suppresses** the effector module. A knockdown mimics an
**inhibitor**. So a target nominated by our gate is only a valid drug target for autoimmunity if the
therapeutic direction is *loss of function*: reduce the gene, reduce inflammation. That holds for a
**positive** regulator of activation. It fails for a **negative** regulator — a brake whose loss of
function *causes* autoimmunity — because inhibiting it would make disease worse.

This was not an axis in the gate or in N10. It surfaced because the 2026-07-08 literature review flagged
that our headline novel nomination, **PTPN2**, is a negative regulator: germline loss of function causes
autoimmunity (IBD, type 1 diabetes, rheumatoid arthritis), and PTPN2 *inhibition* is an oncology strategy
to *unleash* T cells (Manguso et al. 2017, *Nature*; Wiede & Tiganis, reviewed). PTPN2 passed our gate only
at the bare evidence floor (3 of 32 module genes down) and was promoted by N10's genetic tier on strong
autoimmune GWAS support — but that support is direction-agnostic, and for PTPN2 the disease-causing
direction is *loss*, so a knockdown screen selecting inhibitor-mimics points the wrong way.

## 2. Method (literature-identified)

The literature review (`docs/literature/target_discovery_ideas_2026_07_08.md`, §1.4, §4.1) specifies a
direction-of-effect / mechanism-of-action concordance table. The rigorous, non-curated version is
direction-aware eQTL colocalisation — is the disease-risk allele associated with *higher* target
expression (concordant with a knockdown therapeutic) or lower (discordant) — which needs the Open Targets
colocalisation dataset and is the next task (N12).

This is **v1**, using direction signals already on disk, applied to the 214 safe genes:
- **CONCORDANT** if the gene has an approved loss-of-function-mimicking drug (Open Targets `actionType` in
  {INHIBITOR, ANTAGONIST, BINDING AGENT, …}), or if loss of function causes immunodeficiency (IUIS/IEI), so
  the gene *promotes* immunity and an inhibitor reduces it. (The immunodeficiency risk itself is the
  separate, already-reported IEI safety annotation.)
- **DISCORDANT** if the gene is a canonical negative regulator of T-cell activation — a curated,
  individually-auditable set of co-inhibitory receptors and intracellular brakes (phosphatases, E3 ligases,
  SOCS/CISH, DGK, TNFAIP3/A20, Roquin/Regnase, and monogenic Treg-restraint genes; sources cited in the
  script). Loss increases activation, so an inhibitor is anti-therapeutic.
- **UNKNOWN** otherwise: no direction signal on disk. This is where the eQTL-coloc method (N12) is needed.

## 3. Controls, and two v1 errors they caught (step 2/3 discipline)

**Control:** no recovered approved-drug target may be mis-called DISCORDANT (UNKNOWN is an acceptable
coverage gap). First run **FAILED** the control and surfaced two curation errors, both fixed and recorded:

- **STAT3** was wrongly placed in the negative-regulator set. STAT3 is context-dependent, but as a CD4 drug
  target it is pro-Th17/pro-inflammatory and STAT3 inhibitors are anti-inflammatory, so knockdown is
  CONCORDANT. Removed.
- **MLST8** was wrongly called discordant by an agonist-drug heuristic. MLST8 is an mTOR-complex component
  whose inhibitors (rapamycin axis) are immunosuppressive. The heuristic was too noisy and was dropped;
  agonist-drug presence is now reported as an annotation only.

After the fixes the control **PASSES**: of the 6 recovered drugs, 5 CONCORDANT, **0 discordant**, and 1
UNKNOWN. The UNKNOWN is **PPP3R1**: ciclosporin/tacrolimus are annotated by Open Targets to the calcineurin
*complex*, not the regulatory subunit — the same whole-complex annotation issue N8 documented — so v1 has no
per-gene signal for it. That is a coverage gap, not a direction error.

## 4. Result

| verdict | n of 214 safe genes |
|---|---|
| CONCORDANT | 21 |
| DISCORDANT | **2** |
| UNKNOWN | 191 |

**The 2 discordant safe genes are `PTPN2` and `RC3H1`** (Roquin-1, a canonical repressor of T-cell
activation whose loss causes the sanroque autoimmune phenotype). Both are correctly flagged: neither is a
valid inhibitor target for autoimmunity.

**PTPN2 verdict: DISCORDANT.** It was N10 nomination rank 4, Tier 1. It is demoted. Re-examining Tier 1
(genetically supported + LoF-tolerant + tractable): `CD2`, `CD28`, `IL4R` are CONCORDANT (and are recovered
drugs, i.e. validation, not novelty); `PTPN2` is DISCORDANT; `ICAM2` is UNKNOWN.

## 5. Critical review — does it survive literature comparison?

- **PTPN2 and RC3H1 as negative regulators: confirmed by the literature.** PTPN2 (TC-PTP) restrains TCR and
  cytokine (JAK/STAT) signalling; reduced-function variants drive IBD/T1D/RA; its inhibitors boost antitumor
  immunity (Manguso 2017; Wiede & Tiganis). RC3H1/Roquin represses ICOS and Tfh programs; its loss causes
  lupus-like autoimmunity (Vinuesa et al. 2005, sanroque). The two DISCORDANT calls agree with established
  biology.
- **The controls agree with the literature:** every recovered approved immunosuppressant target reads
  CONCORDANT (or UNKNOWN by annotation gap), none discordant.
- **The honest limitation, stated plainly:** 191 of 214 safe genes are UNKNOWN. v1 can adjudicate direction
  only for genes with a drug or a curated annotation — which is *not* the novel candidates, exactly where
  the axis is most needed. So this axis, in v1, correctly removes the known error (PTPN2) but cannot yet
  vet the novel shortlist. That is not a null result; it defines the next task.
- **The curation is incomplete and was already wrong twice** (STAT3, MLST8), caught only by the recovered-
  drug control. A curated set cannot be the final method. The rigorous, data-driven replacement is
  direction-aware eQTL colocalisation (N12), which reads direction from the sign of the risk-allele-to-
  expression effect and needs no gene list.
- **What this does to the thesis.** It strengthens it. PTPN2 is now a clean worked example of *why the
  decision layer matters*: a gene with strong autoimmune genetics, nominated by a naive knockdown-suppression
  ranking, that a direction-of-effect check correctly rejects. "Reversal is not enough, and neither is
  genetics" — genetics is direction-agnostic, and a knockdown screen has a fixed direction, so the two must
  be reconciled explicitly.

## 6. Verdict and next step

Direction of effect is a real, previously missing axis. PTPN2 is a **direction-of-effect false positive**
and is demoted from the headline. The gate remains sound for the recovered drugs (all direction-correct).
The novel shortlist cannot be vetted for direction until **N12: direction-aware eQTL colocalisation**
(Open Targets coloc dataset; risk-allele-to-expression sign), which is the rigorous, non-curated method and
the next task in the loop.

**Do not headline PTPN2 as a novel inhibitor target.** If it appears in the report at all, it appears here,
as the worked example of the direction-of-effect check.
