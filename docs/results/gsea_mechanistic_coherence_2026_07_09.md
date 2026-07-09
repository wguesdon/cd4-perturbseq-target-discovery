# N18 — GSEA mechanistic coherence: is the efficacy axis real suppression, or just magnitude?

Run under RULE #9, confirmatory. It does not find a target and does not resolve direction; it answers the
N16 skeptic lens, which asked whether a top-efficacy perturbation genuinely suppresses the T-cell activation
program or just moves a large number of genes (a magnitude artifact). Script
`scripts/28_gsea_mechanistic_coherence.py`; table `results/tables/gsea_hallmark.csv`.

## 1. Method, and the control that makes it non-circular

Preranked GSEA (gseapy 1.3.0) on the stimulated (Stim48hr) panel-wide DE z-scores against MSigDB Hallmark.
The ranking uses only the 10,282 panel genes, so the panel is the background by construction, never the
genome. The key design is a **magnitude-matched control**: the top 25 safe perturbations by efficacy are
compared with 25 low-efficacy safe perturbations matched on overall effect magnitude (mean |z_l2| 152.7 vs
154.3; mean efficacy 0.981 vs 0.283). If the activation suppression were merely "big perturbations move
everything", it would appear equally in the magnitude-matched control. If it is specific to the efficacy
axis, it will not.

## 2. Result: specific suppression of IL-2/STAT5 and TNF-alpha/NF-kB

| Hallmark set | top-efficacy NES (FDR) | magnitude-matched NES (FDR) | reading |
|---|---|---|---|
| IL-2/STAT5 Signaling | **-1.68 (0.009)** | +1.00 (0.574) | suppressed only in the efficacy set |
| TNF-alpha via NF-kB | **-1.66 (0.009)** | +1.86 (0.001) | suppressed in efficacy, *induced* in control |
| Inflammatory Response | -1.33 (0.104) | +1.41 (0.076) | down in efficacy, up in control (weaker) |
| Adipogenesis (housekeeping) | -1.17 (0.200) | -1.25 (0.106) | non-significant in both (negative control) |

**Controls (RULE #1) all pass:** the positive control (IL-2/STAT5 down in the efficacy set) fires at
FDR 0.009; the housekeeping negative control is non-significant; and the specificity contrast holds —
IL-2/STAT5 and TNF-alpha/NF-kB are suppressed in the top-efficacy set but not in the magnitude-matched
control, which has the same overall effect size. The suppression tracks efficacy, not magnitude.

## 3. Honest caveats, stated plainly

- **It is not "the whole activation program".** Interferon-alpha response is *induced*, not suppressed, in
  both groups (NES +2.45 and +2.14, FDR 0.000). That is a non-specific stress signature of strong knockdowns,
  present regardless of efficacy. So the accurate claim is narrow: the efficacy axis specifically suppresses
  **IL-2/STAT5 and TNF-alpha/NF-kB signalling**, not every activation set.
- **The top-efficacy genes are proliferation and metabolism genes** (PGK1, ENO1, KRR1, PRKAR1A, UBR4). Their
  suppression of activation signalling partly reflects impairing an activated, dividing cell. The
  magnitude-matched control accounts for gross effect size, but not for the fact that these are
  proliferation genes. So this confirms coherent program-level suppression; it does not claim the top genes
  are immune-specific regulators.
- **Inflammatory Response is only weakly down** (NES -1.33, FDR 0.104, not significant at 0.05). The clean,
  significant signals are IL-2/STAT5 and TNF-alpha/NF-kB.

## 4. Critical review — does it survive literature comparison?

- **The direction is biologically right.** IL-2/STAT5 and TNF/NF-kB are the canonical CD4 activation and
  effector-cytokine programs; a perturbation that suppresses the effector module suppressing these is
  coherent with the intended readout, and the magnitude-matched control shows it is not a generic large-effect
  artifact. This is exactly the confirmation the skeptic lens asked for.
- **Circularity is bounded by the control.** The efficacy score is defined on a small effector module, and
  the Hallmark sets (~200 genes) are far larger and mostly disjoint, so coherent NES reflects program-level
  suppression rather than the handful of module genes. More importantly, both groups were selected against
  the same module definition; the low-efficacy control does not suppress IL-2/STAT5, so the contrast is not
  an artifact of the module.
- **What it does not do:** it does not nominate a target, and it does not resolve therapeutic direction
  (N17). A perturbation can coherently suppress IL-2/STAT5 and still be the wrong way round for autoimmunity
  if the gene is a brake. This is coherence of the measurement, not a claim about a drug.

## 5. Verdict

The efficacy axis measures **specific suppression of the IL-2/STAT5 and TNF-alpha/NF-kB activation programs**,
beyond what effect magnitude predicts (positive control FDR 0.009; specificity vs a magnitude-matched set).
The skeptic's "just magnitude" concern is answered for these programs. The honest scope: it is IL-2/STAT5 and
TNF/NF-kB specifically (not interferon, which is induced as stress), and the top-efficacy genes are
proliferation-enriched. This is a confirmatory strengthening of the screen's readout for the report, not a
new target and not a direction call.
