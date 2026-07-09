# N19 — Cross-screen concordance with independent T-cell screens (did more data find a target?)

Run under RULE #9. The discovery push chosen after N16/N17: instead of resolving direction (hard, N17), ask
whether our candidates **replicate** as functional regulators in independent T-cell CRISPR screens. A gene
that is a hit in our Perturb-seq screen and in a separate functional screen is much less likely to be an
artifact of our one dataset. This is the most direct lightweight way "combining datasets" can raise a
candidate's credibility. Script `scripts/27_cross_screen_concordance.py`; table
`results/tables/cross_screen_concordance.csv`.

## 1. The independent screens

Vendored from the source paper's repository (MIT; `resources/external_screens/PROVENANCE.md`):
- **Freimer 2022** (Marson lab): CRISPR screen for regulators of IL2RA, IL2 and CTLA4 marker expression in
  primary human CD4 T cells. Aggregated across the three marker arms to one call per gene.
- **Arce 2025**: fitness screen in Resting_Teff, Stimulated_Teff and Resting_Treg (Stimulated_Teff used as
  the primary condition, matching our stimulated screen).
- **Schmidt & Steinhart 2022 is NOT used here** — it is the held-out validation set (RULE #3).

A gene is an independent-screen "hit" at MAGeCK FDR < 0.05 in either screen. The test universe is genes
tested in both our screen and the independent screen, so there is no collider on our own filters.

## 2. Result: adding this data did not surface a new target

**Coverage is the binding limitation.** These are focused libraries of ~1351 genes. Of our 6,371
DE-tested genes, only 472 are co-tested; of our 214 safe genes, only 21; and of the 7 recovered approved-drug
targets, only 1. The overlap is small because the screens target curated regulator sets that only partly
intersect our panel.

**No enrichment.** Safe genes are not enriched for independent-screen hits beyond a magnitude-matched null:

| screen | safe hits | matched expectation | p |
|---|---|---|---|
| Freimer | 12 / 21 | 11.0 | 0.38 |
| Arce | 6 / 21 | 5.3 | 0.44 |
| either | 14 / 21 | 12.2 | 0.21 |

**Controls.** The negative control (label shuffle) passes (p = 0.73, no spurious enrichment). The positive
control is **underpowered, not confirmatory**: the single recovered-drug target that is co-tested does
replicate, but n = 1 gives binomial p = 0.17, so we cannot claim the drugs replicate here. (Before a merge
fix — Freimer has three marker arms per gene — duplicate rows inflated this to a spurious 3/3, p = 0.001;
that number was wrong and is discarded. RULE #1.)

## 3. The replicated genes (a confidence annotation, not a discovery)

Ten safe, genetically-supported genes are also independent-screen hits:

| gene | our efficacy (rank) | genetics | Freimer | Arce | note |
|---|---|---|---|---|---|
| STAT3 | 0.64 (38) | 0.78 / 8 | CTLA4 arm | suppresses fitness | strongest; both screens; known target, concordant (N17) |
| BATF | 0.88 (7) | 0.58 / 4 | IL2RA arm | suppresses fitness | strong effector TF; both screens; not tractable |
| TBX21 | 0.76 (11) | 0.62 / 5 | CTLA4 arm | – | Th1 master TF |
| PTPRC | 0.11 (2136) | 0.83 / 7 | IL2RA arm | – | CD45; proximal, risky |
| STAT5A | 0.35 (4058) | 0.12 / 2 | IL2RA arm | suppresses fitness | IL2/Treg adjacent |
| PTEN | 0.28 (2099) | 0.34 / 2 | CTLA4 arm | suppresses fitness | tumour suppressor, do-not-headline |
| EGR2, FOXO1 | low | 0.63/10, 0.51/7 | yes | – | tolerance promoters; modality-mismatch (likely discordant) |
| STAT6, IL4R | low | 0.74/2, 0.92/2 | – | Arce | cytokine axis; STAT6 do-not-headline (RULE #6) |

STAT3 and BATF are the only genes strong in our screen and replicated in both independent screens. Both are
already known immune regulators, neither is a novel target, and their therapeutic direction is not resolved
here (STAT3 read concordant but LoF-intolerant; BATF is a transcription factor and untractable). This is a
"these are real regulators" annotation, not a new nomination.

## 4. Critical review — does it survive?

- **The null enrichment is honest, and partly a power problem.** With only 21 safe genes and 1 drug
  co-tested, the test cannot detect a modest effect. So "no enrichment" here does not prove replication fails
  in general; it means these particular focused libraries overlap our candidates too little to decide. A
  genome-wide independent screen would test this properly.
- **Replication is not direction.** Even the genes that replicate carry no resolved therapeutic direction
  (N17). A gene required for stimulated Teff fitness, or an IL2RA regulator, is a real functional node, but
  whether inhibiting it helps or harms in autoimmunity is exactly the question this cannot answer.
- **No new target, consistent with N16.** The replicated set is a subset of the existing watchlist. Adding
  independent-screen data raised confidence in a few genes as genuine regulators; it did not produce a novel,
  vetted nomination.

## 5. Verdict

Combining our screen with two independent T-cell screens **did not find a new drug target**. It provides a
modest, honest replication annotation: STAT3 and BATF are the watchlist genes best supported as real
functional regulators across independent datasets, but both are known, direction-unresolved, and (BATF)
untractable. The main quantitative result is null and underpowered by library overlap. This reinforces the
N16 conclusion rather than overturning it: the honest deliverable remains the audited decision layer and a
credibility-ranked watchlist, not a discovery. The right next dataset for this question would be a
genome-wide independent T-cell screen, where the overlap with our candidates would be large enough to test
replication with power.
