<!-- Synthesised 2026-07-08 by a 4-reader / 152-agent workflow over the source preprint
     (docs/2025.12.23.696273v1.full.pdf, gitignored). 147 load-bearing claims were each
     independently re-checked against this repo's code and committed tables; 16 were dropped
     as misreads of the PDF; 11 CONTRADICT our code. Quotes are from a CC BY-4.0 preprint. -->

# BRIEFING: Zhu, Dann et al. vs. our repo — what has to change tonight

Source: `docs/2025.12.23.696273v1.full.pdf`, 63 pp. Page numbers below are printed page = PDF page.
All numbers below were recomputed against committed tables in `results/tables/` unless marked otherwise. No h5ad layers were read.

Repo state at time of writing: `git status` clean. `results/tables/selection_funnel.csv` at HEAD **already reads 12,779 / 230 / 135**. The "uncommitted fix" from an earlier reader has landed. Only prose is stale.

---

# 1. THINGS THAT CHANGE OUR ANALYSIS

Ordered by blast radius.

## 1.1 We rank on `log_fc`. The authors' default effect estimate is the z-score, and they say why.

**Status: not refuted, but we are off the authors' stated default with no defence in the code.**

Paper, p.37, Methods: *"To account for differences in noise between perturbations with different number of recovered cells across replicates, unless otherwise specified we use the z-score of the DE log-fold change (β̃_j,k) as normalized perturbation effect estimates, where β̃ = β / SE(β). The standard error of the log-fold change is obtained from DESeq2 outputs."*

Paper, p.7, main text: *"increasing cell numbers per perturbation consistently improved differential expression estimates … though estimates stabilized at approximately 200 cells per perturbation."*

Our code: `scripts/04_window_score.py:186-187`
```
frame["efficacy"] = -frame["eff_mean_lfc"]
frame["tolerance_loss"] = -frame["tol_mean_lfc"]
```
`eff_mean_z` is computed at `:124`, stored at `:171`, and never consumed. `scripts/02_risk_kill_reversal.py:199` and `scripts/12_magnitude_matched.py:167` already use `zscore`. The repo is internally inconsistent.

Recomputed from `results/tables/window_score.csv` (n = 6,371):

| quantity | value |
|---|---|
| spearman(`eff_mean_lfc`, `eff_mean_z`) | 0.8664 |
| spearman(&#124;`eff_mean_lfc`&#124;, `n_cells_target`) | **−0.2828** |
| spearman(&#124;`eff_mean_z`&#124;, `n_cells_target`) | −0.1453 |
| safe genes with `n_cells_target` < 200 | **20 of 66** (base rate 9.3%) |
| top 20 by `window_rank` below 200 cells | **9 of 20** |
| spearman(window_score_lfc, window_score_z) | 0.9273 |
| safe-set top-10 overlap under z | **7 of 10** |

Under `efficacy = -eff_mean_z`, the safe-set top ten **loses CARS2, CD3E, IBA57** and **gains CD28, PPP3R1, ZNF649**. IMPDH2 falls 12 → 40. CARS2 falls 10 → 26. Both are cases of a large fold change with almost no statistical support: IMPDH2 `eff_mean_lfc` = −0.598 at `eff_mean_z` = −0.163; CARS2 −0.490 at −0.155.

Our efficacy axis is twice as coupled to recovered cell number as the authors' own effect estimate, and 45% of our top 20 sits below the cell count at which the authors say the estimate has not stabilised. That is a winner's-curse signature.

**Change:** `scripts/04_window_score.py`
- `:173` add `"tol_mean_z": tol["mean_z"]` to the frame dict. `_module_stats` already returns it.
- `:186-187` compute `efficacy = -eff_mean_z`, `tolerance_loss = -tol_mean_z`.
- `:47` and `:51` add the p.37 citation.

**Mandatory / discretionary:** the *sensitivity run* is mandatory. The *primary-measure switch* is a preregistration deviation, so the N6/N9 primaries in `docs/preregistration_n6_2026_07_08.md` and `docs/preregistration_n9_2026_07_08.md` must be reported unchanged, with a dated amendment recording the z-based re-rank. Do not silently swap. Do not ship CARS2, IBA57 or IMPDH2 in a top-ten without disclosing that they are measure-fragile.

Note the paper's "unless otherwise specified" does license `log_fc`. What it does not license is using `log_fc` with no stated reason while the cell-count coupling is visible in our own table.

## 1.2 `selectivity` is built from `n_total_de_genes`. It must be built from `n_downstream`.

**Status: correctness bug. Our own `scripts/12_magnitude_matched.py:57-59` already documents it and already fixed it there. It never propagated to the script that produces the headline table.**

Paper, p.39, Methods: perturbed genes qualify by *"significant trans-effects on at least one gene other than the target gene"*. Paper, p.42: *"To preclude strong on-target knockdown effects from driving cluster assignment, on-target z-scores were masked to zero."* Paper, p.48 and p.50: *"The knock-down effect of each perturbed gene on itself was masked before training"*; *"We masked the effects of knock-down gene j on itself, as it does not reflect a trans-regulatory effect."*

Verified identity, 0 mismatches over all 33,983 rows in **all three** conditions (our notes say "both conditions"; it is all three):
```
n_total_de_genes == n_downstream + ontarget_significant
```

The Stim arm is QC-filtered on `ontarget_significant`, so it always carries the on-target +1. The Rest arm is unfiltered, so it carries the +1 for only 5,683 of 6,320. The subtraction at `scripts/04_window_score.py:185` is asymmetric by construction, one-sidedly inflating selectivity.

Recomputed impact, using `stim_downstream` / `rest_downstream` already present in `results/tables/risk_kill_naive_reversal.csv`:

| | |
|---|---|
| `fail_homeostasis` flips | 51 (50 fail→pass, 1 pass→fail) |
| `safe` | 66 → **67** |
| gained | **TBX21** (T-bet, `window_rank` 20, currently `reject_reason = homeostasis`, selectivity 2.28784 vs gate 2.30259) |
| lost | none |
| spearman(old, new) | 0.9788 |

**Change:** `scripts/04_window_score.py:175` and `:181` source `n_downstream`. Emit both `stim_downstream`/`rest_downstream` and the totals into `window_score.csv`. Mirror in `scripts/02_risk_kill_reversal.py:214-220`, which already carries `n_downstream` alongside; make it primary.

**Do not touch:** the `ACTIVE_DE >= 100` VIF gate at `:162`. It is an activity filter, and the authors' own strong-perturbation filter (p.42) uses total DE genes. Also `scripts/18_n6_selectivity_validation.py:98-99` asserts the shipped formula and will fire. Update in the same commit.

**Mandatory.** It is a correctness fix, it costs nothing, and it only adds a gene.

## 1.3 The Rest arm is never QC-filtered, and 51 genes get a fabricated denominator.

**Status: our construction defect. The paper is silent on per-condition QC, so this is ours to defend.**

`scripts/04_window_score.py:180-182`
```
rest = obs[obs["culture_condition"] == REST].drop_duplicates(gene_col).set_index(gene_col)
frame["rest_de_genes"] = frame["gene_name"].map(rest["n_total_de_genes"])
frame["rest_de_genes"] = frame["rest_de_genes"].fillna(frame["rest_de_genes"].median())
```
`:146` applies `qc_mask` to the Stim arm. `:180` applies nothing to Rest.

Verified: **the imputed value is exactly 3.0**, and it is assigned to the 51 rows with `viability_tier == "unknown"`.

Of the 6,371 QC-passing Stim48hr perturbations:
- 51 have **no Rest row at all** and are median-imputed at 3.0.
- 6,320 have one, of which **802 would fail the Stim QC**: 637 on `ontarget_significant == False`, 107 `neighboring_gene_KD`, 69 `distal_offtarget_flag`, 66 `low_target_gex`.
- Of 370 selectivity-gate passers, **84 (22.7%)** rest on an imputed or QC-failing denominator.

Two directional biases, both in the permissive direction:
1. Rest rows failing `ontarget_significant` have median `rest_de_genes` = 1.0 versus 3.0 for passers. A guide that never knocked its target down at rest is scored as context selectivity.
2. Absence from the Rest arm means the gene failed the authors' ≥3-pseudobulk floor (p.36). That is evidence of **depletion**, not of resting inertness. The 51 include MTOR, RHOA, CENPE, EEF1G, DDOST, YARS1, ANAPC5, NUP133, RPN1, PSME3, DCK, TXNRD1. We impute them as inert.

The same function already refuses to guess for exactly these rows on the viability axis. `scripts/04_window_score.py:199-201`: *"A gene with no resting arm has unknown viability, and NaN >= 0.5 is False, so a plain np.where would silently file missing data as evidence of toxicity."* Line 182 commits the mirror-image sin.

**Change:** `scripts/04_window_score.py`
- `:180` build `rest` from `obs[(obs.culture_condition == REST) & qc_mask(obs)]`, or carry an explicit `rest_qc_pass` boolean.
- `:182` delete the `fillna(median)`. Let `rest_de_genes` be NaN, propagate NaN into `selectivity`, add a third `selectivity_tier` value mirroring the existing three-way `viability_tier`.
- add `rest_qc_pass` and `rest_de_imputed` columns to `window_score.csv`.

Expected: selectivity-gate passers 370 → ~286. Top-20 casualties, all already `safe == False`: LMLN (rank 11, rest `ontarget_significant` False), ABCA3 (rank 17, same), SMIM7 (rank 16, no Rest row, imputed).

**Mandatory.** Median imputation on the denominator of a safety gate is indefensible in review.

## 1.4 REFUTED. "Cytokine signalling targets are close to invisible."

**Refuted by our own committed table, and separately by the paper.**

`README.md:51-53` currently reads: *"The screen stimulates through the TCR with no polarising cytokines. Cytokine signalling targets are therefore close to invisible. `JAK2` ranks 5392, `TYK2` 5618, `S1PR1` 5993, `IL4R` 6047 of 6371."*

The four ranks are correct. They are the `rank` column of `results/tables/risk_kill_naive_reversal.csv`, i.e. a **signed effector-module suppression rank**, not a measure of signal. Verified exact.

| gene | naive rank | `window_rank` | `stim_de_genes` | `safe` |
|---|---|---|---|---|
| JAK2 | 5392 | 6181 | 6 | False |
| TYK2 | 5618 | 3738 | 54 | False |
| S1PR1 | 5993 | 5243 | 3 | False |
| **IL4R** | 6047 | **1785** | **729** | **True** |
| STAT6 | 5916 | **251** | **1128** | **True** |
| STAT5B | **1** | 163 | 1060 | False |
| IL2RB | **4** | 42 | 1070 | False |

IL4R is in our own safe set. STAT6 is in our own safe set. STAT5B is the single strongest effector-suppressing perturbation in the entire screen and IL2RB is fourth. The word "invisible" is false, and a reviewer holding `window_score.csv` and `README.md` finds this in one minute.

The mechanism is also wrong on two counts:

1. **The culture is not cytokine-free.** Paper, p.32, Methods: *"seeded at 1 x 10^6 cells/ml cXVIVO supplemented with 200 IU/ml of recombinant human IL-2"*, re-supplied at every media change including day 12. IL-7 at 5 ng/mL from day 5, but **absent from the day-12 split medium** (p.33), therefore absent from the entire DE measurement window. The `no polarising cytokines` claim survives verbatim; the implied "no cytokines" does not. Exogenous IL-2 is precisely why the IL-2R/STAT5 axis is our loudest axis.
2. **The authors' own explanation is different.** Paper, p.16: *"As with other cytokine perturbations, IL4 knockdown showed relatively weak transcriptional effects, likely due to paracrine signaling between cells in the pooled experimental set-up dampening the perturbation response."* That is a **ligand** effect. Our data agrees: IL4 has 3 DE genes, naive rank 3756. It says nothing about receptors or kinases.
3. **The paper recovers JAK2 and IL4R from the same data.** Paper, p.15-16, Fig 4F: *"Top regulators promoting a Th1 state included: … JAK2 … Top Th2 regulators included: IL4R; STAT6 … and GATA3"*, via signature projection at **Stim8hr**, a condition our pipeline discards entirely.
4. **JAK1, JAK3 and IL2RG are not in the library.** IL6R and IL7R are perturbed but fail QC. Those nodes are invisible by design, not by biology.

**Change, mandatory, three files:**
- `README.md:51-53`
- `docs/ground_truth_immunomodulators.md:91-95` (heading "The assay cannot see cytokine-signalling targets" must go)
- `docs/preregistration_n9_2026_07_08.md:169-170` — **append a dated amendment, do not rewrite a preregistration in place.**

Replacement claim, defensible: *a directional effector-suppression ranking cannot surface cytokine-signalling targets, because their knockdown moves the module the wrong way, not because they produce no signal. IL4R carries 729 DE genes and passes our window gate at rank 1785. Only JAK2 (6 DE genes) and S1PR1 (3) are genuinely near-silent. The culture carries IL-2 at 200 IU/mL throughout (p.32) and no lineage-polarising cytokine (p.18, p.30), so the JAK1/JAK3/STAT5 axis is live and the JAK2/TYK2-coupled axes are not. Zhu, Dann et al. recover JAK2 and IL4R from the same data by projecting perturbation signatures onto a Th2/Th1 signature at Stim8hr (p.15-16). Magnitude ranking and signature projection are not interchangeable. This bounds our method, not the assay.*

## 1.5 REFUTED. Design facts in our prose.

Three separate errors in one sentence, appearing in three places.

`docs/preregistration_n9_2026_07_08.md:20` and `scripts/17_tolerance_is_special.py:26,569`:
> *"a 48 h anti-CD3/CD28 stimulation of **bulk primary human CD4 conventional T cells**"*

Refuted:
- **"bulk"** — the cells are naive. Paper p.3: *"Primary naive CD4+ T cells isolated from four human donors"*. Paper p.32: *"isolated with the EasySep Human Naive CD4+ T Cell Isolation Kit (catalog no. #19555)"*. This **strengthens** our FOXP3 argument: naive isolation depletes thymic Tregs at day 0.
- **"anti-CD3/CD28"** — the activator is a soluble three-antibody tetramer. Paper p.32: *"ImmunoCult human CD3/CD28/CD2 T cell activator (10990, STEMCELL Technologies) at 25 ul/ml"* at day 0; p.33: *"12.5 uL/ml"* at day 12. No beads. No plate-bound.
- **"48 h stimulation"** — Stim48hr is 48 h after **restimulation** of cells expanded for 12 days. Paper p.33: *"At day 12 morning, cells were … split into three populations. The first population (Rest) was kept at a rested state without restimulation and harvested in the afternoon after 8 hours. The second population (Stim8hr) was restimulated … harvested … after 8 hours. The third population (Stim48hr) was restimulated … and harvested after 48 hours."*

Consequence for `scripts/04c_dropout_selectivity.py:3` (*"Resting CD4 T cells are quiescent"*), `scripts/04_window_score.py:190-194`, `scripts/07_organism_safety.py:3-5`: **the Rest arm is not quiescent and is not a within-experiment control for viability.** It is an 8 h no-restim aliquot of a day-12 IL-2 blast pool. Guide dropout accrued during the 12-day expansion shared by all three arms. Median log2(Stim48hr / Rest) `n_cells_target` = +0.000; spearman(rest_cells_ratio, stim n_cells) = 0.94. There is no arm-specific dropout to measure.

This does **not** kill the organism-safety layer. It strengthens it: the screen contains no quiescent-cell arm at all, so a mitochondrial or one-carbon gene has no in-screen safety readout, which is exactly why LOEUF and HPA are needed.

**Mandatory prose fix.** The string feeds no computation. `grep -riE "immunocult|dynabead|bead|plate-bound"` returns 0 hits, so we never asserted a bead format. Do not now write "beads".

## 1.6 REFUTED. `src/cd4_perturbseq/programs.py:17`

Line 17 reads:
```python
        # Regulatory T cell identity and suppressive effectors.
```
directly above FOXP3, IKZF2, IL10, TGFB1, LRRC32. The paper contains **zero** occurrences of "Treg", "regulatory T cell" outside a reference title, and **zero** occurrences of "suppressive". Our own `docs/preregistration_n9_2026_07_08.md:29-30` already forbids this claim, and `scripts/17_tolerance_is_special.py:26-29` already states the correct science. The library module contradicts both.

`programs.py:3-6` docstring also says *"sparing the tolerance program"*.

**Change, mandatory:** delete the comment. Rewrite the docstring. See §5 for the naming.
**Do not rename** `TOLERANCE_GENES` / `tolerance_module()` / the six `tolerance_*.csv` files without keeping aliases: 7 scripts and 5 committed tables depend on them, and the preregistration cites those filenames.

## 1.7 The Rest arm is 25% shallower. `selectivity` is anti-conservative, not conservative.

Paper, p.4, main text: *"mean UMI/cell per condition: Rest = 10,080; Stim8hr = 14,977; Stim48hr = 13,373."* Cells per perturbation are matched (Fig 1C medians 500 / 509 / 494); depth is not. Ratio Stim48hr/Rest = 1.327.

Shallower Rest → fewer Rest DE calls → inflated selectivity → `fail_homeostasis` is too permissive. `grep` for UMI/depth/sequencing across `scripts/` and `src/`: **zero hits**. We never accounted for it.

First-order correction, recomputed (rest_de_genes × 1.327):

| rest scale | `safe` | lost |
|---|---|---|
| ×1.000 (shipped) | 66 | — |
| **×1.327** | **63** | **ZNF649, AP2A1, DDX39B** |
| ×1.5 | 59 | +4 |
| ×2.0 | 53 | +10 |

**ZNF649 is `window_rank` 19 and in the safe set.** It does not survive the honest first-order depth correction.

**Change, mandatory:** add `DEPTH_RATIO = 13373 / 10080` to `scripts/04_window_score.py` with the p.4 citation, print a depth-sensitivity block every run, amend the `MIN_SELECTIVITY` docstring at `:52` to state that selectivity is an **upper bound**. Annotate or drop ZNF649 from any reported top 20.

Do **not** cite the authors' 10% downsampling result (p.6-7) as evidence for a 1.33× effect. That is a tenfold cut and the inference does not follow.

## 1.8 Our QC mask destroys perturbations the authors independently validated.

This is the strongest available evidence for the `ontarget_significant` collider, and it is better than the essentiality counts we currently use.

Of the authors' **9 core regulators of T cell activation** (Fig 3E legend, p.14: CD3D/E/G, CD247, ZAP70, LAT, LCP2, PLCG1, VAV1), **three are absent from our 6,371**:

| gene | why dropped | `n_total_de_genes` @ Stim48hr |
|---|---|---|
| CD3D | `neighboring_gene_KD` | 1,432 |
| LAT | `distal_offtarget_flag` | 3,187 |
| **LCP2** | **`ontarget_significant == False`** | **1,535** |
| MED24 | `neighboring_gene_KD` | — |

MED24 is the most strongly arrayed-validated IL10 positive regulator in Fig 2F (p.10), three stars on both bulk RNA-seq and protein staining. Our mask deletes it.

LCP2 is the money case: an author-published core regulator of T cell activation, 1,535 DE genes, discarded solely because its own transcript was not called significantly knocked down.

The paper states the collider mechanism itself, twice. Paper p.4-5: *"The subset of guides that failed to achieve significant knockdown predominantly targeted genes with very low baseline expression."* Paper p.37: *"Beyond low baseline expression, perturbations without detectable on-target knockdown also showed significantly fewer cells per perturbation (Kolmogorov-Smirnov test, p = 6.57e-12)."*

So `ontarget_significant` is a collider on **two** parents: target expression **and** cell number. Our `src/cd4_perturbseq/de_stats.py:59-61` names only the expression parent. And `src/cd4_perturbseq/magnitude.py:3-6` currently asserts *"Cell count is a viability readout, not a power proxy"* as a general law. Read literally that argues against the collider three files away. It is a **post-selection** observation, true only among rows that already passed the gate.

**Changes, mandatory:**
- `src/cd4_perturbseq/de_stats.py` `qc_mask` docstring: record both collider parents, cite p.4-5 and p.37.
- `src/cd4_perturbseq/magnitude.py:3-6`: scope the sentence to "among perturbations that already passed `ontarget_significant`".
- Add `CORE_REGULATORS = frozenset({"CD3D","CD3E","CD3G","CD247","ZAP70","LAT","LCP2","PLCG1","VAV1"})` to `programs.py`, cite Fig 3E legend p.14, and replace the ad-hoc `{CD3D,CD3E,CD3G,CD52,DHFR}` literal at `scripts/04_window_score.py:382`.
- Lead the collider argument with "3 of 9 author-named core regulators, plus a three-star arrayed-validated positive control, are dropped by our mask". Keep the 3-of-141 nonessential figure as support.

## 1.9 REFUTED. "The paper uses zero druggability **or therapeutic** language."

The druggability half survives. The therapeutic half does not.

Full-text counts over all 63 pages: `druggab*` = 0, `tractab*` = 0, `"therapeutic window"` = 0, `"approved drug"` = 0, `"clinical trial"` = 0, `triage` = 0, `"drug target"` = 0. `drug` = 2, both in reference titles.

But four authorial sentences:
- p.15: *"Furthermore, regulators driving disease-associated or treatment-responsive states could serve as therapeutic targets."*
- p.24: *"translating these findings into actionable therapeutics requires understanding the causal genes, the cell types in which they exert their effects, and the precise nature of their molecular effects."*
- p.29: *"…can be harnessed for insights into human immunology, immunotherapy design, dynamic gene regulatory control of human cells and human disease genetics."*
- p.30 (Limitations): *"…in order to harness these data for immunology and immunotherapy development."*

**Change, mandatory, prose only:** `docs/strategy_2026_07_07.md:91-92` and `docs/claude_tooling_log.md:17`. Replace with the narrow, stronger claim. See §6.
**Leave `docs/research_plan.md:21` alone.** Its three sub-claims all verify.

## 1.10 REFUTED. "The paper recommends no thresholds" on the confidence columns.

Fig 2C legend, p.11: *"Volcano plots of log2 fold change of IL10/IL21 expression relative to NTCs and −log10 FDR for each regulator that has **minimum cross-donor correlation > 0.35** (Suppl. Table 5). Significant regulators (FDR < 10%, |Log2FC| > 1) are colored blue…"*

The column is `donor_correlation_hits_min`, not `_all_min`. On Stim48hr QC-passers with a value (n = 1,106): `hits_min > 0.35` keeps 953 (86.2%); `all_min > 0.35` keeps 28 (2.5%). A 28-point volcano is not Fig 2C.

**Do not adopt as a QC gate.** `donor_correlation_*` is NaN for 82.6% of QC-passing Stim48hr rows, and the NaN pattern is the authors' own eligibility rule, not a quality signal. Gating would cut the universe 6,371 → 953 and stack a second collider on the first.

**Do adopt as free external validation.** 18 of our top 20 by `window_rank` pass `hits_min > 0.35`; 43 of top 50; 84 of top 100. The threshold never enters our ranking. Add `donor_correlation_hits_min` as a reported column on `window_score.csv` and say so.

It does **not** exonerate PPP3R1: `hits_min` = 0.737 while `guide_correlation_signif` = −0.390. Cross-donor reproducibility and cross-guide concordance are orthogonal.

## 1.11 Stale numbers in prose. Code is already right.

`results/tables/selection_funnel.csv` at HEAD: `in sgRNA library,682,141,12779` and `reaches DE_stats,230,135,11526`. Committed. Correct.

Still stale in prose:
- `scripts/02_risk_kill_reversal.py:343` docstring: *"465 of 682 library essentials never reach the DE table"*. Correct figure is **452** (682 − 230).
- `docs/results/risk_kill_2026_07_08.md:30` still says *"only 217 (31.8%)"* although line 14 is an errata correcting it to 230 (33.7%).
- `scripts/11_selection_funnel.py:57,127` mention 13,129 correctly, as the *diagnosed bug*. Leave.

The library, verbatim from the paper and from `data/raw/suppl/sgrna_library_metadata.suppl_table.csv`:
- **12,748 perturbed genes**, **26,504 gRNAs** (p.31 Methods, Fig 1A panel p.7). Row count in the shipped table is exactly 26,504.
- **12,779** unique `designed_target_gene_name` in the shipped table. **31-gene discrepancy with the paper, unexplained.** Carry it as a footnote. Do not adopt either number silently.
- **992** non-targeting guides (3.74%, not the stated 5%). 25,512 targeting guides.
- Guides per gene: 1 → 104 genes, 2 → 12,621, 3 → 50, 4 → 4. **"2 gRNAs per gene" is the mode, not the design.** Write "median 2, range 1 to 4".
- Never write **13,129**. It is an alias-union artefact of two symbol vocabularies that disagree on 350 genes.

**Mandatory prose.** `grep -rn "13,129\|13129\|only 217\|465 of 682" docs/ scripts/ reports/ README.md`.

## 1.12 Prior art we do not cite. This is a novelty risk, not a correctness risk.

Three author-built objects that a reviewer will expect us to know:

1. **Regulator-burden correlation** (p.22, Fig 6A). *"to identify putative core genes, we tested for each gene j whether there is a significant correlation between the knockdown effects of regulators on j and the LoF effects of those regulators on lymphocyte count in the UK Biobank."* This is the paper's own candidate-prioritisation framework. `grep` for it in `scripts/` and `src/`: 0 hits. Cite it. It is on a different axis: they prioritise **downstream measured genes** as trait core genes, we prioritise **perturbations** as druggable targets. Of the 15 named putative core genes rankable by us, median `window_rank` = 3,582 of 6,371 (uniform expectation 3,186), and **0 of 15 pass our safety gate**. Report this as evidence of orthogonality, not as a concession.
2. **`pert2state`** (Fig 4A, p.45-46). Reconstructs an *external* observational DE signature Δ_j ≈ Σ_r β_{r,j} w_r from a linear combination of perturbation effects, yielding regulator weights w_r. 15 train/test splits (5-fold CV × 3 initialisations). On-target effect masked before training. Signatures: Th2/Th1 from Ota 2021 (n = 79 discovery) and Hollbacher 2021 (n = 3 replication); OneK1K aging (782 / 199 donors); a TCR activation signature from ref [13] = Arce et al., used as a **negative control** (p.48). Weights in Suppl. Table 12. **The paper does not say "ElasticNet" anywhere in the pages read. Stop asserting it.** Cite as prior art for signature-projection nomination, which is how they recover JAK2 and IL4R where our magnitude ranking cannot.
3. **Suppl. Tables 9 and 10.** We already downloaded `data/external/gwt_priors/clustering_results_and_annotations.csv` via `scripts/fetch_priors.sh:47` and **never opened it**. Zero references in `src/` or `scripts/`. See §4.

**Mandatory** for the manuscript. Discretionary for tonight's code.

## 1.13 Smaller items, all discretionary tonight, all mandatory before submission

| # | item | file:line |
|---|---|---|
| a | `results/tables/selection_funnel.csv` stage label "tested in Stim48hr" is wrong. Paper p.36 reports 11,333 genes DE-tested at Stim48hr; the released h5ad has 11,281. Rest and Stim8hr are short by 2 each. Relabel to "present in released DE_stats". | `scripts/11_selection_funnel.py:138` |
| b | CD3E is `window_rank` 9, `safe = True`, `tolerance_loss` at the 94.7th percentile, just under the p75 gate. Fig 2F (p.10) shows CD3E knockdown reduces IL-10 **protein**. Disclose, or tighten the tolerance gate and re-run. | `results/tables/window_score.csv` |
| c | PPP3R1 `guide_correlation_signif` = −0.390, p = 1.17e-21, 523 DE genes, 811 cells, both guides on-target-significant. Not low power. The authors' reliability ceiling (p.38) shrinks the maximum *positive* r and cannot rescue a significantly negative one. Disclose. | shortlist table |
| d | `.obs` carries `chunk` = the authors' random 50-gene DESeq2 partition (p.37). We have never used or disclosed it. Permutation ICC on `window_score` = 0.006 (perm p = 0.055), on `efficacy` = 0.017 (perm p < 0.005). All 15 top safe genes sit in 15 distinct chunks. Disclose; a reviewer who reads p.37 will ask. | `src/cd4_perturbseq/de_stats.py` |
| e | `scripts/18_n6_selectivity_validation.py:87-88` asserts *"`low_target_gex` at rest mechanically entails `ontarget_significant == False`"*. **False.** 64 counterexamples genome-wide, 31 at Rest. The 4-row conclusion survives; the stated reason is wrong. | `18:87-88` |
| f | `distal_offtarget_flag` uniquely removes 81 perturbations that pass everything else (6,452 → 6,371). It is undefined anywhere in the paper. Run and report a sensitivity with it relaxed. | `scripts/11_selection_funnel.py` |
| g | Free validations to add: (i) `donor_correlation_hits_min > 0.35`, 18/20 of top 20 pass; (ii) authors' strong-perturbation set (`n_total_de_genes > 75 & n_cells_target > 50`), 20/20 of our top 20 and 91/100 of our top 100 fall inside it; (iii) IL10-regulator sign concordance, 7/7 positive-arm vs 1/5 negative-arm, Fisher p = 0.0101, which proves the 9-gene module is not an IL10 proxy. | new scripts |
| h | `scripts/10_tolerance_is_real.py:38` declares `FDR = 0.10` and never uses it. Dead constant. | `10:38` |

---

# 2. THE DEFINITIONS WE WERE GUESSING AT

## 2.1 The four QC flags

| flag | defined by the paper? | verbatim / recovered rule |
|---|---|---|
| `neighboring_gene_KD` | **YES**, p.37 Methods | *"To estimate off-target prevalence on proximal genes, we annotated each gRNA's distance to the nearest non-target TSS and assessed DE effects on that gene (Suppl. Table 3). Putative off-targets were flagged where we detected significant knock-down of the nearest TSS gene (**DE z-score < −1, FDR < 10%**) (Suppl. Table 5)."* Corroborated: 1,045 genes flagged in ≥1 condition vs the paper's implied 0.60 × 1,782 = 1,069 (p.5). |
| `low_target_gex` | **NO.** Column never named. | Recovered exactly: `low_target_gex == target_baseMean.isna() \| (target_baseMean < 1.0)`. **0 mismatches in all 33,983 rows.** Perfect separation at the boundary: max flagged baseMean 0.999978, min unflagged 1.000201. 5,850 rows have NaN `target_baseMean` (target unmeasured); all flagged. Consistent with p.36 Methods, *"we first removed expression outliers (mean counts > 10,000 or **< 1**…)"*. **Condition-specific**: 358 genes low at Rest but not Stim48hr, 96 the reverse. |
| `ontarget_significant` | **NO.** Column never named. Mechanism documented, definition not. | p.35: *"we calculated the expression of the targeted gene in cells with each guide and compared it to expression in NTC cells with a one-sided T-test, using the Benjamini-Hochberg correction… The knockdown was considered significant if FDR < 10% and t-statistic < 0."* That is **guide-level and preliminary**. The `.obs` column is a perturbation × condition flag and is provably **not** `guide_n_signif_ontarget > 0`: 4,612 rows are True with 0 significant guides, 5,038 rows have ≥1 significant guide and are False. p.37: *"We detected significant on-target knock-down in **73%** of tested condition-target pairs."* Our row-level rate is **62.45%** (26,177-row denominator gives 80.8%; no `.obs` subset gives 73%). **Never quote 73% as ours.** |
| `distal_offtarget_flag` | **NO. Nowhere in 63 pages.** | Only near-mention is p.5, speculative, not a rule: *"potentially reflecting either insufficient sensitivity to detect on-target effects or consistent off-target effects on distal genes."* Tested and rejected: `~ontarget_significant & n_downstream > 0` gives 7,536 rows, not the 433 actually flagged; 296 of the 433 have `ontarget_significant == True`. **We do not know what this flag is.** 433/33,983 rows; 127 at Stim48hr; removes 81 unique perturbations from our Stim48hr set. |

**Purpose.** The paper never presents these four flags as an analysis filter. Its own high-quality-perturbation filter is effect-size based: p.42, *"requiring >75 differentially expressed genes (DEGs) and >50 cells per perturbation. This yielded 3341 perturbations across three conditions (Rest: 1088; Stim8hr: 1136; Stim48hr: 1117)."* That reproduces **exactly** from `.obs` with `n_total_de_genes > 75 & n_cells_target > 50`. **The authors never filter on `ontarget_significant`.** Our mask is ours. Stop writing "standard QC". `docs/results/risk_kill_2026_07_08.md:97` says *"Standard QC was applied first"*. Change to "Our QC mask".

Attrition, Stim48hr, from `selection_funnel.csv`: 11,281 → 11,154 (distal, −127) → 10,324 (neighbour, −830) → 6,388 (`ontarget_significant`, −3,936, **80.2% of all attrition**) → 6,371 (low_target_gex, −17). **The collider is `ontarget_significant` alone.** `low_target_gex` is near-redundant given it. Do not repeat the "18% of tested perturbed genes" figure for `neighboring_gene_KD`; that parenthetical in p.5 attaches to the at-risk denominator, and our observed rate is 9.1% of genes / 7.4% of Stim48hr rows.

## 2.2 `n_total_de_genes` and `n_downstream`

**Neither column is named in the paper.** Both semantics are nonetheless pinned by three independent reproductions of the authors' own headline numbers.

Verified identity, 0 mismatches, 33,983 / 33,983 rows, all three conditions:
```
n_total_de_genes == n_downstream + ontarget_significant
n_total_de_genes == n_up_genes + n_down_genes
```

| paper statistic | page | reproduces on |
|---|---|---|
| *"7,807 (67%) perturbed genes had significant DE (FDR < 10%) effects on at least 3 genes in at least one condition"* | 5 | `n_total_de_genes >= 3` → we get **7,806**. On `n_downstream >= 3` we get 6,928. |
| *"we discovered 2,035,311 significant regulator-to-gene trans-effects"* | 5 | `sum(n_downstream)` → we get **2,035,203**. On `n_total_de_genes` we get 2,056,424. |
| *"median number of trans effects per perturbation with significant knockdown was 2 genes, while the top 5% of perturbations affect > 427 genes (mean = 81.61)"* | 6 | `n_downstream`, pooled, `ontarget_significant` only, n = 21,221 → median **2.00**, p95 **427.00**, mean **81.64**. On `n_total_de_genes`: 3 / 428 / 82.64. |
| *"3341 perturbations (Rest 1088; Stim8hr 1136; Stim48hr 1117)"* | 42 | `n_total_de_genes > 75 & n_cells_target > 50` → **exact**. On `n_downstream`: 3,321 / 1,081 / 1,133 / 1,107. |

Conclusion: `n_downstream` = **trans** effects, on-target excluded. `n_total_de_genes` = trans + cis. The authors use both, with the trans count for trans claims and the total for activity/eligibility filters. So do we, once §1.2 lands.

**Regulator**, verbatim, p.5: *"Throughout the manuscript, we use the term 'regulator' to indicate a perturbed gene with significant trans-effects on at least one measured gene."* By that definition our QC-passing Stim48hr set is 4,720 regulators and 1,651 non-regulators. Adopt the term. Note the paper's own "67%" statistic uses ≥3 **total** DE genes, not ≥3 trans effects. Do not conflate.

## 2.3 DE method and FDR

**Stated. No inference required.**

- Method, p.36: *"We performed differential expression (DE) analysis with the DESeq2 method [143], using the python implementation [144,145]."* Ref [145] = Muzellec et al., PyDESeq2. **The literal string "pydeseq2" does not appear in the text.** Do not quote it as verbatim.
- Model, p.36: negative-binomial GLM, `log μ_{j,p} = β₀ + log10(n_p)·β_{j,n} + Σ_d x_{d,p}β_{j,d} + Σ_k x_{k,p}β_{j,k} + log(s_p)`. Covariates: `log10` of aggregated cell count, donor identity, perturbation identity. `log(library size)` is an **offset**, no coefficient. No condition term.
- Test, p.36: *"we estimated the log2 fold-change … and test for differential expression using the **Wald test**, contrasting expression estimates for each perturbed gene against **non-targeting control pseudobulks**."*
- Fitting, p.36: *"For each culture condition, we fit DE models **independently**, to obtain context-specific perturbation estimates (β^Rest, β^Stim8hr, β^Stim48hr)."*
- Scalability, p.37: perturbed genes randomly partitioned into **groups of 50**, one GLM per group plus NTCs, per condition. `.obs["chunk"]` is that partition id. 681 chunks, exactly 50 per chunk, re-randomised independently per condition.
- Effect estimate, p.37: **z-score = β / SE(β)** is the default.
- **FDR = 0.10.** Stated in main text (p.5: *"FDR < 10%"*) and Methods (p.44: *"significantly differentially expressed (10% FDR)"*; p.41, p.35, p.37, p.39 all use 10%). Our `scripts/04_window_score.py:47` and `scripts/10_tolerance_is_real.py:38` are correct. Cite p.5.
  - Beware two 5% thresholds that are **not** the DE threshold: p.43 cluster annotation (`FDR < 0.05` hypergeometric), p.44 GO BP enrichment (`FDR < 0.05`). And p.25's *"FDR < 0.1"* is a cluster-enrichment FDR, a numeric coincidence. Never cite it as the DE threshold.
- Calibration, p.37: NTC-vs-NTC null. Five random splits, 2 NTC guides assigned to a simulated target group, cell counts restricted to the 5th–95th percentile of the NTC distribution, run through the identical model in all three conditions. Suppl. Fig 4B. **This is a type-I check only. It does not license using unfiltered Rest rows.**
- Gene universe, p.36: outliers removed (mean counts > 10,000 or < 1; % expressing pseudobulks > 99.9% or < 1%), then **top 10,000 HVGs per condition** via `scanpy.pp.highly_variable_genes(flavor='seurat_v3')` on pseudobulk profiles. Released `.var` is **10,282**, i.e. the union. The Rest and Stim48hr panels therefore differ by at most 282 genes (2.8%).
- Pseudobulk, p.36: unit = **donor × condition × guide**. Discard pseudobulks with < 5 cells or total mRNA < 0.5th percentile. *"Only perturbed genes represented by at least 3 pseudobulks across guides and donors were included in differential expression testing (N perturbed genes per condition: Rest = 11,289; Stim8hr = 11,417; Stim48hr = 11,333)."*
- **Released h5ad is short:** Rest 11,287, Stim8hr 11,415, Stim48hr 11,281. Deltas −2 / −2 / **−52**. 56 tested perturbations are absent from the release, 52 of them in our primary arm. **The paper does not say why.** Do not quote 11,281 as "tested".
- Cell QC, p.35: pre-assignment, drop > 20% mito or < 200 genes. Guide assignment by Poisson-Gaussian mixture (Replogle 2022 model, `crispat` implementation), one model per guide per sample. Post-assignment, drop > 5% mito, < 500 genes, `multi_sgRNA`, or no detectable gRNA. Per-guide QC in Suppl. Table 4.
- Guide efficiency, p.35-36: *"Guides for which no significant knock-down was detected in any condition (with at least 10 cells and mean log-normalized expression in NTCs > 0.01) were considered inefficient and excluded from differential expression analysis (**n=869**)."* Cross-condition, applied before DE.

## 2.4 The confidence columns

| column | defined? | what the paper says |
|---|---|---|
| `guide_correlation_signif` | **YES**, p.38 | *"we calculated Pearson correlations between **log fold changes** of guide pairs targeting the same gene. Correlations were computed using only genes showing significant differential expression (10% FDR) in **at least one** of the two guides, restricted to genes with **> 2 total DE genes**… Where computed, cross-guide correlation statistics are reported in Suppl. Table 5."* Guide-level DE was refit with the same DESeq2 framework, each guide as an independent perturbation. **Conflict:** `data/raw/data_sharing_readme.md:101,103` says these are Pearson r of per-gene DE **z-scores**, not log fold changes. Paper and readme disagree. |
| `guide_correlation_all` | **NO.** | Not defined. It is not the Suppl. Table 5 statistic. p.37 eligibility (≥3 biological replicates per guide; both guides testable; > 75 cells and > 30 significant DE genes at 10% FDR) admits only **2,118 guides ≈ 1,059 genes** at Stim48hr, yet `guide_correlation_all` is populated for **8,719** Stim48hr rows, including 559 rows with `n_guides == 1`. Its computation rule is undocumented. Median 0.065 across QC-passers vs 0.743 for `_signif`: it is a transcriptome-wide r swamped by nulls. **Report `_signif`, never `_all`.** |
| `donor_correlation_hits_{mean,min}` | **YES**, p.39 | *"we conducted 6 separate DE analyses… each time subsetting the expression data to include only two donors, representing all possible pairwise combinations of the 4 donors. We then calculated the cross-donor correlation by comparing tests performed on **disjoint donor pairs**… 1. Identify significantly DE genes in each test using a **lenient FDR threshold (20%)** 2. Take the **union** of DE genes from both tests being compared 3. Calculate the **Pearson correlation of log-fold changes** between the two tests for genes in this union 4. Compute the theoretical maximum correlation using the reliability estimate… For each target-condition pair, there are **3 possible comparisons** between disjoint donor pairs."* **Our C(4,2)=6 → 3 disjoint splits reverse-engineering is correct. Our "all genes" gene set is wrong.** The union at 20% FDR is the `_hits` variant. |
| `donor_correlation_all_{mean,min}` | **NO.** | Never defined. Empirically distinct: median `_all_mean` = 0.111 vs `_hits_mean` = 0.731; Pearson between them 0.444; `hits > all` on 94.8% of rows. Presumed all-measured-genes. **State as inference.** |
| eligibility for `donor_correlation_*` | **stated, but does not reproduce.** | p.39: *"restricted to condition-target tests meeting the following criteria in the full dataset: **>75 perturbed cells and >30 significantly differentially expressed genes at 10% FDR**."* The released mask is reproduced with **0/33,983 mismatches** by `n_cells_target > 50 & n_total_de_genes > 30`. The stated 75 leaves 259 counterexamples, all with `n_cells_target` in [51, 75]. Either "perturbed cells" ≠ `n_cells_target`, or the paper misstates the threshold. Do not cite 75 as verified. |
| `single_guide_estimate` | **NO.** Never named. | Exactly `n_guides == 1`. 100% agreement, all three conditions, 3,859 / 33,983 rows. Cause traceable to the p.35 n=869 inefficient-guide exclusion and the p.36 ≥3-pseudobulk floor. **But that story only fits CD3E** (`n_guides` = 1 in all three conditions). COMMD5 has 2 guides at Stim8hr and 1 at Stim48hr, so its Stim48hr status comes from a **condition-specific** pseudobulk filter, not the global 869. |
| `n_guides` | **NO.** | `{1: 3859, 2: 30108, 3: 13, NaN: 3}`. Also condition-specific: 404 genes differ across conditions. |
| reliability / correlation ceiling | **YES**, p.38 | `σ²_error = mean_j(SE(β_j))²`; `σ²_true = var_j(β_j) − σ²_error`; `reliability = σ²_true / (σ²_true + σ²_error)`; *"The correlation ceiling was then calculated as the geometric mean of reliabilities for the two guides being compared."* Reused for cross-donor (p.39, step 4). **Not computable from the released `.obs`.** No SE, no reliability, no ceiling column. It requires per-guide DE SEs, released only in `GWCD4i.DE_stats.by_guide.h5mu`. Reported in Suppl. Table 5. |

**Guide discordance.** The paper's only statement is p.5: *"Perturbation effects were generally consistent between guides targeting the same gene (median cross-guide correlation: Rest = 0.49, Stim48hr = 0.52, Stim8hr = 0.47), with most discrepancies attributable to differential knockdown efficiency between guides… or to one guide producing substantially larger downstream transcriptional effects."* And p.30 Limitations: *"future work would still be required to fully assess potential off-target guide effects with additional guides per target gene."* **No threshold. No policy for a negative r.** PPP3R1 at −0.390 has no authorial remedy to cite. The ceiling would make it worse, not better: reliability ∈ (0,1], so dividing a negative r by a positive number < 1 pushes it further from zero.

**Supplementary CSV.** `data/raw/suppl/DE_stats.suppl_table.csv` (mtime Dec 23 2025) has 33,983 rows × **16** columns. The h5ad `.obs` has 33,983 × **27**. Absent from the CSV: all 14 QC/confidence columns including `distal_offtarget_flag`, `neighboring_gene_KD`, `low_target_gex`, `n_guides`, `single_guide_estimate`, and all six correlation columns. It carries two the h5ad lacks (`offtarget_flag`, `ontarget_effect_category`). It is an **earlier snapshot with partly renamed columns**, and it does not contain the statistics that Fig 2C's legend (p.11) and Methods pp.37-39 attribute to Suppl. Table 5. Our factual statement is correct. Our **label** is sloppy: name the file, do not assert what Suppl. Table 5 contains. Fix `src/cd4_perturbseq/paths.py:26`.

## 2.5 The design, verbatim

| fact | value | source |
|---|---|---|
| donors | 4, healthy, PBMC-enriched leukopaks (STEMCELL #70500) | p.3, p.32 |
| cell type | **naive CD4+ T cells**, EasySep Human Naive CD4+ T Cell Isolation Kit #19555. Not bulk. Not sorted Tconv. Tregs depleted by negative selection at day 0. | p.32, Fig 1A p.7 |
| cells profiled | 22 million | p.3 |
| library | 12,748 perturbed genes, 26,504 gRNAs, "(2 gRNAs/gene)" nominal | Fig 1A panel p.7; Methods p.31 |
| library design | union of all expressed CD4 genes + all TFs from Lambert 2018; guides cross-referenced against hCRISPRiv2 and Dolcetto; *"5% of non-targeting gRNAs from the Dolcetto dataset"* (actual 992 = 3.74%) | p.31 |
| CRISPRi effector | **ZIM3-KRAB-dCas9-BlastR** (plasmid `EF1a-Zim-3-dCas9-P2A-BSD`, pZR316). Not dCas9-KRAB. | p.32 |
| MOI | ~0.2 | p.32 |
| medium | cXVIVO = X-VIVO 15 (Lonza 04-418Q) + 5% heat-inactivated FCS | p.32 |
| **IL-2** | **200 IU/mL, day 0 through harvest**, re-supplied at every media change (d0, 3, 5, 7, 8, 10, 12) | p.32-33 |
| **IL-7** | 5 ng/mL from day 5. **Absent from the day-12 split medium**, therefore absent from the entire DE window. | p.32-33 |
| polarising cytokines | **none.** No IL-12, IL-4, IL-6, TGF-β, IL-23, IFN-γ. Confirmed twice: p.18, *"The population of perturbed CD4+ T cells was not cultured in polarizing conditions"*; p.30, *"Our study focuses on a non-polarized culture condition."* | p.18, p.30 |
| stimulation | **ImmunoCult human CD3/CD28/CD2 T cell activator** (STEMCELL 10990). Day 0: **25 µL/mL**. Day 12 restim: **12.5 µL/mL**. Soluble tetramer. Not beads. Not plate-bound. | p.32-33 |
| selection | blasticidin 10 µg/mL + puromycin 2 µg/mL from day 3 | p.32 |
| timepoints | day-12 morning split into three: **Rest** = no restim, harvested +8 h. **Stim8hr** = restim, +8 h. **Stim48hr** = restim, +48 h. | p.33 |
| assay | probe-based scRNA-seq, 10X Flex; 768 out-of-library gRNA probes as a detection-specificity control (< 0.2% of cells) | Fig 1A p.7; p.33 |
| depth | mean UMI/cell: **Rest 10,080; Stim8hr 14,977; Stim48hr 13,373** | p.4 |
| coverage | ~575 cells per perturbed gene on average; Fig 1C medians 500 / 509 / 494 | p.4, p.7 |
| batches | batch 1 = Rest D1/D2 + Stim8hr D1/D2; batch 2 = Rest D3/D4 + Stim8hr D3/D4; **batch 3 = Stim48hr D1-D4**. Loading 2.5× / 2× / 1.7× overloaded. Condition is confounded with batch for Stim48hr. Authors report *"no substantial batch effects between 10x lanes or experimental batches"* on 395,030 NTC cells by scVI. | p.33-35 |
| on-target KD | *"target gene expression was significantly reduced for 73% of tested guides"* (24,174 detected guides) | p.4-5 |
| cross-platform | 27 of 32 tested targets correlate with arrayed KO screens, **using the Rest condition** | p.5 |

---

# 3. THE RESTING ARM

**Are Rest and Stim DE statistics comparable? As effect estimates, yes. As DE-gene counts, no.**

## What the paper says

- **Fitted separately.** p.36: *"For each culture condition, we fit DE models **independently**, to obtain context-specific perturbation estimates."* Design formula has no condition term. Each perturbation is contrasted against **the same condition's** NTC pseudobulks. Therefore β^Rest and β^Stim48hr are on the same scale, log2 versus within-condition NTC. Nothing in our code ever differences a Rest `log_fc` against a Stim `log_fc`, and nothing should.
- **The authors use Rest at face value.** p.5: *"We compared the DE effect estimated by our perturb-seq screens with DE effects measured in published arrayed knockout screens in CD4+ T cells (**using the Rest condition, to match the cellular context**)."* p.19: Rest gives the best cross-validated fit to the OneK1K aging signature of any arm (R = 0.366) *"with significantly better fit than perturbation effects in stimulated CD4+T conditions or in K562 cells"*. p.20: *"the effects of these regulators on these processes appear to be **specific to the Rest condition** and are not detected post-stimulation."*
- **Rest is not underpowered.** Median cross-guide correlation Rest 0.49 vs Stim48hr 0.47 (p.5). Median `n_cells_target` 543 vs 548. Paired median log2(stim/rest) cells = −0.007. Total DE calls ratio 1.20×. Median `n_total_de_genes` = 2 in every arm.
- **Rest also carries less trait signal.** p.24: autoimmune genes enrich in Stim8hr/Stim48hr regulator clusters *"but not in clusters whose regulators only show correlated effects in Rest condition"*. p.22: regulator-burden correlations show *"we did not find higher than expected correlations in either the rested T cells, K562 cells or in a screen of essential genes in Jurkat cells"*. Both **support the sign** of our selectivity axis.
- **Per-condition QC: the paper does not say.** No per-condition QC is prescribed anywhere. The four flags are nonetheless empirically condition-specific: across 11,102 genes with both rows, Rest and Stim48hr disagree on `low_target_gex` for 454, `ontarget_significant` for 1,074, `neighboring_gene_KD` for 214, `distal_offtarget_flag` for 203.
- **Power/cell-number difference the authors flag:** none between arms. The flagged asymmetry is **sequencing depth** (p.4, above) and **direction of effect** (p.9).

## Why the DE-count ratio is not comparable

Four independent confounds, all sourced:

1. **Depth.** Rest 10,080 vs Stim48hr 13,373 mean UMI/cell (p.4). Shallower arm → fewer Rest DE calls → selectivity inflated. Our gate is anti-conservative.
2. **Panel.** Top 10,000 HVGs selected **per condition** (p.36); released `.var` is 10,282. The Rest and Stim48hr denominators can differ by up to 282 genes.
3. **Floor effect.** p.9: *"For stimulation-responsive cytokines like IL-2 and IL-13, negative regulators (knockdown log2FC > 0) were detected predominantly in the Rest condition, where cytokines exhibit minimal constitutive expression. In contrast, positive regulators (knockdown log2FC < 0) were more readily identified following re-stimulation (e.g., CD3 and LCP2 for IL2; GATA3 for IL13)."* Our showcase examples are exactly that class. `scripts/04_window_score.py:20-25` calls CD3G (1,374 / 5), ITK (2,566 / 2), PLCG1 (2,218 / 3) *"the therapeutic window in its purest form"*. For proximal TCR genes, selectivity is confounded with detectability. **Delete that phrase.**
4. **Arithmetic.** Stim always carries the on-target +1; Rest carries it for only 5,683 of 6,320. §1.2.

Counter-evidence that limits confounds 1 and 3: `rest_de_genes` is direction-agnostic across all 10,282 genes, and the floor cannot explain GATA3 at 888, NSD1 at 1,767, STAT5B at 1,135 resting DE genes. The axis is not vacuous. It is biased, in one direction, in our favour.

## Plus our own construction defects

- Rest is joined **unfiltered** (`04:180`), so 802 of 6,320 denominators come from rows we would reject.
- 51 genes with no Rest row are **median-imputed at exactly 3.0** (`04:182`), and absence from the Rest arm means the gene failed the ≥3-pseudobulk floor, i.e. depletion. We impute depletion as inertness.
- 84 of 370 selectivity-gate passers (22.7%) rest on an imputed or QC-failing denominator.
- The word "quiescent" in `04:190-194`, `04c:3`, `07:3-5` is false. Rest is an 8 h no-restim aliquot of a day-12 IL-2 blast pool (p.33). `rest_cells_ratio` is not an independent viability signal; spearman with Stim `n_cells_target` = 0.94, and all dropout accrued during the shared 12-day expansion.

## Is the N6 demotion for the right reason?

**Yes, and independently confirmed by our own pre-registered null.** `results/tables/reversal_specificity.csv`: *"rejected on homeostasis, top 20"* observed 14 vs magnitude-matched mean 13.625, **p = 0.5296, specific_to_reversal = False**. *"rejected on tolerance, top 20"* observed 12 vs 3.9918, **p = 0.0, specific_to_reversal = True**. The homeostasis axis was already non-specific before we read the paper.

**But state the reason correctly.** Do not write "Rest and Stim DE statistics are not comparable". They are, as within-condition NTC-relative effect estimates, and the authors use Rest as a first-class arm. Write: *the selectivity axis is a ratio of DE-gene counts across two independently fitted models with different sequencing depth (p.4), different per-condition HVG panels (p.36), a documented detectability floor for stimulation-responsive transcripts (p.9), an asymmetric on-target count, and a denominator that is unfiltered and partly imputed. Our own magnitude-matched null already shows it is not reversal-specific (p = 0.5296). We demote it.*

Also state in Methods, to pre-empt the reviewer who reads p.9: **our objective is Stim48hr only. We never compute a Rest-minus-Stim `log_fc` or a sign flip.** The Rest arm enters as two `.obs` scalars.

---

# 4. GENE PROGRAMS

## Do the authors publish modules we should be using instead?

**No. Their "programs" live on the perturbation axis. Ours live on the gene axis. They are not substitutable, and the one authorial object on our axis is circular.**

### What they actually publish

| object | what it is | where |
|---|---|---|
| **111 regulator clusters** | Clusters of **perturbations**, not genes. Consensus Leiden (27 hyperparameter combos × 400 bootstraps = 10,800 clusterings) → co-occurrence distance → HDBSCAN (`min_cluster_size=4, min_samples=1, eom`). Input: the z-score B̃ matrix over 3,341 strong perturbations (`>75 DEGs and >50 cells`), 10,941 HVGs from 13,959, **on-target z-scores masked to zero**. | main text p.11; Methods p.42; **Suppl. Table 9** |
| cluster annotations | Hypergeometric enrichment of the clustered **regulator** genes against CORUM, STRINGDB, KEGG, Reactome (terms < 200 genes), BH, `FDR < 0.05`. Unenriched clusters fell back to *"Gene Ontology analysis, LLM lookup, and manual literature search"*. Names are complex/pathway names: "Integrator complex", "TCR signaling", "Core BAF". | Methods p.43 |
| condition specificity | CV of intra-condition mean pairwise regulator correlations. `CV < 0.5` → `across_condition`. Otherwise condition c enriched iff `ρ_c / mean(M) > 1.2` **and** `ρ_c > 0.2`. Labels `Rest`, `Rest_Stim8hr`, etc. | Methods p.43 |
| **downstream gene sets** | Per cluster **per condition**: *"the union of all genes significantly differentially expressed (10% FDR) by at least one regulator within that cluster in the corresponding condition."* Ranked by upstream regulator count, sign coherence, aggregate rank score. GO BP 2025 enrichment via gseapy, `FDR < 0.05`. | Methods p.44; **Suppl. Table 10** |
| activation regulator clusters | Three clusters, **291 regulators**: early **97**, general **46**, late **148**. Definitions p.14: early = Stim8hr only; late = Stim48hr only; general = both stimulation timepoints but not Rest. | Fig 3E p.13; p.14; Suppl. Table 9 |
| **9 core regulators** | CD3D/E/G, CD247, ZAP70, LAT, LCP2, PLCG1, VAV1. Used as a correlation reference axis. | Fig 3E legend, p.14 |
| 30-cytokine panel | IL32, LTB, VEGFA, CCL3, CCL4, CSF2, IL3, IFNG, IL2, TNF, TNFSF14, CD40LG, IL23A, IL10, CCL5, CXCL8, CSF1, IL24, TGFB1, IL21, IL22, LTA, IL13, IL4, IL5, MYDGF, VEGFB, TNFSF10, IL16, TNFSF9. A **measurement panel**, Rest and Stim8hr only, not a directional program. | Fig 2A y-axis, p.10 |
| Th2/Th1 signature | Bottom-50 / top-50 by z-score of DE log fold change at 1% FDR in the Ota 2021 discovery cohort (n=79); replication Hollbacher 2021 (n=3). | Methods p.48; Suppl. Table 11 |
| aging signature | OneK1K, DE with age across CD4 subsets, 782 discovery / 199 held-out donors. | p.19; Methods p.49 |
| TCR activation signature | *"a TCR activation signature derived from DESeq2 results comparing stimulated and resting Teff cells from [13]"*, used as a **pert2state negative control**. Ref [13] = Arce et al. | Methods p.48 |

### The three reasons we cannot switch

1. **Wrong axis.** Suppl. Table 9 clusters 1,699 perturbation targets (`.obs` rows). Our modules are readout genes (`.var` columns). Overlap of the 61 Fig 3E named regulators with `effector_core()` = **0**. With `tolerance_module()` = **0**.
2. **Circular.** Suppl. Table 10's downstream gene sets are the union of genes DE at 10% FDR **under the same z-score matrix we score against**. Using it as our objective makes the objective a function of the data. That is exactly the failure `scripts/02b_check_circularity.py` exists to bound.
3. **It does not contain what we need.** `data/external/gwt_priors/clustering_results_and_annotations.csv` is already on disk: 112 rows, ids 0-111 (the paper says 111; the released table is off by one). `manual_annotation` value counts: 62 are `"unknown"`. A case-insensitive regex for `treg|regulatory|foxp3|inhibit|toleran|exhaust|effector|checkpoint` returns exactly one hit, cluster 40, *"Protein phosphatase inhibitor 2"*, a false positive on the substring. **There is no effector program and no Treg or co-inhibitory program among the 111.** Fig 3B is explicitly *"a representative set"* and *"a subset of clusters"*, so scope the negative claim to the annotated clusters unless we open the full table.

Also: their clustering input is 10,941 genes from 13,959; our released `.var` is 10,282. **The B̃ matrix they clustered is not the released DE_stats gene space.** We cannot reproduce the 111 clusters and must not claim to.

## Are our sets ad hoc?

**The effector module: no. The co-inhibitory module: yes, and we should say so.**

`scripts/01_build_activation_program.py:42-45` builds the program from two sources:
- `curated` = `data/external/gwt_priors/immune_effector_genes.csv`, fetched by `scripts/fetch_priors.sh:31` from **the authors' own analysis repo**, `emdann/GWT_perturbseq_analysis_2025@master:metadata/immune_effector_genes.csv`. Verified: **117 rows, 117 unique symbols**, categories `Receptor 46 / Cytokine 36 / TF 33 / Others 2`. The docstring's "117" is correct.
- `external` = the Arce 2024 `AAVS1_Teff_Stimulation_vs_Resting` DESeq2 table, `padj < 0.05 & log2FC > 1`. **That is the paper's own ref [13] activation-signature source** (Methods p.48).

`activation_program.csv` = 1,640 genes; `{external: 1523, curated: 77, both: 40}`.
`effector_core()` = the 40 "both" genes minus the tolerance set = **34 genes**. After intersecting with the 10,282 measured genes it is **32** (CCL20 and IL17A are not measured). **Both numbers are right. Our prose says 32 without saying "measured", and the code returns 34. Fix the prose or pin the count in a test.**

So the effector module is *the authors' released curated list intersected with the paper's own external activation contrast*. **Stop calling it "hand-built".**

`tolerance_module()` = 9 genes: CTLA4, FOXP3, IKZF2, IL10, LAG3, LRRC32, PDCD1, TGFB1, TIGIT. `TOLERANCE_GENES` declares 12; ENTPD1, NT5E, HAVCR2 fall out on intersection. **This one is genuinely ours.** No authorial counterpart.

The authors' own cytokine taxonomy **endorses our split**. p.8: *"pro-inflammatory cytokines (IFN-γ, TNF), Th2 cytokines (IL-4/5/13), **regulatory cytokines (IL-10, TGF-β)**, chemokines (CCL3/4/5, CXCL8)."* Our effector core draws from the first, second and fourth bins and excludes the third. That exclusion **is** the therapeutic window. Cite p.8.

## Cost to switch: zero benefit, high cost. Do not.

**Do instead, all cheap and additive:**
1. Cite Fig 3 / Suppl. Table 9 as **regulator clusters**, never as gene programs. Cite Fig 3E / p.14 for the temporal taxonomy. Cite Suppl. Table 10 as per-cluster downstream sets and state the circularity objection in one sentence.
2. Join `clustering_results_and_annotations.csv` onto our shortlist as a **descriptive annotation only**. 45 of the 61 Fig 3E symbols already join `window_score.csv`. It costs one join and no h5ad reads. **We fetched Suppl. Table 9 and never opened it. That is embarrassing in review.**
3. Report the recovery statistic: **12 of the 15 canonical activation regulators named on p.12 survive our QC.** CD3D, LAT, LCP2 do not. State that our mask is stricter than theirs (§1.8).
4. Add a sensitivity in `scripts/17_tolerance_is_special.py` rerunning the module-level null with the effector set replaced by the paper's Fig 2A 30-cytokine panel minus IL10 and TGFB1 (28 genes, all measured). If the separation holds on the authors' own panel, the hand-built module stops being an attack surface.
5. Fix `scripts/01_build_activation_program.py:8`: the 117-gene list is the authors' **repo metadata file**, not the paper's Figure 3 gene programs. Say so.
6. Explain the 10,941 vs 10,282 measured-gene gap before publication.

---

# 5. FOXP3 AND WHAT WE MAY CALL OUR MODULE

## What the paper says about FOXP3, Tregs, and regulatory programs

Full-text counts, all 63 pages:

| term | hits | where |
|---|---|---|
| `FOXP3` | **1** | reference [17] title: *"Umhoefer JM, … FOXP3 expression depends on cell-type-specific cis-regulatory elements and transcription factor circuitry. Immunity. 2025"* |
| `Treg` | **0** | — |
| `"regulatory T"` | **1** | reference [40] title (Hollbacher) |
| `suppressive` | **0** | — |
| `"co-inhibit"` | **0** | — |
| `PDCD1`, `TIGIT`, `IKZF2`, `TGFB1` | **0 each** | — |
| `CTLA4` | 1 | p.25 main text |
| `LAG3` | 1 | p.14 main text |
| `LRRC32` | 1 | p.16 main text |

FOXP3 appears in the paper **only as a rendered axis label** inside Figure 1F (p.7), where it is one of five validation readouts against a published CRISPRi FACS screen (`FOXP3 / CD4+ Stim48hr / Umhoefer 2025`). It never enters prose. **The paper makes no claim, either way, about FOXP3 in these cells.**

What the paper does say, and it all helps us:

- **p.14, main text:** *"Eight hours post-stimulation, perturbations of early and general regulators induced changes of **early activation markers including IL2RA and LAG3**."* The authors call LAG3 an **activation marker**. In-paper support for the word "activation-induced".
- **p.16, main text:** *"whereas STAT6 regulated GATA3 and TGF-β pathway components (SMAD3, SMAD2, TGFBR2, **LRRC32**; FDR = 0.04)."* LRRC32 is a receptor/signal-transduction component of a Th2 polarisation signature. Not a Treg gene. Note their TGF-β set **excludes TGFB1**, so they do not bundle the ligand with GARP.
- **p.25, main text, Fig 7C cluster 80:** *"This STAT3-related network positively regulated autoimmune disease-associated genes encoding pro-inflammatory cytokines (IFNG, CXCL8), signaling proteins (TNFAIP3, TNIP2, STING1, SPSB1), transcriptional regulators (FOSL1, Polycomb components), and **immune receptors (CTLA4, LRBA)**."* CTLA4 is filed under "Receptors", a separate box from "Cytokines". The paper never calls CTLA4 inflammatory. It reports that CTLA4 is **positively co-regulated with IFNG and CXCL8 by a shared regulator network**. That is exactly the co-induction confound we already model, and we already quantify it: `results/tables/tolerance_induction_verdict.csv` gives `coinduction_share_of_scripts10_effect = 0.341`, `T_real = 2.600`, `p_module_level = 0.004975`, positive control `T = 10.48`. Note this is a **Stim8hr** result and our score is Stim48hr.
- **p.8:** IL-10 and TGF-β are *"regulatory cytokines"*. **p.9:** *"IL10 is a potent anti-inflammatory cytokine gene with pleiotropic functions in suppressing excessive immune responses."* These are the only two of our nine that carry authorial "regulatory" framing, and it is **cytokine pharmacology, not lineage identity**.
- **p.41, Methods:** the paper's **only protein-level assay** is intracellular flow for IL-10 and IL-21 after a 6 h Activation Cocktail with Brefeldin A restimulation. **No suppression assay. No CD25/CD127 sort. No Treg phenotyping anywhere in 63 pages.** And one published discordance: p.11, *"A notable exception was observed in the regulation of IL10 by NFKB2, where protein level changes did not mirror the transcriptional phenotype."*
- **p.32:** the cells are naive-isolated. Thymic Tregs are depleted at day 0 by negative selection. **FOXP3 detected at 48 h post-restimulation is therefore activation-induced in conventional T cells by construction, not carried in by a contaminating Treg fraction.** This is the strongest single sentence available to us, and it comes from the authors' Methods.

## What our own data says

`results/tables/tolerance_per_gene_induction.csv`, all nine:

| gene | `induction_log2fc` | `induced` | `mwu_p` |
|---|---|---|---|
| IL10 | +4.734 | True | 1.79e-16 |
| LRRC32 | +3.961 | True | 0.4303 (n.s.) |
| FOXP3 | +3.038 | True | 3.43e-07 |
| LAG3 | +2.777 | True | 9.59e-08 |
| CTLA4 | +2.294 | True | 1.55e-23 |
| PDCD1 | +1.577 | True | 2.41e-17 |
| IKZF2 | +1.195 | True | 5.48e-03 |
| **TIGIT** | **−0.042** | **False** | 8.15e-10 |
| **TGFB1** | **−0.539** | **False** | 0.1288 (n.s.) |

**7 of 9 induced. Two are not.** "Activation-induced" over-claims for TIGIT and TGFB1. "Co-inhibitory" over-claims for IL10, TGFB1, FOXP3, IKZF2, LRRC32.

And: `results/tables/il10_regulator_concordance` (to be committed): against the paper's arrayed-validated IL10 regulators, our `tolerance_loss` is sign-concordant **7/7 on the positive arm** (KDM1A, NFKB2, GATA3, CD3E, ZAP70, LCK, MED12) and **1/5 on the negative arm** (MEN1, SGF29, ATXN7L3, USP22, ELOB). Fisher p = 0.0101. **The module is not an IL10 proxy.** Commit this. It is the cheapest strong validation available against external, protein-level ground truth.

## The name

**Forbidden, permanently:** `tolerance`, `Treg`, `Treg-like`, `regulatory program`, `suppressive`, `regulatory T cell identity`. The paper supports none of these and never uses them of these cells. `src/cd4_perturbseq/programs.py:17` currently asserts two of them.

**Also forbidden without qualification:** the bare phrase "9-gene module". The paper owns "9 core regulators" (CD3D/E/G, CD247, ZAP70, LAT, LCP2, PLCG1, VAV1; Fig 3E legend p.14). In the same manuscript an unqualified "nine-gene module" reads as theirs. Always write **"downstream"**.

**Recommendation: keep the pre-registered name, add two footnotes.** `docs/preregistration_n9_2026_07_08.md:29` already locks `activation-induced co-inhibitory module` and line 121 already promised the global rename. Changing the name now is itself a preregistration deviation. So:

> **activation-induced co-inhibitory module** (nine downstream genes: CTLA4, PDCD1, LAG3, TIGIT, FOXP3, IKZF2, LRRC32, IL10, TGFB1)
> Footnote 1: seven of nine are stimulation-induced in this dataset; TIGIT (−0.04) and TGFB1 (−0.54) are not.
> Footnote 2: only IL10 and TGFB1 carry authorial framing, as *"regulatory cytokines"* (Zhu, Dann et al., p.8). The other seven are our addition and correspond to no program the authors define. LAG3 is called an *"early activation marker"* by the authors (p.14).
> Footnote 3: FOXP3 here is a transcript in naive-isolated, 12-day-expanded, restimulated conventional CD4 T cells (Methods p.32). It is activation-induced without conferring suppressive function (Wang 2007 *Eur J Immunol*; Tran, Ramsey & Shevach 2007 *Blood*; Allan 2007 *Int Immunol*). Cite those, never this paper.

If the naming has to be defended harder, the maximally conservative fallback is **"nine-gene downstream immunoregulatory module"**, which describes the transcripts and claims nothing about cell state, induction, or function. Prefer the prereg name; it is defensible with the footnotes.

**Code changes, mandatory:**
- `src/cd4_perturbseq/programs.py:17` delete the comment. Replace with: *"Activation-induced co-inhibitory and immunoregulatory transcripts. Names describe transcripts, not a cell state. IL10 and TGFB1 are the authors' 'regulatory cytokines' (p.8); LAG3 is an early activation marker (p.14). No suppressive function is claimed or measured. FOXP3 is activation-induced in Tconv, not a Treg-identity call."*
- `programs.py:3-6` module docstring: drop *"sparing the tolerance program"*.
- `programs.py:33` docstring: drop *"tolerance"*.
- Alias `TOLERANCE_GENES = COINHIBITORY_GENES` and `tolerance_module = coinhibitory_module` so the 7 consuming scripts keep running. **Do not rename the six committed `tolerance_*.csv` files.** Add one README line stating that `tolerance_*` is a legacy prefix carrying no claim of suppressive function.
- `README.md:133-134`, `docs/results/risk_kill_2026_07_08.md:160-161`, `docs/handoffs/HANDOFF_02_reviewer_verify_riskkill.md:85`: prose only.
- `scripts/17_tolerance_is_special.py:26,569` already says the right thing but says "bulk" and "anti-CD3/CD28". Fix per §1.5.

**One more thing to run:** `scripts/17_tolerance_is_special.py` has `--n-null` and `--n-modules` but no `--drop-gene`. Add it and rerun the N9 test three times: all 9, minus IL10, minus IL10+TGFB1. IL10 has by far the largest induction (+4.73) and the smallest p. **If the PASS at prereg line 293 does not survive dropping IL10, the N9 result is an IL10 result and must be reported as such.** This is the single most likely reviewer attack on §5 and we have not run it.

---

# 6. NOVELTY CHECK

## Does the paper do druggability, therapeutic window, drug-target recovery, or safety triage?

**No. None of the four. Anywhere, including supplementary.**

Full-text counts over all 63 pages:

| term | hits |
|---|---|
| `druggab*` | **0** |
| `tractab*` | **0** |
| `"therapeutic window"` | **0** |
| `"approved drug"` | **0** |
| `"clinical trial"` | **0** |
| `max_phase`, `"drug target"`, `triage` | **0** |
| `"chemical probe"`, `PROTAC`, `degrader`, `ligandability` | **0** |
| `drug` | 2, **both in reference titles** |

No tractability tiering. No mechanism-of-action annotation. No clinical-phase annotation. No direction-of-effect gate. No approved-drug enrichment. No safety, toxicity, or on-target-liability analysis. No ranked target list.

## What it DOES do, and we must cite

1. **Therapeutic framing, four sentences, all aspirational, none operationalised:**
   - p.15: *"Furthermore, regulators driving disease-associated or treatment-responsive states could serve as therapeutic targets."*
   - p.24: *"Genetic association studies have identified hundreds of loci linked to autoimmune diseases, yet translating these findings into actionable therapeutics requires understanding the causal genes, the cell types in which they exert their effects, and the precise nature of their molecular effects."*
   - p.29: *"Our analyses showcase how this genome-scale perturb-seq map … can be harnessed for insights into human immunology, immunotherapy design, dynamic gene regulatory control of human cells and human disease genetics."*
   - p.30: *"…in order to harness these data for immunology and immunotherapy development."*

2. **Open Targets, twice, and only for genetic evidence.** p.24: *"We annotated a set of genes associated with 14 autoimmune conditions (hereafter autoimmune genes) from Open Targets [48]."* Methods p.50: *"We queried the OpenTargets Platform API (v4) to retrieve genes with genetic association with autoimmune diseases … We included genes with a minimum **genetic evidence score of 0.1** … derived from GWAS associations, gene burden studies, and somatic mutation data, which includes evidence from ClinVar."* Plus three non-immune negative controls: coronary artery disease, macular degeneration, chronic kidney disease. Second citation, p.16: GCSAM / eosinophil counts.
   Our `scripts/15_fetch_open_targets.sh` pulls the **bulk 26.06 release** `molecule`, `drug_mechanism_of_action`, `clinical_indication`, `target`, `disease` parquets. `scripts/16_open_targets_benchmark.py` gates on `actionType` and per-indication `APPROVAL`. **Same database, disjoint endpoints.** Say so explicitly, because a reviewer greps for "Open Targets", hits p.24, and concludes we did not read it.

3. **Regulator-burden correlation**, their own prioritisation method. p.22: *"Specifically, to identify putative core genes, we tested for each gene j whether there is a significant correlation between the knockdown effects of regulators on j and the LoF effects of those regulators on lymphocyte count in the UK Biobank. We refer to the correlation between LoF effects and knockdown effects for each putative core gene as the regulator-burden correlation (Figure 6A)."* Prioritises **downstream measured genes**. We prioritise **perturbations**. Median `window_rank` of their 15 rankable named core genes = 3,582 of 6,371; **0 of 15 pass our safety gate**; 14 of 15 are rejected on homeostasis. Report as orthogonality.

4. **HPA tissue specificity**, descriptively. Methods p.45: same HPA v25.0 consensus table we use, Tau > 0.5 and per-tissue log2FC(nTPM vs median over 51 tissues) > 2.5, used to annotate regulator clusters for Suppl. Fig 13. **Never paired with constraint. Never inverted into a systemic-risk gate. Never excludes lymphoid tissue as on-target.** Our `priors.hpa_tissue_breadth()` counts non-immune tissues with nTPM ≥ 1 and excludes six lymphoid tissues. Different statistic, different purpose. Cite them, differentiate.

## The defensible novelty sentence

Replace `docs/strategy_2026_07_07.md:91-92` and `docs/claude_tooling_log.md:17` with:

> The paper never uses druggability or tractability language. `druggable`, `druggability`, `tractability`, `therapeutic window`, `approved drug`, `clinical trial` appear zero times across 63 pages. It invokes therapeutic framing four times, aspirationally and without analysis: *"regulators driving disease-associated or treatment-responsive states could serve as therapeutic targets"* (p.15), *"translating these findings into actionable therapeutics"* (p.24), *"immunotherapy design"* (p.29), and *"immunology and immunotherapy development"* (p.30). It performs no druggability annotation, no tractability tiering, no therapeutic-window analysis, no drug-target-recovery benchmark, no approved-drug enrichment, and no safety triage. Its Open Targets use is confined to the genetic-evidence score (≥ 0.1, from GWAS, gene burden, somatic mutation and ClinVar; p.24, Methods p.50), disjoint from the mechanism-of-action and clinical-indication tables our benchmark reads. It does prioritise candidate genes, via regulator-burden correlation against UK Biobank LoF burden on lymphocyte count (p.22, Fig 6A), but on the downstream-gene axis, not the perturbation axis. That is the gap we occupy.

**Leave `docs/research_plan.md:21` unchanged.** All three of its sub-claims verify.

## Their limitations, verbatim, in full

Main text, printed p.30, section heading **"Limitations of current study"**. One paragraph, five sentences, quoted complete:

> *"Despite performing large-scale perturb-seq with 2 gRNAs per target across multiple human blood donors, future work would still be required to fully assess potential off-target guide effects with additional guides per target gene and to differentiate donor-specific biology from experimental batch effects. Our analysis uses pseudobulk aggregation to compare mean expression profiles across cells with the same perturbation. This enhances statistical robustness, but potentially masks heterogeneity of response to perturbation across cells, which could be resolved by future distribution-aware analyses. Our study focuses on a non-polarized culture condition; consequently, the regulatory rewiring driven by polarizing cytokines - central to CD4+ T cell plasticity - remains unmapped. Future work extending this framework to include diverse in vitro polarization contexts and complex in vivo signals is essential to generate a comprehensive model of context-dependent immune regulation. Although perturb-seq elucidates causal gene regulatory relationships, future experimental and computational efforts will be required to link the transcription effects of each perturbation to specific immune cell functions in order to harness these data for immunology and immunotherapy development."*

Note the paper writes *"transcription effects"*, not "transcriptional effects".

**Our claim that their roadmap is "all tech-scaling" is wrong.** Four of the five items are scaling asks. The fifth is a **translation ask**, and it explicitly names immunotherapy. It is stated as an **open gap** ("will be required"), not as work performed. It names function-linking, not target selection, therapeutic window, drug-target recovery, or reversal.

**This is a weapon, not a wound. Quote it and answer it.** Our partial answer already exists and is buried: `scripts/02c_validate_efficacy_axis.py` validates the naive suppression score against Schmidt & Steinhart 2022, an independent genome-wide CRISPRi screen with an **IL-2 protein** readout, AUROC 0.702, 95% CI [0.591, 0.814] against 33 significant IL-2-reducing hits. Put it in the Limitations opener. Then concede, honestly, that the window score itself has no functional readout.

**Discussion vision, pp.28-29, which is NOT a limitation:**
- p.28: *"While we employ a relatively simple model here, we view this as a proof-of-concept and a baseline for more sophisticated approaches, for example, models that directly incorporate context-dependent perturbation effects."*
- p.29: scaling to hundreds of millions of cells; and the closing paragraph, *"Now, perturb-seq in human primary cells has potential to map how genetic variation controls cell states, offering new hope that we can systematically link genome sequences to cell programs to human health outcomes."*
- p.27, against a program-first stance: *"While gene program signatures can provide a useful abstraction [3], here, quantification of gene-level effects enabled us to discover novel cytokine regulators (Figure 2) and query regulators of specific disease-associated genes (Figure 7)."*

**Fix `docs/strategy_vs_authors_future_work.md`:** its "five items" list flattens Limitations (p.30) and Discussion vision (pp.28-29) into one, files "tighter pop-gen" as future work when the pop-gen analyses are **completed** (Methods pp.50-51), and **omits the fifth limitation entirely** — the one sentence that names immunotherapy. Add it. Split the list by source.

---

# 7. WHAT WE GOT RIGHT

Stop worrying about these.

| claim | confirmed by |
|---|---|
| FDR = 0.10 is the DE threshold | p.5 main text, p.44 Methods. `04:47` and `10:38` correct. No script uses another threshold against `adj_p_value`. |
| `n_total_de_genes == n_downstream + ontarget_significant` | 0 mismatches, **all three** conditions, 33,983 rows. Independently corroborated by three reproductions of the authors' headline numbers (§2.2). Our "both conditions" understated it. |
| The on-target gene is counted in `n_total_de_genes` and excluded from `n_downstream` | Their 7,807/67% reproduces only on `n_total_de_genes`; their 2,035,311 trans-effects only on `sum(n_downstream)`; their median 2 / p95 427 / mean 81.61 only on `n_downstream`. |
| 6,371 of 11,281 Stim48hr survive our mask | reproduced exactly |
| 51 with no Rest row; 802 Rest-QC-fail; 637 on `ontarget_significant`; 66 on `low_target_gex` | all four reproduced exactly |
| 11,526 unique perturbed genes; 96.18% in all three conditions | paper says 11,527 and *"96% of them tested in all conditions"* (p.5). Off by one. Footnote it. |
| `ontarget_significant` is a collider | The authors say so themselves. p.4-5 (expression) and p.37 (cell number, KS p = 6.57e-12). We reproduce KS D = 0.146, p = 1.4e-25 in Stim48hr among `~low_target_gex`. Strengthened, not weakened. |
| `low_target_gex` is condition-specific | 358 genes low at Rest but not Stim48hr, 96 the reverse. Our 66 Rest-only-flagged rows are expected behaviour, not an anomaly. |
| Rest and Stim were fitted **separately**, per condition, each against its own NTC pseudobulks | p.36. Our code never contrasts `log_fc` across conditions. |
| `donor_correlation_*` is the mean/min Pearson r over **3 disjoint donor-pair splits** from C(4,2)=6 pairwise DE runs | p.39, verbatim. Only the gene set was wrong: it is the union of genes DE at lenient 20% FDR, not all genes. |
| `donor_correlation_*` NaN for 86% of rows | 85.64%. Exact. |
| `guide_correlation_all` present for 77% of rows | 76.63%. Exact. |
| PPP3R1 `guide_correlation_signif` = −0.39 | −0.390491, p = 1.17e-21, `n_guides` = 2, `guide_n_signif_ontarget` = 2. Both guides knocked down the target. Not a low-power artefact. |
| CD3E and COMMD5 are `single_guide_estimate` | both `n_guides` = 1. Also: 461 of our 6,371 (7.2%) are single-guide. That is a designed-in feature of the library (104 genes carry one guide), not an anomaly. |
| PPP3CA 1 DE gene, PPP3CB 3, PPP3R1 523 | exact, on `n_total_de_genes`. On `n_downstream`: 0, 2, 522. Quote the downstream figures if quoting the paper's convention. **Hedge the mechanism.** MYD88 (p.24) is a published counterexample: a top LoF-burden hit with < 3 downstream DE genes despite significant knockdown, and redundancy was not the explanation. |
| Our QC mask drops off-target and neighbour-knockdown rows in the correct direction | p.5, p.37 |
| IL10 and TGFB1 belong on the penalty side, not the objective side | p.8, the authors' own bin: *"regulatory cytokines (IL-10, TGF-β)"* |
| The screen uses **no polarising cytokines** | p.18 and p.30, verbatim. Keep this. Only the implied "no cytokines" dies. |
| The paper contains no druggability, tractability, therapeutic-window, drug-target-recovery, or safety-triage analysis | 0 hits, 63 pages |
| Our shortlist sits inside the authors' own quality envelope, unprompted | 20/20 of top 20 and 91/100 of top 100 by `window_rank` fall inside their strong-perturbation set (`n_total_de_genes > 75 & n_cells_target > 50`, p.42). 18/20 of top 20 pass their `donor_correlation_hits_min > 0.35` (Fig 2C, p.11). Neither threshold enters our ranking. |
| The tolerance axis, not the homeostasis axis, is the reversal-specific one | our own `reversal_specificity.csv`: homeostasis p = 0.5296 (not specific), tolerance p = 0.0 (specific). The paper wounds only the axis our own null already demoted. |
| The 9-gene module is not an IL10 proxy | 7/7 positive-arm vs 1/5 negative-arm concordance against the paper's arrayed-validated IL10 regulators (Fig 2F, p.10). Fisher p = 0.0101. Commit this test. |

---

## Minimum viable commit for tonight

1. `scripts/04_window_score.py`: `n_downstream` for selectivity (`:175`, `:181`); drop `fillna` (`:182`); QC the Rest join (`:180`); add `rest_qc_pass`, `rest_de_imputed`, `tol_mean_z`, `donor_correlation_hits_min`, `below_stabilization` (`n_cells_target < 200`) columns; add `DEPTH_RATIO` sensitivity print. Update the `:20-25` and `:190-194` docstrings.
2. `scripts/18_n6_selectivity_validation.py:98-99`: update the asserts. `:87-88`: fix the false mechanical-entailment comment.
3. `src/cd4_perturbseq/programs.py:17` and `:3-6`: delete the Treg claim.
4. `src/cd4_perturbseq/de_stats.py`: `qc_mask` docstring records the two collider parents, the recovered `low_target_gex` rule, the p.37 `neighboring_gene_KD` definition, and that `distal_offtarget_flag` is undefined in the paper.
5. `src/cd4_perturbseq/magnitude.py:3-6`: scope the "not a power proxy" sentence.
6. Prose: `README.md:51-53`; `docs/ground_truth_immunomodulators.md:91-95`; `docs/preregistration_n9_2026_07_08.md:20` and a dated amendment at `:169-170`; `scripts/17_tolerance_is_special.py:26,569`; `docs/strategy_2026_07_07.md:91-92`; `docs/claude_tooling_log.md:17`; `scripts/02_risk_kill_reversal.py:343`; `docs/results/risk_kill_2026_07_08.md:30`.
7. Re-run 04, 05, 11, 12, 18. Expect `safe` 66 → 67 (TBX21 in), selectivity-gate passers 370 → ~286, LMLN/ABCA3/SMIM7 out of the top 30. Diff the N6/N9 preregistered primaries and record any change as an amendment, not a correction.

Deferred to before submission, not tonight: the z-score sensitivity run, the `--drop-gene` IL10 leave-one-out on script 17, the Suppl. Table 9 annotation join, the IL10-regulator concordance script, the Stim8hr arm, and the 31-gene and 52-row discrepancy footnotes.
