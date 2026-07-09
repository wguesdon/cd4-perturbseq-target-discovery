"""Why does a perturbation with 1 significant DE gene score efficacy_t = 3.5?

Hypotheses to separate:
  H1  VIF is near 1, so CAMERA applied no correction. Would mean the co-regulation premise
      is wrong, because rho_bar is diluted by the mass of null perturbations.
  H2  The competitive statistic divides by the within-perturbation SD of z. Quiet perturbations
      have small SD, so a tiny module shift becomes a large t. CAMERA REWARDS quietness.
  H3  There is no effect-size floor. mean module z of -0.6 is statistically "significant" against
      10,250 background genes but is biologically nothing.

Also: is the safety gate rejecting real drug targets? Locate PPP3CA, MTOR, PIK3CD, IL2RA,
CD3E, CD3G on every axis.
"""

from __future__ import annotations

import pathlib

import sys

import numpy as np
import pandas as pd

# Resolve the checkout from this file, so the script runs from any clone.
REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, f"{REPO}/src")

from cd4_perturbseq import de_stats, priors, programs, score  # noqa: E402

pd.set_option("display.width", 220)
pd.set_option("display.max_columns", 40)

obs = de_stats.read_obs().reset_index(drop=True)
var = de_stats.read_var()
names = var["gene_name"].astype(str).to_numpy()
idx_of = {g: i for i, g in enumerate(names)}

program = programs.load_activation_program()
effector = [g for g in programs.effector_core(program) if g in idx_of]
eff_idx = np.array([idx_of[g] for g in effector])

obs["target_contrast_gene_name"] = obs["target_contrast_gene_name"].astype(str)
keep = (
    (obs["culture_condition"] == "Stim48hr")
    & ~obs["distal_offtarget_flag"].astype(bool)
    & ~obs["neighboring_gene_KD"].astype(bool)
    & obs["ontarget_significant"].astype(bool)
    & ~obs["low_target_gex"].astype(bool)
)
rows = obs.index[keep].to_numpy()
sub = obs.loc[rows].reset_index(drop=True)
print(f"QC-passing Stim48hr: {rows.size}")

# ---- module-only layers: cheap, 32 columns each
z_mod = de_stats.read_layer_columns(eff_idx, layer="zscore")[rows]
lfc_mod = de_stats.read_layer_columns(eff_idx, layer="log_fc")[rows]
padj_mod = de_stats.read_layer_columns(eff_idx, layer="adj_p_value")[rows]
print(f"module blocks: {z_mod.shape}")

# ---- H1: what is rho_bar / VIF actually?
vif, rho = score.inter_gene_vif(z_mod)
print(f"\n[H1] effector module rho_bar = {rho:.4f}   VIF = {vif:.3f}")
print(f"     effective module size = {len(effector)/vif:.1f} of {len(effector)}")
print("     -> VIF near 1 means CAMERA applied essentially NO correction.")

# rho_bar restricted to perturbations that actually did something
active = sub["n_total_de_genes"].to_numpy() >= 100
vif_a, rho_a = score.inter_gene_vif(z_mod[active])
print(f"     among the {active.sum()} perturbations with >=100 DE genes: rho_bar={rho_a:.3f} VIF={vif_a:.2f}")
print("     the co-regulation is real, but diluted to ~0 by the mass of null perturbations.")

# ---- H2: does the within-perturbation SD reward quiet perturbations?
z_all_sd = []
# approximate global SD per perturbation using the module block is wrong; use DE counts as proxy
sub["mean_mod_z"] = np.nanmean(z_mod, axis=1)
sub["mean_mod_lfc"] = np.nanmean(lfc_mod, axis=1)
sub["n_mod_down"] = np.nansum((padj_mod < 0.10) & (lfc_mod < 0), axis=1)
sub["n_mod_sig"] = np.nansum(padj_mod < 0.10, axis=1)

print("\n[H3] the evidence floor problem")
print("     among perturbations with <=5 significant DE genes genome-wide:")
quiet = sub[sub["n_total_de_genes"] <= 5]
print(f"       n = {len(quiet)}")
print(f"       mean module z: min {quiet['mean_mod_z'].min():.2f}  p1 {quiet['mean_mod_z'].quantile(0.01):.2f}")
print(f"       how many have mean module z < -0.5? {(quiet['mean_mod_z'] < -0.5).sum()}")
print(f"       how many module genes are significantly DOWN in these? max = {quiet['n_mod_down'].max()}")

print("\n     module DE-gene count vs total DE-gene count:")
for lo, hi in [(0, 5), (6, 50), (51, 200), (201, 1000), (1001, 10**9)]:
    m = sub[(sub["n_total_de_genes"] >= lo) & (sub["n_total_de_genes"] <= hi)]
    if len(m) == 0:
        continue
    print(f"       total_de {lo:>5}-{hi:<6} n={len(m):5d}  median n_mod_down={m['n_mod_down'].median():.0f}  "
          f"median mean_mod_lfc={m['mean_mod_lfc'].median():+.3f}")

# ---- where do real drug targets sit?
rest = obs[obs["culture_condition"] == "Rest"].drop_duplicates("target_contrast_gene_name").set_index("target_contrast_gene_name")
sub["rest_de"] = sub["target_contrast_gene_name"].map(rest["n_total_de_genes"])

self_idx = np.array([idx_of.get(g, -1) for g in sub["target_contrast_gene_name"]], dtype=np.int64)

print("\n[gate audit] where do approved-drug targets sit on each axis?")
probe = ["PPP3CA", "PPP3CB", "PPP3R1", "MTOR", "PIK3CD", "IL2RA", "CD3E", "CD3G", "CD247",
         "TRAF6", "SCRIB", "KIDINS220", "VAV1", "STAT5B", "NSD1", "CAST", "UPF3A", "CXCR5"]
cols = ["target_contrast_gene_name", "n_total_de_genes", "rest_de", "n_mod_down",
        "mean_mod_z", "mean_mod_lfc", "n_cells_target"]
view = sub[sub["target_contrast_gene_name"].isin(probe)][cols].copy()
view = view.sort_values("n_mod_down", ascending=False)
print(view.to_string(index=False))

p90_stim = sub["n_total_de_genes"].quantile(0.90)
p75_rest = sub["rest_de"].quantile(0.75)
print(f"\n  current gate: stim_de <= {p90_stim:.0f} (p90), rest_de <= {p75_rest:.0f} (p75)")
print("  NOTE: a collateral cap on STIM DE genes penalises efficacy. A real target changes many")
print("  activation genes. The discriminator should be REST disruption, not Stim breadth.")

# ---- how many Schmidt IL-2 hits are quiet?
schmidt = priors.schmidt_cd4_il2_screen()
m = sub.merge(schmidt, left_on="target_contrast_gene_name", right_on="gene_name", how="left")
m["il2_hit"] = (m["il2_neg_fdr"] < 0.05) & (m["il2_lfc"] < 0)
hits = m[m["il2_hit"]]
print(f"\n[held-out] {len(hits)} Schmidt IL-2 hits in the QC set")
print(f"  their median total DE genes: {hits['n_total_de_genes'].median():.0f}  (all perturbations: {sub['n_total_de_genes'].median():.0f})")
print(f"  their median n_mod_down:     {hits['n_mod_down'].median():.0f}")
print(f"  their median rest_de:        {hits['rest_de'].median():.0f}")
print(f"  how many have n_mod_down >= 3? {(hits['n_mod_down'] >= 3).sum()} of {len(hits)}")
print(f"  how many pass rest_de <= {p75_rest:.0f}? {(hits['rest_de'] <= p75_rest).sum()} of {len(hits)}")
print(f"  how many pass BOTH? {((hits['n_mod_down'] >= 3) & (hits['rest_de'] <= p75_rest)).sum()}")
print("\n  the Schmidt hits that would survive an evidence floor + rest gate:")
surv = hits[(hits["n_mod_down"] >= 3) & (hits["rest_de"] <= p75_rest)]
if len(surv):
    print(surv[["target_contrast_gene_name", "n_total_de_genes", "rest_de", "n_mod_down", "mean_mod_lfc", "il2_lfc"]].to_string(index=False))
else:
    print("    NONE")
