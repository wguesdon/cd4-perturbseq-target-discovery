"""Render the three demo figures, in light and dark, as PNG and SVG.

Figure 1  The naive ranking's top 20, and what the safety gate does to it. Emphasis form:
          the rejected are gray, the survivor carries the accent. The story is "almost
          nothing at the top survives", so one accent against gray says it, and five
          categorical hues would bury it.
Figure 2  Twelve knockdowns that an independent lab's protein-level screen proves reduce
          IL-2, all refused by the gate. One series, therefore one hue. The reject reason
          is a direct label, not a second colour channel.
Figure 3  The gate's geometry, with five approved-drug targets placed on it. Pass and reject
          are a genuine status, so they use the reserved status colours. Green against red
          is deutan dE 12.4, which sits just over the floor, so the encoding is never colour
          alone: marker shape and a PASS or REJECT label carry it too.

Usage:
    uv run python scripts/05_figures.py
"""

from __future__ import annotations

import argparse

import matplotlib as mpl
import numpy as np
import pandas as pd

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from cd4_perturbseq import paths, priors  # noqa: E402

THEMES = {
    "light": {
        "surface": "#fcfcfb",
        "ink": "#0b0b0b",
        "ink2": "#52514e",
        "muted": "#898781",
        "grid": "#e1e0d9",
        "axis": "#c3c2b7",
        "series": "#2a78d6",
        "recede": "#d8d7d0",
        "good": "#0ca30c",
        "critical": "#d03b3b",
    },
    "dark": {
        "surface": "#1a1a19",
        "ink": "#ffffff",
        "ink2": "#c3c2b7",
        "muted": "#898781",
        "grid": "#2c2c2a",
        "axis": "#383835",
        "series": "#3987e5",
        "recede": "#3a3a37",
        "good": "#0ca30c",
        "critical": "#d03b3b",
    },
}

REASON_LABEL = {
    "immune_essential": "immunodeficiency",
    "homeostasis": "disrupts rest",
    "tolerance": "kills tolerance",
    "evidence": "no effect",
}

# Hand-placed label offsets. CD3E and CD3G sit almost on top of one another (resting-cell
# ratio 0.396 vs 0.389), as do PPP3R1 and IL4R. Automatic placement collides.
DRUG_OFFSETS = {
    "CD3E": (14, 12),
    "CD3G": (14, -24),
    "IMPDH2": (14, -3),
    "PPP3R1": (14, 10),
    "IL4R": (14, -24),
}


def _pretty_reason(reason: str) -> str:
    """Turn a comma-joined machine reason into readable prose.

    Args:
        reason: The ``reject_reason`` field, e.g. ``"immune_essential,tolerance"``.

    Returns:
        A short human-readable phrase, or an empty string when nothing failed.
    """
    if reason == "-":
        return ""
    return " · ".join(REASON_LABEL.get(part, part) for part in reason.split(","))


def _style(ax, theme: dict[str, str], xgrid: bool = True) -> None:
    """Apply recessive chrome: hairline solid grid, no box, muted ticks."""
    ax.set_facecolor(theme["surface"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(theme["axis"])
        ax.spines[side].set_linewidth(0.8)
    if xgrid:
        ax.grid(axis="x", color=theme["grid"], linewidth=0.6, linestyle="-", zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(colors=theme["muted"], labelsize=10, length=0)


def _save(fig, stem: str, mode: str) -> None:
    """Write PNG and SVG for one figure and mode."""
    paths.ensure_dirs()
    for ext in ("png", "svg"):
        out = paths.FIGURES / f"{stem}_{mode}.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  wrote {stem}_{mode}.png / .svg")


def load() -> pd.DataFrame:
    """Load the window score and merge the held-out Schmidt screen."""
    frame = pd.read_csv(paths.TABLES / "window_score.csv")
    schmidt = priors.schmidt_cd4_il2_screen()
    merged = frame.merge(schmidt, on="gene_name", how="left")
    merged["il2_hit"] = (merged["il2_neg_fdr"] < 0.05) & (merged["il2_lfc"] < 0)
    # Rank 1 must be the STRONGEST suppressor, i.e. the largest value of -mean(module z).
    # Ranking ascending would crown the perturbation that raises the program most.
    merged["naive_rank"] = (-merged["eff_mean_z"]).rank(ascending=False, method="first")
    assert merged.nsmallest(1, "naive_rank")["eff_mean_z"].iloc[0] < 0, "naive rank 1 must suppress"
    return merged


def figure_1(data: pd.DataFrame, mode: str) -> None:
    """The naive top 20, judged by the gate. Emphasis: gray rejected, accent survivor."""
    theme = THEMES[mode]
    top = data.nsmallest(20, "naive_rank").sort_values("naive_rank", ascending=False)

    fig, ax = plt.subplots(figsize=(11, 7.5), facecolor=theme["surface"])
    _style(ax, theme)

    values = -top["eff_mean_z"].to_numpy()
    colors = [theme["good"] if s else theme["recede"] for s in top["safe"]]
    ax.barh(np.arange(len(top)), values, height=0.72, color=colors, zorder=2)

    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(top["gene_name"], fontsize=11, color=theme["ink"])
    ax.set_xlabel("suppression of the inflammatory effector program  (naive score)",
                  fontsize=11, color=theme["ink2"])
    ax.set_xlim(0, values.max() * 2.10)

    for i, value in enumerate(values):
        row = top.iloc[i]
        label = "PASSES THE GATE" if row["safe"] else _pretty_reason(row["reject_reason"])
        color = theme["good"] if row["safe"] else theme["muted"]
        weight = "bold" if row["safe"] else "normal"
        ax.text(value + values.max() * 0.025, i, label,
                va="center", fontsize=9, color=color, fontweight=weight)

    n_safe = int(top["safe"].sum())
    ax.set_title(
        f"Rank perturbations purely by how well they suppress inflammation,\n"
        f"and {20 - n_safe} of the top 20 are targets you must not drug.",
        fontsize=15, color=theme["ink"], loc="left", pad=16, fontweight="bold",
    )
    fig.text(0.005, -0.02,
             "Bars: naive effector-suppression score, 6,371 QC-passing CRISPRi perturbations, CD4+ T cells, 48 h stimulation.\n"
             "Labels: why the safety gate refuses each one. Green = survives the gate.",
             fontsize=9, color=theme["muted"])
    _save(fig, "fig1_naive_ranking_is_toxic", mode)


def figure_2(data: pd.DataFrame, mode: str) -> None:
    """The knockdowns an independent screen proves work, that the gate refuses."""
    theme = THEMES[mode]
    rejected = data[data["il2_hit"] & ~data["safe"] & ~data["fail_evidence"]]
    rejected = rejected.sort_values("il2_lfc")

    fig, ax = plt.subplots(figsize=(11, 7), facecolor=theme["surface"])
    _style(ax, theme)

    y = np.arange(len(rejected))
    ax.barh(y, rejected["il2_lfc"], height=0.62, color=theme["series"], zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(rejected["gene_name"], fontsize=11, color=theme["ink"])
    ax.invert_xaxis()
    ax.set_xlabel("IL-2 production after knockdown, log2 fold change\n"
                  "(Schmidt & Steinhart 2022, independent CRISPRi screen, protein readout)",
                  fontsize=10.5, color=theme["ink2"])
    ax.axvline(0, color=theme["axis"], linewidth=0.8)

    span = abs(rejected["il2_lfc"].min())
    for i, (_, row) in enumerate(rejected.iterrows()):
        ax.text(row["il2_lfc"] - span * 0.03, i, _pretty_reason(row["reject_reason"]),
                va="center", ha="right", fontsize=9.5, color=theme["muted"])
    ax.set_xlim(rejected["il2_lfc"].min() * 1.70, 0.05)

    ax.set_title(
        "These knockdowns really do shut inflammation down.\n"
        "A separate lab proved it. Our gate refuses every one of them.",
        fontsize=15, color=theme["ink"], loc="left", pad=16, fontweight="bold",
    )
    fig.text(0.005, -0.04,
             "The TCR signalosome. Effective, and loss of function of these genes causes human immunodeficiency.\n"
             "The Schmidt screen was held out: it never entered the score, the gate, or any threshold.",
             fontsize=9, color=theme["muted"])
    _save(fig, "fig2_effective_but_rejected", mode)


def figure_3(data: pd.DataFrame, mode: str) -> None:
    """The gate geometry, with approved-drug targets placed on it."""
    theme = THEMES[mode]
    truth = pd.read_csv(paths.GROUND_TRUTH / "immunomodulator_targets.csv")
    positives = truth[truth["include_as_positive"]].set_index("gene_symbol")

    context = data[~data["fail_evidence"] & data["rest_cells_ratio"].notna()]
    drugs = context[context["gene_name"].isin(positives.index)]

    fig, ax = plt.subplots(figsize=(11, 7.5), facecolor=theme["surface"])
    _style(ax, theme, xgrid=False)
    ax.grid(color=theme["grid"], linewidth=0.6, zorder=0)

    ax.scatter(context["rest_cells_ratio"], context["selectivity"], s=14,
               color=theme["recede"], zorder=2, linewidths=0)

    # The horizontal line is a GATE: below it, a perturbation is rejected. The vertical line is
    # a TIER boundary: it sorts survivors into signalling blockade versus antiproliferation and
    # rejects nobody. Drawing them identically would say the wrong thing, so the gate is a solid
    # rule and the tier boundary is dashed.
    ax.axhline(np.log(10), color=theme["ink2"], linewidth=1.2, zorder=1)
    ax.axvline(0.5, color=theme["axis"], linewidth=1.0, linestyle=(0, (6, 4)), zorder=1)

    ymin, ymax = ax.get_ylim()
    ax.text(0.52, ymax - 0.25, "tier boundary\ndepleted at rest  <  |  >  viable at rest",
            fontsize=8.5, color=theme["muted"], ha="left", va="top", linespacing=1.4)
    ax.text(ax.get_xlim()[1] * 0.98, np.log(10) - 0.42,
            "SAFETY GATE  ·  stimulated effect must exceed 10x the resting effect",
            fontsize=9, color=theme["ink2"], ha="right", va="top")

    for _, row in drugs.iterrows():
        passes = bool(row["safe"])
        gene = row["gene_name"]
        ax.scatter(row["rest_cells_ratio"], row["selectivity"],
                   s=190, marker="o" if passes else "X",
                   color=theme["good"] if passes else theme["critical"],
                   edgecolors=theme["surface"], linewidths=2, zorder=4)
        drug = str(positives.loc[gene, "drug_examples"]).split(";")[0]
        tier = "antiproliferative" if row["viability_tier"] == "depleting-at-rest" else "non-depleting"
        verdict = f"PASS · {tier}" if passes else "REJECT · immunodeficiency gene"
        ax.annotate(f"{gene}  ({drug})\n{verdict}",
                    (row["rest_cells_ratio"], row["selectivity"]),
                    textcoords="offset points", xytext=DRUG_OFFSETS.get(gene, (14, -4)),
                    fontsize=9.5, color=theme["ink"], fontweight="bold", zorder=5,
                    linespacing=1.5)

    ax.set_xlabel("resting cells retained, relative to the median perturbation",
                  fontsize=11, color=theme["ink2"])
    ax.set_ylabel("context selectivity\nlog(stimulated effect) - log(resting effect)",
                  fontsize=11, color=theme["ink2"])
    ax.set_xscale("log")
    ax.set_title(
        "Nobody told the score which drugs are well tolerated.\n"
        "It sorted them anyway.",
        fontsize=15, color=theme["ink"], loc="left", pad=16, fontweight="bold",
    )

    legend = [
        Line2D([], [], marker="o", linestyle="", markersize=11, markerfacecolor=theme["good"],
               markeredgecolor=theme["surface"], label="passes the safety gate"),
        Line2D([], [], marker="X", linestyle="", markersize=12, markerfacecolor=theme["critical"],
               markeredgecolor=theme["surface"], label="rejected by the safety gate"),
        Line2D([], [], marker="o", linestyle="", markersize=6, markerfacecolor=theme["recede"],
               markeredgecolor=theme["recede"], label="all other perturbations"),
    ]
    leg = ax.legend(handles=legend, loc="lower left", frameon=False, fontsize=10)
    for text in leg.get_texts():
        text.set_color(theme["ink2"])

    fig.text(0.005, -0.05,
             "PPP3R1 (ciclosporin) and IL4R (dupilumab) are non-depleting signalling blockers. IMPDH2 (mycophenolate)\n"
             "is an antiproliferative and is depleted at rest. CD3E and CD3G are immunodeficiency genes; muromonab was withdrawn.\n"
             "CD3E and CD3G clear the selectivity gate drawn here. They are rejected on the immunodeficiency axis, which is\n"
             "annotated rather than plotted, because this figure has only two spatial dimensions and the gate has four.",
             fontsize=9, color=theme["muted"])
    _save(fig, "fig3_gate_geometry", mode)


def main() -> None:
    """Render every figure in both modes."""
    argparse.ArgumentParser(description=__doc__).parse_args()
    data = load()
    for mode in ("light", "dark"):
        print(f"{mode}:")
        figure_1(data, mode)
        figure_2(data, mode)
        figure_3(data, mode)


if __name__ == "__main__":
    main()
