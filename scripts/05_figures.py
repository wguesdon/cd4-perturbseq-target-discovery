"""Render the four demo figures, in light and dark, as PNG and SVG.

Figure 1  The naive ranking's top 20, and what the safety gate does to it. Emphasis form:
          the rejected are gray, the survivor carries the accent. The story is "almost
          nothing at the top survives", so one accent against gray says it, and five
          categorical hues would bury it.
Figure 2  Nine knockdowns that an independent lab's protein-level screen proves reduce
          IL-2, all refused by the gate. One series, therefore one hue. The reject reason
          is a direct label, not a second colour channel.
Figure 3  The gate's geometry, with five approved-drug targets placed on it. Pass and reject
          are a genuine status, so they use the reserved status colours. Green against red
          is deutan dE 12.4, which sits just over the floor, so the encoding is never colour
          alone: marker shape and a PASS or REJECT label carry it too.
Figure 4  The magnitude-matched control, and the reason the other three can be trusted at all.
          Figure 1 says 18 of 20 are rejected. Figure 4 says a shortlist of the same effect
          magnitude is rejected 15.3 times, so figure 1 is worth 2.7 rejections. The bars are
          observed against a null, so the null recedes; the gap is the finding, so it is
          annotated with its p-value rather than left to the reader.

**Every label in figures 2 and 3 is composed from the gate's own output.** They used to be
hand-written, and on 2026-07-08 they were found asserting that CD3E and CD3G are "rejected on
the immunodeficiency axis". CD3E passes. CD3G is rejected on tolerance. The immunodeficiency axis
had been removed from the gate hours earlier. :func:`_assert_no_dead_axis` and the ``KeyError`` in
:func:`_pretty_reason` exist so that a renamed or removed axis fails the render instead of
printing a confident falsehood at 200 dpi.

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
    "tolerance": "collapses co-inhibitory module",
    "evidence": "no effect",
}
"""Two gate axes. `homeostasis` (selectivity) was demoted to a reported annotation on 2026-07-08
(N6): it lost a pre-registered dominance test to a raw resting-DE threshold, and the association it
rested on appears just as strongly among perturbations the gate never touches. It no longer enters
`reject_reason`, so it is not a label here. `immune_essential` was removed earlier for the same
class of reason (more enriched among approved drugs than among the toxic top 100)."""
"""The gate has THREE reject axes. ``immune_essential`` was a fourth and was removed on
2026-07-08, because the IUIS flag is more enriched among approved drug targets (OR 8.31) than
among the perturbations a naive ranking calls toxic (OR 4.16). :func:`_assert_no_dead_axis`
fails the build if it ever reappears, because these figures once labelled every rejected drug
target an "immunodeficiency gene" from a hardcoded string, and it was on camera."""

TIER_LABEL = {
    "non-depleting": "non-depleting",
    "depleting-at-rest": "depleted at rest",
    "unknown": "viability unknown",
}
"""Viability is a measurement. Whether depletion at rest is *antiproliferation* or *toxicity* is an
interpretation of it, true of IMPDH2 and unearned for CD3E, so the figures state the tier and leave
the interpreting to the footer."""

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
        reason: The ``reject_reason`` field, e.g. ``"homeostasis,tolerance"``.

    Returns:
        A short human-readable phrase, or an empty string when nothing failed.

    Raises:
        KeyError: If the reason names an axis this module does not know about. Silently passing
            an unknown axis through would let a renamed gate axis reach a figure unlabelled.
    """
    if reason == "-":
        return ""
    parts = reason.split(",")
    unknown = [p for p in parts if p not in REASON_LABEL]
    if unknown:
        raise KeyError(f"unknown reject axis {unknown} in {reason!r}; update REASON_LABEL")
    return " · ".join(REASON_LABEL[part] for part in parts)


def _assert_no_dead_axis(data: pd.DataFrame) -> None:
    """Fail loudly if a removed gate axis is still rejecting anything.

    Args:
        data: The window score frame.

    Raises:
        AssertionError: If ``immune_essential`` appears in any ``reject_reason``, or if the frame
            still carries a ``fail_immune_essential`` column.
    """
    if "fail_immune_essential" in data.columns:
        raise AssertionError("fail_immune_essential is back in window_score.csv; the IEI gate was removed")
    offenders = data.loc[data["reject_reason"].str.contains("immune_essential", na=False), "gene_name"]
    if len(offenders):
        raise AssertionError(f"immune_essential is still rejecting: {offenders.tolist()[:5]}")


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
        f"Rank perturbations purely by how well they suppress the activation program,\n"
        f"and the triage layer refuses {20 - n_safe} of the top 20.",
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
    # Counted, not asserted. The previous footer said "the TCR signalosome ... loss of function of
    # these genes causes human immunodeficiency". Two of the genes shown are not IUIS genes at all,
    # and three are transcriptional machinery rejected on homeostasis rather than TCR components.
    n_tol = int(rejected["reject_reason"].str.contains("tolerance").sum())
    n_hom = int(rejected["reject_reason"].str.contains("homeostasis").sum())
    n_iei = int(rejected["is_iei"].sum())
    fig.text(0.005, -0.04,
             f"{len(rejected)} knockdowns, all confirmed to reduce IL-2 protein. "
             f"{n_tol} are refused for collapsing tolerance, {n_hom} for disrupting the resting transcriptome.\n"
             f"{n_iei} of the {len(rejected)} are IUIS immunodeficiency genes. That flag is reported, never gated on: "
             "it is more enriched among approved\ndrug targets than among the perturbations a naive ranking calls toxic. "
             "The Schmidt screen was held out. It never entered\nthe score, the gate, or any threshold.",
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
        # The tier is a measurement. "Antiproliferative" is an interpretation of that measurement,
        # true of IMPDH2 and unearned for CD3E, so the label states the tier and the footer does the
        # interpreting. Likewise the verdict is read off the gate: this line was once a hardcoded
        # "REJECT · immunodeficiency gene", applied to every rejected drug target, naming an axis
        # the gate had already stopped using.
        tier = TIER_LABEL[row["viability_tier"]]
        verdict = f"PASS · {tier}" if passes else f"REJECT · {_pretty_reason(row['reject_reason'])}"
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

    # Composed from the gate's own output. The previous footer was hand-written and asserted that
    # CD3E and CD3G are "rejected on the immunodeficiency axis". CD3E passes, CD3G is rejected on
    # tolerance, and the immunodeficiency axis no longer exists.
    passing = drugs.loc[drugs["safe"], "gene_name"].tolist()
    rejected_here = drugs.loc[~drugs["safe"]]
    reject_phrases = [f"{r['gene_name']} ({_pretty_reason(r['reject_reason'])})" for _, r in rejected_here.iterrows()]
    fig.text(0.005, -0.05,
             f"Passes: {', '.join(passing) if passing else 'none'}.   "
             f"Rejected: {'; '.join(reject_phrases) if reject_phrases else 'none'}.\n"
             "PPP3R1 (ciclosporin) and IL4R (dupilumab) are non-depleting signalling blockers. IMPDH2 (mycophenolate) is\n"
             "depleted at rest, and for that drug the depletion IS the mechanism, not a disqualification. Depletion at rest is\n"
             "reported as a tier and never rejects, because the screen cannot tell antiproliferation from toxicity.\n"
             "Tolerance is annotated rather than plotted: this figure has two spatial dimensions and the gate has three axes.\n"
             "The screen-native axes cannot see cytokine release syndrome, which is why muromonab (anti-CD3) was withdrawn.",
             fontsize=9, color=theme["muted"])
    _save(fig, "fig3_gate_geometry", mode)


def figure_4(mode: str) -> None:
    """The magnitude-matched control. The one panel that carries the whole argument.

    "18 of the naive top 20 are rejected" sounds decisive until you ask what a shortlist of the
    same transcriptome-wide effect magnitude scores. It scores 15.3. The headline is worth 2.7
    rejections. Everything specific to reversal lives in the tolerance axis, and this figure exists
    to make that impossible to miss.

    Encoding: two bars per group, observed against a magnitude-matched null. The null is muted
    because it is the reference, not a series. The gap IS the finding, so it is annotated directly
    with the Monte-Carlo p-value rather than left to the reader to subtract. Colour is never the
    only channel: the p-value and the "not specific" / "specific" verdict carry it too.

    Args:
        mode: ``"light"`` or ``"dark"``.
    """
    theme = THEMES[mode]
    spec = pd.read_csv(paths.TABLES / "reversal_specificity.csv")

    labels, observed, matched, pvals = [], [], [], []
    for pretty, key in (
        ("SELECTIVITY\n(we demoted it)", "fails selectivity annotation"),
        ("CO-INHIBITORY\n(we kept it)", "rejected on tolerance"),
    ):
        row = spec[spec["quantity"].str.startswith(key)].iloc[0]
        labels.append(pretty)
        observed.append(float(row["observed"]))
        matched.append(float(row["matched_mean"]))
        pvals.append(float(row["p_matched_ge_observed"]))

    n_draws = int(spec["n_draws"].iloc[0])
    pool = int(spec["pool"].iloc[0])

    fig, ax = plt.subplots(figsize=(11, 6.6), facecolor=theme["surface"])
    _style(ax, theme, xgrid=False)

    x = np.arange(len(labels))
    width = 0.34
    ax.bar(x - width / 2, observed, width, color=theme["series"],
           edgecolor=theme["surface"], linewidth=1.5, zorder=3, label="naive top 20, ranked by suppression")
    ax.bar(x + width / 2, matched, width, color=theme["recede"],
           edgecolor=theme["surface"], linewidth=1.5, zorder=3,
           label=f"effect-magnitude-matched 20, mean of {n_draws:,} draws")

    # The verdicts sit on a constant baseline rather than floating above each bar. Aligned, they
    # read as one sentence across the figure; ragged, they read as three unrelated annotations.
    verdict_y = 19.2
    for i, (obs, mat, p) in enumerate(zip(observed, matched, pvals, strict=True)):
        ax.text(i - width / 2, obs + 0.35, f"{obs:.0f}", ha="center", va="bottom",
                fontsize=12, fontweight="bold", color=theme["ink"], zorder=4)
        ax.text(i + width / 2, mat + 0.35, f"{mat:.1f}", ha="center", va="bottom",
                fontsize=12, color=theme["muted"], zorder=4)

        specific = p < 0.01
        verdict = "SPECIFIC TO REVERSAL" if specific else "not specific"
        p_text = "P < 0.0001" if p == 0 else f"P = {p:.2f}"
        ax.text(i, verdict_y, f"{p_text}\n{verdict}", ha="center", va="bottom",
                fontsize=10, linespacing=1.4, zorder=4, fontweight="bold" if specific else "normal",
                color=theme["critical"] if specific else theme["muted"])

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11, color=theme["ink2"], linespacing=1.4)
    ax.set_ylabel("perturbations rejected, of 20", fontsize=11, color=theme["ink2"])
    ax.set_ylim(0, 26)
    ax.set_yticks(np.arange(0, 21, 5))

    ax.set_title(
        "We tested our own headline. Most of it was bought by effect size.\n"
        "One axis was not.",
        fontsize=15, color=theme["ink"], loc="left", pad=16, fontweight="bold",
    )

    # Above the verdict band, clear of every bar. ylim 26 against a tallest bar of 18 buys the room.
    leg = ax.legend(loc="upper left", frameon=False, fontsize=10, bbox_to_anchor=(0.0, 1.0),
                    handlelength=1.4, handleheight=0.9)
    for text in leg.get_texts():
        text.set_color(theme["ink2"])

    fig.text(0.005, -0.06,
             f"Both groups are drawn from the same {pool} evidence-passing perturbations, matched on decile of z_L2, the\n"
             "transcriptome-wide effect magnitude, computed off-module and off-target. Homeostasis does most of the rejecting\n"
             "and explains none of it: a shortlist of equal effect size is rejected just as often. Tolerance collapse is what\n"
             "ranking by reversal specifically brings in. It is also the axis a 127-agent adversarial audit attacked hardest.",
             fontsize=9, color=theme["muted"])
    _save(fig, "fig4_magnitude_matched", mode)


def main() -> None:
    """Render every figure in both modes."""
    argparse.ArgumentParser(description=__doc__).parse_args()
    data = load()
    _assert_no_dead_axis(data)
    for mode in ("light", "dark"):
        print(f"{mode}:")
        figure_1(data, mode)
        figure_2(data, mode)
        figure_3(data, mode)
        figure_4(mode)


if __name__ == "__main__":
    main()
