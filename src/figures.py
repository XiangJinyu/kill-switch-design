"""
Publication-quality figures for the kill switch design paper.
Generates 6 main figures + supplementary.
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec
import json
from pathlib import Path

RESULTS_DIR = Path(__file__).parent.parent / "results"
FIGURES_DIR = Path(__file__).parent.parent / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

# Publication style
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "lines.linewidth": 1.2,
    }
)

COLORS = {
    "toxin_antitoxin": "#E74C3C",
    "crispr_single": "#E67E22",
    "crispr_multi": "#3498DB",
    "overlapping_gene": "#2ECC71",
    "auxotrophy": "#9B59B6",
    "integrase_differentiation": "#1ABC9C",
    "combination": "#34495E",
    "nih_threshold": "#C0392B",
}

LABELS = {
    "toxin_antitoxin": "Toxin-Antitoxin",
    "crispr_single": "CRISPR (1 gRNA)",
    "crispr_multi": "CRISPR (multi-gRNA)",
    "overlapping_gene": "Overlapping Gene",
    "auxotrophy": "Synthetic Auxotrophy",
    "integrase_differentiation": "Integrase Diff.",
}


def load(name):
    with open(RESULTS_DIR / f"{name}.json") as f:
        return json.load(f)


def fig1_architecture_comparison():
    """Figure 1: Single-layer kill switch comparison."""
    data = load("exp1_single_layer_baselines")

    fig, axes = plt.subplots(1, 2, figsize=(7, 3.2))

    # Panel A: Escape rates (bar chart)
    ax = axes[0]
    names = list(data.keys())
    rates = [data[n]["escape_rate"] for n in names]
    colors = [COLORS.get(n, "#666") for n in names]
    labels = [LABELS.get(n, n) for n in names]

    y_pos = np.arange(len(names))
    ax.barh(
        y_pos,
        [np.log10(r) for r in rates],
        color=colors,
        height=0.6,
        edgecolor="white",
        linewidth=0.5,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("log$_{10}$(Escape Rate per Cell per Generation)")
    ax.axvline(
        x=-8, color=COLORS["nih_threshold"], linestyle="--", linewidth=1, alpha=0.8
    )
    ax.text(
        -7.8,
        len(names) - 0.5,
        "NIH\nThreshold\n(10$^{-8}$)",
        fontsize=7,
        color=COLORS["nih_threshold"],
        va="top",
    )
    ax.set_xlim(-10.5, 0)
    ax.set_title("A. Single-Layer Escape Rates", fontweight="bold", loc="left")
    ax.invert_yaxis()

    # Panel B: Escape dynamics over time
    ax = axes[1]
    for name in [
        "crispr_single",
        "toxin_antitoxin",
        "overlapping_gene",
        "crispr_multi",
        "auxotrophy",
    ]:
        ts = data[name]["deterministic_timeseries"]
        gen = np.array(ts["generations"])
        frac = np.array(ts["escape_fraction"])
        ax.semilogy(
            gen,
            np.clip(frac, 1e-15, 1),
            color=COLORS[name],
            label=LABELS[name],
            linewidth=1.2,
        )
    ax.axhline(y=0.5, color="gray", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.text(5, 0.55, "50% escape", fontsize=7, color="gray")
    ax.set_xlabel("Generations")
    ax.set_ylabel("Escape Fraction")
    ax.set_ylim(1e-12, 2)
    ax.set_xlim(0, 500)
    ax.legend(loc="lower right", framealpha=0.9, edgecolor="none")
    ax.set_title("B. Escape Dynamics Over Time", fontweight="bold", loc="left")

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig1_architecture_comparison.pdf")
    fig.savefig(FIGURES_DIR / "fig1_architecture_comparison.png")
    plt.close()
    print("  Saved fig1_architecture_comparison")


def fig2_combination_heatmap():
    """Figure 2: Pairwise combination heatmap."""
    data = load("exp2_multilayer_combinations")

    architectures = [
        "toxin_antitoxin",
        "crispr_multi",
        "overlapping_gene",
        "auxotrophy",
    ]
    n = len(architectures)

    # Build matrix
    matrix = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            if i == j:
                # Single layer
                bl = load("exp1_single_layer_baselines")
                matrix[i, j] = bl[architectures[i]]["log10_escape_rate"]
            elif i < j:
                key = f"{architectures[i]}+{architectures[j]}_rho=low"
                if key in data["pairwise"]:
                    matrix[i, j] = data["pairwise"][key]["log10_combined"]
                    matrix[j, i] = matrix[i, j]

    fig, ax = plt.subplots(figsize=(5, 4.5))
    labels = [LABELS.get(a, a) for a in architectures]

    im = ax.imshow(matrix, cmap="RdYlGn", vmin=-20, vmax=0, aspect="auto")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)

    # Add text annotations
    for i in range(n):
        for j in range(n):
            if not np.isnan(matrix[i, j]):
                val = matrix[i, j]
                color = "white" if val < -12 else "black"
                ax.text(
                    j,
                    i,
                    f"{val:.1f}",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color=color,
                    fontweight="bold",
                )

    # NIH threshold contour
    for i in range(n):
        for j in range(n):
            if not np.isnan(matrix[i, j]) and matrix[i, j] <= -8:
                rect = Rectangle(
                    (j - 0.5, i - 0.5),
                    1,
                    1,
                    fill=False,
                    edgecolor=COLORS["nih_threshold"],
                    linewidth=2,
                    linestyle="--",
                )
                ax.add_patch(rect)

    cb = plt.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("log$_{10}$(Combined Escape Rate)")

    ax.set_title(
        "Pairwise Combination Escape Rates ($\\rho$ = 0.01)", fontweight="bold"
    )
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig2_combination_heatmap.pdf")
    fig.savefig(FIGURES_DIR / "fig2_combination_heatmap.png")
    plt.close()
    print("  Saved fig2_combination_heatmap")


def fig3_correlation_sensitivity():
    """Figure 3: How correlation between layers affects escape rate."""
    data = load("exp3_correlation_sensitivity")

    fig, axes = plt.subplots(1, 2, figsize=(7, 3))

    # Panel A: Escape rate vs correlation
    ax = axes[0]
    rho = np.array(data["rho_values"])
    rates = np.array(data["log10_rates"])

    # Separate zero from positive rho for log scale
    mask = rho > 0
    ax.plot(rho[mask], rates[mask], color=COLORS["crispr_multi"], linewidth=1.5)
    ax.axhline(
        y=-8,
        color=COLORS["nih_threshold"],
        linestyle="--",
        linewidth=1,
        label="NIH threshold (10$^{-8}$)",
    )

    if data["nih_threshold_rho"]:
        ax.axvline(
            x=data["nih_threshold_rho"], color="gray", linestyle=":", linewidth=0.8
        )
        ax.text(
            data["nih_threshold_rho"] * 1.2,
            -6,
            f"$\\rho_{{crit}}$ = {data['nih_threshold_rho']:.3f}",
            fontsize=7,
            color="gray",
        )

    ax.set_xscale("log")
    ax.set_xlabel("Correlation ($\\rho$) Between Kill Switch Layers")
    ax.set_ylabel("log$_{10}$(Combined Escape Rate)")
    ax.legend(fontsize=7, loc="lower right")
    ax.set_title("A. Escape Rate vs. Layer Correlation", fontweight="bold", loc="left")
    ax.set_ylim(-18, -4)

    # Panel B: Comparison of independent vs correlated for different combos
    ax = axes[1]
    combo_data = load("exp2_multilayer_combinations")
    combos = [
        "toxin_antitoxin+crispr_multi",
        "crispr_multi+overlapping_gene",
        "crispr_multi+auxotrophy",
        "overlapping_gene+auxotrophy",
    ]
    x = np.arange(len(combos))
    width = 0.35

    independent_rates = []
    correlated_rates = []
    for c in combos:
        ind_key = f"{c}_rho=independent"
        cor_key = f"{c}_rho=medium"
        if ind_key in combo_data["pairwise"]:
            independent_rates.append(combo_data["pairwise"][ind_key]["log10_combined"])
        else:
            independent_rates.append(-20)
        if cor_key in combo_data["pairwise"]:
            correlated_rates.append(combo_data["pairwise"][cor_key]["log10_combined"])
        else:
            correlated_rates.append(-20)

    bars1 = ax.bar(
        x - width / 2,
        independent_rates,
        width,
        label="Independent ($\\rho$=0)",
        color=COLORS["crispr_multi"],
        alpha=0.8,
    )
    bars2 = ax.bar(
        x + width / 2,
        correlated_rates,
        width,
        label="Correlated ($\\rho$=0.1)",
        color=COLORS["toxin_antitoxin"],
        alpha=0.8,
    )
    ax.axhline(y=-8, color=COLORS["nih_threshold"], linestyle="--", linewidth=1)

    combo_labels = ["TA + CRISPR", "CRISPR + OG", "CRISPR + Aux", "OG + Aux"]
    ax.set_xticks(x)
    ax.set_xticklabels(combo_labels, rotation=30, ha="right", fontsize=7)
    ax.set_ylabel("log$_{10}$(Escape Rate)")
    ax.legend(fontsize=7, loc="upper right")
    ax.set_title("B. Effect of Correlation", fontweight="bold", loc="left")

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig3_correlation_sensitivity.pdf")
    fig.savefig(FIGURES_DIR / "fig3_correlation_sensitivity.png")
    plt.close()
    print("  Saved fig3_correlation_sensitivity")


def fig4_optimal_design():
    """Figure 4: Optimal design - Pareto front of escape rate vs fitness cost."""
    data = load("exp4_optimal_design")
    designs = data["designs"]

    fig, ax = plt.subplots(figsize=(5.5, 4))

    for d in designs:
        n = d["n_layers"]
        marker = ["o", "s", "^", "D"][n - 1]
        color = ["#3498DB", "#2ECC71", "#E74C3C", "#9B59B6"][n - 1]
        alpha = 0.9 if d["meets_NIH"] else 0.3

        ax.scatter(
            d["fitness_cost"] * 100,
            d["log10_rate"],
            s=40,
            marker=marker,
            color=color,
            alpha=alpha,
            edgecolors="white",
            linewidth=0.3,
        )

    # Add NIH threshold line
    ax.axhline(y=-8, color=COLORS["nih_threshold"], linestyle="--", linewidth=1.2)
    ax.text(
        1, -7.5, "NIH Threshold (10$^{-8}$)", fontsize=7, color=COLORS["nih_threshold"]
    )

    # Highlight top design
    best = designs[0]
    ax.scatter(
        best["fitness_cost"] * 100,
        best["log10_rate"],
        s=120,
        marker="*",
        color="gold",
        edgecolors="black",
        linewidth=0.8,
        zorder=10,
    )
    ax.annotate(
        f"Best: {'+'.join(best['layers'][:2])}+...\n"
        f"Rate: 10$^{{{best['log10_rate']:.0f}}}$",
        xy=(best["fitness_cost"] * 100, best["log10_rate"]),
        xytext=(best["fitness_cost"] * 100 + 2, best["log10_rate"] + 3),
        fontsize=7,
        arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
    )

    # Legend for n_layers
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="#3498DB",
            markersize=6,
            label="1 layer",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor="#2ECC71",
            markersize=6,
            label="2 layers",
        ),
        Line2D(
            [0],
            [0],
            marker="^",
            color="w",
            markerfacecolor="#E74C3C",
            markersize=6,
            label="3 layers",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            color="w",
            markerfacecolor="#9B59B6",
            markersize=6,
            label="4 layers",
        ),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=7, framealpha=0.9)

    ax.set_xlabel("Fitness Cost (%)")
    ax.set_ylabel("log$_{10}$(Combined Escape Rate)")
    ax.set_title("Pareto Front: Escape Rate vs. Fitness Cost", fontweight="bold")
    ax.set_ylim(-28, 0)

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig4_optimal_design.pdf")
    fig.savefig(FIGURES_DIR / "fig4_optimal_design.png")
    plt.close()
    print("  Saved fig4_optimal_design")


def fig5_validation():
    """Figure 5: Validation against experimental data."""
    data = load("exp5_validation")

    fig, axes = plt.subplots(1, 2, figsize=(7, 3.2))

    # Panel A: Predicted vs Observed escape rates
    ax = axes[0]
    rate_results = [r for r in data["results"] if "log10_observed" in r]

    obs = [r["log10_observed"] for r in rate_results]
    pred = [r["log10_predicted"] for r in rate_results]
    names = [r["source"].split(" - ")[1] for r in rate_results]

    ax.scatter(
        obs,
        pred,
        s=60,
        c=COLORS["crispr_multi"],
        edgecolors="white",
        linewidth=0.5,
        zorder=5,
    )

    # Perfect prediction line
    lims = [-10, -3]
    ax.plot(lims, lims, "k--", linewidth=0.8, alpha=0.5, label="Perfect prediction")
    # 1 OOM bounds
    ax.fill_between(
        lims,
        [l - 1 for l in lims],
        [l + 1 for l in lims],
        alpha=0.1,
        color="gray",
        label="Within 1 OOM",
    )

    for i, name in enumerate(names):
        ax.annotate(
            name,
            (obs[i], pred[i]),
            fontsize=6,
            textcoords="offset points",
            xytext=(5, 5),
            ha="left",
        )

    ax.set_xlabel("log$_{10}$(Observed Escape Rate)")
    ax.set_ylabel("log$_{10}$(Predicted Escape Rate)")
    ax.set_xlim(-10, -3)
    ax.set_ylim(-10, -3)
    ax.set_aspect("equal")
    ax.legend(fontsize=7, loc="upper left")
    ax.set_title("A. Escape Rate Calibration", fontweight="bold", loc="left")

    # Panel B: Time to escape comparison
    ax = axes[1]
    time_results = [r for r in data["results"] if "observed_generations" in r]
    if time_results:
        obs_t = [r["observed_generations"] for r in time_results]
        pred_t = [r["predicted_generations"] for r in time_results]
        labels_t = ["Unprotected\n(+ile)", "Entangled\n(-ile)"]

        x = np.arange(len(time_results))
        width = 0.35
        ax.bar(
            x - width / 2,
            obs_t,
            width,
            label="Observed",
            color=COLORS["overlapping_gene"],
            alpha=0.8,
        )
        ax.bar(
            x + width / 2,
            pred_t,
            width,
            label="Predicted",
            color=COLORS["crispr_multi"],
            alpha=0.8,
        )
        ax.set_xticks(x)
        ax.set_xticklabels(labels_t, fontsize=8)
        ax.set_ylabel("Generations to 10% Escape")
        ax.legend(fontsize=7)
        ax.set_title("B. Time to Escape Validation", fontweight="bold", loc="left")

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig5_validation.pdf")
    fig.savefig(FIGURES_DIR / "fig5_validation.png")
    plt.close()
    print("  Saved fig5_validation")


def fig6_ablation_sensitivity():
    """Figure 6: Ablation and sensitivity analysis."""
    ablation = load("exp6_ablation")
    sensitivity = load("exp7_sensitivity")

    fig, axes = plt.subplots(1, 2, figsize=(7, 3.2))

    # Panel A: Ablation
    ax = axes[0]
    full_rate = ablation["full_system"]["log10_rate"]
    names = [a["removed"] for a in ablation["ablations"]]
    rates = [a["log10_rate"] for a in ablation["ablations"]]
    labels = [LABELS.get(n, n) for n in names]
    colors = [COLORS.get(n, "#666") for n in names]

    y_pos = np.arange(len(names) + 1)
    all_labels = ["Full System"] + [f"Remove\n{l}" for l in labels]
    all_rates = [full_rate] + rates
    all_colors = [COLORS["combination"]] + colors

    ax.barh(
        y_pos, all_rates, color=all_colors, height=0.6, edgecolor="white", linewidth=0.5
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(all_labels, fontsize=7)
    ax.axvline(x=-8, color=COLORS["nih_threshold"], linestyle="--", linewidth=1)
    ax.set_xlabel("log$_{10}$(Escape Rate)")
    ax.set_title("A. Ablation Study", fontweight="bold", loc="left")
    ax.set_xlim(-20, 0)
    ax.invert_yaxis()

    # Panel B: Sensitivity
    ax = axes[1]
    for name in ["crispr_multi", "overlapping_gene", "auxotrophy"]:
        if name in sensitivity:
            factors = np.array(sensitivity[name]["factors"])
            rates_s = np.array(sensitivity[name]["log10_rates"])
            ax.plot(
                np.log10(factors),
                rates_s,
                color=COLORS[name],
                label=f"{LABELS[name]} (S={sensitivity[name]['local_sensitivity']:.2f})",
                linewidth=1.2,
            )

    ax.axhline(
        y=-8, color=COLORS["nih_threshold"], linestyle="--", linewidth=1, alpha=0.7
    )
    ax.axvline(x=0, color="gray", linestyle=":", linewidth=0.5)
    ax.set_xlabel("log$_{10}$(Factor Change in Component Rate)")
    ax.set_ylabel("log$_{10}$(Combined Escape Rate)")
    ax.legend(fontsize=6, loc="lower right")
    ax.set_title("B. Sensitivity Analysis", fontweight="bold", loc="left")
    ax.set_ylim(-25, -5)

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig6_ablation_sensitivity.pdf")
    fig.savefig(FIGURES_DIR / "fig6_ablation_sensitivity.png")
    plt.close()
    print("  Saved fig6_ablation_sensitivity")


def generate_all():
    """Generate all figures."""
    print("Generating publication-quality figures...")
    fig1_architecture_comparison()
    fig2_combination_heatmap()
    fig3_correlation_sensitivity()
    fig4_optimal_design()
    fig5_validation()
    fig6_ablation_sensitivity()
    print(f"All figures saved to {FIGURES_DIR}/")


if __name__ == "__main__":
    generate_all()
