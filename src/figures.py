"""
Publication-quality figures for the kill switch design paper.
Generates 8 figures (6 main + 2 supplementary) with confidence intervals.
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import json
from pathlib import Path

RESULTS_DIR = Path(__file__).parent.parent / "results"
FIGURES_DIR = Path(__file__).parent.parent / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

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
    "nih": "#C0392B",
}
LABELS = {
    "toxin_antitoxin": "Toxin-Antitoxin",
    "crispr_single": "CRISPR (1 gRNA)",
    "crispr_multi": "CRISPR (multi-gRNA)",
    "overlapping_gene": "Overlapping Gene",
    "auxotrophy": "Synth. Auxotrophy",
    "integrase_differentiation": "Integrase Diff.",
}
MUT_TYPES = [
    "point_mutation",
    "small_indel",
    "is_element",
    "large_deletion",
    "recombination",
]
MUT_COLORS = ["#3498DB", "#E74C3C", "#F39C12", "#2ECC71", "#9B59B6"]
MUT_LABELS = [
    "Point Mutation",
    "Small Indel",
    "IS Element",
    "Large Deletion",
    "Recombination",
]


def load(name):
    with open(RESULTS_DIR / f"{name}.json") as f:
        return json.load(f)


def fig1_architecture_overview():
    """Figure 1: Architecture comparison + escape dynamics with CI."""
    data = load("exp1_single_layer")
    spec = load("exp9_mutation_spectrum")

    fig = plt.figure(figsize=(7.5, 7))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.35)

    # Panel A: Escape rates bar chart
    ax = fig.add_subplot(gs[0, 0])
    names = list(data.keys())
    rates = [data[n]["escape_rate"] for n in names]
    colors = [COLORS.get(n, "#666") for n in names]
    labels = [LABELS.get(n, n) for n in names]
    y = np.arange(len(names))
    ax.barh(
        y,
        [np.log10(r) for r in rates],
        color=colors,
        height=0.6,
        edgecolor="white",
        linewidth=0.5,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(x=-8, color=COLORS["nih"], linestyle="--", linewidth=1)
    ax.text(-7.8, -0.3, "NIH (10$^{-8}$)", fontsize=6, color=COLORS["nih"])
    ax.set_xlabel("log$_{10}$(Escape Rate)")
    ax.set_xlim(-10, -3)
    ax.set_title("A", fontweight="bold", loc="left")
    ax.invert_yaxis()

    # Panel B: Mutation spectrum stacked bar
    ax = fig.add_subplot(gs[0, 1])
    bottom = np.zeros(len(names))
    for mi, mt in enumerate(MUT_TYPES):
        fracs = [spec[n]["fractions"][mt] for n in names]
        ax.barh(
            y,
            fracs,
            left=bottom,
            height=0.6,
            color=MUT_COLORS[mi],
            label=MUT_LABELS[mi],
            edgecolor="white",
            linewidth=0.3,
        )
        bottom += fracs
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Fraction of Escape Events")
    ax.legend(fontsize=6, loc="lower right", framealpha=0.9)
    ax.set_title("B", fontweight="bold", loc="left")
    ax.invert_yaxis()

    # Panel C: Deterministic escape dynamics
    ax = fig.add_subplot(gs[1, 0])
    for name in [
        "crispr_single",
        "toxin_antitoxin",
        "overlapping_gene",
        "crispr_multi",
        "auxotrophy",
    ]:
        det = data[name]["deterministic"]
        gen = np.array(det["generations"])
        frac = np.array(det["escape_fraction"])
        ax.semilogy(
            gen, np.clip(frac, 1e-15, 1), color=COLORS[name], label=LABELS[name]
        )
    ax.axhline(y=0.5, color="gray", linestyle=":", linewidth=0.6)
    ax.set_xlabel("Generations")
    ax.set_ylabel("Escape Fraction")
    ax.set_ylim(1e-10, 2)
    ax.set_xlim(0, 1000)
    ax.legend(fontsize=6, loc="lower right")
    ax.set_title("C", fontweight="bold", loc="left")

    # Panel D: Stochastic with CI bands
    ax = fig.add_subplot(gs[1, 1])
    for name in [
        "crispr_single",
        "toxin_antitoxin",
        "overlapping_gene",
        "crispr_multi",
    ]:
        st = data[name]["stochastic"]
        gen = np.array(st["generations"])
        med = np.array(st["median_fraction"])
        p5 = np.array(st["pct5"])
        p95 = np.array(st["pct95"])
        c = COLORS[name]
        ax.plot(gen, med, color=c, label=LABELS[name])
        ax.fill_between(gen, p5, p95, alpha=0.15, color=c)
    ax.axhline(y=0.5, color="gray", linestyle=":", linewidth=0.6)
    ax.set_xlabel("Generations")
    ax.set_ylabel("Escape Fraction (median, 90% CI)")
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlim(0, 1000)
    ax.legend(fontsize=6, loc="center right")
    ax.set_title("D", fontweight="bold", loc="left")

    fig.savefig(FIGURES_DIR / "fig1_architecture_overview.pdf")
    fig.savefig(FIGURES_DIR / "fig1_architecture_overview.png")
    plt.close()
    print("  Saved fig1")


def fig2_combination_heatmap():
    """Figure 2: Pairwise combination heatmap at rho=0.01."""
    data = load("exp2_combinations")
    bl = load("exp1_single_layer")

    archs = [
        "toxin_antitoxin",
        "crispr_multi",
        "overlapping_gene",
        "auxotrophy",
        "integrase_differentiation",
    ]
    n = len(archs)
    matrix = np.full((n, n), np.nan)

    for i in range(n):
        matrix[i, i] = bl[archs[i]]["log10_escape_rate"]
    for combo_name, combo_data in data["pairwise"].items():
        parts = combo_name.split("+")
        if len(parts) == 2:
            i = archs.index(parts[0])
            j = archs.index(parts[1])
            r01 = combo_data.get("rho=0.01", {})
            if r01:
                matrix[i, j] = r01["log10_rate"]
                matrix[j, i] = r01["log10_rate"]

    fig, ax = plt.subplots(figsize=(5.5, 5))
    labels = [LABELS.get(a, a) for a in archs]
    im = ax.imshow(matrix, cmap="RdYlGn", vmin=-18, vmax=-4, aspect="auto")

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
                    fontsize=7,
                    color=color,
                    fontweight="bold",
                )
                if val <= -8:
                    ax.add_patch(
                        Rectangle(
                            (j - 0.5, i - 0.5),
                            1,
                            1,
                            fill=False,
                            edgecolor=COLORS["nih"],
                            linewidth=1.5,
                            linestyle="--",
                        )
                    )

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_yticklabels(labels, fontsize=7)
    cb = plt.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("log$_{10}$(Combined Escape Rate)")
    ax.set_title("Pairwise Combinations ($\\rho$ = 0.01)", fontweight="bold")

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig2_combination_heatmap.pdf")
    fig.savefig(FIGURES_DIR / "fig2_combination_heatmap.png")
    plt.close()
    print("  Saved fig2")


def fig3_correlation():
    """Figure 3: Correlation sensitivity with critical thresholds."""
    data = load("exp3_correlation")

    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.2))

    # Panel A: All three combos
    ax = axes[0]
    style = {
        "crispr_multi+overlapping_gene": ("-", COLORS["crispr_multi"]),
        "crispr_multi+auxotrophy": ("--", COLORS["auxotrophy"]),
        "triple_best": ("-.", COLORS["overlapping_gene"]),
    }
    for name, (ls, c) in style.items():
        d = data[name]
        rho = np.array(d["rho"])
        rates = np.array(d["log10_rate"])
        mask = rho > 0
        ax.plot(
            rho[mask],
            rates[mask],
            linestyle=ls,
            color=c,
            label=name.replace("+", " + ").replace("_", " "),
        )
        if d["nih_threshold_rho"]:
            ax.axvline(
                x=d["nih_threshold_rho"],
                color=c,
                linestyle=":",
                linewidth=0.6,
                alpha=0.5,
            )

    ax.axhline(
        y=-8, color=COLORS["nih"], linestyle="--", linewidth=1, label="NIH threshold"
    )
    ax.set_xscale("log")
    ax.set_xlabel("Inter-layer Correlation ($\\rho$)")
    ax.set_ylabel("log$_{10}$(Escape Rate)")
    ax.legend(fontsize=5.5, loc="lower right")
    ax.set_ylim(-28, -4)
    ax.set_title("A", fontweight="bold", loc="left")

    # Panel B: Phase diagram (rho vs OG factor)
    ax = axes[1]
    sens = load("exp7_sensitivity")
    phase = np.array(sens["phase_diagram"]["log10_rates"])
    rho_2d = np.array(sens["phase_diagram"]["rho_values"])
    og_2d = np.array(sens["phase_diagram"]["og_factors"])
    im = ax.contourf(
        np.log10(og_2d),
        np.log10(rho_2d),
        phase,
        levels=np.arange(-25, -4, 1),
        cmap="RdYlGn",
        extend="both",
    )
    cs = ax.contour(
        np.log10(og_2d),
        np.log10(rho_2d),
        phase,
        levels=[-8],
        colors=[COLORS["nih"]],
        linewidths=2,
    )
    ax.clabel(cs, fmt="NIH: $10^{-8}$", fontsize=6)
    ax.set_xlabel("log$_{10}$(OG Rate Factor)")
    ax.set_ylabel("log$_{10}$($\\rho$)")
    cb = plt.colorbar(im, ax=ax)
    cb.set_label("log$_{10}$(Rate)")
    ax.set_title("B", fontweight="bold", loc="left")

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig3_correlation.pdf")
    fig.savefig(FIGURES_DIR / "fig3_correlation.png")
    plt.close()
    print("  Saved fig3")


def fig4_pareto():
    """Figure 4: Pareto front of escape rate vs fitness cost."""
    data = load("exp4_pareto")
    designs = data["all_designs"]
    pareto = data["pareto_front"]

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    layer_colors = {
        1: "#3498DB",
        2: "#2ECC71",
        3: "#E74C3C",
        4: "#9B59B6",
        5: "#F39C12",
    }
    markers = {1: "o", 2: "s", 3: "^", 4: "D", 5: "p"}

    for d in designs:
        n = d["n_layers"]
        alpha = 0.9 if d["meets_NIH"] else 0.25
        ax.scatter(
            d["fitness_cost"] * 100,
            d["log10_rate"],
            s=35,
            marker=markers.get(n, "o"),
            color=layer_colors.get(n, "#666"),
            alpha=alpha,
            edgecolors="white",
            linewidth=0.3,
            zorder=3,
        )

    # Pareto front line
    pareto_sorted = sorted(pareto, key=lambda x: x["fitness_cost"])
    ax.plot(
        [p["fitness_cost"] * 100 for p in pareto_sorted],
        [p["log10_rate"] for p in pareto_sorted],
        "k-",
        linewidth=0.8,
        alpha=0.5,
        zorder=2,
    )

    ax.axhline(y=-8, color=COLORS["nih"], linestyle="--", linewidth=1.2)
    ax.text(1, -7.3, "NIH (10$^{-8}$)", fontsize=6, color=COLORS["nih"])

    # Highlight best triple (practical recommendation)
    best3 = [
        d
        for d in designs
        if d["n_layers"] == 3
        and "crispr_multi" in d["layers"]
        and "overlapping_gene" in d["layers"]
        and "auxotrophy" in d["layers"]
    ]
    if best3:
        b = best3[0]
        ax.scatter(
            b["fitness_cost"] * 100,
            b["log10_rate"],
            s=150,
            marker="*",
            color="gold",
            edgecolors="black",
            linewidth=1,
            zorder=10,
        )
        ax.annotate(
            f"Recommended\n({b['log10_rate']:.0f}, {b['fitness_cost'] * 100:.0f}%)",
            xy=(b["fitness_cost"] * 100, b["log10_rate"]),
            xytext=(b["fitness_cost"] * 100 + 3, b["log10_rate"] + 4),
            fontsize=6,
            arrowprops=dict(arrowstyle="->", lw=0.6),
        )

    handles = [
        Line2D(
            [0],
            [0],
            marker=markers[n],
            color="w",
            markerfacecolor=layer_colors[n],
            markersize=6,
            label=f"{n} layer{'s' if n > 1 else ''}",
        )
        for n in range(1, 6)
    ]
    ax.legend(handles=handles, fontsize=6, loc="upper right")
    ax.set_xlabel("Fitness Cost (%)")
    ax.set_ylabel("log$_{10}$(Escape Rate)")
    ax.set_title(
        "Pareto Front: Escape Rate vs. Fitness Cost ($\\rho$ = 0.01)",
        fontweight="bold",
        fontsize=9,
    )
    ax.set_ylim(-32, 0)

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig4_pareto.pdf")
    fig.savefig(FIGURES_DIR / "fig4_pareto.png")
    plt.close()
    print("  Saved fig4")


def fig5_validation():
    """Figure 5: Cross-validation results."""
    data = load("exp5_validation")

    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.5))

    # Panel A: Predicted vs Observed
    ax = axes[0]
    cv = data["cross_validation"]
    obs = [r["observed"] for r in cv]
    pred = [r["predicted_raw"] for r in cv]
    studies = [r["study"] for r in cv]

    for i, r in enumerate(cv):
        c = "#3498DB" if r["study"] == "rottinghaus" else "#E74C3C"
        m = "o" if r["study"] == "rottinghaus" else "s"
        ax.scatter(
            np.log10(r["observed"]),
            np.log10(r["predicted_raw"]),
            s=60,
            c=c,
            marker=m,
            edgecolors="white",
            linewidth=0.5,
            zorder=5,
        )
        ax.annotate(
            r["name"].replace("_", "\n"),
            (np.log10(r["observed"]), np.log10(r["predicted_raw"])),
            fontsize=5,
            textcoords="offset points",
            xytext=(5, 5),
        )

    lims = [-10, -3]
    ax.plot(lims, lims, "k--", linewidth=0.6, alpha=0.5)
    ax.fill_between(
        lims, [l - 1 for l in lims], [l + 1 for l in lims], alpha=0.08, color="gray"
    )
    ax.set_xlabel("log$_{10}$(Observed)")
    ax.set_ylabel("log$_{10}$(Predicted)")
    ax.set_xlim(-10, -3)
    ax.set_ylim(-10, -3)
    ax.set_aspect("equal")
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="#3498DB",
            markersize=6,
            label="Rottinghaus 2022",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor="#E74C3C",
            markersize=6,
            label="Chlebek 2023",
        ),
    ]
    ax.legend(handles=handles, fontsize=6)
    ax.set_title("A", fontweight="bold", loc="left")
    ax.text(
        -9.5, -4, f"Mean |log err| = {data['mean_abs_log_error_raw']:.3f}", fontsize=7
    )

    # Panel B: Time to escape
    ax = axes[1]
    tv = data["time_validation"]
    x = np.arange(len(tv))
    obs_t = [t["observed_gen"] for t in tv]
    pred_t = [t["stochastic_median"] for t in tv]
    pred_lo = [t["stochastic_ci95"][0] for t in tv]
    pred_hi = [t["stochastic_ci95"][1] for t in tv]
    yerr = [
        [p - lo for p, lo in zip(pred_t, pred_lo)],
        [hi - p for p, hi in zip(pred_t, pred_hi)],
    ]

    w = 0.35
    ax.bar(
        x - w / 2,
        obs_t,
        w,
        label="Observed",
        color=COLORS["overlapping_gene"],
        alpha=0.8,
    )
    ax.bar(
        x + w / 2,
        pred_t,
        w,
        label="Predicted (stochastic)",
        color=COLORS["crispr_multi"],
        alpha=0.8,
        yerr=yerr,
        capsize=3,
        error_kw={"linewidth": 0.8},
    )
    ax.set_xticks(x)
    ax.set_xticklabels(["TA (+ile)\n30 gen", "OG (-ile)\n>130 gen"], fontsize=7)
    ax.set_ylabel("Generations to 10% Escape")
    ax.legend(fontsize=6)
    ax.set_title("B", fontweight="bold", loc="left")

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig5_validation.pdf")
    fig.savefig(FIGURES_DIR / "fig5_validation.png")
    plt.close()
    print("  Saved fig5")


def fig6_ablation_sensitivity():
    """Figure 6: Ablation + sensitivity analysis."""
    ablation = load("exp6_ablation")
    sens = load("exp7_sensitivity")

    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.5))

    # Panel A: Ablation
    ax = axes[0]
    full = ablation["full_system"]
    abl = ablation["ablations"]
    all_names = ["Full System"] + [
        f"Remove\n{LABELS.get(a['removed'], a['removed'])}" for a in abl
    ]
    all_rates = [full["log10_rate"]] + [a["log10_rate"] for a in abl]
    all_colors = ["#34495E"] + [COLORS.get(a["removed"], "#666") for a in abl]

    y = np.arange(len(all_names))
    ax.barh(
        y, all_rates, color=all_colors, height=0.6, edgecolor="white", linewidth=0.5
    )
    ax.axvline(x=-8, color=COLORS["nih"], linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(all_names, fontsize=7)
    ax.set_xlabel("log$_{10}$(Escape Rate)")
    ax.set_xlim(-20, 0)
    ax.invert_yaxis()

    # Add fold-increase annotations
    for i, a in enumerate(abl):
        ax.text(
            a["log10_rate"] + 0.3,
            i + 1,
            f"{a['fold_increase']:.0e}x",
            fontsize=6,
            va="center",
            color="black",
        )

    ax.set_title("A", fontweight="bold", loc="left")

    # Panel B: Sensitivity
    ax = axes[1]
    for name in ["crispr_multi", "overlapping_gene", "auxotrophy"]:
        d = sens[name]
        ax.plot(
            np.log10(d["factors"]),
            d["log10_rates"],
            color=COLORS[name],
            label=f"{LABELS[name]} (S={d['local_sensitivity']:.2f})",
        )

    ax.axhline(y=-8, color=COLORS["nih"], linestyle="--", linewidth=1, alpha=0.7)
    ax.axvline(x=0, color="gray", linestyle=":", linewidth=0.5)
    ax.set_xlabel("log$_{10}$(Factor Change)")
    ax.set_ylabel("log$_{10}$(Combined Rate)")
    ax.legend(fontsize=6, loc="lower right")
    ax.set_ylim(-25, -5)
    ax.set_title("B", fontweight="bold", loc="left")

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig6_ablation_sensitivity.pdf")
    fig.savefig(FIGURES_DIR / "fig6_ablation_sensitivity.png")
    plt.close()
    print("  Saved fig6")


def fig_s1_convergence():
    """Supplementary: Convergence analysis."""
    data = load("exp8_convergence")

    fig, ax = plt.subplots(figsize=(4, 3))
    ns = [d["n_replicates"] for d in data]
    means = [d["mean_t50"] for d in data]
    ses = [d["se_t50"] for d in data]

    ax.errorbar(
        ns,
        means,
        yerr=[2 * s for s in ses],
        fmt="o-",
        color=COLORS["crispr_multi"],
        markersize=4,
        capsize=3,
        linewidth=0.8,
        label="Mean t$_{50}$ (95% CI)",
    )
    ax.set_xscale("log")
    ax.set_xlabel("Number of Replicates")
    ax.set_ylabel("Time to 50% Escape (gen)")
    ax.legend(fontsize=7)
    ax.set_title("Stochastic Simulation Convergence", fontsize=9, fontweight="bold")

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig_s1_convergence.pdf")
    fig.savefig(FIGURES_DIR / "fig_s1_convergence.png")
    plt.close()
    print("  Saved fig_s1_convergence")


def generate_all():
    print("Generating publication-quality figures...")
    fig1_architecture_overview()
    fig2_combination_heatmap()
    fig3_correlation()
    fig4_pareto()
    fig5_validation()
    fig6_ablation_sensitivity()
    fig_s1_convergence()
    print(f"All figures saved to {FIGURES_DIR}/")


if __name__ == "__main__":
    generate_all()
