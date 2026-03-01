"""
Supplementary figures for iteration 2 experiments.
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import json
from pathlib import Path

RESULTS_DIR = Path(__file__).parent.parent / "results"
FIGURES_DIR = Path(__file__).parent.parent / "figures"

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
    }
)

COLORS = ["#3498DB", "#E74C3C", "#2ECC71", "#9B59B6", "#F39C12"]


def load(name):
    with open(RESULTS_DIR / f"{name}.json") as f:
        return json.load(f)


def fig_s2_deployment():
    """Supplementary: Long-term deployment P(escape) vs generations."""
    data = load("exp12_deployment")

    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.2))

    # Panel A: P(>1% escape) vs generations
    ax = axes[0]
    labels = {
        "single_crispr_multi": "CRISPR multi (1 layer)",
        "double_crispr_aux": "CRISPR + Auxotrophy (2 layers)",
        "triple_best": "CRISPR + OG + Aux (3 layers)",
    }
    for i, (config, dlist) in enumerate(data.items()):
        gens = [d["generations"] for d in dlist]
        p1 = [d["p_escape_1pct"] for d in dlist]
        p50 = [d["p_escape_50pct"] for d in dlist]
        ax.plot(
            gens,
            p1,
            "o-",
            color=COLORS[i],
            label=labels.get(config, config),
            markersize=4,
            linewidth=1.2,
        )

    ax.axhline(y=0.01, color="gray", linestyle=":", linewidth=0.6)
    ax.set_xlabel("Deployment Duration (generations)")
    ax.set_ylabel("P(Escape > 1% of population)")
    ax.legend(fontsize=6, loc="upper left")
    ax.set_xscale("log")
    ax.set_ylim(-0.02, 1.05)
    ax.set_title(
        "A. Probability of Containment Breach",
        fontweight="bold",
        loc="left",
        fontsize=9,
    )

    # Panel B: Mean escape fraction vs generations
    ax = axes[1]
    for i, (config, dlist) in enumerate(data.items()):
        gens = [d["generations"] for d in dlist]
        frac = [d["mean_final_frac"] for d in dlist]
        ax.semilogy(
            gens,
            [max(f, 1e-12) for f in frac],
            "s-",
            color=COLORS[i],
            label=labels.get(config, config),
            markersize=4,
            linewidth=1.2,
        )

    ax.set_xlabel("Deployment Duration (generations)")
    ax.set_ylabel("Mean Escape Fraction")
    ax.legend(fontsize=6, loc="lower right")
    ax.set_xscale("log")
    ax.set_ylim(1e-12, 2)
    ax.set_title("B. Mean Escape Fraction", fontweight="bold", loc="left", fontsize=9)

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig_s2_deployment.pdf")
    fig.savefig(FIGURES_DIR / "fig_s2_deployment.png")
    plt.close()
    print("  Saved fig_s2_deployment")


def fig_s3_environment():
    """Supplementary: Environment comparison."""
    data = load("exp11_environment")

    fig, ax = plt.subplots(figsize=(6, 4))

    configs = list(data.keys())
    envs = list(data[configs[0]].keys())
    env_labels = [data[configs[0]][e]["description"] for e in envs]

    x = np.arange(len(envs))
    width = 0.25

    for i, config in enumerate(configs):
        t50s = [data[config][e]["t50_median"] for e in envs]
        bars = ax.bar(
            x + i * width,
            t50s,
            width,
            label=config.replace("_", " "),
            color=COLORS[i],
            alpha=0.8,
        )

    ax.set_xticks(x + width)
    ax.set_xticklabels(env_labels, rotation=30, ha="right", fontsize=6)
    ax.set_ylabel("Time to 50% Escape (generations)")
    ax.legend(fontsize=6)
    ax.set_title(
        "Kill Switch Performance Across Growth Environments",
        fontweight="bold",
        fontsize=9,
    )

    # Cap the y-axis and add ">" for never-escaped
    ax.set_ylim(0, 1100)
    for i, config in enumerate(configs):
        for j, e in enumerate(envs):
            t50 = data[config][e]["t50_median"]
            if t50 >= 1000:
                ax.text(
                    j + i * width, 1020, ">1000", ha="center", fontsize=5, rotation=90
                )

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig_s3_environment.pdf")
    fig.savefig(FIGURES_DIR / "fig_s3_environment.png")
    plt.close()
    print("  Saved fig_s3_environment")


def fig_s4_correlation_models():
    """Supplementary: Alternative correlation model comparison."""
    data = load("exp13_correlation_models")
    sweeps = data["sweeps"]

    fig, ax = plt.subplots(figsize=(5, 3.5))

    rho_vals = []
    models = {
        "model1_pairwise": [],
        "model2_inflated": [],
        "model3_max": [],
        "model4_copula": [],
    }
    for key in sorted(sweeps.keys(), key=lambda k: sweeps[k]["rho"]):
        val = sweeps[key]
        if val["rho"] > 0:
            rho_vals.append(val["rho"])
            for m in models:
                models[m].append(val[m])

    rho_arr = np.array(rho_vals)
    labels = {
        "model1_pairwise": "Pairwise additive (this work)",
        "model2_inflated": "Inflated individual rates",
        "model3_max": "Max(indep., worst-pair)",
        "model4_copula": "Copula-inspired",
    }
    styles = ["-", "--", "-.", ":"]

    for i, (m, vals) in enumerate(models.items()):
        ax.plot(
            rho_arr,
            vals,
            linestyle=styles[i],
            color=COLORS[i],
            label=labels[m],
            linewidth=1.2,
        )

    ax.axhline(
        y=-8, color="#C0392B", linestyle="--", linewidth=1, label="NIH threshold"
    )
    ax.set_xscale("log")
    ax.set_xlabel("Inter-layer Correlation ($\\rho$)")
    ax.set_ylabel("log$_{10}$(Combined Escape Rate)")
    ax.legend(fontsize=5.5, loc="lower right")
    ax.set_ylim(-28, -4)
    ax.set_title(
        "Sensitivity to Correlation Model Assumption\n(Triple: CRISPR + OG + Aux)",
        fontweight="bold",
        fontsize=9,
    )

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig_s4_correlation_models.pdf")
    fig.savefig(FIGURES_DIR / "fig_s4_correlation_models.png")
    plt.close()
    print("  Saved fig_s4_correlation_models")


def generate_all_supp():
    print("Generating supplementary figures...")
    fig_s2_deployment()
    fig_s3_environment()
    fig_s4_correlation_models()
    print("Done.")


if __name__ == "__main__":
    generate_all_supp()
