"""
Iteration 2 experiments: Robustness checks and supplementary analyses.
These address anticipated reviewer concerns.
"""

import numpy as np
import json
import time
from pathlib import Path
from .models import (
    MutationSpectrum,
    MUTATION_SPECTRA,
    KillSwitchLayer,
    MultiLayerKillSwitch,
    stochastic_passage_simulation,
    run_replicate_simulations,
    time_to_escape_analytical,
)
from .experiments import make_layer, _save, bootstrap_ci, N_REPLICATES

RESULTS_DIR = Path(__file__).parent.parent / "results"


# =========================================================================
# Exp 10: Population size (N0) sensitivity
# =========================================================================
def exp10_N0_sensitivity(n_reps=2000):
    """How does effective population size affect escape dynamics?"""
    t0 = time.time()
    print("=" * 70)
    print(f"Experiment 10: Population size sensitivity ({n_reps} reps)")
    print("=" * 70)

    N0_values = [1e5, 1e6, 1e7, 1e8, 1e9]
    K_values = [1e7, 1e8, 1e9, 1e10, 1e11]  # K = N0 * 100

    # Test on TA (well-calibrated) and best triple
    configs = {
        "toxin_antitoxin": MultiLayerKillSwitch(layers=[make_layer("toxin_antitoxin")]),
        "triple_best": MultiLayerKillSwitch(
            layers=[
                make_layer("crispr_multi"),
                make_layer("overlapping_gene"),
                make_layer("auxotrophy"),
            ],
            genomic_correlation=0.01,
        ),
    }

    results = {}
    for config_name, ks in configs.items():
        data = []
        for N0, K in zip(N0_values, K_values):
            stats = run_replicate_simulations(
                n_replicates=n_reps,
                N0=int(N0),
                generations=1000,
                kill_switch=ks,
                K=K,
                dilution_factor=100,
                passage_gens=6.6,
                seed=42,
            )
            t10_lo, t10_med, t10_hi = bootstrap_ci(
                stats["escape_times_10"], stat_fn=np.median
            )
            t50_lo, t50_med, t50_hi = bootstrap_ci(
                stats["escape_times_50"], stat_fn=np.median
            )

            data.append(
                {
                    "N0": float(N0),
                    "K": float(K),
                    "t10_median": t10_med,
                    "t10_ci95": [t10_lo, t10_hi],
                    "t50_median": t50_med,
                    "t50_ci95": [t50_lo, t50_hi],
                    "final_frac_mean": float(np.mean(stats["final_fractions"])),
                }
            )
            print(f"  {config_name} N0={N0:.0e}: t10={t10_med:.0f}, t50={t50_med:.0f}")

        results[config_name] = data

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")
    _save("exp10_N0_sensitivity", results)
    return results


# =========================================================================
# Exp 11: Multi-environment sensitivity (growth condition variation)
# =========================================================================
def exp11_environment(n_reps=2000):
    """Escape dynamics under different growth conditions."""
    t0 = time.time()
    print("\n" + "=" * 70)
    print(f"Experiment 11: Multi-environment sensitivity ({n_reps} reps)")
    print("=" * 70)

    # Rottinghaus 2022 showed strong correlation between generation time and escape
    # We model this by varying the passage generations (slower growth = more gens per passage)
    # and dilution factors
    environments = {
        "rich_lab": {
            "passage_gens": 6.6,
            "dilution": 100,
            "K": 1e9,
            "description": "LB, 37C, 1:100 dilution",
        },
        "minimal_lab": {
            "passage_gens": 10.0,
            "dilution": 50,
            "K": 5e8,
            "description": "M9 minimal, 37C, 1:50 dilution",
        },
        "gut_aerobic": {
            "passage_gens": 8.0,
            "dilution": 10,
            "K": 1e8,
            "description": "Gut-like, aerobic niche, 1:10 dilution",
        },
        "gut_anaerobic": {
            "passage_gens": 12.0,
            "dilution": 5,
            "K": 5e7,
            "description": "Gut-like, anaerobic, 1:5 dilution",
        },
        "soil": {
            "passage_gens": 20.0,
            "dilution": 2,
            "K": 1e7,
            "description": "Soil, slow growth, minimal dilution",
        },
    }

    configs = {
        "toxin_antitoxin": MultiLayerKillSwitch(layers=[make_layer("toxin_antitoxin")]),
        "crispr_multi": MultiLayerKillSwitch(layers=[make_layer("crispr_multi")]),
        "triple_best": MultiLayerKillSwitch(
            layers=[
                make_layer("crispr_multi"),
                make_layer("overlapping_gene"),
                make_layer("auxotrophy"),
            ],
            genomic_correlation=0.01,
        ),
    }

    results = {}
    for config_name, ks in configs.items():
        env_data = {}
        for env_name, env in environments.items():
            N0 = int(env["K"] / env["dilution"])
            stats = run_replicate_simulations(
                n_replicates=n_reps,
                N0=N0,
                generations=1000,
                kill_switch=ks,
                K=env["K"],
                dilution_factor=env["dilution"],
                passage_gens=env["passage_gens"],
                seed=42,
            )
            t50_lo, t50_med, t50_hi = bootstrap_ci(
                stats["escape_times_50"], stat_fn=np.median
            )

            env_data[env_name] = {
                "description": env["description"],
                "N0": N0,
                "K": env["K"],
                "passage_gens": env["passage_gens"],
                "dilution": env["dilution"],
                "t50_median": t50_med,
                "t50_ci95": [t50_lo, t50_hi],
                "final_frac_mean": float(np.mean(stats["final_fractions"])),
            }
            print(f"  {config_name}/{env_name}: t50={t50_med:.0f} gen")

        results[config_name] = env_data

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")
    _save("exp11_environment", results)
    return results


# =========================================================================
# Exp 12: Long-term deployment simulation
# =========================================================================
def exp12_deployment(n_reps=5000):
    """Simulate realistic deployment scenarios over extended time."""
    t0 = time.time()
    print("\n" + "=" * 70)
    print(f"Experiment 12: Long-term deployment ({n_reps} reps)")
    print("=" * 70)

    generations_range = [100, 200, 500, 1000, 2000, 5000]

    configs = {
        "single_crispr_multi": MultiLayerKillSwitch(
            layers=[make_layer("crispr_multi")]
        ),
        "double_crispr_aux": MultiLayerKillSwitch(
            layers=[make_layer("crispr_multi"), make_layer("auxotrophy")],
            genomic_correlation=0.01,
        ),
        "triple_best": MultiLayerKillSwitch(
            layers=[
                make_layer("crispr_multi"),
                make_layer("overlapping_gene"),
                make_layer("auxotrophy"),
            ],
            genomic_correlation=0.01,
        ),
    }

    results = {}
    for config_name, ks in configs.items():
        data = []
        for gens in generations_range:
            stats = run_replicate_simulations(
                n_replicates=n_reps,
                N0=int(1e7),
                generations=gens,
                kill_switch=ks,
                K=1e9,
                dilution_factor=100,
                passage_gens=6.6,
                seed=42,
            )
            # Probability of any escape (>1% of population)
            p_escape_1pct = float(np.mean(stats["final_fractions"] > 0.01))
            p_escape_50pct = float(np.mean(stats["final_fractions"] > 0.5))

            data.append(
                {
                    "generations": gens,
                    "p_escape_1pct": p_escape_1pct,
                    "p_escape_50pct": p_escape_50pct,
                    "mean_final_frac": float(np.mean(stats["final_fractions"])),
                    "median_final_frac": float(np.median(stats["final_fractions"])),
                }
            )
            print(
                f"  {config_name}/{gens}gen: P(>1%escape)={p_escape_1pct:.3f}, "
                f"P(>50%)={p_escape_50pct:.3f}"
            )

        results[config_name] = data

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")
    _save("exp12_deployment", results)
    return results


# =========================================================================
# Exp 13: Alternative correlation models
# =========================================================================
def exp13_correlation_models():
    """Test sensitivity of results to the correlation model assumption."""
    t0 = time.time()
    print("\n" + "=" * 70)
    print("Experiment 13: Alternative correlation models")
    print("=" * 70)

    layers = [
        make_layer("crispr_multi"),
        make_layer("overlapping_gene"),
        make_layer("auxotrophy"),
    ]
    rates = [l.escape_rate for l in layers]

    rho_range = np.logspace(-4, 0, 100)
    results = {}

    for rho in rho_range:
        rho_f = float(rho)

        # Model 1: Current (pairwise additive correction)
        ks1 = MultiLayerKillSwitch(layers=layers, genomic_correlation=rho_f)
        rate1 = ks1.combined_escape_rate()

        # Model 2: Pure multiplicative with inflated individual rates
        # Each layer's effective rate is mu_i * (1 + rho * n_other_layers)
        inflated = [r * (1 + rho_f * (len(rates) - 1)) for r in rates]
        rate2 = float(np.prod(inflated))

        # Model 3: Maximum of (independent, worst-pair-correlated)
        # P = max(prod(mu_i), rho * max(mu_i))
        rate3 = max(np.prod(rates), rho_f * max(rates))

        # Model 4: Copula-inspired: P = prod(mu_i)^(1-rho) * max(mu_i)^rho
        rate4 = float(np.prod(rates) ** (1 - rho_f) * max(rates) ** rho_f)

        results[f"rho={rho_f:.6f}"] = {
            "rho": rho_f,
            "model1_pairwise": float(np.log10(rate1)) if rate1 > 0 else -30,
            "model2_inflated": float(np.log10(rate2)) if rate2 > 0 else -30,
            "model3_max": float(np.log10(rate3)) if rate3 > 0 else -30,
            "model4_copula": float(np.log10(rate4)) if rate4 > 0 else -30,
        }

    # Summarize: at which rho does each model lose NIH compliance?
    summary = {}
    for model in ["model1_pairwise", "model2_inflated", "model3_max", "model4_copula"]:
        threshold_rho = None
        for key, val in results.items():
            if val[model] > -8 and threshold_rho is None and val["rho"] > 0:
                threshold_rho = val["rho"]
        summary[model] = {"nih_threshold_rho": threshold_rho}
        print(f"  {model}: NIH lost at rho = {threshold_rho}")

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")
    _save("exp13_correlation_models", {"sweeps": results, "summary": summary})
    return results


# =========================================================================
# Run all iteration 2 experiments
# =========================================================================
def run_iter2():
    t_start = time.time()
    print("=" * 70)
    print("Iteration 2: Robustness and Supplementary Analyses")
    print("=" * 70)

    exp10_N0_sensitivity()
    exp11_environment()
    exp12_deployment()
    exp13_correlation_models()

    total = time.time() - t_start
    print(f"\nIteration 2 complete. Total: {total:.1f}s ({total / 60:.1f} min)")


if __name__ == "__main__":
    run_iter2()
