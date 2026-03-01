"""
Main experiment runner: executes all analyses for the paper.

Experiments:
1. Single-layer baseline characterization
2. Multi-layer combination analysis
3. Correlation sensitivity analysis
4. Optimal design search
5. Validation against published data
6. Ablation studies
"""

import numpy as np
import json
import os
from pathlib import Path
from .models import (
    KillSwitchLayer,
    MultiLayerKillSwitch,
    deterministic_escape_dynamics,
    stochastic_escape_simulation,
    analytical_escape_probability,
    time_to_escape,
    sweep_escape_rates,
)
from .parameters import (
    ESCAPE_RATES,
    FITNESS_ADVANTAGE_ESCAPER,
    DEFAULT_SEED,
    DEFAULT_N0,
    DEFAULT_GENERATIONS,
    DEFAULT_REPLICATES,
    MU_MAX,
    K,
    SERIAL_PASSAGE,
    CORRELATION_FACTORS,
)

RESULTS_DIR = Path(__file__).parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)


def _save(name, data):
    path = RESULTS_DIR / f"{name}.json"

    # Convert numpy types for JSON serialization
    def convert(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if hasattr(obj, "item") and callable(obj.item):
            return obj.item()
        if isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [convert(v) for v in obj]
        return obj

    with open(path, "w") as f:
        json.dump(convert(data), f, indent=2)
    print(f"  Saved: {path}")


# =========================================================================
# Experiment 1: Single-layer baseline characterization
# =========================================================================
def exp1_single_layer_baselines():
    """Characterize each kill switch architecture independently."""
    print("=" * 60)
    print("Experiment 1: Single-layer baseline characterization")
    print("=" * 60)

    results = {}
    for name, params in ESCAPE_RATES.items():
        layer = KillSwitchLayer(
            name=name,
            escape_rate=params["rate"],
            target_size_bp=params["target_size_bp"],
        )
        ks = MultiLayerKillSwitch(layers=[layer])

        # Deterministic dynamics
        det = deterministic_escape_dynamics(
            N0=DEFAULT_N0,
            generations=DEFAULT_GENERATIONS,
            kill_switch=ks,
            mu_max=MU_MAX,
            K=K,
        )

        # Stochastic replicates
        escape_times = []
        final_fractions = []
        rng = np.random.default_rng(DEFAULT_SEED)
        for rep in range(DEFAULT_REPLICATES):
            stoch = stochastic_escape_simulation(
                N0=int(DEFAULT_N0),
                generations=DEFAULT_GENERATIONS,
                kill_switch=ks,
                mu_max=MU_MAX,
                K=K,
                dilution_factor=100,
                passage_interval_gen=6.6,
                rng=rng,
            )
            final_fractions.append(float(stoch["escape_fraction"][-1]))
            # Find generation where escape > 50%
            above_50 = np.where(stoch["escape_fraction"] > 0.5)[0]
            if len(above_50) > 0:
                escape_times.append(float(stoch["generations"][above_50[0]]))
            else:
                escape_times.append(float(DEFAULT_GENERATIONS))

        # Analytical
        p_escape = analytical_escape_probability(
            params["rate"],
            DEFAULT_GENERATIONS,
            DEFAULT_N0,
            FITNESS_ADVANTAGE_ESCAPER,
        )
        t_escape_analytical = time_to_escape(
            params["rate"],
            DEFAULT_N0,
            FITNESS_ADVANTAGE_ESCAPER,
        )

        results[name] = {
            "escape_rate": params["rate"],
            "log10_escape_rate": float(np.log10(params["rate"])),
            "mechanism": params["mechanism"],
            "deterministic_final_escape_fraction": float(det["escape_fraction"][-1]),
            "stochastic_median_final_fraction": float(np.median(final_fractions)),
            "stochastic_mean_escape_time_gen": float(np.mean(escape_times)),
            "stochastic_std_escape_time_gen": float(np.std(escape_times)),
            "analytical_escape_prob": float(p_escape),
            "analytical_time_to_50pct": float(t_escape_analytical),
            "deterministic_timeseries": {
                "generations": det["generations"].tolist(),
                "escape_fraction": det["escape_fraction"].tolist(),
            },
        }
        print(
            f"  {name}: escape_rate={params['rate']:.1e}, "
            f"t_50%={np.mean(escape_times):.1f} gen, "
            f"analytical_t_50%={t_escape_analytical:.1f} gen"
        )

    _save("exp1_single_layer_baselines", results)
    return results


# =========================================================================
# Experiment 2: Multi-layer combination analysis
# =========================================================================
def exp2_multilayer_combinations():
    """Test all pairwise and triple combinations of kill switch layers."""
    print("\n" + "=" * 60)
    print("Experiment 2: Multi-layer combination analysis")
    print("=" * 60)

    architectures = [
        ("toxin_antitoxin", ESCAPE_RATES["toxin_antitoxin"]),
        ("crispr_multi", ESCAPE_RATES["crispr_multi"]),
        ("overlapping_gene", ESCAPE_RATES["overlapping_gene"]),
        ("auxotrophy", ESCAPE_RATES["auxotrophy"]),
    ]

    results = {"pairwise": {}, "triple": {}}
    rng = np.random.default_rng(DEFAULT_SEED)

    # Pairwise combinations
    for i in range(len(architectures)):
        for j in range(i + 1, len(architectures)):
            name_i, params_i = architectures[i]
            name_j, params_j = architectures[j]
            combo_name = f"{name_i}+{name_j}"

            layer_i = KillSwitchLayer(
                name=name_i,
                escape_rate=params_i["rate"],
                target_size_bp=params_i["target_size_bp"],
            )
            layer_j = KillSwitchLayer(
                name=name_j,
                escape_rate=params_j["rate"],
                target_size_bp=params_j["target_size_bp"],
            )

            # Test with different correlation levels
            for rho_label, rho in [
                ("independent", 0.0),
                ("low", 0.01),
                ("medium", 0.1),
                ("high", 0.3),
            ]:
                ks = MultiLayerKillSwitch(layers=[layer_i, layer_j], correlation=rho)
                combined_rate = ks.combined_escape_rate()

                det = deterministic_escape_dynamics(
                    N0=DEFAULT_N0,
                    generations=DEFAULT_GENERATIONS,
                    kill_switch=ks,
                    mu_max=MU_MAX,
                    K=K,
                )

                # Stochastic
                final_fracs = []
                esc_times = []
                for rep in range(DEFAULT_REPLICATES):
                    stoch = stochastic_escape_simulation(
                        N0=int(DEFAULT_N0),
                        generations=DEFAULT_GENERATIONS,
                        kill_switch=ks,
                        mu_max=MU_MAX,
                        K=K,
                        dilution_factor=100,
                        passage_interval_gen=6.6,
                        rng=rng,
                    )
                    final_fracs.append(float(stoch["escape_fraction"][-1]))
                    above = np.where(stoch["escape_fraction"] > 0.5)[0]
                    esc_times.append(
                        float(stoch["generations"][above[0]])
                        if len(above) > 0
                        else float(DEFAULT_GENERATIONS)
                    )

                key = f"{combo_name}_rho={rho_label}"
                results["pairwise"][key] = {
                    "layers": [name_i, name_j],
                    "correlation": rho,
                    "combined_escape_rate": float(combined_rate),
                    "log10_combined": float(np.log10(combined_rate))
                    if combined_rate > 0
                    else -30,
                    "independent_product": float(params_i["rate"] * params_j["rate"]),
                    "det_final_escape": float(det["escape_fraction"][-1]),
                    "stoch_median_final": float(np.median(final_fracs)),
                    "stoch_mean_t50": float(np.mean(esc_times)),
                    "meets_NIH_threshold": combined_rate < 1e-8,
                }
                if rho_label == "low":
                    print(
                        f"  {combo_name} (rho={rho}): rate={combined_rate:.2e}, "
                        f"log10={np.log10(combined_rate):.1f}, "
                        f"NIH={'PASS' if combined_rate < 1e-8 else 'FAIL'}"
                    )

    # Triple combinations (all 4 choose 3 = 4 triples)
    from itertools import combinations

    for combo in combinations(range(len(architectures)), 3):
        names = [architectures[k][0] for k in combo]
        combo_name = "+".join(names)
        layers = [
            KillSwitchLayer(
                name=architectures[k][0],
                escape_rate=architectures[k][1]["rate"],
                target_size_bp=architectures[k][1]["target_size_bp"],
            )
            for k in combo
        ]

        for rho_label, rho in [("independent", 0.0), ("low", 0.01)]:
            ks = MultiLayerKillSwitch(layers=layers, correlation=rho)
            combined_rate = ks.combined_escape_rate()
            t_esc = time_to_escape(combined_rate, DEFAULT_N0, FITNESS_ADVANTAGE_ESCAPER)

            key = f"{combo_name}_rho={rho_label}"
            results["triple"][key] = {
                "layers": names,
                "correlation": rho,
                "combined_escape_rate": float(combined_rate),
                "log10_combined": float(np.log10(combined_rate))
                if combined_rate > 0
                else -30,
                "time_to_50pct": float(t_esc),
                "meets_NIH_threshold": combined_rate < 1e-8,
            }
            if rho_label == "low":
                print(
                    f"  {combo_name} (rho={rho}): rate={combined_rate:.2e}, "
                    f"NIH={'PASS' if combined_rate < 1e-8 else 'FAIL'}"
                )

    _save("exp2_multilayer_combinations", results)
    return results


# =========================================================================
# Experiment 3: Correlation sensitivity analysis
# =========================================================================
def exp3_correlation_sensitivity():
    """How does correlation between layers affect combined escape rate?"""
    print("\n" + "=" * 60)
    print("Experiment 3: Correlation sensitivity analysis")
    print("=" * 60)

    rho_range = np.logspace(-4, 0, 50)
    rho_range = np.concatenate([[0], rho_range])

    # Test the most promising combination
    layer_configs = [
        {
            "name": "crispr_multi",
            "escape_rate": ESCAPE_RATES["crispr_multi"]["rate"],
            "target_size_bp": ESCAPE_RATES["crispr_multi"]["target_size_bp"],
        },
        {
            "name": "overlapping_gene",
            "escape_rate": ESCAPE_RATES["overlapping_gene"]["rate"],
            "target_size_bp": ESCAPE_RATES["overlapping_gene"]["target_size_bp"],
        },
    ]

    results = {
        "rho_values": [],
        "combined_rates": [],
        "log10_rates": [],
        "time_to_escape": [],
        "nih_threshold_rho": None,
    }

    for rho in rho_range:
        layers = [KillSwitchLayer(**cfg) for cfg in layer_configs]
        ks = MultiLayerKillSwitch(layers=layers, correlation=float(rho))
        rate = ks.combined_escape_rate()
        t_esc = time_to_escape(rate, DEFAULT_N0, FITNESS_ADVANTAGE_ESCAPER)

        results["rho_values"].append(float(rho))
        results["combined_rates"].append(float(rate))
        results["log10_rates"].append(float(np.log10(rate)) if rate > 0 else -30)
        results["time_to_escape"].append(float(t_esc))

        # Find where it crosses NIH threshold
        if rate > 1e-8 and results["nih_threshold_rho"] is None:
            results["nih_threshold_rho"] = float(rho)

    print(f"  NIH threshold crossed at rho={results['nih_threshold_rho']}")
    print(f"  At rho=0: rate={results['combined_rates'][0]:.2e}")
    print(f"  At rho=0.01: rate={results['combined_rates'][5]:.2e}")
    print(f"  At rho=1.0: rate={results['combined_rates'][-1]:.2e}")

    _save("exp3_correlation_sensitivity", results)
    return results


# =========================================================================
# Experiment 4: Optimal design search
# =========================================================================
def exp4_optimal_design():
    """Find the optimal combination that minimizes escape rate
    while keeping fitness cost below a threshold."""
    print("\n" + "=" * 60)
    print("Experiment 4: Optimal design search")
    print("=" * 60)

    from itertools import combinations

    all_layers = {
        name: KillSwitchLayer(
            name=name,
            escape_rate=params["rate"],
            target_size_bp=params["target_size_bp"],
            fitness_cost=0.03 if name == "auxotrophy" else 0.05,
        )
        for name, params in ESCAPE_RATES.items()
        if name not in ("crispr_single",)  # exclude weak single-gRNA
    }

    max_fitness_cost = 0.20  # 20% maximum acceptable growth penalty
    max_layers = 4
    rho = 0.01  # realistic low correlation

    results = []
    for n in range(1, max_layers + 1):
        for combo in combinations(all_layers.keys(), n):
            layers = [all_layers[name] for name in combo]
            ks = MultiLayerKillSwitch(layers=layers, correlation=rho)
            cost = ks.total_fitness_cost()

            if cost > max_fitness_cost:
                continue

            rate = ks.combined_escape_rate()
            t_esc = time_to_escape(rate, DEFAULT_N0, FITNESS_ADVANTAGE_ESCAPER)

            results.append(
                {
                    "layers": list(combo),
                    "n_layers": n,
                    "fitness_cost": float(cost),
                    "combined_escape_rate": float(rate),
                    "log10_rate": float(np.log10(rate)) if rate > 0 else -30,
                    "time_to_50pct_gen": float(t_esc),
                    "meets_NIH": rate < 1e-8,
                    "correlation": rho,
                }
            )

    # Sort by escape rate
    results.sort(key=lambda x: x["combined_escape_rate"])

    # Print top 10
    print(f"  Total valid combinations: {len(results)}")
    print(
        f"  Combinations meeting NIH (<10^-8): {sum(1 for r in results if r['meets_NIH'])}"
    )
    print("\n  Top 10 designs:")
    for i, r in enumerate(results[:10]):
        print(
            f"  {i + 1}. {'+'.join(r['layers'])}: rate={r['combined_escape_rate']:.2e}, "
            f"cost={r['fitness_cost']:.0%}, t50={r['time_to_50pct_gen']:.0f} gen"
        )

    _save(
        "exp4_optimal_design",
        {"designs": results, "max_fitness_cost": max_fitness_cost, "correlation": rho},
    )
    return results


# =========================================================================
# Experiment 5: Validation against published experimental data
# =========================================================================
def exp5_validation():
    """Compare model predictions against published experimental results."""
    print("\n" + "=" * 60)
    print("Experiment 5: Validation against published data")
    print("=" * 60)

    # Validation dataset from literature
    validation_data = [
        {
            "source": "Rottinghaus 2022 - CRISPR single gRNA",
            "observed_escape_rate": 3e-5,
            "observed_escape_rate_range": (1e-5, 1e-4),
            "condition": "single gRNA, E. coli Nissle",
        },
        {
            "source": "Rottinghaus 2022 - CRISPR optimized 2-gRNA",
            "observed_escape_rate": 2.5e-9,
            "observed_escape_rate_range": (1e-9, 1e-8),
            "condition": "2-gRNA + SOS knockout, E. coli Nissle",
        },
        {
            "source": "Chlebek 2023 - Toxin-antitoxin baseline",
            "observed_escape_rate": 1e-6,
            "observed_escape_rate_range": (5e-7, 5e-6),
            "condition": "RelE toxin, P. protegens",
        },
        {
            "source": "Chlebek 2023 - Overlapping gene",
            "observed_escape_rate": 1.4e-7,
            "observed_escape_rate_range": (5e-8, 5e-7),
            "condition": "ilvA/relE entanglement, -ile, P. protegens",
        },
        {
            "source": "Chlebek 2023 - Auxotrophy",
            "observed_escape_rate": 1e-9,
            "observed_escape_rate_range": (1e-10, 1e-8),
            "condition": "delta-ilvA delta-TdcB, P. protegens",
        },
        {
            "source": "Chlebek 2023 - Time to 10% escape (no entanglement, +ile)",
            "observed_value": 30,
            "metric": "generations_to_10pct_escape",
            "condition": "relE, +ile, serial passage 1:100",
        },
        {
            "source": "Chlebek 2023 - Time to 10% escape (entangled, -ile)",
            "observed_value": 130,
            "metric": "generations_to_10pct_escape_lower_bound",
            "condition": "ilvA/relE entangled, -ile, serial passage 1:100",
        },
    ]

    results = []
    rng = np.random.default_rng(DEFAULT_SEED)

    for vd in validation_data:
        if "observed_escape_rate" in vd:
            # Model the corresponding architecture
            if "single gRNA" in vd["source"]:
                layer = KillSwitchLayer(
                    name="crispr_single",
                    escape_rate=ESCAPE_RATES["crispr_single"]["rate"],
                    target_size_bp=50,
                )
            elif "2-gRNA" in vd["source"]:
                layer = KillSwitchLayer(
                    name="crispr_multi",
                    escape_rate=ESCAPE_RATES["crispr_multi"]["rate"],
                    target_size_bp=150,
                )
            elif "baseline" in vd["source"]:
                layer = KillSwitchLayer(
                    name="toxin_antitoxin",
                    escape_rate=ESCAPE_RATES["toxin_antitoxin"]["rate"],
                    target_size_bp=300,
                )
            elif "Overlapping" in vd["source"]:
                layer = KillSwitchLayer(
                    name="overlapping_gene",
                    escape_rate=ESCAPE_RATES["overlapping_gene"]["rate"],
                    target_size_bp=288,
                )
            elif "Auxotrophy" in vd["source"]:
                layer = KillSwitchLayer(
                    name="auxotrophy",
                    escape_rate=ESCAPE_RATES["auxotrophy"]["rate"],
                    target_size_bp=0,
                    is_element_susceptible=False,
                )
            else:
                continue

            ks = MultiLayerKillSwitch(layers=[layer])
            predicted = ks.combined_escape_rate()

            observed = vd["observed_escape_rate"]
            log_ratio = np.log10(predicted / observed)

            results.append(
                {
                    "source": vd["source"],
                    "observed": float(observed),
                    "predicted": float(predicted),
                    "log10_observed": float(np.log10(observed)),
                    "log10_predicted": float(np.log10(predicted)),
                    "log10_ratio": float(log_ratio),
                    "within_order_of_magnitude": abs(log_ratio) < 1.0,
                    "condition": vd["condition"],
                }
            )
            status = "OK" if abs(log_ratio) < 1.0 else "MISMATCH"
            print(
                f"  [{status}] {vd['source']}: obs={observed:.1e}, pred={predicted:.1e}, "
                f"ratio={10**log_ratio:.2f}x"
            )

        elif "metric" in vd and "generations" in vd["metric"]:
            # Validate time-to-escape predictions
            if "no entanglement" in vd["source"]:
                layer = KillSwitchLayer(
                    name="toxin_antitoxin",
                    escape_rate=ESCAPE_RATES["toxin_antitoxin"]["rate"],
                    target_size_bp=300,
                )
            else:
                layer = KillSwitchLayer(
                    name="overlapping_gene",
                    escape_rate=ESCAPE_RATES["overlapping_gene"]["rate"],
                    target_size_bp=288,
                )

            ks = MultiLayerKillSwitch(layers=[layer])
            # Simulate to find time to 10% escape
            esc_times = []
            for rep in range(DEFAULT_REPLICATES):
                stoch = stochastic_escape_simulation(
                    N0=int(1e6),
                    generations=500,
                    kill_switch=ks,
                    mu_max=MU_MAX,
                    K=K,
                    dilution_factor=100,
                    passage_interval_gen=6.6,
                    rng=rng,
                )
                above_10 = np.where(stoch["escape_fraction"] > 0.1)[0]
                if len(above_10) > 0:
                    esc_times.append(float(stoch["generations"][above_10[0]]))
                else:
                    esc_times.append(500.0)

            predicted_t = float(np.median(esc_times))
            observed_t = vd["observed_value"]

            results.append(
                {
                    "source": vd["source"],
                    "observed_generations": observed_t,
                    "predicted_generations": predicted_t,
                    "ratio": predicted_t / observed_t,
                    "condition": vd["condition"],
                }
            )
            print(
                f"  {vd['source']}: obs={observed_t} gen, pred={predicted_t:.0f} gen, "
                f"ratio={predicted_t / observed_t:.2f}x"
            )

    # Overall validation metrics
    escape_rate_results = [r for r in results if "log10_ratio" in r]
    mean_abs_log_error = None
    all_within_oom = None
    if escape_rate_results:
        log_ratios = [abs(r["log10_ratio"]) for r in escape_rate_results]
        mean_abs_log_error = float(np.mean(log_ratios))
        all_within_oom = all(
            r["within_order_of_magnitude"] for r in escape_rate_results
        )
        print(f"\n  Mean |log10(pred/obs)| = {mean_abs_log_error:.3f}")
        print(f"  All within 1 order of magnitude: {all_within_oom}")

    _save(
        "exp5_validation",
        {
            "results": results,
            "mean_abs_log_error": mean_abs_log_error,
            "all_within_oom": all_within_oom,
        },
    )
    return results


# =========================================================================
# Experiment 6: Ablation study
# =========================================================================
def exp6_ablation():
    """Ablation: remove each component to measure its contribution."""
    print("\n" + "=" * 60)
    print("Experiment 6: Ablation study")
    print("=" * 60)

    # Use the best triple combination
    full_layers = [
        KillSwitchLayer(
            name="crispr_multi",
            escape_rate=ESCAPE_RATES["crispr_multi"]["rate"],
            target_size_bp=ESCAPE_RATES["crispr_multi"]["target_size_bp"],
        ),
        KillSwitchLayer(
            name="overlapping_gene",
            escape_rate=ESCAPE_RATES["overlapping_gene"]["rate"],
            target_size_bp=ESCAPE_RATES["overlapping_gene"]["target_size_bp"],
        ),
        KillSwitchLayer(
            name="auxotrophy",
            escape_rate=ESCAPE_RATES["auxotrophy"]["rate"],
            target_size_bp=0,
            is_element_susceptible=False,
            fitness_cost=0.03,
        ),
    ]

    rho = 0.01
    full_ks = MultiLayerKillSwitch(layers=full_layers, correlation=rho)
    full_rate = full_ks.combined_escape_rate()
    full_cost = full_ks.total_fitness_cost()

    results = {
        "full_system": {
            "layers": [l.name for l in full_layers],
            "escape_rate": float(full_rate),
            "log10_rate": float(np.log10(full_rate)) if full_rate > 0 else -30,
            "fitness_cost": float(full_cost),
        },
        "ablations": [],
    }

    print(f"  Full system: rate={full_rate:.2e}, cost={full_cost:.0%}")

    for i, removed in enumerate(full_layers):
        remaining = [l for j, l in enumerate(full_layers) if j != i]
        ks = MultiLayerKillSwitch(layers=remaining, correlation=rho)
        rate = ks.combined_escape_rate()
        cost = ks.total_fitness_cost()
        contribution = np.log10(rate / full_rate) if full_rate > 0 and rate > 0 else 0

        ablation = {
            "removed": removed.name,
            "remaining": [l.name for l in remaining],
            "escape_rate": float(rate),
            "log10_rate": float(np.log10(rate)) if rate > 0 else -30,
            "fitness_cost": float(cost),
            "log10_rate_increase": float(contribution),
            "fold_increase": float(rate / full_rate) if full_rate > 0 else float("inf"),
        }
        results["ablations"].append(ablation)
        print(
            f"  Remove {removed.name}: rate={rate:.2e} ({rate / full_rate:.0f}x worse), "
            f"cost={cost:.0%}"
        )

    _save("exp6_ablation", results)
    return results


# =========================================================================
# Experiment 7: Parameter sensitivity analysis
# =========================================================================
def exp7_sensitivity():
    """Sensitivity analysis: vary each parameter and measure effect on escape rate."""
    print("\n" + "=" * 60)
    print("Experiment 7: Parameter sensitivity analysis")
    print("=" * 60)

    # Base configuration: best triple
    base_rates = {
        "crispr_multi": ESCAPE_RATES["crispr_multi"]["rate"],
        "overlapping_gene": ESCAPE_RATES["overlapping_gene"]["rate"],
        "auxotrophy": ESCAPE_RATES["auxotrophy"]["rate"],
    }
    rho = 0.01

    # Vary each parameter by factors of 0.01x to 100x
    factors = np.logspace(-2, 2, 50)
    results = {}

    for vary_name in base_rates:
        rates_vs_factor = []
        for f in factors:
            layers = []
            for name, rate in base_rates.items():
                r = rate * f if name == vary_name else rate
                layers.append(
                    KillSwitchLayer(
                        name=name,
                        escape_rate=r,
                        target_size_bp=ESCAPE_RATES[name]["target_size_bp"],
                        is_element_susceptible=(name != "auxotrophy"),
                        fitness_cost=0.03 if name == "auxotrophy" else 0.05,
                    )
                )
            ks = MultiLayerKillSwitch(layers=layers, correlation=rho)
            rates_vs_factor.append(float(ks.combined_escape_rate()))

        results[vary_name] = {
            "factors": factors.tolist(),
            "combined_rates": rates_vs_factor,
            "log10_rates": [
                float(np.log10(r)) if r > 0 else -30 for r in rates_vs_factor
            ],
        }
        # Compute local sensitivity: d(log rate) / d(log factor) at factor=1
        idx_base = len(factors) // 2
        if idx_base > 0:
            dlr = np.log10(rates_vs_factor[idx_base + 1]) - np.log10(
                rates_vs_factor[idx_base - 1]
            )
            dlf = np.log10(factors[idx_base + 1]) - np.log10(factors[idx_base - 1])
            sensitivity = dlr / dlf
        else:
            sensitivity = 0
        results[vary_name]["local_sensitivity"] = float(sensitivity)
        print(f"  Sensitivity to {vary_name}: {sensitivity:.3f}")

    # Also vary correlation
    rho_range = np.logspace(-4, 0, 50)
    rho_rates = []
    for r in rho_range:
        layers = [
            KillSwitchLayer(
                name=n,
                escape_rate=base_rates[n],
                target_size_bp=ESCAPE_RATES[n]["target_size_bp"],
                is_element_susceptible=(n != "auxotrophy"),
                fitness_cost=0.03 if n == "auxotrophy" else 0.05,
            )
            for n in base_rates
        ]
        ks = MultiLayerKillSwitch(layers=layers, correlation=float(r))
        rho_rates.append(float(ks.combined_escape_rate()))

    results["correlation"] = {
        "rho_values": rho_range.tolist(),
        "combined_rates": rho_rates,
        "log10_rates": [float(np.log10(r)) if r > 0 else -30 for r in rho_rates],
    }

    _save("exp7_sensitivity", results)
    return results


# =========================================================================
# Run all experiments
# =========================================================================
def run_all():
    """Run all experiments and save results."""
    print("Kill Switch Design Framework - Running All Experiments")
    print("=" * 60)
    exp1_single_layer_baselines()
    exp2_multilayer_combinations()
    exp3_correlation_sensitivity()
    exp4_optimal_design()
    exp5_validation()
    exp6_ablation()
    exp7_sensitivity()
    print("\n" + "=" * 60)
    print("All experiments complete. Results saved to results/")
    print("=" * 60)


if __name__ == "__main__":
    run_all()
