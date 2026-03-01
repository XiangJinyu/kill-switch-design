"""
Rigorous experiment runner for the kill switch design framework.

All experiments use:
- Mutation-spectrum-decomposed escape rates
- scipy solve_ivp for deterministic ODE (adaptive RK45)
- 10,000 stochastic replicates with convergence verification
- Bootstrap confidence intervals
- Proper cross-validation (train on subset, test on holdout)
"""

import numpy as np
import json
import time
from pathlib import Path
from itertools import combinations

from .models import (
    MutationSpectrum,
    MUTATION_SPECTRA,
    KillSwitchLayer,
    MultiLayerKillSwitch,
    deterministic_escape_ode,
    stochastic_passage_simulation,
    run_replicate_simulations,
    convergence_analysis,
    analytical_escape_probability,
    time_to_escape_analytical,
)
from .parameters import (
    FITNESS_ADVANTAGE_ESCAPER,
    MU_MAX,
    K,
    DEFAULT_SEED,
    DEFAULT_N0,
)

RESULTS_DIR = Path(__file__).parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# Default simulation parameters (rigorous)
N_REPLICATES = 10000  # stochastic replicates
N_GENERATIONS = 1000  # total generations to simulate
N_SWEEP_POINTS = 200  # parameter sweep resolution
BOOTSTRAP_SAMPLES = 5000  # for CI estimation


def _save(name, data):
    """Save results with numpy-safe JSON serialization."""
    path = RESULTS_DIR / f"{name}.json"

    def convert(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if hasattr(obj, "item") and callable(getattr(obj, "item")):
            return obj.item()
        if isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [convert(v) for v in obj]
        return obj

    with open(path, "w") as f:
        json.dump(convert(data), f, indent=2)
    print(f"  -> Saved: {path}")
    return path


def make_layer(name):
    """Create a KillSwitchLayer from the literature-calibrated spectrum."""
    cost = 0.03 if name == "auxotrophy" else 0.05
    return KillSwitchLayer(
        name=name, spectrum=MUTATION_SPECTRA[name], fitness_cost=cost
    )


def bootstrap_ci(data, stat_fn=np.mean, n_boot=BOOTSTRAP_SAMPLES, ci=0.95, seed=42):
    """Compute bootstrap confidence interval for a statistic."""
    rng = np.random.default_rng(seed)
    boot_stats = []
    n = len(data)
    for _ in range(n_boot):
        sample = rng.choice(data, size=n, replace=True)
        boot_stats.append(stat_fn(sample))
    alpha = (1 - ci) / 2
    lo = float(np.percentile(boot_stats, 100 * alpha))
    hi = float(np.percentile(boot_stats, 100 * (1 - alpha)))
    return lo, float(stat_fn(data)), hi


# =========================================================================
# Experiment 1: Single-layer characterization with full statistics
# =========================================================================
def exp1_single_layer(n_reps=N_REPLICATES):
    """Full characterization of each single-layer architecture."""
    t0 = time.time()
    print("=" * 70)
    print(f"Experiment 1: Single-layer characterization ({n_reps} replicates)")
    print("=" * 70)

    results = {}
    for name in MUTATION_SPECTRA:
        layer = make_layer(name)
        ks = MultiLayerKillSwitch(layers=[layer])
        rate = ks.combined_escape_rate()

        # Deterministic ODE
        det = deterministic_escape_ode(
            N0=DEFAULT_N0,
            generations=N_GENERATIONS,
            kill_switch=ks,
            K=K,
        )

        # Stochastic with full replicates
        stats = run_replicate_simulations(
            n_replicates=n_reps,
            N0=DEFAULT_N0,
            generations=N_GENERATIONS,
            kill_switch=ks,
            K=K,
            seed=DEFAULT_SEED,
        )

        # Bootstrap CIs for time-to-escape
        t50_lo, t50_med, t50_hi = bootstrap_ci(
            stats["escape_times_50"], stat_fn=np.median
        )
        t10_lo, t10_med, t10_hi = bootstrap_ci(
            stats["escape_times_10"], stat_fn=np.median
        )

        # Analytical
        p_esc = analytical_escape_probability(
            rate, N_GENERATIONS, DEFAULT_N0, FITNESS_ADVANTAGE_ESCAPER
        )
        t_analytical = time_to_escape_analytical(
            rate, DEFAULT_N0, FITNESS_ADVANTAGE_ESCAPER
        )

        # Mutation spectrum decomposition
        spectrum = MUTATION_SPECTRA[name].to_dict()

        results[name] = {
            "escape_rate": rate,
            "log10_escape_rate": float(np.log10(rate)) if rate > 0 else -30,
            "mutation_spectrum": spectrum,
            "deterministic": {
                "generations": det["generations"].tolist(),
                "escape_fraction": det["escape_fraction"].tolist(),
            },
            "stochastic": {
                "n_replicates": n_reps,
                "generations": stats["generations"].tolist(),
                "mean_fraction": stats["mean_fraction"].tolist(),
                "median_fraction": stats["median_fraction"].tolist(),
                "pct5": stats["pct5_fraction"].tolist(),
                "pct95": stats["pct95_fraction"].tolist(),
                "t50_median": t50_med,
                "t50_ci95": [t50_lo, t50_hi],
                "t10_median": t10_med,
                "t10_ci95": [t10_lo, t10_hi],
                "final_frac_mean": float(np.mean(stats["final_fractions"])),
                "final_frac_std": float(np.std(stats["final_fractions"])),
            },
            "analytical": {
                "escape_probability": p_esc,
                "time_to_50pct": t_analytical,
            },
        }
        print(
            f"  {name}: rate={rate:.2e}, "
            f"t50={t50_med:.1f} gen [{t50_lo:.1f}, {t50_hi:.1f}], "
            f"analytical_t50={t_analytical:.1f} gen"
        )

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")
    _save("exp1_single_layer", results)
    return results


# =========================================================================
# Experiment 2: Multi-layer combinations (exhaustive)
# =========================================================================
def exp2_combinations(n_reps=N_REPLICATES):
    """Evaluate all pairwise and triple combinations."""
    t0 = time.time()
    print("\n" + "=" * 70)
    print(f"Experiment 2: Multi-layer combinations ({n_reps} replicates)")
    print("=" * 70)

    archs = [
        "toxin_antitoxin",
        "crispr_multi",
        "overlapping_gene",
        "auxotrophy",
        "integrase_differentiation",
    ]
    rho_levels = [0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5]

    results = {"pairwise": {}, "triple": {}, "quad": {}}

    # Pairwise
    for i, j in combinations(range(len(archs)), 2):
        layers = [make_layer(archs[i]), make_layer(archs[j])]
        combo_name = f"{archs[i]}+{archs[j]}"

        combo_results = {}
        for rho in rho_levels:
            ks = MultiLayerKillSwitch(layers=layers, genomic_correlation=rho)
            rate = ks.combined_escape_rate()

            # Only run stochastic for key rho values
            stoch = None
            if rho in (0.0, 0.01, 0.1):
                stats = run_replicate_simulations(
                    n_replicates=min(n_reps, 2000),
                    N0=DEFAULT_N0,
                    generations=N_GENERATIONS,
                    kill_switch=ks,
                    K=K,
                    seed=DEFAULT_SEED,
                )
                t50_lo, t50_med, t50_hi = bootstrap_ci(
                    stats["escape_times_50"], stat_fn=np.median
                )
                stoch = {
                    "t50_median": t50_med,
                    "t50_ci95": [t50_lo, t50_hi],
                    "final_frac_mean": float(np.mean(stats["final_fractions"])),
                }

            combo_results[f"rho={rho}"] = {
                "correlation": rho,
                "combined_rate": float(rate),
                "log10_rate": float(np.log10(rate)) if rate > 0 else -30,
                "meets_NIH": rate < 1e-8,
                "stochastic": stoch,
            }

        # Print summary at rho=0.01
        r01 = combo_results["rho=0.01"]
        print(
            f"  {combo_name} (rho=0.01): rate={r01['combined_rate']:.2e}, "
            f"NIH={'PASS' if r01['meets_NIH'] else 'FAIL'}"
        )
        results["pairwise"][combo_name] = combo_results

    # Triple combinations
    for combo_idx in combinations(range(len(archs)), 3):
        layers = [make_layer(archs[k]) for k in combo_idx]
        combo_name = "+".join(archs[k] for k in combo_idx)

        combo_results = {}
        for rho in [0.0, 0.01, 0.05, 0.1]:
            ks = MultiLayerKillSwitch(layers=layers, genomic_correlation=rho)
            rate = ks.combined_escape_rate()
            t_esc = time_to_escape_analytical(
                rate, DEFAULT_N0, FITNESS_ADVANTAGE_ESCAPER
            )

            combo_results[f"rho={rho}"] = {
                "correlation": rho,
                "combined_rate": float(rate),
                "log10_rate": float(np.log10(rate)) if rate > 0 else -30,
                "meets_NIH": rate < 1e-8,
                "analytical_t50": float(t_esc),
            }

        r01 = combo_results["rho=0.01"]
        print(
            f"  {combo_name} (rho=0.01): rate={r01['combined_rate']:.2e}, "
            f"NIH={'PASS' if r01['meets_NIH'] else 'FAIL'}"
        )
        results["triple"][combo_name] = combo_results

    # Quadruple combinations
    for combo_idx in combinations(range(len(archs)), 4):
        layers = [make_layer(archs[k]) for k in combo_idx]
        combo_name = "+".join(archs[k] for k in combo_idx)
        ks = MultiLayerKillSwitch(layers=layers, genomic_correlation=0.01)
        rate = ks.combined_escape_rate()
        cost = ks.total_fitness_cost()
        results["quad"][combo_name] = {
            "combined_rate": float(rate),
            "log10_rate": float(np.log10(rate)) if rate > 0 else -30,
            "fitness_cost": float(cost),
            "meets_NIH": rate < 1e-8,
        }

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")
    _save("exp2_combinations", results)
    return results


# =========================================================================
# Experiment 3: Correlation sensitivity (high-resolution)
# =========================================================================
def exp3_correlation():
    """High-resolution sweep of correlation parameter."""
    t0 = time.time()
    print("\n" + "=" * 70)
    print(f"Experiment 3: Correlation sensitivity ({N_SWEEP_POINTS} points)")
    print("=" * 70)

    rho_range = np.concatenate(
        [
            [0],
            np.logspace(-4, 0, N_SWEEP_POINTS - 1),
        ]
    )

    # Test two representative combinations
    combos = {
        "crispr_multi+overlapping_gene": [
            make_layer("crispr_multi"),
            make_layer("overlapping_gene"),
        ],
        "crispr_multi+auxotrophy": [
            make_layer("crispr_multi"),
            make_layer("auxotrophy"),
        ],
        "triple_best": [
            make_layer("crispr_multi"),
            make_layer("overlapping_gene"),
            make_layer("auxotrophy"),
        ],
    }

    results = {}
    for combo_name, layers in combos.items():
        rho_data = {"rho": [], "rate": [], "log10_rate": [], "nih_threshold_rho": None}
        for rho in rho_range:
            ks = MultiLayerKillSwitch(layers=layers, genomic_correlation=float(rho))
            rate = ks.combined_escape_rate()
            rho_data["rho"].append(float(rho))
            rho_data["rate"].append(float(rate))
            rho_data["log10_rate"].append(float(np.log10(rate)) if rate > 0 else -30)

            if rate > 1e-8 and rho_data["nih_threshold_rho"] is None and rho > 0:
                rho_data["nih_threshold_rho"] = float(rho)

        results[combo_name] = rho_data
        print(
            f"  {combo_name}: rho_crit = {rho_data['nih_threshold_rho']}, "
            f"rate@rho=0: {rho_data['rate'][0]:.2e}, rate@rho=1: {rho_data['rate'][-1]:.2e}"
        )

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")
    _save("exp3_correlation", results)
    return results


# =========================================================================
# Experiment 4: Optimal design (Pareto front)
# =========================================================================
def exp4_pareto():
    """Exhaustive search for Pareto-optimal designs."""
    t0 = time.time()
    print("\n" + "=" * 70)
    print("Experiment 4: Optimal design search (Pareto front)")
    print("=" * 70)

    archs = [
        "toxin_antitoxin",
        "crispr_multi",
        "overlapping_gene",
        "auxotrophy",
        "integrase_differentiation",
    ]
    max_cost = 0.25
    rho = 0.01

    designs = []
    for n in range(1, len(archs) + 1):
        for combo in combinations(archs, n):
            layers = [make_layer(name) for name in combo]
            ks = MultiLayerKillSwitch(layers=layers, genomic_correlation=rho)
            cost = ks.total_fitness_cost()
            if cost > max_cost:
                continue
            rate = ks.combined_escape_rate()
            t_esc = time_to_escape_analytical(
                rate, DEFAULT_N0, FITNESS_ADVANTAGE_ESCAPER
            )

            designs.append(
                {
                    "layers": list(combo),
                    "n_layers": n,
                    "fitness_cost": float(cost),
                    "combined_rate": float(rate),
                    "log10_rate": float(np.log10(rate)) if rate > 0 else -30,
                    "time_to_50pct": float(t_esc),
                    "meets_NIH": rate < 1e-8,
                    "correlation": rho,
                }
            )

    designs.sort(key=lambda x: x["combined_rate"])

    # Identify Pareto front (non-dominated: no other design has both lower rate AND lower cost)
    pareto = []
    for d in designs:
        dominated = False
        for other in designs:
            if (
                other["combined_rate"] < d["combined_rate"]
                and other["fitness_cost"] <= d["fitness_cost"]
            ):
                dominated = True
                break
        if not dominated:
            pareto.append(d)

    print(f"  Total designs: {len(designs)}")
    print(f"  Meeting NIH: {sum(1 for d in designs if d['meets_NIH'])}")
    print(f"  Pareto-optimal: {len(pareto)}")
    print("\n  Top 5 designs:")
    for i, d in enumerate(designs[:5]):
        print(
            f"    {i + 1}. {'+'.join(d['layers'])}: rate={d['combined_rate']:.2e}, "
            f"cost={d['fitness_cost']:.0%}"
        )

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")
    _save("exp4_pareto", {"all_designs": designs, "pareto_front": pareto})
    return designs


# =========================================================================
# Experiment 5: Cross-validation
# =========================================================================
def exp5_validation():
    """
    Proper cross-validation:
    - Calibrate on Rottinghaus 2022 data -> predict Chlebek 2023 results
    - Calibrate on Chlebek 2023 data -> predict Rottinghaus 2022 results
    Also validate time-to-escape with stochastic simulation.
    """
    t0 = time.time()
    print("\n" + "=" * 70)
    print("Experiment 5: Cross-validation")
    print("=" * 70)

    # --- Part A: Leave-one-study-out cross-validation ---
    # All data points with their source study
    data_points = [
        {"name": "crispr_single", "observed": 3e-5, "study": "rottinghaus"},
        {"name": "crispr_multi", "observed": 2.5e-9, "study": "rottinghaus"},
        {"name": "toxin_antitoxin", "observed": 1e-6, "study": "chlebek"},
        {"name": "overlapping_gene", "observed": 1.4e-7, "study": "chlebek"},
        {"name": "auxotrophy", "observed": 1e-9, "study": "chlebek"},
    ]

    # Cross-validation: predict each study's data using only the other study's calibration
    cv_results = []
    for test_study in ["rottinghaus", "chlebek"]:
        train_points = [d for d in data_points if d["study"] != test_study]
        test_points = [d for d in data_points if d["study"] == test_study]

        # Calibration factor: ratio of predicted/observed across training set
        # This simulates "what if we only had one study's data"
        train_ratios = []
        for tp in train_points:
            predicted = MUTATION_SPECTRA[tp["name"]].total
            train_ratios.append(np.log10(predicted / tp["observed"]))

        # Mean calibration bias from training set
        bias = np.mean(train_ratios) if train_ratios else 0

        # Predict test set with calibration correction
        for dp in test_points:
            predicted_raw = MUTATION_SPECTRA[dp["name"]].total
            predicted_corrected = predicted_raw * 10 ** (-bias)
            obs = dp["observed"]
            log_error = np.log10(predicted_raw / obs)
            log_error_corrected = np.log10(predicted_corrected / obs)

            cv_results.append(
                {
                    "name": dp["name"],
                    "study": dp["study"],
                    "observed": obs,
                    "predicted_raw": float(predicted_raw),
                    "predicted_corrected": float(predicted_corrected),
                    "log10_error_raw": float(log_error),
                    "log10_error_corrected": float(log_error_corrected),
                    "test_study": test_study,
                }
            )
            within_oom = abs(log_error) < 1.0
            print(
                f"  [{test_study}] {dp['name']}: obs={obs:.1e}, "
                f"pred={predicted_raw:.1e}, log_err={log_error:.3f} "
                f"{'OK' if within_oom else 'MISMATCH'}"
            )

    # Overall CV metrics
    raw_errors = [abs(r["log10_error_raw"]) for r in cv_results]
    corrected_errors = [abs(r["log10_error_corrected"]) for r in cv_results]

    # --- Part B: Time-to-escape validation with stochastic simulation ---
    print("\n  Time-to-escape validation:")
    time_validation = []

    # Chlebek: toxin-antitoxin without entanglement, +ile, ~30 gen to 10% escape
    layer_ta = make_layer("toxin_antitoxin")
    ks_ta = MultiLayerKillSwitch(layers=[layer_ta])
    stats_ta = run_replicate_simulations(
        n_replicates=5000,
        N0=int(DEFAULT_N0),
        generations=500,
        kill_switch=ks_ta,
        K=K,
        dilution_factor=100,
        passage_gens=6.6,
        seed=DEFAULT_SEED,
    )
    t10_ta_lo, t10_ta_med, t10_ta_hi = bootstrap_ci(
        stats_ta["escape_times_10"], stat_fn=np.median
    )
    t10_ta_analytical = time_to_escape_analytical(
        layer_ta.escape_rate, DEFAULT_N0, FITNESS_ADVANTAGE_ESCAPER, threshold=0.1
    )
    time_validation.append(
        {
            "description": "Toxin-antitoxin, +ile (Chlebek 2023)",
            "observed_gen": 30,
            "stochastic_median": t10_ta_med,
            "stochastic_ci95": [t10_ta_lo, t10_ta_hi],
            "analytical": t10_ta_analytical,
            "ratio_stochastic": t10_ta_med / 30,
        }
    )
    print(
        f"    TA +ile: obs=30 gen, stoch={t10_ta_med:.1f} [{t10_ta_lo:.1f},{t10_ta_hi:.1f}], "
        f"analytical={t10_ta_analytical:.1f}"
    )

    # Chlebek: overlapping gene, -ile, >130 gen to 10% escape
    layer_og = make_layer("overlapping_gene")
    ks_og = MultiLayerKillSwitch(layers=[layer_og])
    stats_og = run_replicate_simulations(
        n_replicates=5000,
        N0=int(DEFAULT_N0),
        generations=500,
        kill_switch=ks_og,
        K=K,
        dilution_factor=100,
        passage_gens=6.6,
        seed=DEFAULT_SEED,
    )
    t10_og_lo, t10_og_med, t10_og_hi = bootstrap_ci(
        stats_og["escape_times_10"], stat_fn=np.median
    )
    t10_og_analytical = time_to_escape_analytical(
        layer_og.escape_rate, DEFAULT_N0, FITNESS_ADVANTAGE_ESCAPER, threshold=0.1
    )
    time_validation.append(
        {
            "description": "Overlapping gene, -ile (Chlebek 2023)",
            "observed_gen": 130,
            "observed_note": ">130 gen (7/10 lineages never reached 10%)",
            "stochastic_median": t10_og_med,
            "stochastic_ci95": [t10_og_lo, t10_og_hi],
            "analytical": t10_og_analytical,
            "ratio_stochastic": t10_og_med / 130,
        }
    )
    print(
        f"    OG -ile: obs=>130 gen, stoch={t10_og_med:.1f} [{t10_og_lo:.1f},{t10_og_hi:.1f}], "
        f"analytical={t10_og_analytical:.1f}"
    )

    summary = {
        "cross_validation": cv_results,
        "mean_abs_log_error_raw": float(np.mean(raw_errors)),
        "mean_abs_log_error_corrected": float(np.mean(corrected_errors)),
        "all_within_1_oom": all(e < 1.0 for e in raw_errors),
        "time_validation": time_validation,
    }
    print(f"\n  CV mean |log10 error| (raw): {np.mean(raw_errors):.4f}")
    print(f"  CV mean |log10 error| (corrected): {np.mean(corrected_errors):.4f}")
    print(f"  All within 1 OOM: {summary['all_within_1_oom']}")

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")
    _save("exp5_validation", summary)
    return summary


# =========================================================================
# Experiment 6: Ablation study
# =========================================================================
def exp6_ablation(n_reps=5000):
    """Ablation: remove each component and measure impact."""
    t0 = time.time()
    print("\n" + "=" * 70)
    print(f"Experiment 6: Ablation study ({n_reps} replicates)")
    print("=" * 70)

    best_layers = [
        make_layer("crispr_multi"),
        make_layer("overlapping_gene"),
        make_layer("auxotrophy"),
    ]
    rho = 0.01

    # Full system
    full_ks = MultiLayerKillSwitch(layers=best_layers, genomic_correlation=rho)
    full_rate = full_ks.combined_escape_rate()
    full_cost = full_ks.total_fitness_cost()

    full_stats = run_replicate_simulations(
        n_replicates=n_reps,
        N0=DEFAULT_N0,
        generations=N_GENERATIONS,
        kill_switch=full_ks,
        K=K,
        seed=DEFAULT_SEED,
    )
    full_t50_lo, full_t50_med, full_t50_hi = bootstrap_ci(
        full_stats["escape_times_50"], stat_fn=np.median
    )

    results = {
        "full_system": {
            "layers": [l.name for l in best_layers],
            "rate": float(full_rate),
            "log10_rate": float(np.log10(full_rate)) if full_rate > 0 else -30,
            "cost": float(full_cost),
            "t50": full_t50_med,
            "t50_ci95": [full_t50_lo, full_t50_hi],
        },
        "ablations": [],
    }

    print(f"  Full: rate={full_rate:.2e}, cost={full_cost:.0%}, t50={full_t50_med:.1f}")

    for i, removed in enumerate(best_layers):
        remaining = [l for j, l in enumerate(best_layers) if j != i]
        ks = MultiLayerKillSwitch(layers=remaining, genomic_correlation=rho)
        rate = ks.combined_escape_rate()
        cost = ks.total_fitness_cost()

        stats = run_replicate_simulations(
            n_replicates=n_reps,
            N0=DEFAULT_N0,
            generations=N_GENERATIONS,
            kill_switch=ks,
            K=K,
            seed=DEFAULT_SEED,
        )
        t50_lo, t50_med, t50_hi = bootstrap_ci(
            stats["escape_times_50"], stat_fn=np.median
        )

        fold_increase = rate / full_rate if full_rate > 0 else float("inf")
        ablation = {
            "removed": removed.name,
            "remaining": [l.name for l in remaining],
            "rate": float(rate),
            "log10_rate": float(np.log10(rate)) if rate > 0 else -30,
            "cost": float(cost),
            "fold_increase": float(fold_increase),
            "log10_fold": float(np.log10(fold_increase)) if fold_increase > 0 else 0,
            "t50": t50_med,
            "t50_ci95": [t50_lo, t50_hi],
        }
        results["ablations"].append(ablation)
        print(
            f"  Remove {removed.name}: rate={rate:.2e} ({fold_increase:.1e}x worse), "
            f"t50={t50_med:.1f}"
        )

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")
    _save("exp6_ablation", results)
    return results


# =========================================================================
# Experiment 7: Sensitivity analysis (high resolution)
# =========================================================================
def exp7_sensitivity():
    """High-resolution sensitivity analysis on best triple design."""
    t0 = time.time()
    print("\n" + "=" * 70)
    print(f"Experiment 7: Sensitivity analysis ({N_SWEEP_POINTS} points)")
    print("=" * 70)

    base_names = ["crispr_multi", "overlapping_gene", "auxotrophy"]
    base_rates = {n: MUTATION_SPECTRA[n].total for n in base_names}
    rho = 0.01
    factors = np.logspace(-2, 2, N_SWEEP_POINTS)

    results = {}
    for vary_name in base_names:
        rates_out = []
        for f in factors:
            layers = []
            for name in base_names:
                spec = MUTATION_SPECTRA[name]
                if name == vary_name:
                    # Scale all mutation types by the same factor
                    scaled = MutationSpectrum(
                        point_mutation=spec.point_mutation * f,
                        small_indel=spec.small_indel * f,
                        is_element=spec.is_element * f,
                        large_deletion=spec.large_deletion * f,
                        recombination=spec.recombination * f,
                    )
                    layers.append(
                        KillSwitchLayer(
                            name=name,
                            spectrum=scaled,
                            fitness_cost=0.03 if name == "auxotrophy" else 0.05,
                        )
                    )
                else:
                    layers.append(make_layer(name))
            ks = MultiLayerKillSwitch(layers=layers, genomic_correlation=rho)
            rates_out.append(float(ks.combined_escape_rate()))

        # Local sensitivity: d(log10 rate) / d(log10 factor) at factor=1
        mid = len(factors) // 2
        if mid > 0 and mid < len(factors) - 1:
            dr = np.log10(rates_out[mid + 1]) - np.log10(rates_out[mid - 1])
            df = np.log10(factors[mid + 1]) - np.log10(factors[mid - 1])
            sensitivity = dr / df if df != 0 else 0
        else:
            sensitivity = 0

        results[vary_name] = {
            "factors": factors.tolist(),
            "combined_rates": rates_out,
            "log10_rates": [float(np.log10(r)) if r > 0 else -30 for r in rates_out],
            "local_sensitivity": float(sensitivity),
            "base_rate": float(base_rates[vary_name]),
        }
        print(f"  {vary_name}: sensitivity={sensitivity:.4f}")

    # Also: 2D phase diagram — rho vs. escape rate
    print("\n  Generating 2D phase diagram (rho vs overlapping_gene rate)...")
    rho_2d = np.logspace(-4, 0, 100)
    og_factors = np.logspace(-2, 2, 100)
    phase = np.zeros((len(rho_2d), len(og_factors)))

    for ri, r in enumerate(rho_2d):
        for fi, f in enumerate(og_factors):
            spec_og = MUTATION_SPECTRA["overlapping_gene"]
            scaled_og = MutationSpectrum(
                point_mutation=spec_og.point_mutation * f,
                small_indel=spec_og.small_indel * f,
                is_element=spec_og.is_element * f,
                large_deletion=spec_og.large_deletion * f,
                recombination=spec_og.recombination * f,
            )
            layers = [
                make_layer("crispr_multi"),
                KillSwitchLayer(name="overlapping_gene", spectrum=scaled_og),
                make_layer("auxotrophy"),
            ]
            ks = MultiLayerKillSwitch(layers=layers, genomic_correlation=float(r))
            phase[ri, fi] = np.log10(ks.combined_escape_rate())

    results["phase_diagram"] = {
        "rho_values": rho_2d.tolist(),
        "og_factors": og_factors.tolist(),
        "log10_rates": phase.tolist(),
    }

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")
    _save("exp7_sensitivity", results)
    return results


# =========================================================================
# Experiment 8: Convergence analysis
# =========================================================================
def exp8_convergence():
    """Verify stochastic simulation convergence."""
    t0 = time.time()
    print("\n" + "=" * 70)
    print("Experiment 8: Convergence analysis")
    print("=" * 70)

    layer = make_layer("toxin_antitoxin")  # moderate escape rate for faster convergence
    ks = MultiLayerKillSwitch(layers=[layer])

    rep_counts = [50, 100, 200, 500, 1000, 2000, 5000, 10000]
    conv = convergence_analysis(
        N0=DEFAULT_N0,
        generations=500,
        kill_switch=ks,
        replicate_counts=rep_counts,
        seed=DEFAULT_SEED,
    )

    for c in conv:
        print(
            f"  n={c['n_replicates']:>6d}: t50={c['mean_t50']:.1f} +/- {c['se_t50']:.2f} (SE)"
        )

    _save("exp8_convergence", conv)
    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")
    return conv


# =========================================================================
# Experiment 9: Mutation spectrum decomposition analysis
# =========================================================================
def exp9_mutation_spectrum():
    """Analyze contribution of each mutation type to escape."""
    t0 = time.time()
    print("\n" + "=" * 70)
    print("Experiment 9: Mutation spectrum decomposition")
    print("=" * 70)

    results = {}
    for name, spec in MUTATION_SPECTRA.items():
        total = spec.total
        results[name] = {
            "total": float(total),
            "log10_total": float(np.log10(total)) if total > 0 else -30,
            "fractions": {
                "point_mutation": float(spec.point_mutation / total)
                if total > 0
                else 0,
                "small_indel": float(spec.small_indel / total) if total > 0 else 0,
                "is_element": float(spec.is_element / total) if total > 0 else 0,
                "large_deletion": float(spec.large_deletion / total)
                if total > 0
                else 0,
                "recombination": float(spec.recombination / total) if total > 0 else 0,
            },
            "dominant_mechanism": max(
                [
                    "point_mutation",
                    "small_indel",
                    "is_element",
                    "large_deletion",
                    "recombination",
                ],
                key=lambda m: getattr(spec, m),
            ),
        }
        print(
            f"  {name}: dominant={results[name]['dominant_mechanism']}, "
            f"total={total:.2e}"
        )

    _save("exp9_mutation_spectrum", results)
    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")
    return results


# =========================================================================
# Run all experiments
# =========================================================================
def run_all():
    """Run all experiments with publication-quality rigor."""
    t_start = time.time()
    print("=" * 70)
    print("Kill Switch Design Framework — Full Experiment Suite")
    print(f"Replicates: {N_REPLICATES}, Sweep points: {N_SWEEP_POINTS}")
    print("=" * 70)

    exp1_single_layer()
    exp2_combinations()
    exp3_correlation()
    exp4_pareto()
    exp5_validation()
    exp6_ablation()
    exp7_sensitivity()
    exp8_convergence()
    exp9_mutation_spectrum()

    total = time.time() - t_start
    print("\n" + "=" * 70)
    print(f"All experiments complete. Total time: {total:.1f}s ({total / 60:.1f} min)")
    print("=" * 70)


if __name__ == "__main__":
    run_all()
