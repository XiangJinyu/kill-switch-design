"""
Iteration 3 experiments: Out-of-sample validation + theoretical derivation support.
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
    run_replicate_simulations,
)
from .experiments import make_layer, _save, bootstrap_ci

RESULTS_DIR = Path(__file__).parent.parent / "results"


# =========================================================================
# Exp 14: Out-of-sample validation
# =========================================================================
def exp14_out_of_sample():
    """
    Predict escape rates for systems NOT used in calibration:
    1. Foo 2025: Entangled STALEMATE system in E. coli Nissle
    2. Hayashi 2024: thyA auxotrophy + Cas9 device in B. thetaiotaomicron
    """
    t0 = time.time()
    print("=" * 70)
    print("Experiment 14: Out-of-sample validation")
    print("=" * 70)

    predictions = []

    # =====================================================================
    # Foo 2025: STALEMATE system
    # =====================================================================
    print("\n  --- Foo et al. 2025 (STALEMATE, E. coli Nissle) ---")

    # Prediction 1: Non-entangled endonuclease (pEndo)
    # This is an inducible killing mechanism similar to our CRISPR single-gRNA
    # (single-component inducible kill, dominated by IS element inactivation)
    # Foo observed: ~10^-5 escape, dominated by IS911 (94.9% of escapers)
    # Our model: crispr_single has escape rate 2.9e-5 (IS + recombination dominated)
    # But the endonuclease is a single gene, more like toxin-antitoxin
    # Use TA model (single toxic protein, IS-susceptible, ~300bp target)
    pred1_rate = MUTATION_SPECTRA["toxin_antitoxin"].total
    obs1_rate = 1e-5  # Foo 2025, pEndo
    predictions.append(
        {
            "source": "Foo 2025 - Non-entangled endonuclease (pEndo)",
            "observed": obs1_rate,
            "predicted": float(pred1_rate),
            "log10_obs": float(np.log10(obs1_rate)),
            "log10_pred": float(np.log10(pred1_rate)),
            "log10_error": float(np.log10(pred1_rate / obs1_rate)),
            "model_used": "toxin_antitoxin",
            "reasoning": "Single inducible killing protein, IS-susceptible, similar to TA",
        }
    )
    print(
        f"  pEndo: obs={obs1_rate:.0e}, pred={pred1_rate:.1e}, "
        f"log_err={np.log10(pred1_rate / obs1_rate):.2f}"
    )

    # Prediction 2: Entangled (pEnt-533, initial)
    # Overlap reduces IS element escape by redirecting insertions to non-overlap region
    # Foo observed: ~10^-7 for pEnt-533 (reoptimized, IS hotspot removed)
    # Our model: overlapping_gene = 1.01e-7
    pred2_rate = MUTATION_SPECTRA["overlapping_gene"].total
    obs2_rate = 1e-7  # Foo 2025, pEnt-533 reoptimized
    predictions.append(
        {
            "source": "Foo 2025 - Entangled (pEnt-533 reop)",
            "observed": obs2_rate,
            "predicted": float(pred2_rate),
            "log10_obs": float(np.log10(obs2_rate)),
            "log10_pred": float(np.log10(pred2_rate)),
            "log10_error": float(np.log10(pred2_rate / obs2_rate)),
            "model_used": "overlapping_gene",
            "reasoning": "Gene entanglement with IS hotspot removal, ~28% overlap",
        }
    )
    print(
        f"  pEnt-533 reop: obs={obs2_rate:.0e}, pred={pred2_rate:.1e}, "
        f"log_err={np.log10(pred2_rate / obs2_rate):.2f}"
    )

    # Prediction 3: Entangled + strong ColE9 (pEnt-533, ColE9 0.47)
    # Additional population-level policing: ColE9 kills cells that lose Im9
    # Foo observed: ~10^-10
    # Our model: overlapping_gene + auxotrophy-like mechanism (ColE9 = essential)
    # This is a 2-layer system: entanglement + toxin policing
    og_layer = make_layer("overlapping_gene")
    # ColE9 policing acts like auxotrophy: losing the protection is lethal
    # But escape rate for this is harder to predict from first principles
    # Use our pairwise model: OG + Auxotrophy
    aux_layer = make_layer("auxotrophy")
    ks_2layer = MultiLayerKillSwitch(
        layers=[og_layer, aux_layer],
        genomic_correlation=0.05,  # on same plasmid
    )
    pred3_rate = ks_2layer.combined_escape_rate()
    obs3_rate = 1e-10  # Foo 2025, pEnt-533 (0.47)
    predictions.append(
        {
            "source": "Foo 2025 - Entangled + strong ColE9 (pEnt-533, 0.47)",
            "observed": obs3_rate,
            "predicted": float(pred3_rate),
            "log10_obs": float(np.log10(obs3_rate)),
            "log10_pred": float(np.log10(pred3_rate)),
            "log10_error": float(np.log10(pred3_rate / obs3_rate)),
            "model_used": "overlapping_gene + auxotrophy (rho=0.05)",
            "reasoning": "Entanglement + ColE9 policing modeled as OG + Aux 2-layer",
        }
    )
    print(
        f"  pEnt-533+ColE9: obs={obs3_rate:.0e}, pred={pred3_rate:.1e}, "
        f"log_err={np.log10(pred3_rate / obs3_rate):.2f}"
    )

    # Prediction 4: Foo in vivo (pEnt-533 in mouse gut)
    # Foo observed: 10^-6 to 10^-7 in vivo (~100x worse than in vitro)
    # Our model: OG = 1.01e-7 in vitro; applying 100x in vivo penalty
    obs4_rate = 3e-7  # geometric mean of 10^-6 to 10^-7
    pred4_rate_invitro = MUTATION_SPECTRA["overlapping_gene"].total
    # In vivo: stress-induced transposition increases IS rate
    # Foo data: ~100x higher escape in vivo vs in vitro
    # We don't model this directly, but note the in vitro prediction
    predictions.append(
        {
            "source": "Foo 2025 - Entangled in vivo (mouse gut)",
            "observed": obs4_rate,
            "predicted": float(pred4_rate_invitro),
            "log10_obs": float(np.log10(obs4_rate)),
            "log10_pred": float(np.log10(pred4_rate_invitro)),
            "log10_error": float(np.log10(pred4_rate_invitro / obs4_rate)),
            "model_used": "overlapping_gene (in vitro calibration)",
            "reasoning": "In vitro model applied to in vivo; expect underprediction",
            "note": "In vivo rates ~100x higher due to stress-induced transposition",
        }
    )
    print(
        f"  In vivo: obs={obs4_rate:.0e}, pred(in vitro)={pred4_rate_invitro:.1e}, "
        f"log_err={np.log10(pred4_rate_invitro / obs4_rate):.2f}"
    )

    # =====================================================================
    # Hayashi 2024: B. thetaiotaomicron
    # =====================================================================
    print("\n  --- Hayashi et al. 2024 (B. thetaiotaomicron) ---")

    # Prediction 5: thyA auxotrophy alone
    # Hayashi: 10^4-fold CFU reduction in 6 days without thymidine
    # This corresponds to escape frequency ~10^-4 (fraction surviving)
    # Our auxotrophy model: 1e-9 — but this is for complete metabolic bypass
    # thyA auxotrophy is less strict (thymidine can be scavenged from dead cells)
    # Better to use a relaxed version
    obs5_rate = 1e-4  # fraction surviving at day 6
    pred5_rate = MUTATION_SPECTRA["auxotrophy"].total
    predictions.append(
        {
            "source": "Hayashi 2024 - thyA auxotrophy alone (day 6 survival)",
            "observed": obs5_rate,
            "predicted": float(pred5_rate),
            "log10_obs": float(np.log10(obs5_rate)),
            "log10_pred": float(np.log10(pred5_rate)),
            "log10_error": float(np.log10(pred5_rate / obs5_rate)),
            "model_used": "auxotrophy",
            "reasoning": "thyA auxotrophy; our model underpredicts because it assumes strict dependence",
            "note": "Different organism (Bacteroides), different auxotrophy mechanism",
        }
    )
    print(
        f"  thyA auxotrophy: obs={obs5_rate:.0e}, pred={pred5_rate:.1e}, "
        f"log_err={np.log10(pred5_rate / obs5_rate):.2f}"
    )

    # Prediction 6: thyA auxotrophy + Cas9 device
    # Hayashi: 10^6-fold reduction in 6 days
    # This is a 2-layer system: auxotrophy + CRISPR
    obs6_rate = 1e-6
    crispr_layer = make_layer("crispr_multi")  # Cas9 targeting, multi-component
    aux_layer2 = make_layer("auxotrophy")
    ks_hayashi = MultiLayerKillSwitch(
        layers=[crispr_layer, aux_layer2], genomic_correlation=0.01
    )
    pred6_rate = ks_hayashi.combined_escape_rate()
    predictions.append(
        {
            "source": "Hayashi 2024 - thyA auxotrophy + Cas9 device (day 6 survival)",
            "observed": obs6_rate,
            "predicted": float(pred6_rate),
            "log10_obs": float(np.log10(obs6_rate)),
            "log10_pred": float(np.log10(pred6_rate)),
            "log10_error": float(np.log10(pred6_rate / obs6_rate)),
            "model_used": "auxotrophy + crispr_multi (rho=0.01)",
            "reasoning": "2-layer system; different organism but analogous mechanisms",
        }
    )
    print(
        f"  thyA+Cas9: obs={obs6_rate:.0e}, pred={pred6_rate:.1e}, "
        f"log_err={np.log10(pred6_rate / obs6_rate):.2f}"
    )

    # =====================================================================
    # Summary statistics
    # =====================================================================
    # Only use predictions where the comparison is fair (in vitro, same mechanism class)
    fair_predictions = [
        p
        for p in predictions
        if "in vivo" not in p["source"].lower() and "day 6 survival" not in p["source"]
    ]
    all_log_errors = [abs(p["log10_error"]) for p in predictions]
    fair_log_errors = [abs(p["log10_error"]) for p in fair_predictions]

    summary = {
        "predictions": predictions,
        "n_total": len(predictions),
        "n_fair": len(fair_predictions),
        "all_mean_abs_log_error": float(np.mean(all_log_errors)),
        "fair_mean_abs_log_error": float(np.mean(fair_log_errors))
        if fair_log_errors
        else None,
        "all_within_2_oom": all(e < 2.0 for e in all_log_errors),
        "fair_within_1_oom": all(e < 1.0 for e in fair_log_errors)
        if fair_log_errors
        else None,
    }

    print(
        f"\n  All predictions ({len(predictions)}): mean |log err| = {np.mean(all_log_errors):.2f}"
    )
    if fair_log_errors:
        print(
            f"  Fair comparisons ({len(fair_predictions)}): mean |log err| = {np.mean(fair_log_errors):.2f}"
        )
    print(f"  All within 2 OOM: {summary['all_within_2_oom']}")

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")
    _save("exp14_out_of_sample", summary)
    return summary


def run_iter3():
    t_start = time.time()
    print("=" * 70)
    print("Iteration 3: Out-of-sample validation")
    print("=" * 70)
    exp14_out_of_sample()
    total = time.time() - t_start
    print(f"\nIteration 3 complete. Total: {total:.1f}s")


if __name__ == "__main__":
    run_iter3()
