"""
Core mathematical models for kill switch escape dynamics.

Three levels of modeling:
1. Deterministic ODE: mean-field population dynamics
2. Stochastic (Gillespie-like): finite population effects
3. Analytical: closed-form expressions for escape probability
"""

import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field


@dataclass
class KillSwitchLayer:
    """A single kill switch mechanism."""

    name: str
    escape_rate: float  # per cell per generation (mu_e), empirically calibrated
    fitness_cost: float = 0.05  # fractional growth penalty when active
    target_size_bp: int = 300  # size of the vulnerability target
    is_element_susceptible: bool = True
    empirically_calibrated: bool = (
        True  # if True, escape_rate already includes IS effects
    )

    def effective_escape_rate(self, growth_rate: float = 2.0, is_rate: float = 1e-5):
        """
        Total escape rate. If empirically calibrated, use base rate directly
        (it already includes all mutation sources measured in experiments).
        Only add IS element contribution for hypothetical/uncalibrated designs.
        """
        if self.empirically_calibrated:
            return self.escape_rate
        base = self.escape_rate
        if self.is_element_susceptible:
            base += is_rate * (self.target_size_bp / 1000)
        return base


@dataclass
class MultiLayerKillSwitch:
    """Combined multi-layer kill switch system."""

    layers: list[KillSwitchLayer] = field(default_factory=list)
    correlation: float = 0.01  # pairwise correlation between layers

    @property
    def n_layers(self):
        return len(self.layers)

    def total_fitness_cost(self) -> float:
        """Cumulative fitness cost of all layers."""
        return sum(l.fitness_cost for l in self.layers)

    def combined_escape_rate(self, growth_rate: float = 2.0):
        """
        Combined escape rate accounting for correlation.

        For n independent layers with escape rates p_i:
            P_escape = prod(p_i)  (if independent)

        With pairwise correlation rho:
            P_escape = prod(p_i) * (1 + rho * sum_{i<j} 1/(p_i * p_j) * prod(p_k))

        Simplified: we use the formula
            P_combined ≈ prod(p_i) + rho * max(p_i) * prod_{j≠argmax} p_j

        This captures that correlation allows the highest-rate layer to
        "drag along" escape of other layers.
        """
        rates = [l.effective_escape_rate(growth_rate) for l in self.layers]
        if not rates:
            return 1.0

        independent = np.prod(rates)
        if len(rates) == 1:
            return rates[0]

        # Correlation correction: pairs of layers can co-escape
        correction = 0.0
        for i in range(len(rates)):
            for j in range(i + 1, len(rates)):
                # Probability that layers i,j co-escape due to correlation
                pair = max(rates[i], rates[j]) * self.correlation
                # Combined with independent escape of remaining layers
                remaining = [rates[k] for k in range(len(rates)) if k != i and k != j]
                correction += pair * np.prod(remaining) if remaining else pair

        return independent + correction


def deterministic_escape_dynamics(
    N0: float,
    generations: int,
    kill_switch: MultiLayerKillSwitch,
    mu_max: float = 2.0,
    K: float = 1e9,
    fitness_advantage: float = 40.0,
    dt_gen: float = 0.1,
) -> dict:
    """
    Deterministic ODE model of escape dynamics during serial passage.

    Tracks two populations:
    - N_w: wild-type (kill switch functional) cells
    - N_e: escape mutant cells

    dN_w/dt = mu_w * N_w * (1 - N_total/K) - mu_escape * N_w
    dN_e/dt = mu_e * N_e * (1 - N_total/K) + mu_escape * N_w

    where:
    - mu_w = mu_max * (1 - fitness_cost)
    - mu_e = mu_max (escapers have no burden)
    - mu_escape = combined escape rate per generation * growth rate
    """
    fitness_cost = kill_switch.total_fitness_cost()
    escape_rate = kill_switch.combined_escape_rate(mu_max)

    mu_w = mu_max * (1.0 - fitness_cost)
    mu_e = mu_max  # escapers regain full growth

    # Convert escape rate per generation to per hour
    gen_time = np.log(2) / mu_max
    escape_per_hour = escape_rate / gen_time

    steps = int(generations / dt_gen)
    t = np.zeros(steps)
    Nw = np.zeros(steps)
    Ne = np.zeros(steps)

    Nw[0] = N0 * (1 - escape_rate)
    Ne[0] = N0 * escape_rate

    for i in range(1, steps):
        gen = i * dt_gen
        t[i] = gen

        Ntot = Nw[i - 1] + Ne[i - 1]
        growth_factor = max(0, 1 - Ntot / K)

        dNw = mu_w * Nw[i - 1] * growth_factor - escape_per_hour * Nw[i - 1] * gen_time
        dNe = mu_e * Ne[i - 1] * growth_factor + escape_per_hour * Nw[i - 1] * gen_time

        Nw[i] = max(0, Nw[i - 1] + dNw * dt_gen * gen_time)
        Ne[i] = max(0, Ne[i - 1] + dNe * dt_gen * gen_time)

    escape_fraction = Ne / (Nw + Ne + 1e-30)

    return {
        "generations": t,
        "N_wildtype": Nw,
        "N_escape": Ne,
        "escape_fraction": escape_fraction,
        "parameters": {
            "N0": N0,
            "mu_w": mu_w,
            "mu_e": mu_e,
            "escape_rate": escape_rate,
            "fitness_cost": fitness_cost,
            "K": K,
        },
    }


def stochastic_escape_simulation(
    N0: int,
    generations: int,
    kill_switch: MultiLayerKillSwitch,
    mu_max: float = 2.0,
    K: float = 1e9,
    dilution_factor: int = 100,
    passage_interval_gen: float = 6.6,
    rng: np.random.Generator | None = None,
) -> dict:
    """
    Stochastic simulation of escape dynamics with serial passage.

    Each passage:
    1. Population grows to carrying capacity
    2. Mutations occur stochastically (binomial sampling)
    3. Population is diluted

    Returns time series of escape fraction.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    escape_rate = kill_switch.combined_escape_rate(mu_max)
    fitness_cost = kill_switch.total_fitness_cost()

    n_passages = int(generations / passage_interval_gen)
    gen_record = []
    Nw_record = []
    Ne_record = []

    Nw = int(N0 * (1 - escape_rate))
    Ne = int(N0 * escape_rate)

    for p in range(n_passages):
        gen = p * passage_interval_gen
        gen_record.append(gen)
        Nw_record.append(Nw)
        Ne_record.append(Ne)

        # Growth phase: both populations grow to fill capacity
        Ntot = Nw + Ne
        if Ntot == 0:
            break

        # Number of generations in this passage
        n_doublings = passage_interval_gen

        # Mutations during growth: new escapers from Nw
        # Expected new escapers = Nw * escape_rate * n_doublings
        # (each cell has escape_rate chance per generation, over n_doublings generations)
        expected_new_escapers = Nw * escape_rate * n_doublings
        if expected_new_escapers < 100:
            new_escapers = rng.poisson(expected_new_escapers)
        else:
            new_escapers = int(
                rng.normal(expected_new_escapers, np.sqrt(expected_new_escapers))
            )
            new_escapers = max(0, new_escapers)

        # Growth: both populations expand
        growth_w = (1 - fitness_cost) ** n_doublings
        growth_e = 1.0  # no fitness cost for escapers

        # After growth
        Nw_grown = Nw * (2**n_doublings) * (1 - fitness_cost) ** n_doublings
        Ne_grown = (Ne + new_escapers) * (2**n_doublings)

        # Cap at carrying capacity
        Ntot_grown = Nw_grown + Ne_grown
        if Ntot_grown > K:
            scale = K / Ntot_grown
            Nw_grown *= scale
            Ne_grown *= scale

        # Dilution: sample N/dilution_factor cells
        Ntot_after = Nw_grown + Ne_grown
        if Ntot_after == 0:
            break

        n_sample = int(Ntot_after / dilution_factor)
        if n_sample == 0:
            n_sample = 1

        frac_w = Nw_grown / Ntot_after
        # Binomial sampling for dilution
        Nw = rng.binomial(n_sample, frac_w)
        Ne = n_sample - Nw

    gen_record = np.array(gen_record)
    Nw_record = np.array(Nw_record, dtype=float)
    Ne_record = np.array(Ne_record, dtype=float)
    escape_fraction = Ne_record / (Nw_record + Ne_record + 1e-30)

    return {
        "generations": gen_record,
        "N_wildtype": Nw_record,
        "N_escape": Ne_record,
        "escape_fraction": escape_fraction,
        "n_passages": n_passages,
    }


def analytical_escape_probability(
    escape_rate,
    generations,
    N,
    fitness_advantage=40.0,
):
    """
    Analytical approximation for probability that at least one escaper
    establishes in the population.

    P(escape) ≈ 1 - exp(-N * escape_rate * generations * p_establish)

    where p_establish = 1 - 1/fitness_advantage (Haldane's approximation
    for fixation probability of a beneficial mutation).
    """
    p_establish = 1 - 1 / fitness_advantage if fitness_advantage > 1 else 0.5
    expected_established = N * escape_rate * generations * p_establish
    return 1 - np.exp(-expected_established)


def time_to_escape(
    escape_rate,
    N,
    fitness_advantage=40.0,
    threshold=0.5,
):
    """
    Expected number of generations until escape fraction reaches threshold.

    From mutation-selection balance, the escape fraction grows as:
    f(t) ≈ f0 * exp(s * t)

    where s = selective advantage and f0 = mu_e (mutation-selection equilibrium).

    Solving for t when f(t) = threshold:
    t = ln(threshold / f0) / s
    """
    s = np.log(fitness_advantage) / np.log(2)  # selective advantage per generation
    f0 = escape_rate  # initial escape fraction at mutation-selection balance
    if f0 <= 0 or s <= 0:
        return np.inf
    t = np.log(threshold / f0) / s
    return max(0, t)


def sweep_escape_rates(
    layers_config: list[dict],
    correlation_range: np.ndarray,
    N0: float = 1e6,
    generations: int = 500,
):
    """
    Sweep over correlation values to see how non-independence
    affects combined escape dynamics.
    """
    results = []
    for rho in correlation_range:
        layers = [KillSwitchLayer(**cfg) for cfg in layers_config]
        ks = MultiLayerKillSwitch(layers=layers, correlation=rho)
        rate = ks.combined_escape_rate()
        t_esc = time_to_escape(rate, N0)
        results.append(
            {
                "correlation": rho,
                "combined_escape_rate": rate,
                "time_to_50pct_escape": t_esc,
                "log10_escape_rate": np.log10(rate) if rate > 0 else -np.inf,
            }
        )
    return results
