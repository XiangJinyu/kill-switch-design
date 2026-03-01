"""
Core mathematical models for kill switch escape dynamics.

Four levels of modeling:
1. Mutation-spectrum-resolved: decomposed per-mutation-type escape rates
2. Deterministic ODE with proper selective sweeps (scipy solve_ivp)
3. Stochastic serial passage (tau-leaping with exact passage bottleneck)
4. Analytical closed-form approximations for rapid parameter sweeps
"""

import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field


# =============================================================================
# Mutation spectrum decomposition
# =============================================================================
@dataclass
class MutationSpectrum:
    """
    Decomposed mutation rates by mechanism.
    All rates are per cell per generation.
    """

    point_mutation: float = 0.0  # SNPs that inactivate the kill switch component
    small_indel: float = 0.0  # <50bp insertions/deletions (frameshifts, etc.)
    is_element: float = 0.0  # IS element transposition into the target
    large_deletion: float = 0.0  # >50bp deletions removing regulatory/coding regions
    recombination: float = 0.0  # homologous recombination between tandem repeats

    @property
    def total(self):
        return (
            self.point_mutation
            + self.small_indel
            + self.is_element
            + self.large_deletion
            + self.recombination
        )

    def to_dict(self):
        return {
            "point_mutation": self.point_mutation,
            "small_indel": self.small_indel,
            "is_element": self.is_element,
            "large_deletion": self.large_deletion,
            "recombination": self.recombination,
            "total": self.total,
        }


# Literature-calibrated mutation spectra for each architecture
MUTATION_SPECTRA = {
    "toxin_antitoxin": MutationSpectrum(
        # Chlebek 2023: RelE toxin, total ~1e-6
        # Dominant: promoter/regulator mutations + small deletions in toxin
        point_mutation=2e-7,  # missense/nonsense in toxin active site
        small_indel=3e-7,  # frameshift in toxin gene
        is_element=1e-7,  # IS insertion into toxin gene (~300bp target)
        large_deletion=2e-7,  # deletion of promoter/regulator region
        recombination=1e-7,  # tandem repeat slippage if present
    ),
    "crispr_single": MutationSpectrum(
        # Rottinghaus 2022: single gRNA, total ~3e-5
        # Dominant: 25bp deletions in Ptet promoter (tandem repeat recombination)
        point_mutation=1e-6,  # promoter/cas9 coding mutations
        small_indel=2e-6,  # small promoter deletions
        is_element=1e-6,  # IS into cas9 or gRNA cassette
        large_deletion=5e-6,  # regulatory region deletions
        recombination=2e-5,  # tandem-repeat-mediated 25bp deletion (DOMINANT)
    ),
    "crispr_multi": MutationSpectrum(
        # Rottinghaus 2022: multi-gRNA + SOS knockout, total ~2.5e-9
        # SOS knockout eliminates most error-prone repair
        # Multiple gRNAs require multiple simultaneous hits
        point_mutation=5e-10,
        small_indel=5e-10,
        is_element=5e-10,
        large_deletion=5e-10,
        recombination=5e-10,  # greatly reduced by eliminating tandem repeats
    ),
    "overlapping_gene": MutationSpectrum(
        # Chlebek 2023: ilvA/relE entanglement, total ~1e-7
        # Large deletions eliminated by ilvA essentiality
        # Remaining: small deletions in relE, regulator mutations
        point_mutation=1e-8,  # reduced by codon overlap
        small_indel=3e-8,  # small indels in relE still possible
        is_element=1e-8,  # IS into overlap region
        large_deletion=1e-9,  # strongly suppressed by essential gene coupling
        recombination=5e-8,  # regulator (rhaR/rhaS) mutations
    ),
    "auxotrophy": MutationSpectrum(
        # Chlebek 2023: delta-ilvA delta-TdcB, total ~1e-9
        # Requires gain-of-function: new metabolic pathway or HGT
        point_mutation=1e-10,  # reversion extremely unlikely with double knockout
        small_indel=0,  # cannot gain function through deletion
        is_element=1e-10,  # IS-mediated promoter capture of cryptic gene
        large_deletion=0,  # cannot gain function through deletion
        recombination=8e-10,  # HGT acquisition from environmental DNA
    ),
    "integrase_differentiation": MutationSpectrum(
        # Williams & Murray 2022: Bxb1 integrase circuit, total ~1e-6
        # Dominant: integrase inversion errors (recombine in wrong orientation)
        point_mutation=5e-8,
        small_indel=5e-8,
        is_element=1e-7,
        large_deletion=1e-7,
        recombination=7e-7,  # integrase inversion error (DOMINANT)
    ),
}


# =============================================================================
# Kill switch layer and multi-layer system
# =============================================================================
@dataclass
class KillSwitchLayer:
    """A single kill switch mechanism with decomposed mutation spectrum."""

    name: str
    spectrum: MutationSpectrum
    fitness_cost: float = 0.05

    @property
    def escape_rate(self):
        return self.spectrum.total

    def shared_mutation_classes(self, other):
        """
        Identify mutation classes that could simultaneously affect both layers.
        Returns the correlation contribution from shared vulnerabilities.
        """
        shared = 0.0
        # Large deletions can span adjacent genomic loci
        if self.spectrum.large_deletion > 0 and other.spectrum.large_deletion > 0:
            shared += (
                min(self.spectrum.large_deletion, other.spectrum.large_deletion) * 0.1
            )
        # IS elements: same IS family can hit multiple targets
        if self.spectrum.is_element > 0 and other.spectrum.is_element > 0:
            shared += min(self.spectrum.is_element, other.spectrum.is_element) * 0.05
        return shared


@dataclass
class MultiLayerKillSwitch:
    """Combined multi-layer kill switch system with mechanistic correlation."""

    layers: list[KillSwitchLayer] = field(default_factory=list)
    genomic_correlation: float = 0.01  # base correlation from genomic context
    shared_plasmid: bool = False  # if True, adds correlation from co-location

    @property
    def n_layers(self):
        return len(self.layers)

    def total_fitness_cost(self):
        return sum(l.fitness_cost for l in self.layers)

    def effective_correlation(self, i, j):
        """Compute effective correlation between layers i and j."""
        rho = self.genomic_correlation
        if self.shared_plasmid:
            rho += 0.15  # plasmid loss eliminates both
        # Add mechanistic correlation from shared mutation classes
        rho += self.layers[i].shared_mutation_classes(self.layers[j])
        return min(rho, 1.0)

    def combined_escape_rate(self):
        """
        Combined escape rate with pairwise mechanistic correlations.

        P_combined = prod(p_i) + sum_{i<j} rho_ij * max(p_i, p_j) * prod_{k!=i,j} p_k

        The second term accounts for correlated co-escape: a single mutational
        event (e.g., large deletion, IS element) that simultaneously
        inactivates layers i and j, with remaining layers escaping independently.
        """
        rates = [l.escape_rate for l in self.layers]
        if not rates:
            return 1.0
        if len(rates) == 1:
            return rates[0]

        independent = np.prod(rates)

        correction = 0.0
        for i in range(len(rates)):
            for j in range(i + 1, len(rates)):
                rho = self.effective_correlation(i, j)
                pair_rate = max(rates[i], rates[j]) * rho
                remaining = [rates[k] for k in range(len(rates)) if k != i and k != j]
                correction += pair_rate * np.prod(remaining) if remaining else pair_rate

        return float(independent + correction)


# =============================================================================
# Deterministic ODE model (scipy solve_ivp for proper integration)
# =============================================================================
def deterministic_escape_ode(
    N0,
    generations,
    kill_switch,
    mu_max=2.0,
    K=1e9,
    dt_gen=0.01,
):
    """
    Deterministic ODE model using scipy's solve_ivp (RK45 adaptive stepping).

    State: [N_w, N_e] — wild-type and escape mutant populations.

    Uses generation-time units (not hours).
    """
    fitness_cost = kill_switch.total_fitness_cost()
    escape_rate = kill_switch.combined_escape_rate()

    # Growth rates in units of doublings per generation
    s = fitness_cost  # selective disadvantage of wild-type
    growth_w = np.log(2) * (1 - s)  # wild-type growth (per generation time)
    growth_e = np.log(2)  # escape mutant growth (per generation time)

    def rhs(t, y):
        Nw, Ne = y
        Ntot = Nw + Ne
        cap = max(0.0, 1.0 - Ntot / K)
        dNw = growth_w * Nw * cap - escape_rate * Nw
        dNe = growth_e * Ne * cap + escape_rate * Nw
        return [dNw, dNe]

    y0 = [N0 * (1 - escape_rate), max(1.0, N0 * escape_rate)]
    t_span = (0, generations)
    t_eval = np.linspace(0, generations, int(generations / dt_gen) + 1)

    sol = solve_ivp(
        rhs,
        t_span,
        y0,
        method="RK45",
        t_eval=t_eval,
        rtol=1e-10,
        atol=1e-12,
        max_step=1.0,
    )

    Nw = np.maximum(sol.y[0], 0)
    Ne = np.maximum(sol.y[1], 0)
    frac = Ne / (Nw + Ne + 1e-30)

    return {
        "generations": sol.t,
        "N_wildtype": Nw,
        "N_escape": Ne,
        "escape_fraction": frac,
        "parameters": {
            "N0": N0,
            "escape_rate": escape_rate,
            "fitness_cost": fitness_cost,
            "K": K,
            "growth_w": growth_w,
            "growth_e": growth_e,
        },
    }


# =============================================================================
# Stochastic serial passage with selective sweep modeling
# =============================================================================
def stochastic_passage_simulation(
    N0,
    generations,
    kill_switch,
    K=1e9,
    dilution_factor=100,
    passage_gens=6.6,
    fitness_advantage=40.0,
    rng=None,
):
    """
    Stochastic serial passage simulation with proper selective sweep dynamics.

    Each passage:
    1. Mutations: Poisson-distributed new escapers per generation during growth
    2. Growth with competition: wild-type and escapers grow at different rates
       within carrying capacity. Growth is modeled across sub-generations.
    3. Bottleneck: binomial sampling at dilution.

    The key improvement over the previous model: within each passage, we
    simulate multiple sub-generations to properly capture the selective
    advantage of escapers.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    escape_rate = kill_switch.combined_escape_rate()
    fitness_cost = kill_switch.total_fitness_cost()

    n_passages = int(generations / passage_gens)
    gens_per_passage = passage_gens

    gen_record = np.zeros(n_passages + 1)
    Nw_record = np.zeros(n_passages + 1)
    Ne_record = np.zeros(n_passages + 1)

    Nw = float(N0)
    Ne = 0.0

    gen_record[0] = 0
    Nw_record[0] = Nw
    Ne_record[0] = Ne

    for p in range(n_passages):
        if Nw + Ne <= 0:
            gen_record[p + 1] = (p + 1) * passage_gens
            Nw_record[p + 1] = 0
            Ne_record[p + 1] = 0
            continue

        # ---- Growth phase ----
        # Standard serial passage: cells grow from post-dilution density
        # back to carrying capacity K. The number of generations within
        # a passage is determined by ln(K / N_after_dilution) / ln(2).
        #
        # During growth, mutations accumulate proportionally to the total
        # cell divisions. We model this accurately per-generation.

        Ntot_start = Nw + Ne
        # Actual generations of growth: how many doublings to reach K
        actual_gens = np.log2(K / max(Ntot_start, 1))
        actual_gens = max(1, min(actual_gens, 30))  # cap at 30 doublings

        n_subgens = int(actual_gens)
        for _ in range(n_subgens):
            # 1. Mutations from wild-type pool
            new_mut = 0
            if Nw > 0:
                lam = Nw * escape_rate
                if 0 < lam < 1e6:
                    new_mut = rng.poisson(lam)
                elif lam >= 1e6:
                    new_mut = max(0, int(rng.normal(lam, np.sqrt(lam))))

            Ne += new_mut
            Nw -= new_mut
            if Nw < 0:
                Nw = 0

            # 2. Growth: one doubling
            #    Wild-type fitness = (1 - fitness_cost) relative to escape
            Nw *= 2.0 * (1.0 - fitness_cost)
            Ne *= 2.0

            # 3. Cap at K
            Ntot = Nw + Ne
            if Ntot > K:
                Nw *= K / Ntot
                Ne *= K / Ntot

        # ---- Bottleneck: dilution to K / dilution_factor ----
        Ntot = Nw + Ne
        if Ntot <= 0:
            Nw, Ne = 0.0, 0.0
        else:
            n_sample = max(1, int(K / dilution_factor))
            frac_w = max(0.0, min(1.0, Nw / Ntot))
            nw_sampled = rng.binomial(n_sample, frac_w)
            Nw = float(nw_sampled)
            Ne = float(n_sample - nw_sampled)

        gen_record[p + 1] = (p + 1) * passage_gens
        Nw_record[p + 1] = Nw
        Ne_record[p + 1] = Ne

    frac = Ne_record / (Nw_record + Ne_record + 1e-30)
    return {
        "generations": gen_record,
        "N_wildtype": Nw_record,
        "N_escape": Ne_record,
        "escape_fraction": frac,
        "n_passages": n_passages,
    }


def run_replicate_simulations(
    n_replicates,
    N0,
    generations,
    kill_switch,
    K=1e9,
    dilution_factor=100,
    passage_gens=6.6,
    fitness_advantage=40.0,
    seed=42,
):
    """
    Run n_replicates stochastic simulations and aggregate statistics.
    Returns per-generation statistics: mean, median, std, 5th/95th percentiles.
    """
    rng = np.random.default_rng(seed)
    all_fracs = []
    escape_times_50 = []
    escape_times_10 = []
    result = None

    for rep in range(n_replicates):
        result = stochastic_passage_simulation(
            N0=N0,
            generations=generations,
            kill_switch=kill_switch,
            K=K,
            dilution_factor=dilution_factor,
            passage_gens=passage_gens,
            fitness_advantage=fitness_advantage,
            rng=rng,
        )
        frac = result["escape_fraction"]
        all_fracs.append(frac)

        # Time to 50% escape
        above = np.where(frac > 0.5)[0]
        t50 = (
            float(result["generations"][above[0]])
            if len(above) > 0
            else float(generations)
        )
        escape_times_50.append(t50)

        # Time to 10% escape
        above10 = np.where(frac > 0.1)[0]
        t10 = (
            float(result["generations"][above10[0]])
            if len(above10) > 0
            else float(generations)
        )
        escape_times_10.append(t10)

    all_fracs = np.array(all_fracs)  # shape: (n_replicates, n_timepoints)
    assert result is not None, "No replicates were run"
    gens = result["generations"]

    return {
        "generations": gens,
        "mean_fraction": np.mean(all_fracs, axis=0),
        "median_fraction": np.median(all_fracs, axis=0),
        "std_fraction": np.std(all_fracs, axis=0),
        "pct5_fraction": np.percentile(all_fracs, 5, axis=0),
        "pct95_fraction": np.percentile(all_fracs, 95, axis=0),
        "pct25_fraction": np.percentile(all_fracs, 25, axis=0),
        "pct75_fraction": np.percentile(all_fracs, 75, axis=0),
        "escape_times_50": np.array(escape_times_50),
        "escape_times_10": np.array(escape_times_10),
        "final_fractions": all_fracs[:, -1],
        "n_replicates": n_replicates,
    }


# =============================================================================
# Analytical approximations
# =============================================================================
def analytical_escape_probability(escape_rate, generations, N, fitness_advantage=40.0):
    """
    Probability that at least one escape mutant establishes (Haldane 1927).

    P_establish per mutant = (s) / (1 + s/2) for branching process (Haldane)
    where s = (w_e - w_w) / w_w is the selective advantage.

    Total: P(>=1 established) = 1 - exp(-N * mu * T * p_establish)
    """
    s = (fitness_advantage - 1) / 1.0  # relative fitness advantage
    p_fix = s / (1 + s / 2)  # Haldane's approximation
    expected = N * escape_rate * generations * p_fix
    return float(1 - np.exp(-min(expected, 500)))  # cap to avoid overflow


def time_to_escape_analytical(escape_rate, N, fitness_cost=0.05, threshold=0.5):
    """
    Analytical estimate: generations until escape fraction reaches threshold.

    Two phases:
    1. Waiting: expected generations until first escaper appears.
       Rate of appearance = N * escape_rate per generation.
       Each mutant survives drift with probability ~ s (selective advantage).
       For serial passage, s = fitness_cost / (1 - fitness_cost).
       Expected wait = 1 / (N * escape_rate * s)

    2. Selective sweep: single escaper grows exponentially relative to
       wild-type. Relative growth advantage per generation = 1/(1 - fitness_cost).
       Time for escape fraction to go from 1/N to threshold:
       t_sweep = ln(threshold * N) / ln(1/(1 - fitness_cost))

    Total: t_wait + t_sweep
    """
    # Per-generation selective advantage of escaper
    s = fitness_cost / (1 - fitness_cost)  # e.g., 0.05/0.95 = 0.0526

    # Probability a new mutant survives stochastic loss (Haldane)
    # For small s: p_survive ≈ 2*s (branching process result)
    p_survive = min(2 * s, 1.0)

    # Rate of established mutations per generation
    rate_established = N * escape_rate * p_survive
    if rate_established <= 0:
        return np.inf

    # Phase 1: waiting time
    t_wait = 1.0 / rate_established

    # Phase 2: sweep time (logistic growth of escape fraction)
    # Relative growth rate per generation: r = ln(1/(1-fitness_cost))
    r = -np.log(1 - fitness_cost)  # e.g., ln(1/0.95) = 0.0513
    if r > 0 and N > 0:
        t_sweep = np.log(threshold * N) / r
    else:
        t_sweep = np.inf

    return float(max(0, t_wait + t_sweep))


# =============================================================================
# Convergence analysis
# =============================================================================
def convergence_analysis(
    N0,
    generations,
    kill_switch,
    replicate_counts=(100, 500, 1000, 2000, 5000, 10000),
    seed=42,
):
    """
    Check convergence of stochastic estimates with increasing replicates.
    Returns statistics at each replicate count for convergence plotting.
    """
    results = []
    for n in replicate_counts:
        stats = run_replicate_simulations(
            n_replicates=n,
            N0=N0,
            generations=generations,
            kill_switch=kill_switch,
            seed=seed,
        )
        results.append(
            {
                "n_replicates": n,
                "mean_t50": float(np.mean(stats["escape_times_50"])),
                "std_t50": float(np.std(stats["escape_times_50"])),
                "se_t50": float(np.std(stats["escape_times_50"]) / np.sqrt(n)),
                "mean_final_frac": float(np.mean(stats["final_fractions"])),
                "std_final_frac": float(np.std(stats["final_fractions"])),
                "se_final_frac": float(np.std(stats["final_fractions"]) / np.sqrt(n)),
                "median_t50": float(np.median(stats["escape_times_50"])),
            }
        )
    return results
