"""
Biologically calibrated parameters for kill switch modeling.
All values sourced from published experimental data.
"""

import numpy as np

# =============================================================================
# Growth parameters
# =============================================================================
MU_MAX = 2.0  # max growth rate (h^-1), E. coli in rich media, Williams 2022
MU_POOR = 0.6  # growth rate in minimal media (h^-1), estimated
K = 1e9  # carrying capacity (cells), Williams 2022
DOUBLING_RICH = np.log(2) / MU_MAX  # ~0.35 h
DOUBLING_POOR = np.log(2) / MU_POOR  # ~1.16 h

# =============================================================================
# Mutation rates (per cell per generation)
# =============================================================================
# Base mutation rates for different escape mechanisms
MU_POINT = 1e-9  # per-bp point mutation rate in E. coli
GENOME_BP = 4.6e6  # E. coli genome size (bp)

# Kill switch-specific escape rates (per cell per generation)
ESCAPE_RATES = {
    "toxin_antitoxin": {
        "description": "Type II toxin-antitoxin (e.g., CcdB/CcdA, RelE/RelB)",
        "rate": 1e-6,  # Chlebek 2023, baseline
        "mechanism": "point_mutation_or_small_deletion",
        "target_size_bp": 300,  # typical toxin gene ~300bp
    },
    "crispr_single": {
        "description": "Single gRNA CRISPR-Cas9 kill switch",
        "rate": 3e-5,  # Rottinghaus 2022, single gRNA
        "mechanism": "promoter_deletion_tandem_repeat",
        "target_size_bp": 50,  # Ptet promoter region
    },
    "crispr_multi": {
        "description": "Multi-gRNA CRISPR + SOS knockout",
        "rate": 2.5e-9,  # Rottinghaus 2022, optimized < 10^-8.6
        "mechanism": "multiple_independent_mutations",
        "target_size_bp": 150,  # multiple promoter regions
    },
    "overlapping_gene": {
        "description": "Gene entanglement (overlapping reading frames)",
        "rate": 1e-7,  # Chlebek 2023, with ilvA selection
        "mechanism": "small_deletion_or_regulator_mutation",
        "target_size_bp": 288,  # relE gene in +1 frame
    },
    "auxotrophy": {
        "description": "Synthetic auxotrophy (essential gene deletion)",
        "rate": 1e-9,  # Chlebek 2023, auxotrophy bypass
        "mechanism": "metabolic_pathway_acquisition",
        "target_size_bp": 0,  # requires gain-of-function
    },
    "integrase_differentiation": {
        "description": "Integrase-mediated terminal differentiation",
        "rate": 1e-6,  # Williams 2022, integrase error rate
        "mechanism": "integrase_inversion_error",
        "target_size_bp": 100,  # attB/attP sites
    },
}

# =============================================================================
# Fitness parameters
# =============================================================================
FITNESS_COST_KILL_SWITCH = 0.05  # 5% growth penalty per kill switch layer
FITNESS_ADVANTAGE_ESCAPER = 40.0  # Chlebek 2023: 37-51x competitive advantage
BURDEN_PENALTY = {
    "low": 0.11,  # 11% growth penalty, Williams 2022 (10 uM IPTG)
    "mid": 0.31,  # 31% growth penalty, Williams 2022 (20 uM IPTG)
    "high": 0.57,  # 57% growth penalty, Williams 2022 (50 uM IPTG)
}

# =============================================================================
# Experimental conditions
# =============================================================================
SERIAL_PASSAGE = {
    "dilution_100x": {
        "dilution_factor": 100,
        "interval_h": 24,
        "generations_per_passage": 6.6,  # Chlebek 2023
    },
    "dilution_250x": {
        "dilution_factor": 250,
        "interval_h": 24,
        "generations_per_passage": 8.0,  # Rottinghaus 2022
    },
    "dilution_50x": {
        "dilution_factor": 50,
        "interval_h": 8,
        "generations_per_passage": 5.6,  # Williams 2022
    },
}

# =============================================================================
# IS element parameters (major source of instability)
# =============================================================================
IS_ELEMENT_INSERTION_RATE = 1e-5  # per gene per generation (literature consensus)
IS_ELEMENTS_IN_ECOLI = 40  # approximate number in E. coli K-12

# =============================================================================
# Non-independence factors for combination strategies
# =============================================================================
# When combining two kill switches, escape is not fully independent because:
# 1. Shared regulatory elements (e.g., same promoter system)
# 2. IS elements can disrupt multiple cassettes
# 3. Large deletions can span multiple components
# 4. Metabolic burden from multiple switches reduces growth → stronger selection for escapers
CORRELATION_FACTORS = {
    "shared_promoter": 0.3,  # 30% correlation if same promoter system
    "same_plasmid": 0.5,  # 50% correlation if on same plasmid
    "shared_IS_target": 0.1,  # 10% correlation from IS element insertion
    "independent_genomic": 0.01,  # 1% residual correlation (large deletions, etc.)
}

# =============================================================================
# Simulation defaults
# =============================================================================
DEFAULT_SEED = 42
DEFAULT_N0 = 1e7  # initial population size (post-dilution from K/100)
DEFAULT_GENERATIONS = 1000
DEFAULT_REPLICATES = 100
