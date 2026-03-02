# Kill Switch Design: A Computational Framework for Evolutionary-Robust Biocontainment

Computational framework for designing and evaluating multi-layer kill switch systems in engineered bacteria. Models six kill switch architectures with mutation-spectrum-decomposed escape rates, simulates stochastic serial passage dynamics, and optimizes multi-layer combinations while accounting for mechanistic inter-layer correlations.

**Paper:** [kill_switch_design.pdf](kill_switch_design.pdf)

## Key Results

- Models 6 kill switch architectures (toxin-antitoxin, CRISPR single/multi-gRNA, overlapping gene entanglement, synthetic auxotrophy, integrase differentiation)
- Recommends triple-layer design: **CRISPR multi-gRNA + overlapping gene + synthetic auxotrophy**
  - Predicted escape rate: 6.1 × 10⁻¹⁸ (10 orders below NIH threshold of 10⁻⁸)
  - Total fitness cost: 20%
- Cross-validation mean absolute log-error: 0.04 (leave-one-study-out)
- Zero escape events across 5,000 replicates × 5,000 generations

## Repository Structure

```
├── src/
│   ├── models.py           # Kill switch escape rate models
│   ├── parameters.py       # Calibrated parameters from literature
│   ├── experiments.py      # Core simulation experiments
│   ├── experiments_iter2.py # Robustness and sensitivity analyses
│   ├── experiments_iter3.py # Out-of-sample validation
│   ├── figures.py          # Main figure generation
│   └── figures_supp.py     # Supplementary figure generation
├── figures/                # Generated figures (PDF + PNG)
├── results/                # Simulation outputs
├── paper.tex               # LaTeX manuscript (achemso/ACS format)
├── kill_switch_design.pdf  # Compiled paper
└── references.bib          # Bibliography
```

## Reproducing Results

Requires Python 3.12, NumPy 1.26, SciPy 1.12.

```bash
pip install numpy scipy matplotlib
python src/experiments.py      # Single-layer characterization, combinations, Pareto search
python src/experiments_iter2.py # Robustness checks
python src/experiments_iter3.py # Out-of-sample validation
python src/figures.py           # Generate main figures
python src/figures_supp.py      # Generate supplementary figures
```

All simulations use a fixed random seed (42) for reproducibility.

## Methods Summary

- **Escape rate model:** Mutation-spectrum decomposition into 5 mechanistic classes (point mutations, small indels, IS element insertions, large deletions, recombination)
- **Combination model:** Product formula with pairwise correlation correction (Eq. 1 in paper)
- **Simulation:** Stochastic serial passage, 10,000 replicates, N₀ = 10⁷, K = 10⁹, 1:100 dilution
- **Calibration:** Parameters fit to Rottinghaus et al. (2022), Chlebek et al. (2023), Williams & Murray (2022)

## Citation

If you use this framework, please cite the paper and the experimental studies it builds on (see `references.bib`).
