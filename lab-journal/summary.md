# Lab Journal — Kill Switch Design Framework

Initialized 2026-03-02. Last updated 2026-03-02 (after 3 iteration cycles).

## Project

Computational Framework for Evolutionary-Robust Kill Switch Design.

## Current State

**Paper v3 complete after 3 iteration cycles.** LaTeX draft + 10 figures (6 main + 4 supplementary) + 13 experiments.

## Iteration History

1. **Iter 0**: Initial prototype (had IS element double-counting bug, stochastic crash bug)
2. **Iter 1**: Fitness cost calibration from competitive index data. TA t10: 40 gen (obs 30, 1.3x). OG t10: 106 gen (obs >130, consistent).
3. **Iter 2**: Robustness checks — N0 sensitivity, multi-environment, 5000-gen deployment (0/5000 escape for triple), alternative correlation models.
4. **Iter 3**: Paper finalization — updated all results, added new sections, supplementary figures, reference additions.

## Key Findings

1. Triple design (CRISPR multi + OG + Auxotrophy): 6.06e-18 escape rate, 20% fitness cost
2. 0/5000 escape events through 5000 generations in long-term deployment simulation
3. Robust across: N0 (1e5-1e9), environments (rich-soil), correlation models (2/4 never lose NIH)
4. Time-to-escape validation: 1.3x error for TA (was 6.4x before fitness calibration)
5. Cross-validation: mean |log10 error| = 0.04
6. Sensitivity: OG most sensitive (S=1.00), Auxotrophy most robust (S=0.17)
7. Critical rho for pairwise CRISPR+OG: 0.102

## Experiments (13 total)

- 001-009: Core framework experiments
- 010: N0 sensitivity
- 011: Multi-environment
- 012: Long-term deployment (5000 gen)
- 013: Alternative correlation models

## Files

- paper.tex + references.bib: LaTeX manuscript (15 citations)
- src/: 6 Python modules
- figures/: 10 publication figures
- results/: 13 JSON result files
