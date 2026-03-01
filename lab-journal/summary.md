# Lab Journal — Kill Switch Design Framework

Initialized 2026-03-02. Last updated 2026-03-02.

## Project

Computational Framework for Evolutionary-Robust Kill Switch Design.

## Current State

**Paper v1 complete.** LaTeX draft + 7 publication figures + 9 experiments with converged stochastic simulations.

## Key Findings

1. Triple-layer design (CRISPR multi + OG + Auxotrophy) achieves 6.06e-18 escape rate at 13% fitness cost
2. Inter-layer correlation rho_crit = 0.102 for pairwise CRISPR+OG combination
3. Triple design maintains NIH compliance across all tested correlations (rho=0 to 1)
4. Mutation spectrum: recombination dominates CRISPR single and OG; point mutations dominate optimized CRISPR multi
5. Cross-validation: mean |log10 error| = 0.04 across Rottinghaus and Chlebek datasets
6. Convergence verified: SE(t50) = 0.02 at 10,000 replicates

## Key Bug Fixes

- IS element double-counting (Phase 1): was inflating rates by ~1000x for CRISPR multi
- Stochastic population crash (Phase 2): growth model wasn't recovering population to K after dilution

## Experiments

- 001: Literature extraction (completed)
- 002: Model overhaul (completed)

## Open Questions

- Time-to-escape still overpredicts by 6-8x vs Chlebek data
- Correlation parameter needs experimental measurement
- Need in vivo validation of triple-layer design
