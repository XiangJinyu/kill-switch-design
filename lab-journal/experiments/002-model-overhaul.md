---
id: "002"
slug: model-overhaul
type: optimization
status: completed
created: 2026-03-02
concluded: 2026-03-02
depends_on: ["001"]
conclusion_type: quantitative
conclusion: "Fixed critical stochastic model bug (population crash), added mutation spectrum decomposition, validated with 10k replicates (SE=0.02)"
tags: [model, stochastic, mutation-spectrum, convergence]
commit: ""
---

# 002: Model Overhaul — Mutation Spectrum + Fixed Stochastic Dynamics

## Question

Can we improve the model to (a) decompose escape by mutation type, (b) fix the stochastic simulation population crash, and (c) properly validate convergence?

## Method

1. Replaced lumped escape rates with 5-component mutation spectra per architecture
2. Rewrote stochastic passage: cells grow to K, then dilute (not the reverse)
3. Ran convergence analysis from 50 to 10,000 replicates
4. Cross-validated against Rottinghaus and Chlebek data separately

## Evidence

- Population now maintains stable K/dilution_factor after each passage
- t50 hierarchy: crispr_single (165) < TA (238) < OG (284) < crispr_multi (403) < auxotrophy (792)
- Convergence: SE(t50) = 0.02 at 10k replicates
- Cross-validation: mean |log10 error| = 0.04, all within 1 OOM

## Interpretation

The critical bug was that the old model's growth dynamics couldn't sustain population through serial passage — cells were declining each round instead of recovering to K. The new model correctly simulates exponential growth back to carrying capacity, producing biologically realistic escape dynamics.

## Next

Generate LaTeX paper (003), improve analytical formula accuracy.
