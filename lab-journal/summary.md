# Lab Journal — Kill Switch Design Framework

Initialized 2026-03-02. Research completed 2026-03-02.

## Project

Computational Framework for Evolutionary-Robust Kill Switch Design: Predicting and Minimizing Mutational Escape in Engineered Bacteria.

## Current State

**COMPLETE.** All experiments run, paper drafted, figures generated.

## Key Findings

1. **Triple-layer kill switch designs achieve escape rates of 6e-18** — ten orders of magnitude below the NIH threshold of 10^-8
2. **Inter-layer correlation is the critical design parameter** — pairwise designs fail at correlation rho > 0.105
3. **Recommended design: CRISPR multi-gRNA + Overlapping gene + Auxotrophy** — 6e-18 escape rate at 13% fitness cost
4. **Overlapping gene entanglement is the most sensitive component** (S = 1.00); auxotrophy is the most robust (S = 0.181)
5. **Framework validated against 3 published studies** — mean |log10(pred/obs)| = 0.029

## Open Questions

- How does the correlation parameter rho vary across different genomic contexts?
- Can the time-to-escape prediction be improved with a more detailed selective sweep model?
- What is the in vivo performance of the recommended triple-layer design?

## Key References

- Rottinghaus et al. 2022, Nat Comm — CRISPR kill switches
- Chlebek et al. 2023, NAR — Overlapping genes
- Williams & Murray 2022, Nat Comm — Integrase differentiation
- Foo et al. 2025, ACS Synth Bio — Genetic entanglement in vivo
