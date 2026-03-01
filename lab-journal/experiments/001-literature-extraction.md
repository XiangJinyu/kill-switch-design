---
id: "001"
slug: literature-extraction
type: exploration
status: completed
created: 2026-03-02
concluded: 2026-03-02
depends_on: []
conclusion_type: quantitative
conclusion: "Extracted key parameters from 3 papers: escape rates 10^-4 to 10^-8, mutation rates ~10^-6/gen, fitness advantages 37-51x for escapers"
tags: [literature, parameters, kill-switch]
commit: ""
---

# 001: Literature Deep-Dive and Parameter Extraction

## Question

What quantitative parameters from published kill switch experiments are needed to calibrate our computational model?

## Method

Deep reading of 3 key papers: Rottinghaus 2022 (CRISPR kill switch), Chlebek 2023 (overlapping genes), Williams & Murray 2022 (differentiation circuits).

## Evidence

### Key Parameters for Modeling

| Parameter                                    | Value                             | Source           |
| -------------------------------------------- | --------------------------------- | ---------------- |
| Baseline mutation rate (toxin inactivation)  | ~10^-6 per cell per gen           | Chlebek 2023     |
| Entangled gene escape rate                   | ~10^-7 per cell per gen           | Chlebek 2023     |
| Auxotrophy bypass rate                       | ~10^-9 per cell per gen           | Chlebek 2023     |
| CRISPR kill switch escape (single gRNA)      | 10^-4 to 10^-5                    | Rottinghaus 2022 |
| CRISPR kill switch escape (optimized 2-gRNA) | < 10^-8.6                         | Rottinghaus 2022 |
| CRISPR kill switch escape (gut-like)         | 10^-3.3 to 10^-4.2                | Rottinghaus 2022 |
| Fitness advantage of escapers                | 37-51x                            | Chlebek 2023     |
| Growth penalty from burden (50uM IPTG, 1x)   | ~57%                              | Williams 2022    |
| Generations per day (serial passage 1:100)   | ~6.6                              | Chlebek 2023     |
| Generations per day (serial passage 1:250)   | ~8                                | Rottinghaus 2022 |
| NIH recommended escape frequency             | < 10^-8                           | Rottinghaus 2022 |
| Dominant mutation type for promoter circuits | 25bp deletions in tandem repeats  | Rottinghaus 2022 |
| Dominant mutation type for overlapping genes | Large deletions (if no selection) | Chlebek 2023     |
| Time to 10% escape (unprotected, +ile)       | ~30 gen                           | Chlebek 2023     |
| Time to 10% escape (entangled, -ile)         | >130 gen (7/10 never)             | Chlebek 2023     |
| Carrying capacity (modeling)                 | 10^9 cells                        | Williams 2022    |
| Burden mutation rate (modeling)              | 10^-6 x mu                        | Williams 2022    |
| Plasmid loss rate (modeling)                 | 10^-4 x mu                        | Williams 2022    |

### Kill Switch Architecture Comparison

| Architecture                          | Escape Rate    | Mechanism                              | Key Weakness                   |
| ------------------------------------- | -------------- | -------------------------------------- | ------------------------------ |
| Toxin-antitoxin (single)              | ~10^-5         | Single point mutation                  | Single point of failure        |
| CRISPR (single gRNA)                  | 10^-4 to 10^-5 | Promoter deletion                      | Tandem repeat recombination    |
| CRISPR (multi-gRNA + SOS knockout)    | < 10^-8.6      | Multiple simultaneous mutations needed | Complex construction           |
| Overlapping genes                     | ~10^-7         | Reduces large deletion class           | Small deletions still possible |
| Synthetic auxotrophy                  | ~10^-9         | Requires new metabolic pathway         | Needs continuous supply        |
| Differentiation circuit (2x terminal) | Best long-term | Decouples burden from escape           | Integrase errors               |

## Interpretation

The literature reveals a clear hierarchy of biocontainment strategies with quantitatively characterized escape rates. Crucially, no published work has built a unified computational framework that:

1. Models all architectures under a common mathematical framework
2. Predicts combination escape rates accounting for non-independence
3. Optimizes multi-layer designs

This is the gap our paper fills.

## Next

Build the computational framework: ODE models for intracellular dynamics (002), stochastic population model (003), combination framework (004).
