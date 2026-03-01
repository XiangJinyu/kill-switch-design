---
id: "003"
slug: fitness-calibration
type: optimization
status: completed
created: 2026-03-02
concluded: 2026-03-02
depends_on: ["002"]
conclusion_type: quantitative
conclusion: "Calibrated fitness costs from competitive index data. TA t10: 40 gen (obs 30, ratio 1.3x). OG t10: 106 gen (obs >130, consistent)."
tags: [fitness, calibration, time-to-escape, validation]
commit: ""
---

# 003: Fitness Cost Calibration from Competitive Index Data

## Question

Can we fix the 6-8x time-to-escape overestimation by using experimentally measured competitive fitness advantages instead of a generic 5% cost?

## Method

- Chlebek 2023 measured 37-51x competitive advantage for TA escapers over 13 generations
- Per-generation advantage: 40^(1/13) = 1.328, implying cost ≈ 0.25 for TA
- Different architectures have different burdens: CRISPR (low, induced system), OG (moderate), Integrase (moderate), Auxotrophy (minimal)
- Assigned architecture-specific costs: TA=0.25, CRISPR=0.05, OG=0.12, Aux=0.03, Integrase=0.15

## Evidence

| Architecture     | Cost | Predicted t10 (gen) | Observed t10 (gen) | Ratio      |
| ---------------- | ---- | ------------------- | ------------------ | ---------- |
| Toxin-antitoxin  | 25%  | 40                  | 30                 | 1.3x       |
| Overlapping gene | 12%  | 106                 | >130 (7/10 never)  | consistent |
| CRISPR multi     | 5%   | 350                 | stable >224        | consistent |

Previous model with uniform 5% cost: TA t10 = 191 gen (6.4x overestimate)

## Interpretation

The dominant source of time-to-escape overestimation was the fitness cost parameter, not the mutation rate or the simulation mechanism. Architecture-specific calibration from competition assay data brings predictions within 1.3x of experimental values for TA, and within the observed range for OG and CRISPR.

The key insight: the "fitness cost" in serial passage is NOT just the metabolic cost of maintaining the circuit. It includes leaky expression of toxic proteins, metabolic diversion to antitoxin production, and any growth-rate effects of carrying extra genetic material. These combined effects are captured by the competitive index measurement.

## Next

Iteration 2: systematic robustness check (004), then literature comparison (005).
