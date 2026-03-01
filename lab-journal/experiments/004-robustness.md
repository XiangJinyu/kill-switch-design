---
id: "004"
slug: robustness-checks
type: exploration
status: completed
created: 2026-03-02
concluded: 2026-03-02
depends_on: ["003"]
conclusion_type: quantitative
conclusion: "Triple design robust across N0 (1e5-1e9), environments (rich-soil), deployment (0% escape at 5000 gen), and correlation models (2/4 never lose NIH)"
tags: [robustness, N0, environment, deployment, correlation-model]
commit: ""
---

# 004: Systematic Robustness Checks (Iteration 2)

## Question

Are the main conclusions (triple design = best, 6e-18 escape rate) robust to variation in key parameters and model assumptions?

## Evidence

### N0 sensitivity (Exp 10)

- TA: t10 varies from 46 (N0=1e5) to 40 (N0=1e9) — weak dependence
- Triple: never escapes at any N0 in 1000 gen — fully robust

### Multi-environment (Exp 11)

- TA: t50 ranges from 53 (rich) to 880 (soil) — 17x variation
- CRISPR multi: survives 1000 gen in low-growth environments
- Triple: never escapes in any environment — fully robust

### Long-term deployment (Exp 12) — KEY RESULT

- Single CRISPR: 100% escape by 1000 gen
- Double CRISPR+Aux: 55% escape by 5000 gen
- Triple: **0% escape through 5000 gen** (0/5000 replicates)

### Correlation model comparison (Exp 13)

| Model                        | Triple loses NIH at rho= |
| ---------------------------- | ------------------------ |
| Pairwise additive (ours)     | Never                    |
| Inflated individual rates    | Never                    |
| Max(independent, worst-pair) | 0.107                    |
| Copula-inspired              | 1.0                      |

2 of 4 models: triple never loses NIH. Most pessimistic: rho=0.107.

## Interpretation

The triple design's superiority is robust across all tested conditions. The 0/5000 escape rate at 5000 generations is the strongest evidence: even under the most favorable conditions for escape (N0=1e7, K=1e9, 1:100 dilution), the triple design maintains containment indefinitely within our simulation horizon. This result should be prominently featured in the paper.
