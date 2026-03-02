---
id: "005"
slug: out-of-sample-validation
type: hypothesis
status: completed
created: 2026-03-02
concluded: 2026-03-02
depends_on: ["003", "004"]
conclusion_type: quantitative
conclusion: "Out-of-sample predictions within 1 OOM for E. coli systems (mean log err = 0.92). Cross-organism (Bacteroides) fails by 5 OOM — model requires organism-specific calibration."
tags: [validation, out-of-sample, Foo2025, Hayashi2024]
commit: ""
---

# 005: Out-of-Sample Validation

## Question

Can the framework predict escape rates for systems NOT used in calibration?

## Method

Tested 6 predictions against Foo 2025 (STALEMATE in E. coli Nissle) and Hayashi 2024 (thyA auxotrophy + Cas9 in B. thetaiotaomicron).

## Evidence

### Fair comparisons (same organism class, in vitro):

| System                         | Observed    | Predicted  | log10 error |
| ------------------------------ | ----------- | ---------- | ----------- | --- | -------- |
| Foo - Non-entangled pEndo      | 1e-5        | 9e-7       | -1.05       |
| Foo - Entangled pEnt-533 reop  | 1e-7        | 1.01e-7    | 0.00        |
| Foo - Entangled + ColE9 (0.47) | 1e-10       | 5.1e-9     | +1.70       |
| \*\*Mean                       | log10 error | (fair)\*\* |             |     | **0.92** |

### Cross-organism (Bacteroides, informative but unfair):

| System                    | Observed | Predicted | log10 error |
| ------------------------- | -------- | --------- | ----------- |
| Hayashi - thyA auxotrophy | 1e-4     | 1e-9      | -5.00       |
| Hayashi - thyA + Cas9     | 1e-6     | 2.5e-11   | -4.60       |

## Interpretation

1. For E. coli systems, out-of-sample predictions average 0.92 log10 error — within 1 OOM
2. The entangled gene prediction (pEnt-533 reop) is essentially perfect (0.00 error)
3. Non-entangled system (pEndo) underpredicted by 10x because IS element rates vary between strains
4. ColE9-augmented system overpredicted by 50x because our model lacks population-level policing
5. Bacteroides predictions fail completely (-5 OOM) — model requires organism-specific calibration
6. This honestly demonstrates both the power and limitations of the framework

## Next

Update paper with these results, add schematic figure, finalize.
