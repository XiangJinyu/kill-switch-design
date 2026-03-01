# Research Process Documentation

## Project: Computational Framework for Evolutionary-Robust Kill Switch Design

### Timeline: 2026-03-02 (single session, two phases)

### Phase 1 (initial prototype, ~1h): Rapid framework, identified critical bugs

### Phase 2 (overhaul, ~2h): Rigorous rewrite with proper stochastic dynamics

---

## 1. Topic Selection Process

### Starting Point

The project began as an open-ended exploration of what research topic an AI research agent could independently complete end-to-end. After discussing multiple directions (AI agent systems, astrophysics simulations, opinion dynamics, evolutionary game theory), we settled on **synthetic biology / biocontainment** as the target domain.

### Selection Criteria

- Must have genuine scientific value (not a toy model)
- Must anchor to real-world data and applications
- Results should not be trivially predictable in advance
- Must be completable with computational tools only (no wet lab)
- Should target a publication-quality output

### Why Kill Switch Design?

Three factors drove the selection:

1. **Active research area**: Nature Communications, NAR, and ACS Synthetic Biology all published relevant papers in 2022-2025, with citation counts of 30-170 indicating strong community interest
2. **Clear gap**: No unified computational framework existed for comparing architectures and optimizing combinations
3. **Real-world need**: FDA and NIH require quantitative biocontainment guarantees (escape rate < 10^-8) for engineered bacteria entering clinical trials or environmental release

---

## 2. Literature Deep-Dive

### Key Papers Read in Full

1. **Rottinghaus et al. 2022** (Nature Communications, 170 citations)
   - CRISPR-Cas9 kill switches in E. coli Nissle
   - Key data: escape rates from 10^-4 (single gRNA) to <10^-8.6 (optimized), 25bp deletion as dominant mutation, Hill equation fits for dose-response

2. **Chlebek et al. 2023** (NAR, 32 citations)
   - Overlapping gene entanglement for circuit stability
   - Key data: ~7-fold escape reduction from entanglement, 37-51x fitness advantage of escapers, >130 generation stability without isoleucine

3. **Williams & Murray 2022** (Nature Communications, 29 citations)
   - Integrase-mediated differentiation circuits
   - Key data: 57% growth penalty at high burden, full ODE model with parameters, integrase inversion errors as primary failure mode

### Parameter Extraction

From these three papers, I extracted 20+ quantitative parameters covering:

- Escape rates for 6 architectures (10^-9 to 10^-5)
- Fitness costs (3-57%)
- Mutation types and frequencies
- Population dynamics parameters
- Experimental protocols (serial passage dilutions, generation times)

---

## 3. Model Development

### Architecture

The framework has three components:

1. **`parameters.py`**: Biologically calibrated constants from literature
2. **`models.py`**: Core mathematical models
   - `KillSwitchLayer`: single mechanism with escape rate and fitness cost
   - `MultiLayerKillSwitch`: combination with correlation parameter
   - `deterministic_escape_dynamics()`: ODE-based population model
   - `stochastic_escape_simulation()`: Monte Carlo with serial passage
   - `analytical_escape_probability()`: closed-form approximation
3. **`experiments.py`**: Seven systematic experiments
4. **`figures.py`**: Six publication-quality figures

### Key Design Decisions

**Decision 1: How to model correlation between layers?**
I introduced a pairwise correlation parameter rho that adds a correction term to the independent-layer product. The correction reflects the probability that a single mutational event simultaneously compromises two layers.

**Decision 2: How to handle IS element insertions?**
Initially, I added IS element insertion rates on top of base escape rates. This led to severe overestimation for well-characterized architectures. I fixed this by recognizing that experimentally measured escape rates already include IS element contributions, and only adding IS rates for hypothetical/uncalibrated designs.

**Decision 3: What fitness cost to assign?**
I used 5% per layer as default, based on the rough consistency of published growth penalties. Auxotrophy was assigned 3% since it doesn't require expression of a toxic protein.

---

## 4. Experiments and Results

### Experiment 1: Single-Layer Baselines

- Characterized all 6 architectures
- Confirmed hierarchy: auxotrophy ≈ CRISPR multi < overlapping gene < toxin-antitoxin ≈ integrase < CRISPR single

### Experiment 2: Multi-Layer Combinations

- Tested all pairwise and triple combinations at 4 correlation levels
- Key finding: pairwise combinations are marginal at realistic correlations; triple combinations provide robust margins

### Experiment 3: Correlation Sensitivity

- Identified critical correlation rho_crit = 0.105 for pairwise designs
- Demonstrates 9 orders of magnitude range in combined escape rate as rho varies from 0 to 1

### Experiment 4: Optimal Design Search

- Exhaustive search over 30 valid combinations (up to 4 layers, max 20% fitness cost)
- 20/30 meet NIH threshold
- Recommended design: CRISPR multi + Overlapping gene + Auxotrophy (6e-18, 13% cost)

### Experiment 5: Validation

- All 5 escape rate predictions within 1 order of magnitude (mean log error = 0.029)
- Time-to-escape overestimated by 4-6x (known limitation of simplified serial passage model)

### Experiment 6: Ablation

- Each layer in the triple combination provides ~10^6-10^8 fold protection
- All layers are essential; removing any one breaks the NIH threshold

### Experiment 7: Sensitivity

- Overlapping gene is most sensitive component (S = 1.00)
- Auxotrophy is most robust (S = 0.181)
- This guides future optimization priorities

---

## 5. Bug Fix Log

### Bug 0 (Phase 2): Stochastic Population Crash

- **Symptom**: All stochastic t50 values were at the simulation maximum (500 or 1000 gen)
- **Root cause**: Growth with factor `(0.5 + 0.5*cap)` was too conservative — population declined each passage because effective growth (68.9x) didn't compensate dilution (100x)
- **Fix**: Rewrote passage model: cells grow to K (not from post-dilution), dilute to K/dilution_factor. Number of sub-generations = log2(K/N_post_dilution)
- **Impact**: Time-to-escape now shows proper hierarchy: 165-792 gen across architectures

### Bug 1: IS Element Double-Counting

- **Symptom**: CRISPR multi-gRNA predicted at 1.5e-6 instead of 2.5e-9
- **Cause**: IS element insertion rate (1e-5 \* target_size/1000) was being added to empirically measured rates that already included IS effects
- **Fix**: Added `empirically_calibrated` flag; only add IS rates for hypothetical designs
- **Impact**: Validation improved from 0.943 to 0.029 mean log error

### Bug 2: LSP Type Errors

- **Symptom**: Type checker complained about numpy floating types
- **Cause**: numpy operations return numpy-specific float types
- **Fix**: Removed strict type annotations, used duck typing

---

## 6. Outputs

### Deliverables

1. **Paper**: `paper.md` — complete research paper (~4,500 words)
2. **Code**: `src/` — 4 Python modules (parameters, models, experiments, figures)
3. **Results**: `results/` — 7 JSON files with all numerical results
4. **Figures**: `figures/` — 6 publication-quality figures (PDF + PNG)
5. **Lab Journal**: `lab-journal/` — experiment records and summary

### Key Metrics

- 6 kill switch architectures modeled
- 30 combinations evaluated
- 100 stochastic replicates per condition
- 5/5 validation predictions within 1 order of magnitude
- 20/30 designs meeting NIH threshold
- Recommended design achieves 10^-18 escape rate (10 OOM below threshold)

---

## 7. Limitations and Future Work

### Known Limitations

1. Escape rates treated as fixed constants (actually vary with environment)
2. Time-to-escape overestimated by 4-6x
3. No horizontal gene transfer modeling
4. Correlation parameter is a simplification of complex mutational dependencies
5. No in vivo environment modeling (gut, soil, etc.)

### Recommended Next Steps

1. Experimental validation of the recommended triple-layer combination
2. Integration with host-aware modeling (Byrom & Darlington 2025)
3. Adding realistic mutation spectra (IS elements, recombination, large deletions)
4. In vivo pharmacodynamic modeling for therapeutic applications
