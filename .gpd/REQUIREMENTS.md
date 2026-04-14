# Requirements: Classification of 4D N=1 Quiver Gauge Theories

**Defined:** 2026-04-14
**Core Research Question:** What is the complete classification of 4D N=1 asymptotically free quiver gauge theories (1-node and 2-node) that admit a large N limit, with their IR superconformal data?

## Primary Requirements

### Enumeration and Classification

- [ ] **ENUM-01**: Verify single-node classification (19 theories) with matter content, R-charges, a/N^2, c/N^2
- [ ] **ENUM-02**: Verify two-node classification (135 universality classes, 4560 theories) from quivers.db is complete and correct

### a-Maximization Results

- [ ] **AMAX-01**: Tabulate R-charges for all 135 two-node universality classes (76 exact + 59 numerical)
- [ ] **AMAX-02**: Compute c/N^2 and a/c ratio for all classes (derive from existing a/N^2 and Tr R, Tr R^3)

### Conformal Window Analysis

- [ ] **CONF-01**: Identify all morphologies with N_rank2=0 nodes and compute x = N_bif/2
- [ ] **CONF-02**: Classify each N_rank2=0 node as inside (3/2 < x < 3) or outside conformal window
- [ ] **CONF-03**: Verify SQCD limit: single SU(N) node reproduces 3/2 < N_f/N_c < 3

### Distribution Analysis

- [ ] **DIST-01**: Generate histogram of a/N^2 values across all universality classes
- [ ] **DIST-02**: Generate histogram of a/c ratios across all universality classes
- [ ] **DIST-03**: Analyze distribution by gauge pair type (SU-SU, SU-SO, SU-Sp, SO-SO, SO-Sp, Sp-Sp)

### Comparison

- [ ] **COMP-01**: Compare single-node results against arXiv:2510.19136 for overlapping theories

### Validation

- [ ] **VALD-01**: Verify a/c > 0 for all classified theories
- [ ] **VALD-02**: Verify anomaly-free conditions hold for all theories in the database
- [ ] **VALD-03**: Check a-theorem consistency (a_UV >= a_IR) where applicable

### Paper

- [ ] **PAPER-01**: Write JHEP paper with classification tables, distribution plots, conformal window analysis
- [ ] **PAPER-02**: Include comparison section with arXiv:2510.19136

## Follow-up Requirements

### Module 4 Extension

- **MOD4-01**: Implement boundary analysis (B_a < 0) to filter for non-trivial IR fixed points
- **MOD4-02**: Apply IR fixed point filter to all 135 classes and report which survive

### Extended Classification

- **EXTD-01**: Extend enumeration to k=3 node quivers
- **EXTD-02**: Obtain exact symbolic results for the 59 numerical-only classes

## Out of Scope

| Topic | Reason |
|-------|--------|
| Module 4 IR fixed point filter | Deferred to future paper; requires sub-quiver a-maximization not yet implemented |
| k >= 3 node quivers | Combinatorial explosion; separate investigation |
| Conformal window for rank-2 nodes | No simple SQCD-like criterion; requires different analysis |
| Exact symbolic solutions for all classes | 59 classes timed out at 30s; numerical values sufficient for classification |

## Accuracy and Validation Criteria

| Requirement | Accuracy Target | Validation Method |
|-------------|-----------------|-------------------|
| AMAX-01 | Numerical: 6 significant figures | Cross-check exact vs numerical for 76 classes |
| AMAX-02 | Same precision as a/N^2 | Derived consistently from Tr R, Tr R^3 |
| CONF-01 | Exact at large N | x = N_bif/2 is exact in the large-N limit |
| COMP-01 | Exact agreement on overlapping theories | Direct comparison with 2510.19136 tables |
| VALD-01 | Exact (a/c > 0) | Automated check across database |

## Contract Coverage

| Requirement | Decisive Output / Deliverable | Anchor / Benchmark | Prior Inputs | False Progress To Reject |
|-------------|-------------------------------|-------------------|--------------|--------------------------|
| ENUM-01 | Single-node classification table | arXiv:2510.19136 | a_maximization_large_N.py | Incomplete table missing known theories |
| ENUM-02 | Two-node classification table | Internal consistency | quivers.db | Table without R-charges or central charges |
| DIST-01/02 | Distribution histograms | None (new result) | quivers.db | Plots without physics interpretation |
| COMP-01 | Comparison section | arXiv:2510.19136 | Single-node results | Qualitative agreement only |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| ENUM-01 | Phase 1: Single-node summary | Planned |
| ENUM-02 | Phase 2: Two-node summary | Planned |
| AMAX-01 | Phase 2: Two-node summary | Planned |
| AMAX-02 | Phase 2: Two-node summary | Planned |
| CONF-01 | Phase 3: Conformal window | Planned |
| CONF-02 | Phase 3: Conformal window | Planned |
| CONF-03 | Phase 3: Conformal window | Planned |
| DIST-01 | Phase 4: Distribution analysis | Planned |
| DIST-02 | Phase 4: Distribution analysis | Planned |
| DIST-03 | Phase 4: Distribution analysis | Planned |
| COMP-01 | Phase 5: Comparison | Planned |
| VALD-01 | Phase 5: Comparison | Planned |
| VALD-02 | Phase 5: Comparison | Planned |
| VALD-03 | Phase 5: Comparison | Planned |
| PAPER-01 | Phase 6: Paper writing | Planned |
| PAPER-02 | Phase 6: Paper writing | Planned |

**Coverage:**

- Primary requirements: 16 total
- Mapped to phases: 16
- Unmapped: 0

---

_Requirements defined: 2026-04-14_
_Last updated: 2026-04-14 after initialization_
