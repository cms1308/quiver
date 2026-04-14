# Roadmap: Classification of 4D N=1 Quiver Gauge Theories

## Overview

Systematic classification of all 4D N=1 asymptotically free quiver gauge theories (1-node and 2-node) that admit a large N limit. The project summarizes existing computational results (Modules 1-3), performs conformal window and distribution analysis, compares with literature, and produces a JHEP paper.

## Phases

- [ ] **Phase 1: Single-node summary** - Compile and verify 67 single-node theories (52 SU + 6 SO + 9 Sp) with R-charges and central charges
- [ ] **Phase 2: Two-node summary** - Compile 135 universality classes with full superconformal data (R, a/N^2, c/N^2, a/c)
- [ ] **Phase 3: Conformal window analysis** - SQCD conformal window for N_rank2=0 nodes
- [ ] **Phase 4: Distribution analysis** - Central charge and a/c ratio distributions with physics interpretation
- [ ] **Phase 5: Validation and comparison** - Internal consistency checks and comparison with arXiv:2510.19136
- [ ] **Phase 6: Paper writing** - JHEP paper with all results

## Phase Details

### Phase 1: Single-node summary

**Goal:** Compile the complete single-node classification table (67 theories: 52 SU + 6 SO + 9 Sp) with matter content, N_f bound, R-charges, a/N^2, c/N^2, a/c for all SU(N), SO(N), Sp(N) theories with rank-2 matter, including Veneziano-limit and conformal manifold theories.
**Depends on:** Nothing (first phase; uses existing `a_maximization_large_N.py`)
**Requirements:** ENUM-01
**Contract coverage:**
- Claim: claim-complete-classification (single-node component)
- Anchor: ref-2510-19136 (comparison target)
**Success Criteria** (what must be TRUE):

1. Table of all 67 single-node theories (52 SU + 6 SO + 9 Sp) with gauge type, matter content, N_f bound, R-charges, a/N^2, c/N^2, a/c
2. Results formatted for direct comparison with arXiv:2510.19136 (35 non-Veneziano theories match; 31 SU Veneziano-limit theories are new)

Plans:

- [ ] 01-01: [TBD — created during /gpd:plan-phase]

### Phase 2: Two-node summary

**Goal:** Compile all 135 two-node universality classes with complete superconformal data including c/N^2 and a/c ratios.
**Depends on:** Nothing (uses existing quivers.db)
**Requirements:** ENUM-02, AMAX-01, AMAX-02
**Contract coverage:**
- Claim: claim-complete-classification (two-node component)
- Deliverable: deliv-classification-tables
**Success Criteria** (what must be TRUE):

1. All 135 classes have a/N^2, c/N^2, and a/c values
2. Data organized by gauge pair type (SU-SU, SU-SO, SU-Sp, SO-SO, SO-Sp, Sp-Sp)
3. Summary statistics per gauge pair type

Plans:

- [ ] 02-01: [TBD — created during /gpd:plan-phase]

### Phase 3: Conformal window analysis

**Goal:** Classify all N_rank2=0 nodes by SQCD conformal window criterion (3/2 < x < 3 where x = N_bif/2).
**Depends on:** Phase 2 (needs morphology data)
**Requirements:** CONF-01, CONF-02, CONF-03
**Contract coverage:**
- Claim: claim-conformal-window
- Deliverable: deliv-conformal-window-table
- Test: test-sqcd-limit
**Success Criteria** (what must be TRUE):

1. All morphologies with N_rank2=0 nodes identified
2. Each classified as inside/outside conformal window with x value
3. SQCD single-node limit verified

Plans:

- [ ] 03-01: [TBD — created during /gpd:plan-phase]

### Phase 4: Distribution analysis

**Goal:** Analyze distribution of central charges and a/c ratios across all classified theories, broken down by gauge pair type.
**Depends on:** Phase 1 (single-node data), Phase 2 (two-node data)
**Requirements:** DIST-01, DIST-02, DIST-03
**Contract coverage:**
- Claim: claim-central-charge-structure
- Deliverable: deliv-distribution-plots
- Forbidden proxy: fp-tables-without-structure (this phase provides the physics content)
**Success Criteria** (what must be TRUE):

1. Histograms of a/N^2 and a/c generated
2. Distributions broken down by gauge pair type
3. Notable features identified and discussed (clustering, bounds, outliers)

Plans:

- [ ] 04-01: [TBD — created during /gpd:plan-phase]

### Phase 5: Validation and comparison

**Goal:** Run internal consistency checks and compare single-node results against arXiv:2510.19136.
**Depends on:** Phases 1-4 (all data must be ready)
**Requirements:** COMP-01, VALD-01, VALD-02, VALD-03
**Contract coverage:**
- Test: test-single-node-comparison, test-internal-consistency
- Anchor: ref-2510-19136
**Success Criteria** (what must be TRUE):

1. All theories satisfy a/c > 0
2. Anomaly-free conditions verified programmatically
3. Single-node results match arXiv:2510.19136 on overlapping theories
4. Any discrepancies identified and explained

Plans:

- [ ] 05-01: [TBD — created during /gpd:plan-phase]

### Phase 6: Paper writing

**Goal:** Write JHEP paper presenting the complete classification with tables, distribution plots, conformal window analysis, and literature comparison.
**Depends on:** Phases 1-5 (all results and validation complete)
**Requirements:** PAPER-01, PAPER-02
**Contract coverage:**
- Deliverable: deliv-paper
- All claims, deliverables, and acceptance tests
**Success Criteria** (what must be TRUE):

1. Complete draft with introduction, methods, results, discussion, conclusions
2. Classification tables (single-node and two-node)
3. Distribution plots with interpretation
4. Conformal window analysis section
5. Comparison with arXiv:2510.19136

Plans:

- [ ] 06-01: [TBD — created during /gpd:plan-phase]

## Progress

| Phase | Plans Complete | Status | Completed |
|-------|---------------|--------|-----------|
| 1. Single-node summary | 0/TBD | Not started | - |
| 2. Two-node summary | 0/TBD | Not started | - |
| 3. Conformal window | 0/TBD | Not started | - |
| 4. Distribution analysis | 0/TBD | Not started | - |
| 5. Validation & comparison | 0/TBD | Not started | - |
| 6. Paper writing | 0/TBD | Not started | - |
