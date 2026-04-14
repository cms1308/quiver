# Classification of 4D N=1 Quiver Gauge Theories

## What This Is

Classification of all 4D N=1 asymptotically free quiver gauge theories (1-node and 2-node) that admit a large N limit, with their IR superconformal R-charges, central charges, and conformal window analysis. The deliverable is a JHEP paper presenting classification tables, distribution analysis, and comparison with existing literature.

## Core Research Question

What is the complete classification of 4D N=1 asymptotically free quiver gauge theories (1-node and 2-node) that admit a large N limit, with their IR superconformal data (R-charges, central charges, a/c ratios)?

## Scoping Contract Summary

### Contract Coverage

- Complete classification tables (1-node + 2-node): All universality classes with R-charges, a/N^2, c/N^2, a/c
- Central charge distribution analysis: Histograms of a/N^2 and a/c ratios revealing structural patterns
- Conformal window analysis: N_rank2=0 nodes classified by SQCD criterion
- Single-node benchmark: Comparison with arXiv:2510.19136
- False progress to reject: Raw tables without distribution analysis or conformal window filtering

### User Guidance To Preserve

- **User-stated observables:** Distribution of central charges and a/c ratios
- **User-stated deliverables:** JHEP paper with classification tables, distribution plots, conformal window tables
- **Must-have references / prior outputs:** arXiv:2510.19136 (single-node classification), quivers.db (4560 theories)
- **Stop / rethink conditions:** Single-node results disagree with arXiv:2510.19136

### Scope Boundaries

**In scope**

- Single-node classification (19 theories) as comparison section
- Two-node classification (135 universality classes, 4560 theories) from Modules 1-3
- R-charges and central charges (a/N^2, c/N^2) for all classes
- Conformal window analysis for N_rank2=0 nodes using SQCD criterion
- Distribution analysis of central charges and a/c ratios
- Comparison with arXiv:2510.19136 for single-node results

**Out of scope**

- Module 4 boundary analysis (B_a < 0 IR fixed point filter)
- k >= 3 node quivers
- Exact symbolic results for the 59 numerical-only classes
- Conformal window analysis for nodes with rank-2 matter

### Active Anchor Registry

- **ref-2510-19136**: arXiv:2510.19136 (single-node classification)
  - Why it matters: Contains partial single-node classification to compare against
  - Carry forward: planning, execution, verification, writing
  - Required action: read, compare, cite

### Carry-Forward Inputs

- `quivers.db` — pre-built database with 4560 theories in 135 universality classes
- `a_maximization_large_N.py` — large-N a-maximization solver

### Skeptical Review

- **Weakest anchor:** Comparison with 2510.19136 is partial — they do not have the full single-node classification
- **Unvalidated assumptions:** Equal-rank assumption (all nodes share same N); large N leading order is sufficient
- **Competing explanation:** None identified
- **Disconfirming observation:** Single-node results disagree with 2510.19136 on overlapping theories
- **False progress to reject:** Raw classification tables without distribution analysis or conformal window filtering

### Open Contract Questions

- How to present the 135 classes compactly in JHEP format (tables vs appendix)
- Whether distribution of a/c ratios reveals interesting structure worth highlighting

## Research Questions

### Answered

- [x] What are the UV asymptotic freedom conditions for quiver gauge theories at large N? — Modules 1-2 established
- [x] How to enumerate all anomaly-free 2-node quiver gauge theories? — Module 2 implemented, 4560 theories found
- [x] What are the IR superconformal R-charges from a-maximization? — Module 3 implemented, 135 universality classes computed
- [x] What is the conformal window criterion for N_rank2=0 nodes? — SQCD analysis: x = N_bif/2, window 3/2 < x < 3

### Active

- [ ] What is the distribution of central charges (a/N^2) and a/c ratios across all classified theories?
- [ ] Which morphologies have N_rank2=0 nodes outside the conformal window?
- [ ] How do the single-node results compare with arXiv:2510.19136?
- [ ] What structural patterns emerge from the classification?

### Out of Scope

- Module 4 IR fixed point boundary analysis — deferred to future work
- k >= 3 node quivers — combinatorial explosion, separate investigation
- Conformal window for rank-2 nodes — no simple SQCD-like criterion available

## Research Context

### Physical System

4D N=1 supersymmetric gauge theories with product gauge groups (quiver gauge theories). Each node carries SU(N), SO(N), or Sp(N) gauge symmetry. Matter consists of bifundamentals between nodes plus rank-2 tensors (adjoint, symmetric, antisymmetric) at individual nodes.

### Theoretical Framework

4D N=1 superconformal field theory. Key tools: NSVZ exact beta function, Intriligator-Wecht a-maximization, 't Hooft anomaly matching, large N expansion.

### Key Parameters and Scales

| Parameter | Symbol | Regime | Notes |
|-----------|--------|--------|-------|
| Rank | N | Large N limit | All nodes share same N |
| Flavor count | N_f | O(1) | Fundamental flavors, subleading at large N |
| Beta function coeff | b_0 | > 0 | UV asymptotic freedom |
| Central charge | a/N^2 | O(1) | Leading large-N behavior |
| Conformal window param | x = N_bif/2 | 3/2 < x < 3 | For N_rank2=0 nodes only |

### Known Results

- Seiberg's SQCD conformal window: 3/2 N < N_f < 3N for SU(N) — Seiberg (1994)
- a-maximization: Intriligator-Wecht (2003), hep-th/0304128
- Boundary analysis for IR fixed points: Beem et al., hep-th/0502049
- Single-node classification (partial): arXiv:2510.19136

### What Is New

Complete systematic classification of all 1-node and 2-node N=1 AF quiver gauge theories at large N, with full superconformal data. Previous work classified single-node theories partially; no systematic two-node classification exists.

### Target Venue

JHEP — natural venue for 4D N=1 SCFT classification results

### Computational Environment

- Local workstation with Python (numpy, scipy, sympy) and Mathematica (wolframscript)
- SQLite database for theory storage and querying
- Mathematica NSolve for numerical a-maximization

## Notation and Conventions

See `.gpd/CONVENTIONS.md` for all notation and sign conventions.
See `.gpd/NOTATION_GLOSSARY.md` for symbol definitions.

## Unit System

Natural units (hbar = c = 1). Dimensionless quantities throughout — R-charges, central charges normalized by N^2.

## Requirements

See `.gpd/REQUIREMENTS.md` for the detailed requirements specification.

Key requirement categories: ENUM (enumeration), AMAX (a-maximization), CONF (conformal window), DIST (distribution analysis), COMP (comparison), PAPER (paper writing)

## Key References

- arXiv:2510.19136 — single-node classification benchmark (must compare)

## Constraints

- **Computational:** Mathematica NSolve for a-maximization; 30s timeout per theory for exact symbolic solve
- **Data:** quivers.db already built with 4560 theories; 59 classes numerical-only
- **Scope:** Module 4 (IR fixed point filter) explicitly deferred

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Skip Module 4 for this paper | Focus on Modules 1-3 results first | Decided |
| Numerical values sufficient | 59 classes without exact symbolic results are acceptable | Decided |
| SQCD conformal window only for N_rank2=0 | No simple criterion for rank-2 nodes | Decided |
| Target JHEP | Standard venue for SCFT classification | Decided |

---

_Last updated: 2026-04-14 after initialization_
