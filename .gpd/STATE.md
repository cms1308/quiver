# Research State

## Project Reference

See: .gpd/PROJECT.md (updated 2026-04-14)

**Core research question:** What is the complete classification of 4D N=1 asymptotically free quiver gauge theories (1-node and 2-node) that admit a large N limit, with their IR superconformal data?
**Current focus:** Phase 2 — Two-node summary

## Current Position

**Current Phase:** 2
**Current Phase Name:** Two-node summary
**Total Phases:** 6
**Current Plan:** 1
**Total Plans in Phase:** 1
**Status:** Plan 02-01 complete — ready for verification
**Last Activity:** 2026-04-14
**Last Activity Description:** Two-node classification: 326 classes, 27 checks pass

**Progress:** [██████████] 100%

## Active Calculations

None.

## Intermediate Results

- 67 single-node theories with exact symbolic R-charges and central charges
- SU universality: effective_T = n_adj + n_tensor/2 determines universality class
- Type II: a/N² = c/N² = (27/128)·dim(G)/N², R_rank2 = 1/2
- Type III: a/N² = c/N² = (1/4)·dim(G)/N², R_rank2 = 2/3
- Type I SO/Sp: a = c = 0 at leading order
- 326 two-node universality classes with complete superconformal data
- Three-tier: 285 unitary + 10 Type-I + 31 non-unitary
- a/c range [0.500, 1.000] for unitary (mean 0.939, median 0.956)
- Non-unitary fraction: SU-SO 14.8%, SU-SU 8.2%, SU-Sp 10.8%, SO/Sp pairs 0%

## Open Questions

- Why does a/c never exceed 1.0 for two-node unitary theories?
- Why is non-unitary fraction zero for orthosymplectic-only pairs?
- Class 31: which value (exact or numerical) is correct for c?
- How to present the 326 classes compactly in JHEP format (tables vs appendix)

## Performance Metrics

| Label | Duration | Tasks | Files |
|-------|----------|-------|-------|
| 01-01 | ~12 min | 2 | 6 |
| 02-01 | ~5 min | 2 | 4 |

## Accumulated Context

### Decisions

- [Init]: Skip Module 4 for this paper — focus on Modules 1-3 results first
- [Init]: Numerical values sufficient — 59 classes without exact symbolic results are acceptable
- [Init]: SQCD conformal window only for N_rank2=0 nodes
- [Init]: Target JHEP
- [01-01]: SU universality class determined by (effective_T, |delta|) not (n_rank2, n_adj)
- [01-01]: Type I Veneziano theories (SU with A or S alone) give negative leading-order a/N² — subleading artifact
- [01-01]: S+2A+3Abar included in classification despite absence from paper Table 2
- [02-01]: Use exact symbolic c when available (181 classes), numerical otherwise (145)
- [02-01]: Class 31 discrepancy documented, not treated as error
- [02-01]: 7537 classified theories (220 unclassified excluded from counts)

### Active Approximations

None yet.

**Convention Lock:**

No conventions locked yet.

### Propagated Uncertainties

None yet.

### Pending Todos

None yet.

### Blockers/Concerns

None yet.

## Session Continuity

**Last session:** 2026-04-14
**Stopped at:** Phase 2 execution complete, awaiting verification
**Resume file:** .gpd/phases/02-two-node-summary/02-01-SUMMARY.md
