# Research State

## Project Reference

See: .gpd/PROJECT.md (updated 2026-04-14)

**Core research question:** What is the complete classification of 4D N=1 asymptotically free quiver gauge theories (1-node and 2-node) that admit a large N limit, with their IR superconformal data?
**Current focus:** Phase 3 — Conformal window analysis

## Current Position

**Current Phase:** 3
**Current Phase Name:** Conformal window analysis
**Total Phases:** 6
**Current Plan:** 1
**Total Plans in Phase:** 1
**Status:** Plan 03-01 complete — ready for verification
**Last Activity:** 2026-04-15
**Last Activity Description:** Conformal window: 237 nodes classified, 79 inside, 7 validation checks pass

**Progress:** [██████████] 100%

## Active Calculations

None.

## Intermediate Results

- 67 single-node theories with exact symbolic R-charges and central charges
- SU universality: effective_T = n_adj + n_tensor/2 determines universality class
- Type II: a/N² = c/N² = (27/128)·dim(G)/N², R_rank2 = 1/2
- Type III: a/N² = c/N² = (1/4)·dim(G)/N², R_rank2 = 2/3
- Type I SO/Sp: a = c = 0 at leading order
- 243 two-node universality classes: 242 with superconformal data + 1 merged below-window class (478 theories, a=c=R=None)
- Non-unitary fraction: SU-SO 14.8%, SU-SU 8.2%, SU-Sp 10.8%, SO/Sp pairs 0%
- 237 N_rank2=0 nodes classified: 79 inside (33%), 130 below (55%), 18 lower boundary, 10 upper boundary
- SO-SO highest inside fraction (56%), SU-Sp lowest (23%)
- 7 morphologies with both nodes inside conformal window out of 38 with both N_rank2=0
- Rank multipliers increase inside fraction: 39% for rank!=(1,1) vs 28% for rank(1,1)

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
| 03-01 | ~8 min | 2 | 4 |

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
- [03-01]: Strict inequality 3/2 < x < 3 for inside classification; boundary values at x=3/2 and x=3 reported separately
- [03-01]: Per-node classification (not per-morphology); general formula x = N_bif * T_lead * m_other/m_self with rank multipliers
- [03-01]: Below-window nodes (x < 3/2) → a, c, R set to None; 84 classes and 571 theories affected

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

**Last session:** 2026-04-15
**Stopped at:** Phase 3 plan 01 complete, awaiting verification
**Resume file:** .gpd/phases/03-conformal-window-analysis/03-01-SUMMARY.md
