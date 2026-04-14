# Research State

## Project Reference

See: .gpd/PROJECT.md (updated 2026-04-14)

**Core research question:** What is the complete classification of 4D N=1 asymptotically free quiver gauge theories (1-node and 2-node) that admit a large N limit, with their IR superconformal data?
**Current focus:** Phase 1 — Single-node summary

## Current Position

**Current Phase:** 1
**Current Phase Name:** Single-node summary
**Total Phases:** 6
**Current Plan:** 1
**Total Plans in Phase:** 1
**Status:** Phase 1 complete — verified
**Last Activity:** 2026-04-14
**Last Activity Description:** Single-node classification: 67 theories, 477 validations pass

**Progress:** [██████████] 100%

## Active Calculations

None.

## Intermediate Results

- 67 single-node theories with exact symbolic R-charges and central charges
- SU universality: effective_T = n_adj + n_tensor/2 determines universality class
- Type II: a/N² = c/N² = (27/128)·dim(G)/N², R_rank2 = 1/2
- Type III: a/N² = c/N² = (1/4)·dim(G)/N², R_rank2 = 2/3
- Type I SO/Sp: a = c = 0 at leading order

## Open Questions

- How to present the 135 classes compactly in JHEP format (tables vs appendix)
- Whether distribution of a/c ratios reveals interesting structure worth highlighting

## Performance Metrics

| Label | Duration | Tasks | Files |
|-------|----------|-------|-------|
| 01-01 | ~12 min | 2 | 6 |

## Accumulated Context

### Decisions

- [Init]: Skip Module 4 for this paper — focus on Modules 1-3 results first
- [Init]: Numerical values sufficient — 59 classes without exact symbolic results are acceptable
- [Init]: SQCD conformal window only for N_rank2=0 nodes
- [Init]: Target JHEP
- [01-01]: SU universality class determined by (effective_T, |delta|) not (n_rank2, n_adj)
- [01-01]: Type I Veneziano theories (SU with A or S alone) give negative leading-order a/N² — subleading artifact
- [01-01]: S+2A+3Abar included in classification despite absence from paper Table 2

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
**Stopped at:** Phase 1 execution complete, awaiting verification
**Resume file:** .gpd/phases/01-single-node-summary/01-01-SUMMARY.md
