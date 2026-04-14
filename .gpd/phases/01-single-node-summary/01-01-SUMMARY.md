---
phase: 01-single-node-summary
plan: 01
depth: full
one-liner: "Complete single-node classification: 67 theories (52 SU + 6 SO + 9 Sp) with exact symbolic R-charges and central charges, validated against universal formulas and arXiv:2510.19136"
subsystem: [computation, validation, classification]
tags: [a-maximization, large-N, quiver-gauge-theory, central-charges, R-charges]

plan_contract_ref: .gpd/phases/01-single-node-summary/01-01-PLAN.md#/contract

contract_results:
  claims:
    claim-single-node-complete:
      status: passed
      summary: "Complete classification of all 67 single-node theories with exact symbolic R-charges, a/N^2, c/N^2, a/c. SU universality governed by effective_T = n_adj + n_tensor/2."
      linked_ids: [deliv-su-table, deliv-so-table, deliv-sp-table, deliv-raw-data, deliv-computation-script, test-theory-count, test-type-ii-universal, test-type-iii-universal, test-type-i-vanish, test-anchor-comparison, test-ac-ratio, test-anomaly-constraint, ref-2510-19136, ref-universal-formulas, ref-a-max-code]
  deliverables:
    deliv-su-table:
      status: passed
      path: results/tables/su_single_node.md
      summary: "52 SU theories with matter, delta, N_f bound, R-charges, a/N^2, c/N^2, a/c"
    deliv-so-table:
      status: passed
      path: results/tables/so_single_node.md
      summary: "6 SO theories with all required columns"
    deliv-sp-table:
      status: passed
      path: results/tables/sp_single_node.md
      summary: "9 Sp theories with all required columns"
    deliv-raw-data:
      status: passed
      path: results/single_node_classification.json
      summary: "67-entry JSON with gauge_type, matter, R_charges, a_over_N2, c_over_N2, all exact symbolic"
    deliv-computation-script:
      status: passed
      path: results/single_node_tables.py
      summary: "Self-contained script: enumerate, compute, validate, format. 477 automated checks, all pass."
  acceptance_tests:
    test-theory-count:
      status: passed
      summary: "52 SU + 6 SO + 9 Sp = 67 total, matching expected counts"
      linked_ids: [claim-single-node-complete, deliv-raw-data]
    test-type-ii-universal:
      status: passed
      summary: "All SO/Sp Type II and SU 2adj theories give a/N^2 = c/N^2 = (27/128)*dim(G)/N^2 with R=1/2 exactly. SU mixed-tensor Type II theories form their own universality classes (see Key Insight)."
      linked_ids: [claim-single-node-complete, deliv-raw-data, ref-universal-formulas]
    test-type-iii-universal:
      status: passed
      summary: "All SO/Sp Type III and SU 3adj theories give a/N^2 = c/N^2 = (1/4)*dim(G)/N^2 with R=2/3 exactly. SU mixed-tensor theories with eff_T=3 also give 1/4."
      linked_ids: [claim-single-node-complete, deliv-raw-data, ref-universal-formulas]
    test-type-i-vanish:
      status: passed
      summary: "All SO/Sp Type I and SU non-Veneziano eff_T=1 theories give a=c=0, R=0. SU Veneziano Type I (A, S) give negative leading-order a (subleading artifact)."
      linked_ids: [claim-single-node-complete, deliv-raw-data]
    test-anchor-comparison:
      status: passed
      summary: "All 6 SO theories match Table 27, all 9 Sp theories match Table 35, SU adj-only theories match Table 2. Extra theory S+2A+3Abar (conformal manifold) is present in our classification but absent from paper Table 2."
      linked_ids: [claim-single-node-complete, deliv-su-table, deliv-so-table, deliv-sp-table, ref-2510-19136]
    test-ac-ratio:
      status: passed
      summary: "a/c = 1 exactly for all non-Type-I theories with eff_T = 2 or 3 and delta = 0. Veneziano theories with delta != 0 have a/c != 1 (physical, not a bug)."
      linked_ids: [claim-single-node-complete, deliv-raw-data]
    test-anomaly-constraint:
      status: passed
      summary: "Anomaly constraint T(adj) + sum T_i(R_i-1) = 0 verified exactly (symbolic zero) for all 67 theories"
      linked_ids: [claim-single-node-complete, deliv-raw-data]
  references:
    ref-2510-19136:
      status: completed
      completed_actions: [compare, cite]
      missing_actions: []
      summary: "Compared SO Table 27 (6 theories), Sp Table 35 (9 theories), and SU adj-only theories from Table 2. All overlapping values match exactly."
    ref-universal-formulas:
      status: completed
      completed_actions: [compare]
      missing_actions: []
      summary: "Type II (27/128)*dim(G) and Type III (1/4)*dim(G) verified for all applicable theories"
    ref-a-max-code:
      status: completed
      completed_actions: [use]
      missing_actions: []
      summary: "a_maximization_large_N.py used as computation engine; no modifications needed"
  forbidden_proxies:
    fp-incomplete-table:
      status: rejected
      notes: "All 67 theories present including 31 Veneziano-limit SU theories"
    fp-numerical-only:
      status: rejected
      notes: "All R-charges and central charges are exact sympy expressions (Rational or algebraic)"
    fp-no-validation:
      status: rejected
      notes: "477 automated validation checks run, covering universal formulas, anomaly constraints, universality, bounds, and anchor comparison"
  uncertainty_markers:
    weakest_anchors:
      - "S+2A+3Abar (conformal manifold) has no independent benchmark"
      - "31 Veneziano SU theories lack direct comparison (paper Table 45 not extracted)"
    unvalidated_assumptions:
      - "Leading-order large N sufficient for Type I Veneziano theories (a~O(N) not O(N^2))"

provides:
  equations:
    - label: "universal-type-II"
      expression: "a/N^2 = c/N^2 = (27/128) * dim(G)/N^2 for eff_T = 2"
      confidence: HIGH
    - label: "universal-type-III"
      expression: "a/N^2 = c/N^2 = (1/4) * dim(G)/N^2 for eff_T = 3"
      confidence: HIGH
    - label: "su-universality"
      expression: "SU universality class = (effective_T, |delta|) where effective_T = n_adj + n_tensor/2"
      confidence: HIGH
  parameters: []
  approximations:
    - "Large N leading order: a/N^2, c/N^2, R as O(1) quantities"
  figures: []

completed: true
---

# Plan 01-01 Summary: Single-Node Classification

## Performance

| Metric | Value |
|--------|-------|
| Tasks completed | 2/2 |
| Total checks | 477 |
| Checks passed | 477 |
| Checks failed | 0 |
| Theories computed | 67 |
| Computation time | ~2 min (sympy exact a-maximization) |

## Key Results

### Theory Count [CONFIDENCE: HIGH]

| Gauge type | Count | Type I | Type II | Type III | Higher |
|------------|-------|--------|---------|----------|--------|
| SU(N) | 52 | 3 | 12 | 14 | 23 |
| SO(N) | 6 | 2 | 3 | 1 | 0 |
| Sp(N) | 9 | 2 | 3 | 4 | 0 |
| **Total** | **67** | 7 | 18 | 19 | 23 |

### Key Insight: SU Universality Classes [CONFIDENCE: HIGH]

For SO/Sp, all rank-2 representations have identical leading-order Dynkin index (T_lead = 1), so universality depends only on n_rank2. For SU, adj has T_lead = 1 while S/Sbar/A/Abar have T_lead = 1/2. The true universality parameter is:

**effective_T = n_adj + (n_S + n_Sbar + n_A + n_Abar) / 2**

All SU theories with the same (effective_T, |chiral_excess_coefficient|) give identical leading-order a/N^2, c/N^2, and R-charges. There are 21 non-Veneziano SU theories falling into 3 universality classes:

| effective_T | a/N^2 | c/N^2 | a/c | # non-Veneziano theories |
|-------------|-------|-------|-----|--------------------------|
| 1 | 0 | 0 | n/a | 4 (adj, A+Abar, S+Sbar, Sbar+A) |
| 2 | 27/128 | 27/128 | 1 | 10 (2adj, A+A+Abar+Abar, adj+A+Abar, ...) |
| 3 | 1/4 | 1/4 | 1 | 7 (3adj, adj+2A+2Abar, 3A+3Abar, ...) |

### Veneziano Type I Theories

2 SU Type I theories (A alone, S alone) have nonzero chiral excess (|delta| = N), forcing O(N) fundamentals. At leading order these give negative a/N^2 ~ -0.065, indicating the true central charge is O(N) not O(N^2). The leading-order R-charges are slightly negative (R ~ -0.15), a subleading artifact. These require subleading analysis (deferred).

### SO/Sp Universality [CONFIDENCE: HIGH]

Confirmed: within each (gauge_type, n_rank2) class, all theories give identical a/N^2, c/N^2, R-charges. This follows from T_rep_lead = 1 for all SO/Sp representations.

| Gauge | n_rank2 | a/N^2 | R |
|-------|---------|-------|---|
| SO | 1 | 0 | 0 |
| SO | 2 | 27/256 | 1/2 |
| SO | 3 | 1/8 | 2/3 |
| Sp | 1 | 0 | 0 |
| Sp | 2 | 27/64 | 1/2 |
| Sp | 3 | 1/2 | 2/3 |

## Conventions

| Convention | Value |
|------------|-------|
| Gauge groups | SU(N), SO(N), Sp(N) = USp(2N) |
| Dynkin index | T(fund_SU) = 1/2, T(V_SO) = 1, T(f_Sp) = 1/2 |
| Representation notation | adj for all; S/Sbar/A/Abar for SU; S for SO symmetric; A for Sp antisymmetric |
| Large N scaling | a/N^2, c/N^2 as O(1) leading-order quantities |
| Type classification | By n_rank2 (matching paper), with SU universality by effective_T |

## Task Commits

| Task | Hash | Description |
|------|------|-------------|
| 1 | 4220345 | Enumerate 67 theories, compute exact data, format tables + JSON |
| 2 | 2571c5a | Comprehensive validation: 477 checks, 0 failures |

## Files Written

| File | Description |
|------|-------------|
| results/single_node_tables.py | Computation + validation script (self-contained) |
| results/single_node_classification.json | Machine-readable JSON (67 entries) |
| results/tables/su_single_node.md | SU table (52 theories) |
| results/tables/so_single_node.md | SO table (6 theories) |
| results/tables/sp_single_node.md | Sp table (9 theories) |

## Validations Performed

1. **Theory count**: 52 + 6 + 9 = 67 (exact)
2. **SO/Sp Type II**: a = c = (27/128) * dim(G)/N^2, R = 1/2 (all 6 theories)
3. **SO/Sp Type III**: a = c = (1/4) * dim(G)/N^2, R = 2/3 (all 5 theories)
4. **SU adj-only**: 2adj -> 27/128, 3adj -> 1/4 (matches universal formula)
5. **Type I vanishing**: a = c = 0 for all SO/Sp Type I and SU non-Veneziano eff_T=1
6. **Anomaly-free R**: T(adj) + sum T_i(R_i - 1) = 0 verified symbolically for all 67
7. **SO/Sp universality**: identical values within each (gauge_type, n_rank2) class
8. **SU universality**: identical values within each (effective_T, |delta|) class
9. **a/c bounds**: 1/2 <= a/c <= 5/4 for all non-Type-I theories
10. **R-charge range**: 0 <= R <= 2 for all (Type I Veneziano exceptions documented)
11. **Anchor comparison**: All SO/Sp and SU adj-only match arXiv:2510.19136

## Deviations

None. All tasks completed as planned.

## Issues

- S+2A+3Abar (conformal manifold) present in our classification but absent from arXiv:2510.19136 Table 2. May be excluded due to finite-N unitarity violation. Flagged for verification in later phases.
- 31 Veneziano-limit SU theories lack direct comparison target until paper Table 45 is extracted (deferred to Phase 5).

## Open Questions

- Whether the effective_T universality for SU constitutes a novel observation worth highlighting in the paper
- Distribution of a/c ratios across Veneziano theories (relevant for Phase 5 analysis)

## Next Phase Readiness

This plan provides:
- Complete single-node baseline for comparison in Phase 2 (two-node classification)
- JSON data file importable by downstream analysis scripts
- Validated computation methodology (a_maximization_large_N.py confirmed working)
