---
phase: 01-single-node-summary
verified: 2026-04-14T20:00:00Z
status: gaps_found
score: 4/6 contract targets verified
consistency_score: 10/12 physics checks passed
independently_confirmed: 6/12 checks independently confirmed
confidence: medium
gaps:
  - subject_kind: acceptance_test
    subject_id: test-type-ii-universal
    expectation: "All n_rank2=2 theories have a/N^2 = c/N^2 = (27/128)*dim(G)/N^2"
    expected_check: "Exact symbolic match for all Type II theories"
    status: failed
    category: limiting_case
    reason: >
      The acceptance test is overly broad. SU Type II theories split into three
      universality classes based on effective_T (= total T_lead of rank-2 matter / T_adj_lead):
      (1) effective_T=1 (2adj, adj+S+Sbar, etc.): a/N^2=27/128, R=1/2 -- matches universal formula.
      (2) effective_T=1/2 (A+Abar, S+Sbar, Sbar+A): a=c=0, R=0 -- does NOT match universal formula.
      (3) Veneziano theories (2A, 2S, S+A, adj+A, adj+S): non-universal values with fundamentals.
      The universal formula (27/128)*dim_lead holds ONLY for SO, Sp, and SU theories where all
      rank-2 fields are adjoint-type with effective_T = T_adj_lead (case 1 above).
      The computed results are physically correct -- the acceptance test specification is too coarse.
    computation_evidence: >
      SU A+Abar: anomaly gives R_A + R_Abar = 0, maximized at R=0, a=c=0 (not 27/128).
      SU 2A (Veneziano): a/N^2 = -3/2 + 19*sqrt(19)/48 = 0.2254 (not 27/128=0.2109).
      SO/Sp Type II: all correctly give 27/256 and 27/64 respectively.
    artifacts:
      - path: results/single_node_classification.json
        issue: "SU Type II theories have three distinct a/N^2 values, not one universal value"
    missing:
      - "Refine test-type-ii-universal to specify it applies only to SO/Sp and SU effective_T=1 theories"
    severity: minor
  - subject_kind: acceptance_test
    subject_id: test-type-i-vanish
    expectation: "All n_rank2=1 theories give a/N^2 = c/N^2 = 0 and R_rank2 = 0"
    expected_check: "All Type I theories give vanishing leading-order central charges"
    status: failed
    category: limiting_case
    reason: >
      The test holds for non-Veneziano Type I (SU adj, SO S, SO adj, Sp A, Sp adj) which give
      a=c=0 and R=0. However, SU Veneziano Type I theories (SU+A and SU+S) have non-trivial
      R-charges and NEGATIVE a/N^2 = -0.0655. The negative central charge indicates these
      theories likely do not have a non-trivial IR fixed point at leading order in N, which is
      consistent with the expected Type I behavior (trivial IR) but through a different mechanism
      than simple vanishing.
    computation_evidence: >
      SU + A: R(A) = -3 + sqrt(73)/3 = -0.152, a/N^2 = -105/16 + 73*sqrt(73)/96 = -0.0655.
      SU + S: identical values by SU S <-> A+charge-conjugation duality.
      SO adj, SO S, Sp adj, Sp A: all give a=c=0, R=0 correctly.
      SU adj: a=c=0, R=0 correctly.
    artifacts:
      - path: results/single_node_classification.json
        issue: "SU+A and SU+S Type I theories have a/N^2 < 0 and R < 0 (unphysical)"
    missing:
      - "Refine test-type-i-vanish to exclude SU Veneziano Type I or note they are unphysical"
      - "Flag theories with a < 0 or R < 0 as 'no interacting IR SCFT'"
    severity: minor
comparison_verdicts:
  - subject_kind: claim
    subject_id: claim-single-node-complete
    reference_id: ref-2510-19136
    comparison_kind: benchmark
    verdict: pass
    metric: "Theory count and central charge values for SO/Sp"
    threshold: "exact match"
  - subject_kind: claim
    subject_id: claim-single-node-complete
    reference_id: ref-universal-formulas
    comparison_kind: benchmark
    verdict: partial
    metric: "Type II and III universal formulas"
    threshold: "exact match where applicable"
suggested_contract_checks:
  - check: "Flag theories with a/N^2 < 0 as potentially lacking interacting IR SCFT"
    reason: "Negative central charge a is unphysical for interacting SCFTs; 2 theories affected"
    suggested_subject_kind: acceptance_test
    suggested_subject_id: test-positivity
    evidence_path: "results/single_node_classification.json"
---

# Phase 1 Verification: Single-Node Summary

**Phase goal:** Compile the complete single-node classification table (67 theories) with matter content, N_f bound, R-charges, a/N^2, c/N^2, a/c.

**Verified:** 2026-04-14
**Status:** gaps_found (minor -- results are correct, two acceptance tests are overly broad)
**Confidence:** MEDIUM (all key results independently confirmed; gaps are in contract spec, not physics)

## Contract Coverage

| ID | Kind | Status | Confidence | Evidence |
|----|------|--------|------------|---------|
| claim-single-node-complete | claim | VERIFIED | INDEPENDENTLY CONFIRMED | 67 theories enumerated, spot-checked computationally |
| deliv-su-table | deliverable | VERIFIED | INDEPENDENTLY CONFIRMED | 52 theories, all columns present |
| deliv-so-table | deliverable | VERIFIED | INDEPENDENTLY CONFIRMED | 6 theories, all columns present |
| deliv-sp-table | deliverable | VERIFIED | INDEPENDENTLY CONFIRMED | 9 theories, all columns present |
| deliv-raw-data | deliverable | VERIFIED | INDEPENDENTLY CONFIRMED | 67 entries, all required fields |
| deliv-computation-script | deliverable | VERIFIED | INDEPENDENTLY CONFIRMED | Script runs, produces correct output |
| test-theory-count | acceptance_test | VERIFIED | INDEPENDENTLY CONFIRMED | 52+6+9=67 confirmed |
| test-type-ii-universal | acceptance_test | PARTIAL | INDEPENDENTLY CONFIRMED | Holds for SO/Sp and SU effective_T=1; fails for SU Veneziano and SU effective_T=1/2 |
| test-type-iii-universal | acceptance_test | VERIFIED | INDEPENDENTLY CONFIRMED | SU 3adj: 1/4, SO 3adj: 1/8, Sp all Type III: 1/2 |
| test-type-i-vanish | acceptance_test | PARTIAL | INDEPENDENTLY CONFIRMED | Holds for 5/7 Type I theories; SU Veneziano Type I give a<0 |
| test-anchor-comparison | acceptance_test | VERIFIED | STRUCTURALLY PRESENT | SO(6) and Sp(9) match anchor paper structure; SU non-Ven(20) match; detailed R-charge comparison requires paper reading |
| test-ac-ratio | acceptance_test | VERIFIED | INDEPENDENTLY CONFIRMED | a/c=1 for all non-Veneziano theories with a!=0 (28 theories) |
| ref-2510-19136 | reference | PARTIAL | STRUCTURALLY PRESENT | Count and structure match; detailed equation-level comparison not yet done |
| ref-universal-formulas | reference | VERIFIED | INDEPENDENTLY CONFIRMED | Independently derived and matched |

## Required Artifacts

| Artifact | Status | Details |
|----------|--------|---------|
| results/single_node_tables.py | VERIFIED | 320+ lines, runs correctly, produces all outputs |
| results/single_node_classification.json | VERIFIED | 67 entries with all required fields |
| results/tables/su_single_node.md | VERIFIED | 52 theories, all columns present |
| results/tables/so_single_node.md | VERIFIED | 6 theories, all columns present |
| results/tables/sp_single_node.md | VERIFIED | 9 theories, all columns present |

## Computational Verification Details

### Spot-Check Results

| Expression | Test Point | Computed | Expected | Match |
|------------|-----------|----------|----------|-------|
| SU 2adj a/N^2 | R=1/2 | 27/128 | 27/128 (Type II universal) | PASS |
| SU 3adj a/N^2 | R=2/3 | 1/4 | 1/4 (Type III universal) | PASS |
| SO 2adj a/N^2 | R=1/2 | 27/256 | 27/256 = (27/128)*(1/2) | PASS |
| SO 3adj a/N^2 | R=2/3 | 1/8 | 1/8 = (1/4)*(1/2) | PASS |
| Sp 2adj a/N^2 | R=1/2 | 27/64 | 27/64 = (27/128)*2 | PASS |
| Sp 3adj a/N^2 | R=2/3 | 1/2 | 1/2 = (1/4)*2 | PASS |
| SU adj+A a/N^2 | R_adj=R_A=sqrt(89)/51+5/17 | 89*sqrt(89)/9248 + 591/4624 = 0.2186 | Independent a-max | PASS |
| SU adj+A c/N^2 | same | 53*sqrt(89)/4624 + 287/2312 = 0.2323 | Independent computation | PASS |
| SU adj+A a/c | same | 179/110 - 4*sqrt(89)/55 = 0.9412 | Independent computation | PASS |
| SU A+Abar a/N^2 | R_A=R_Abar=0 | 0 | 0 (anomaly forces R=0) | PASS |

### Limiting Cases Re-Derived

**Type II universal formula (R = 1/2):**

For n_adj adjoint fields at leading order in N:
1. Anomaly-free: T_adj + n_adj * T_adj * (R-1) = 0 => 1 + n_adj*(R-1) = 0
2. For n_adj = 2: R = 1/2
3. TrR^3 / N^2 = 1 + 2*(1/2-1)^3 = 1 - 1/4 = 3/4
4. TrR / N^2 = 1 + 2*(1/2-1) = 0
5. a/N^2 = (3/32)(3*3/4 - 0) = (3/32)(9/4) = 27/128
6. c/N^2 = (1/32)(9*3/4 - 0) = (1/32)(27/4) = 27/128
7. Scaling to SO/Sp: multiply by dim_group_lead = 1/2 (SO) or 2 (Sp)

Result: a = c = 27/128 * dim_lead. INDEPENDENTLY CONFIRMED.

**Type III universal formula (R = 2/3):**

For n_adj = 3: anomaly gives 1 + 3*(R-1) = 0 => R = 2/3.
1. TrR^3 / N^2 = 1 + 3*(-1/3)^3 = 1 - 1/9 = 8/9
2. TrR / N^2 = 1 + 3*(-1/3) = 0
3. a/N^2 = (3/32)(3*8/9 - 0) = (3/32)(8/3) = 1/4
4. c/N^2 = (1/32)(9*8/9 - 0) = 1/4

Result: a = c = 1/4 * dim_lead. INDEPENDENTLY CONFIRMED.

**SU adj+A (Veneziano, non-trivial):**

Independent two-variable a-maximization over R_adj, R_A with:
- Anomaly: R_adj + R_A/2 + R_Qbar/2 = 1, RQ = 2 - 2*R_adj - R_A
- dim_leads: adj=1, A=1/2, antifund=1
- Found 4 critical points; the physical maximum (negative-definite Hessian, H11=-3.83, H22=-1.18)
  gives R_adj = R_A = sqrt(89)/51 + 5/17, matching JSON exactly.

### Anomaly Constraint Verification

Verified Tr[R G^2] = T_adj + sum_i T_i(R_i - 1) = 0 for 12 sample theories:

| Theory | Anomaly Sum | Status |
|--------|------------|--------|
| SU 2adj | 0 | PASS |
| SU 3adj | 0 | PASS |
| SU A+Abar | 0 | PASS |
| SU adj+A+Abar | 0 | PASS |
| SO 2adj | 0 | PASS |
| SO 3adj | 0 | PASS |
| SO 2S | 0 | PASS |
| SO adj+S | 0 | PASS |
| Sp 2adj | 0 | PASS |
| Sp 3adj | 0 | PASS |
| Sp 2A | 0 | PASS |
| Sp adj+A | 0 | PASS |

All anomaly constraints satisfied. INDEPENDENTLY CONFIRMED.

### Dimensional Analysis

All quantities are dimensionless in natural units:
- a/N^2, c/N^2: pure numbers (rationals or algebraic)
- R-charges: dimensionless (range confirmed: all between -0.15 and 2.0)
- dim_lead = dim(rep)/N^2: dimensionless
- T_lead = T(rep)/N: dimensionless

Status: CONSISTENT. INDEPENDENTLY CONFIRMED.

## Physics Consistency Summary

| Check | Status | Confidence | Notes |
|-------|--------|------------|-------|
| 5.1 Dimensional analysis | CONSISTENT | INDEPENDENTLY CONFIRMED | All quantities dimensionless in natural units |
| 5.2 Numerical spot-check | PASS | INDEPENDENTLY CONFIRMED | 10 spot-checks, all match |
| 5.3 Limiting cases | PASS | INDEPENDENTLY CONFIRMED | Type II, III universals re-derived from scratch |
| 5.4 Cross-check | PASS | INDEPENDENTLY CONFIRMED | 2-variable a-max for SU adj+A confirms JSON |
| 5.5 Intermediate spot-check | N/A | N/A | Single-step computations, no long derivation chain |
| 5.6 Symmetry | PASS | INDEPENDENTLY CONFIRMED | a/c=1 for all non-Veneziano (28 theories); S<->A duality confirmed |
| 5.7 Conservation (anomaly) | PASS | INDEPENDENTLY CONFIRMED | 12/12 sample theories satisfy Tr[R G^2]=0 |
| 5.8 Math consistency | PASS | STRUCTURALLY PRESENT | Hessian check confirms maximum (not saddle); index structure correct |
| 5.10 Literature agreement | PASS | STRUCTURALLY PRESENT | SO(6) and Sp(9) theory counts and universality match arXiv:2510.19136 |
| 5.11 Plausibility | WARNING | INDEPENDENTLY CONFIRMED | 2 theories (SU+A, SU+S) have a/N^2 < 0, indicating no IR SCFT |
| 5.13 Thermodynamic consistency | N/A | N/A | Not applicable (SCFT, not thermal) |
| 5.14 Spectral/analytic | N/A | N/A | Not applicable |

## Forbidden Proxy Audit

No forbidden proxies defined in contract.

## Comparison Verdict Ledger

| Subject ID | Comparison Kind | Verdict | Threshold | Notes |
|-----------|----------------|---------|-----------|-------|
| claim-single-node-complete | benchmark (ref-2510-19136) | pass | Theory count and central charge structure | SO(6), Sp(9) match; SU(20 non-Ven) match |
| claim-single-node-complete | benchmark (ref-universal-formulas) | partial | Exact match where applicable | Holds for SO/Sp/SU-adj; SU Veneziano/zero-excess are different |

## Discrepancies Found

| Severity | Location | Evidence | Root Cause | Fix |
|----------|----------|----------|------------|-----|
| minor | test-type-ii-universal | SU A+Abar gives a=0, not 27/128; SU 2A gives 0.225, not 0.211 | Test spec too broad: SU Type II has 3 universality classes | Restrict test to SO/Sp and SU effective_T=1 |
| minor | test-type-i-vanish | SU+A gives a=-0.065, R=-0.15 | SU Veneziano Type I have fundamentals that drive non-trivial (unphysical) R-charges | Restrict test to non-Veneziano Type I |
| info | Positivity | SU+A and SU+S have a/N^2 < 0 | These Veneziano theories likely have no interacting IR SCFT | Flag in tables; consistent with Type I triviality |

## Suggested Contract Checks

| Check | Reason | Suggested ID | Evidence |
|-------|--------|-------------|----------|
| Positivity check: a/N^2 >= 0 for all claimed SCFTs | 2 theories violate; may not have interacting IR fixed points | test-positivity | results/single_node_classification.json |

## Anti-Patterns Found

| Pattern | Severity | File | Physics Impact |
|---------|----------|------|---------------|
| ASSERT_CONVENTION present | INFO | single_node_tables.py | Conventions declared in file header |
| No TODO/FIXME found | INFO | All artifacts | Code is complete, no stubs |

## Expert Verification Required

- **Detailed comparison with arXiv:2510.19136 Tables 2, 27, 35:** Requires reading the specific tables in the paper to verify R-charge expressions match term-by-term. Currently verified at the level of counting and universality structure. Domain: superconformal field theory.
- **Physical interpretation of SU Veneziano Type I theories with a < 0:** These theories formally give negative central charge at leading order. Needs expert assessment of whether this indicates (a) no interacting IR fixed point, (b) the large-N limit breaks down, or (c) a different IR phase. Domain: large-N gauge theory.

## Confidence Assessment

Overall confidence is MEDIUM because:
1. All key computations (universal formulas, anomaly constraints, a-maximization) independently confirmed via Python/SymPy execution
2. Theory count (67 = 52+6+9) confirmed
3. All required artifacts exist with correct structure and content
4. Two acceptance tests fail due to overly broad specifications, not physics errors
5. Detailed equation-level comparison with arXiv:2510.19136 not yet performed (hybrid test)
6. Two theories have a < 0 which needs physical interpretation

The computed results themselves appear correct -- the gaps are in the contract specification being too coarse about SU Veneziano behavior, not in the actual physics.

## Gaps Summary

Two minor gaps, both stemming from the contract's acceptance tests being overly broad for SU theories with fundamentals (Veneziano limit):

1. **test-type-ii-universal**: The universal formula a = c = (27/128)*dim_lead applies to SO, Sp, and SU theories where effective_T = 1 (all fields are adjoint-type). SU theories with mixed rank-2 types (A+Abar etc.) or Veneziano fundamentals have different universality classes. Fix: refine the test scope.

2. **test-type-i-vanish**: The vanishing a = c = 0 at leading order holds for non-Veneziano Type I theories. SU Veneziano Type I (SU+A, SU+S) give non-zero but negative a, consistent with absence of interacting IR SCFT but through a different mechanism. Fix: refine the test to non-Veneziano or add positivity check.

Neither gap indicates an error in the computed results.
