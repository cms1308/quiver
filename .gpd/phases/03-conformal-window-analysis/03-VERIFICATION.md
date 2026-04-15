---
phase: 03-conformal-window-analysis
verified: 2026-04-15T12:00:00Z
status: passed
score: 3/3 contract targets verified
consistency_score: 8/8 physics checks passed
independently_confirmed: 7/8 checks independently confirmed
confidence: high
---

# Phase 3 Verification: Conformal Window Analysis

**Phase goal:** Classify all N_rank2=0 nodes by SQCD conformal window criterion (3/2 < x < 3 where x_a = N_bif * T_bifund_lead(ga,gb)[a] * m_b / m_a).

**Status: PASSED** | Confidence: HIGH | Score: 3/3 contract targets verified

## Contract Coverage

| ID | Kind | Status | Confidence | Evidence |
|---|---|---|---|---|
| claim-conformal-window | claim | VERIFIED | INDEPENDENTLY CONFIRMED | x computed for all 237 qualifying nodes from 199 morphologies; formula verified from first principles |
| deliv-conformal-window-json | deliverable | VERIFIED | INDEPENDENTLY CONFIRMED | results/conformal_window_classification.json exists, 237 entries, all required fields present |
| deliv-conformal-window-table | deliverable | VERIFIED | INDEPENDENTLY CONFIRMED | results/tables/conformal_window.md exists with distribution analysis and physics interpretation |
| deliv-conformal-window-script | deliverable | VERIFIED | INDEPENDENTLY CONFIRMED | results/conformal_window.py runs successfully, produces correct output |
| test-sqcd-limit | acceptance_test | VERIFIED | INDEPENDENTLY CONFIRMED | SU-SU N_bif=3->x=3/2, N_bif=4->x=2, N_bif=6->x=3 all confirmed |
| test-af-consistency | acceptance_test | VERIFIED | INDEPENDENTLY CONFIRMED | 0 entries with x > 3; 10 entries at x=3 boundary |
| test-su-su-symmetry | acceptance_test | VERIFIED | INDEPENDENTLY CONFIRMED | 12 SU-SU morphologies with both nodes and equal ranks verified x_0=x_1 |

## Required Artifacts

| Artifact | Expected | Status | Details |
|---|---|---|---|
| results/conformal_window.py | Script computing x for all N_rank2=0 nodes | VERIFIED | 544 lines, self-contained, runs cleanly |
| results/conformal_window_classification.json | JSON with per-node classification | VERIFIED | 237 entries with morph_id, node_index, gauge_pair, gauge_type, x, x_float, classification, N_bif, T_lead, rank_self, rank_other |
| results/tables/conformal_window.md | Markdown table with distribution and interpretation | VERIFIED | Summary, distribution by gauge pair, rank sector comparison, physics interpretation section |

## Computational Verification Details

### 5.1 Dimensional Analysis

All quantities are dimensionless:
- x = N_bif * T_lead * m_other / m_self: [integer] * [1] * [1] / [1] = [dimensionless] -- CONSISTENT
- b_0 = m_a * N * (3 - x): [1] * [1] * [1] = [dimensionless] -- CONSISTENT
- Classification thresholds 3/2 and 3 are dimensionless -- CONSISTENT

**Status:** CONSISTENT | **Confidence:** INDEPENDENTLY CONFIRMED

### 5.2 Numerical Spot-Checks

Evaluated x = N_bif * T_lead * m_other / m_self at 8 test points:

| gauge_pair | node | N_bif | rank | T_lead | Expected x | Got x | Match |
|---|---|---|---|---|---|---|---|
| SU-SU | 0 | 4 | (1,1) | 1/2 | 2 | 2 | PASS |
| SU-SU | 0 | 3 | (1,1) | 1/2 | 3/2 | 3/2 | PASS |
| SU-SU | 0 | 6 | (1,1) | 1/2 | 3 | 3 | PASS |
| SO-SO | 0 | 2 | (1,1) | 1 | 2 | 2 | PASS |
| SU-Sp | 0 | 1 | (1,1) | 1 | 1 | 1 | PASS |
| SU-Sp | 1 | 1 | (1,1) | 1/2 | 1/2 | 1/2 | PASS |
| SO-Sp | 0 | 1 | (1,1) | 2 | 2 | 2 | PASS |
| SO-Sp | 1 | 1 | (1,1) | 1/2 | 1/2 | 1/2 | PASS |

Additional spot-checks for asymmetric rank cases:
- SU-Sp node 0, N_bif=2, rank(1,1): x = 2*1*1/1 = 2 -- PASS
- SO-Sp node 0, N_bif=1, rank(1,1): x = 1*2*1/1 = 2 -- PASS
- Sp-Sp node 0, N_bif=2, rank(1,1): x = 2*1*1/1 = 2 -- PASS

**Status:** ALL PASS | **Confidence:** INDEPENDENTLY CONFIRMED

### 5.3 Limiting Cases (SQCD Limit)

The SQCD limit requires that the general two-node formula reduces to the known single-node conformal window for each gauge type.

**SU(N_c) SQCD:** b_0 = 3N_c - N_f. Window: 3N_c/2 < N_f < 3N_c => 3/2 < N_f/N_c < 3.
For SU-SU quiver at large N: b_0 = N(3 - N_bif/2), so x = N_bif/2.
- N_bif=3: x = 3/2 (lower boundary) -- matches Seiberg threshold
- N_bif=4: x = 2 (inside) -- matches
- N_bif=6: x = 3 (upper boundary, b_0=0) -- matches AF boundary

**SO(N_c) SQCD:** b_0 = 3(N_c-2) - N_f ~ 3N_c - N_f at large N. Same universal x.

**Sp(N_c) SQCD:** b_0 = 3(N_c+1) - N_f ~ 3N_c - N_f at large N. Same universal x.

All three gauge types reduce to the same 3/2 < x < 3 window at large N, as expected from the universal form b_0 ~ m_a*N*(3 - x).

**Status:** LIMITS_VERIFIED | **Confidence:** INDEPENDENTLY CONFIRMED

### 5.4 Cross-Check: T_bifund_lead Source Consistency

All 6 canonical T_bifund_lead values in conformal_window.py were cross-checked against the source of truth in a_maximization_large_N.py (lines 116-137). All match exactly:

| gauge_pair | conformal_window.py | a_maximization_large_N.py | Match |
|---|---|---|---|
| SU-SU | (1/2, 1/2) | (1/2, 1/2) | PASS |
| SU-SO | (1/2, 1) | (1/2, 1) | PASS |
| SU-Sp | (1, 1/2) | (1, 1/2) | PASS |
| SO-SO | (1, 1) | (1, 1) | PASS |
| SO-Sp | (2, 1/2) | (2, 1/2) | PASS |
| Sp-Sp | (1, 1) | (1, 1) | PASS |

Additionally, T_bifund_lead was verified from first principles:
T_lead_a = dim_lead_coeff(fund_b) * T(fund_a) where dim_lead_coeff(SU)=1, dim_lead_coeff(SO)=1, dim_lead_coeff(Sp)=2 and T(fund_SU)=1/2, T(V_SO)=1, T(f_Sp)=1/2.

**Status:** CROSS_CHECK_PASS | **Confidence:** INDEPENDENTLY CONFIRMED

### 5.6 Symmetry Check: SU-SU with Equal Ranks

For SU-SU quivers with rank_0 = rank_1, the two nodes are physically equivalent, so x_0 must equal x_1. The script verified 12 such morphologies; independent DB query confirms all satisfy x_0 = x_1.

**Status:** VERIFIED | **Confidence:** INDEPENDENTLY CONFIRMED

### 5.7 Conservation / Consistency: AF Constraint

x < 3 is required for asymptotic freedom (b_0 > 0). All 237 entries satisfy x <= 3, with 10 at the marginal boundary x=3 and 0 above. This is consistent since the DB only contains AF theories.

**Status:** VERIFIED | **Confidence:** INDEPENDENTLY CONFIRMED

### 5.8 Mathematical Consistency

Classification logic independently verified:
- x < 3/2 -> "below" (confining)
- x = 3/2 -> "at_lower_boundary" (free magnetic phase boundary)
- 3/2 < x < 3 -> "inside" (conformal window)
- x = 3 -> "at_upper_boundary" (AF boundary)
- x > 3 -> "above" (not AF)

All 7 distinct x values (1/4, 1/2, 1, 3/2, 2, 5/2, 3) correctly classified.

Total entries: 237 = 130 (below) + 18 (at_lower_boundary) + 79 (inside) + 10 (at_upper_boundary) + 0 (above). Sum = 237. CONSISTENT.

**Status:** CONSISTENT | **Confidence:** INDEPENDENTLY CONFIRMED

### 5.10 Agreement with Literature

The conformal window criterion 3/2 < x < 3 correctly implements:
- Seiberg (hep-th/9411149): SU(N_c) SQCD with 3N_c/2 < N_f < 3N_c
- Intriligator-Seiberg (hep-th/9503179): SO(N_c) SQCD analog
- Intriligator-Pouliot (hep-th/9505006): Sp(N_c) SQCD analog

At large N, all three reduce to the universal form 3/2 < x < 3, which is what the script implements.

**Status:** AGREES | **Confidence:** STRUCTURALLY PRESENT (full finite-N comparison not attempted -- large-N limit only)

## Forbidden Proxy Audit

| Proxy ID | Status | Evidence |
|---|---|---|
| fp-tables-without-physics | REJECTED | Markdown includes Physics Interpretation section with: (1) explicit statement that classification is necessary but not sufficient, (2) gauge pair comparison with fraction inside window, (3) rank multiplier effect discussion, (4) SQCD literature connection |

## Physics Consistency Summary

| Check | Status | Confidence | Notes |
|---|---|---|---|
| 5.1 Dimensional analysis | CONSISTENT | INDEPENDENTLY CONFIRMED | All quantities dimensionless |
| 5.2 Numerical spot-check | PASS (11/11) | INDEPENDENTLY CONFIRMED | 8 primary + 3 additional test points |
| 5.3 Limiting cases (SQCD) | LIMITS_VERIFIED | INDEPENDENTLY CONFIRMED | All 3 gauge types reduce to 3/2 < x < 3 at large N |
| 5.4 Cross-check (T_bifund_lead) | PASS | INDEPENDENTLY CONFIRMED | 6/6 values match source module |
| 5.6 Symmetry (SU-SU equal ranks) | VERIFIED | INDEPENDENTLY CONFIRMED | 12 morphologies verified x_0 = x_1 |
| 5.7 Conservation (AF consistency) | VERIFIED | INDEPENDENTLY CONFIRMED | 0 entries with x > 3 |
| 5.8 Math consistency | CONSISTENT | INDEPENDENTLY CONFIRMED | Classification logic correct, totals sum correctly |
| 5.10 Literature agreement | AGREES | STRUCTURALLY PRESENT | Large-N SQCD window matches Seiberg et al. |

**Overall physics assessment:** SOUND -- All checks pass, 7/8 independently confirmed.

## Confidence Assessment

**HIGH confidence.** The conformal window classification is a straightforward application of a well-established physics criterion (SQCD conformal window) to the enumerated two-node quiver theories. The x formula was verified from first principles using Dynkin index decomposition, cross-checked against the source module (a_maximization_large_N.py), and spot-checked at 11 test parameter sets. The SQCD limiting cases were independently re-derived. All classification logic is exact (using Python Fraction arithmetic), eliminating floating-point concerns. The script runs cleanly and produces self-consistent output.

The one check at STRUCTURALLY PRESENT (literature agreement) is because the Seiberg/Intriligator-Seiberg/Intriligator-Pouliot conformal windows were verified at large N only, not at finite N. This is appropriate since the entire project works at large N leading order.
