---
phase: 03-conformal-window-analysis
plan: 01
depth: full
one-liner: "Classified 237 N_rank2=0 nodes by SQCD conformal window: 79 inside (33%), 130 below (55%), 28 at boundaries, with T_bifund_lead independently verified from Dynkin indices"
subsystem: analysis
tags: [conformal-window, SQCD, a-maximization, large-N, quiver]

requires:
  - phase: 02-two-node-summary
    provides: "morphology_class table with 471 morphologies, T_bifund_lead values"
provides:
  - "Conformal window classification for all 237 N_rank2=0 nodes across 199 morphologies"
  - "JSON classification data (results/conformal_window_classification.json)"
  - "Markdown table with distribution analysis (results/tables/conformal_window.md)"
  - "T_bifund_lead independently verified from first-principles Dynkin index computation"
  - "SQCD limit verified for SU, SO, Sp gauge types at large N"
affects: [paper-writing, distribution-analysis]

methods:
  added: [conformal-window-criterion, per-node-SQCD-analogy]
  patterns: [exact-rational-arithmetic, per-node-classification]

key-files:
  created:
    - results/conformal_window.py
    - results/conformal_window_classification.json
    - results/tables/conformal_window.md

key-decisions:
  - "Strict inequality 3/2 < x < 3 for 'inside' classification; boundary values reported separately"
  - "Per-node classification (not per-morphology): each N_rank2=0 node gets independent x value"
  - "Exact rational arithmetic (fractions.Fraction) for all x computations -- no floating-point in classification logic"

conventions:
  - "T(fund_SU) = 1/2, T(V_SO) = 1, T(f_Sp) = 1/2"
  - "gauge_pair ordering: node 0 = first type, node 1 = second type (e.g., SU-SO: node 0 is SU)"
  - "x_a = N_bif * T_lead_a * m_other / m_self (general formula with rank multipliers)"

plan_contract_ref: ".gpd/phases/03-conformal-window-analysis/03-01-PLAN.md#/contract"
contract_results:
  claims:
    claim-conformal-window:
      status: passed
      summary: "All 237 N_rank2=0 nodes classified by conformal window criterion using general formula x = N_bif * T_lead * m_other/m_self. Distribution analysis by gauge pair shows SO-SO has highest inside fraction (56%), SU-Sp lowest (23%)."
      linked_ids: [deliv-conformal-window-json, deliv-conformal-window-table, deliv-conformal-window-script, test-sqcd-limit, test-af-consistency, test-su-su-symmetry, ref-seiberg-sqcd, ref-phase2-data, ref-t-bifund-lead]
  deliverables:
    deliv-conformal-window-json:
      status: passed
      path: "results/conformal_window_classification.json"
      summary: "237 entries with morph_id, node_index, gauge_type, x (exact rational), classification. All required fields present."
      linked_ids: [claim-conformal-window, test-sqcd-limit, test-af-consistency, test-su-su-symmetry]
    deliv-conformal-window-table:
      status: passed
      path: "results/tables/conformal_window.md"
      summary: "Markdown with summary, distribution by gauge pair, rank sector comparison, most common x values, full classification table, and physics interpretation including necessary-condition caveat and SQCD literature connection."
      linked_ids: [claim-conformal-window, test-sqcd-limit, test-af-consistency]
    deliv-conformal-window-script:
      status: passed
      path: "results/conformal_window.py"
      summary: "Self-contained Python script with T_bifund_lead lookup, x computation, classification, 7 automated validation checks, JSON and markdown output generation."
      linked_ids: [claim-conformal-window]
  acceptance_tests:
    test-sqcd-limit:
      status: passed
      summary: "SU-SU N_bif=4 gives x=2 (inside), N_bif=3 gives x=3/2 (boundary), N_bif=6 gives x=3 (upper boundary). All match SQCD conformal window 3N_c/2 < N_f < 3N_c with x=N_f/N_c."
      linked_ids: [claim-conformal-window, deliv-conformal-window-json, ref-seiberg-sqcd]
    test-af-consistency:
      status: passed
      summary: "No x > 3 for any entry (0 'above' classifications). 10 entries at x=3 (upper boundary, marginal). All strictly AF nodes have x < 3."
      linked_ids: [claim-conformal-window, deliv-conformal-window-json, ref-phase2-data]
    test-su-su-symmetry:
      status: passed
      summary: "12 SU-SU morphologies with both nodes N_rank2=0 and equal rank multipliers all have x_0 = x_1."
      linked_ids: [claim-conformal-window, deliv-conformal-window-json]
  references:
    ref-seiberg-sqcd:
      status: completed
      completed_actions: [cite, compare]
      missing_actions: []
      summary: "SQCD conformal window 3/2 < N_f/N_c < 3 used as the classification criterion and verified via limiting case checks."
    ref-phase2-data:
      status: completed
      completed_actions: [read, use]
      missing_actions: []
      summary: "morphology_class table queried for all 199 morphologies with N_rank2=0 nodes. Count verified."
    ref-t-bifund-lead:
      status: completed
      completed_actions: [use]
      missing_actions: []
      summary: "T_bifund_lead values used directly from a_maximization_large_N.py and independently verified from Dynkin indices."
  forbidden_proxies:
    fp-tables-without-structure:
      status: rejected
      notes: "Output includes distribution analysis by gauge pair, rank sector comparison, physics interpretation with necessary-condition caveat, and SQCD literature connection -- not just raw tables."
    fp-wrong-x-formula:
      status: rejected
      notes: "Used general formula x = N_bif * T_lead * m_other/m_self with gauge-pair-dependent T_bifund_lead, not x = N_bif/2."
  uncertainty_markers:
    weakest_anchors:
      - "Per-node SQCD analogy ignores coupled RG effects -- classification is necessary but not sufficient for IR fixed point"
    unvalidated_assumptions:
      - "Inter-node coupling effects on conformal window boundaries not assessed (Module 4 deferred)"
    competing_explanations: []
    disconfirming_observations: []

comparison_verdicts:
  - subject_id: test-sqcd-limit
    subject_kind: acceptance_test
    subject_role: decisive
    reference_id: ref-seiberg-sqcd
    comparison_kind: benchmark
    metric: exact_match
    threshold: "exact"
    verdict: pass
    recommended_action: "None -- SQCD limit exactly reproduced"
    notes: "x = N_bif/2 for SU-SU with equal ranks matches Seiberg window with x = N_f/N_c"

duration: 8min
completed: 2026-04-15
---

# Phase 3, Plan 01: Conformal Window Classification Summary

**Classified 237 N_rank2=0 nodes by SQCD conformal window: 79 inside (33%), 130 below (55%), 28 at boundaries, with T_bifund_lead independently verified from Dynkin indices**

## Performance

- **Duration:** ~8 min
- **Started:** 2026-04-15
- **Completed:** 2026-04-15
- **Tasks:** 2
- **Files created:** 3

## Key Results

- **237 qualifying nodes** from 199 morphologies classified: 79 inside (33.3%), 130 below (54.9%), 18 at lower boundary (7.6%), 10 at upper boundary (4.2%), 0 above (AF consistency) [CONFIDENCE: HIGH]
- **SO-SO has highest inside fraction** (5/9 = 56%), **SU-Sp has lowest** (11/48 = 23%) [CONFIDENCE: HIGH]
- **Only 7 distinct x values appear:** 1/4, 1/2, 1, 3/2, 2, 5/2, 3 -- all rational with small denominators [CONFIDENCE: HIGH]
- **Both-node analysis:** of 38 morphologies with both nodes N_rank2=0, only 7 (18%) have both inside the window [CONFIDENCE: HIGH]
- **Rank multipliers matter:** rank != (1,1) nodes show 46/117 (39%) inside vs 33/120 (28%) for rank(1,1) [CONFIDENCE: HIGH]
- **T_bifund_lead verified** from first-principles Dynkin indices for all 6 gauge pairs [CONFIDENCE: HIGH]

## Task Commits

1. **Task 1: Compute conformal window classification** - `0ddc7bf` (compute)
2. **Task 2: Cross-validate T_bifund_lead and SQCD limits** - `cbcd6d8` (validate)

## Files Created

- `results/conformal_window.py` - Self-contained script: queries DB, computes x, classifies, validates, outputs JSON + markdown
- `results/conformal_window_classification.json` - 237 entries with exact rational x values and classification
- `results/tables/conformal_window.md` - Full markdown with distribution analysis and physics interpretation

## Next Phase Readiness

- Conformal window classification data available for paper writing (distribution plots, tables)
- JSON format enables programmatic analysis in downstream phases
- Physics interpretation section ready for incorporation into manuscript

## Contract Coverage

- **claim-conformal-window:** passed -- all 237 nodes classified with distribution analysis
- **deliv-conformal-window-json:** passed -- 237 entries, all required fields
- **deliv-conformal-window-table:** passed -- includes distribution analysis and physics interpretation
- **deliv-conformal-window-script:** passed -- self-contained with 7 validation checks
- **test-sqcd-limit:** passed -- exact match with Seiberg window
- **test-af-consistency:** passed -- no x > 3
- **test-su-su-symmetry:** passed -- 12 morphologies verified
- **ref-seiberg-sqcd:** completed (cite, compare)
- **ref-phase2-data:** completed (read, use)
- **ref-t-bifund-lead:** completed (use)
- **fp-tables-without-structure:** rejected (distribution analysis included)
- **fp-wrong-x-formula:** rejected (general formula used)

## Equations Derived

**Eq. (03.1): Conformal window parameter**

$$
x_a = N_{\mathrm{bif}} \cdot T_{\mathrm{lead},a}(G_a, G_b) \cdot \frac{m_b}{m_a}
$$

**Eq. (03.2): Beta function in terms of x**

$$
b_0^{(a)} = m_a N (3 - x_a)
$$

**Eq. (03.3): Conformal window criterion**

$$
\frac{3}{2} < x_a < 3
$$

## Validations Completed

- **AF consistency:** x < 3 for all 237 entries (0 above, 10 at boundary) -- PASS
- **SQCD limit:** SU-SU N_bif=4 gives x=2, N_bif=3 gives x=3/2, N_bif=6 gives x=3 -- PASS
- **SU-SU symmetry:** x_0 = x_1 for all 12 equal-rank morphologies with both nodes N_rank2=0 -- PASS
- **Morphology count:** 199 morphologies (matches research prediction) -- PASS
- **T_bifund_lead first principles:** dim_lead(fund_b) * T(fund_a) matches code for all 6 pairs -- PASS
- **SQCD large-N limit:** SU, SO, Sp all reduce to 3/2 < x < 3 -- PASS
- **Totals consistency:** 79 + 130 + 18 + 10 + 0 = 237 -- PASS

## Decisions & Deviations

None -- plan executed exactly as written.

## Open Questions

- Do the 7 morphologies with both nodes inside the conformal window form a physically interesting subset?
- How do conformal window classifications correlate with a/c ratios from Phase 2?
- Would Module 4 boundary analysis modify any of these classifications?

## Approximations Used

| Approximation | Valid When | Error Estimate | Breaks Down At |
|---|---|---|---|
| Large-N leading order | N >> 1 | O(1/N) corrections to b_0 | Finite N |
| Per-node SQCD analogy | Nodes approximately independent | Ignores coupled RG effects | Strong inter-node coupling |

## Self-Check: PASSED

- [x] results/conformal_window.py exists
- [x] results/conformal_window_classification.json exists (237 entries)
- [x] results/tables/conformal_window.md exists (with all required sections)
- [x] Commit 0ddc7bf exists
- [x] Commit cbcd6d8 exists
- [x] All 7 validation checks pass
- [x] Contract: every claim, deliverable, test, reference, and forbidden proxy ID has explicit status

---

_Phase: 03-conformal-window-analysis_
_Completed: 2026-04-15_
