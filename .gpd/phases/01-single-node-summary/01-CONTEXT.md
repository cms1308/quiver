# Phase 1: Single-node summary - Context

**Gathered:** 2026-04-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Compile the complete single-node classification table for ALL single-node 4D N=1 AF theories (SU(N), SO(N), Sp(N)) with rank-2 matter that admit a large N limit. Include both non-Veneziano (integer N_rank2, O(1) fundamentals) and Veneziano-limit theories (O(N) fundamentals from anomaly cancellation). Include conformal manifold theories at b_0=0.

Total: 67 theories (52 SU + 6 SO + 9 Sp), not 19 as initially estimated.

Requirements: ENUM-01

</domain>

<contract_coverage>
## Contract Coverage

- **Claim/deliverable:** Complete single-node classification table with gauge type, matter content, N_f bound, R-charges, a/N^2, c/N^2, a/c for all 67 theories
- **Acceptance signal:** Table covers all 67 theories; SO/Sp counts match arXiv:2510.19136 (6 SO, 9 Sp); SU count is 52 (paper has 20 non-Veneziano + unknown Veneziano in their Table 45)
- **False progress to reject:** Table with only non-Veneziano theories (missing 31 SU Veneziano-limit theories)

</contract_coverage>

<user_guidance>
## User Guidance To Preserve

- **User-stated observables:** Matter content, a/N^2, c/N^2, a/c, R-charges, N_f bound for every theory
- **User-stated deliverables:** Three tables (one per gauge type: SU, SO, Sp) with all columns above
- **Must-have references / prior outputs:** arXiv:2510.19136 Tables 2, 27, 35 (comparison target); existing `a_maximization_large_N.py` (computation engine)
- **Stop / rethink conditions:** Any SO or Sp theory count disagrees with paper (6 and 9 respectively); any Type II/III central charge disagrees with universal formulas a=c=(27/128)dim(G) or a=c=(1/4)dim(G)

</user_guidance>

<decisions>
## Methodological Decisions

### Scope of classification

- Include ALL single-node theories: non-Veneziano (O(1) fundamentals) AND Veneziano-limit (O(N) fundamentals from anomaly cancellation)
- Include theories at b_0=0 with conformal manifolds (e.g., 3adj, adj+S+S_bar+A+A_bar, S+2A+3A_bar+8Q_bar)
- Total count: 52 SU(N) + 6 SO(N) + 9 Sp(N) = 67

### Type I theories (N_rank2=1)

- At leading order in large N: a/N^2 = c/N^2 = 0, R_rank2 = 0
- Ignore subleading O(N)/O(1) corrections for now
- Present leading-order values as-is; subleading analysis deferred

### Notation

- Use our convention: `adj` for adjoint of all gauge types (not paper's `S` for Sp adjoint or `A` for SO adjoint)
- SO(N): adj = antisymmetric, S = symmetric
- Sp(N): adj = symmetric, A = antisymmetric
- SU(N): adj, S, S_bar, A, A_bar, fund, antifund

### Table format

- One table per gauge type (SU, SO, Sp)
- Columns: matter content, N_f bound (AF condition), R-charges (exact symbolic), a/N^2, c/N^2, a/c
- For SU: include chiral excess delta to show which theories need O(N) vs O(1) fundamentals

### Comparison with arXiv:2510.19136

- Paper classifies 20 SU + 6 SO + 9 Sp = 35 non-Veneziano theories
- We have one extra SU theory: S+2A+3A_bar (delta=-8, b_0=0 conformal manifold) not in their Table 2
- Our SO and Sp tables match exactly
- Paper's Veneziano-limit theories are in their Table 45 (not extracted yet)

### Agent's Discretion

- Table sorting order (by a/N^2, by N_rank2 type, or by matter content complexity)
- Whether to group by Type I/II/III following the paper's classification
- Formatting of exact symbolic R-charges (sympy expressions vs simplified radicals)

</decisions>

<assumptions>
## Physical Assumptions

- Large N leading order sufficient for classification: a/N^2, c/N^2, R-charges computed at O(1) in 1/N | Type I theories have a/N^2=0 at this order, so subleading needed for those eventually
- Equal-rank assumption (all nodes share same N): built into a_maximization_large_N.py | Multi-rank theories would be a different classification
- Anomaly-free R-symmetry uniquely determines IR fixed point via a-maximization: standard Intriligator-Wecht result | Breaks if accidental symmetries or unitarity violations at finite N
- Theories at b_0=0 included as conformal manifold theories: they don't flow but are valid SCFTs | Paper marks these with * and includes them

</assumptions>

<limiting_cases>
## Expected Limiting Behaviors

- N_rank2=2 (Type II): All theories must give a/N^2 = c/N^2 = (27/128)·(dim(G)/N^2) with R_rank2=1/2
- N_rank2=3 (Type III): All theories must give a/N^2 = c/N^2 = (1/4)·(dim(G)/N^2) with R_rank2=2/3
- N_rank2=1 (Type I): All theories must give R_rank2=0 and a/N^2=c/N^2=0 at leading order
- a/c = 1 at leading order for all Type II and III theories (universal result from paper)
- SO/Sp theory counts must be 6 and 9 (already verified)

</limiting_cases>

<anchor_registry>
## Active Anchor Registry

- **arXiv:2510.19136 Tables 2, 27, 35**
  - Why it matters: Primary comparison target for non-Veneziano theories
  - Carry forward: execution, verification, writing
  - Required action: compare, cite

- **arXiv:2510.19136 Table 45**
  - Why it matters: Contains Veneziano-limit theories to compare against our 31 SU Veneziano theories
  - Carry forward: execution, verification
  - Required action: read, compare

- **a_maximization_large_N.py**
  - Why it matters: Computation engine for all R-charges and central charges
  - Carry forward: execution
  - Required action: use (fix _dedup_symmetries bug for n_nodes=1)

- **Universal central charge formulas**
  - Why it matters: Type II: a=c=(27/128)dim(G), Type III: a=c=(1/4)dim(G)
  - Carry forward: verification
  - Required action: verify all computed values match

</anchor_registry>

<skeptical_review>
## Skeptical Review

- **Weakest anchor:** The extra SU theory S+2A+3A_bar — paper excludes it, possibly for good reason (unitarity violation at finite N?)
- **Unvalidated assumptions:** Veneziano-limit theories at large N leading order — R-charges may depend on N_f/N ratio, making a single "leading order" value less meaningful
- **Competing explanation:** Paper may exclude S+2A+3A_bar because it checked finite-N superconformal index and found unitarity violation
- **Disconfirming check:** If any Type II theory gives a/N^2 != 27/128 * dim(G)/N^2, something is wrong in the code
- **False progress to reject:** A table of only 19-20 theories (missing Veneziano-limit theories the user explicitly requested)

</skeptical_review>

<deferred>
## Deferred Ideas

- Subleading a-maximization for Type I theories (a ~ O(N), depends on N_f) — future phase or Phase 1 extension
- Double-scaling limit analysis (N >> N_f >> 1) for universal Type I R-charges — paper Section 2.2
- Comparison with paper's Table 45 (Veneziano-limit classification) — should be done during Phase 5 validation
- Type I/II/III classification labels — could adopt paper's terminology in the JHEP paper

</deferred>

---

_Phase: 01-single-node-summary_
_Context gathered: 2026-04-14_
