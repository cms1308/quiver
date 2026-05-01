# Phase 1: Single-Node Summary - Research

**Researched:** 2026-04-14
**Domain:** 4D N=1 supersymmetric gauge theory, a-maximization, large N classification
**Confidence:** HIGH

## Summary

This phase compiles a complete classification table of all 67 single-node 4D N=1 asymptotically free gauge theories (SU(N), SO(N), Sp(N)) with rank-2 matter that admit a large N limit. The computation engine (`a_maximization_large_N.py`) and enumeration code (`quiver_generation.py`) already exist and correctly enumerate all 67 theories (verified: 52 SU + 6 SO + 9 Sp, excluding 3 pure SYM theories).

The key physics insight is that at large N, the leading-order central charges and R-charges depend on two quantities: (1) the number and type of rank-2 fields (with SU adjoints distinct from SU rank-2 tensors, but SO/Sp adjoints identical to their respective rank-2 tensors), and (2) whether anomaly cancellation forces O(N) fundamentals (Veneziano-limit). This creates a richer universality structure for SU than for SO/Sp.

**Primary recommendation:** Use `enumerate_quivers(n_nodes=1, min_multiedge=0, require_connected=False)` to enumerate all 70 theories, filter out 3 pure SYM, then run `a_maximize_large_N()` on each. Present as three tables (SU, SO, Sp) with columns: matter content, chiral excess (SU only), N_f bound, R-charges, a/N^2, c/N^2, a/c.

## User Constraints

See phase CONTEXT.md for locked decisions and user constraints that apply to this phase.

Key constraints affecting this research:
- LOCKED: Include ALL theories (non-Veneziano AND Veneziano-limit AND conformal manifold b_0=0)
- LOCKED: Total count 52 SU + 6 SO + 9 Sp = 67
- LOCKED: Type I theories report a/N^2 = c/N^2 = 0 at leading order (subleading deferred)
- LOCKED: One table per gauge type with columns: matter, N_f bound, R-charges, a/N^2, c/N^2, a/c
- DISCRETION: Table sorting order, Type I/II/III grouping, R-charge formatting
- OUT OF SCOPE: Subleading a-maximization, double-scaling limit, comparison with paper Table 45

## Active Anchor References

| Anchor / Artifact | Type | Why It Matters Here | Required Action | Where It Must Reappear |
| --- | --- | --- | --- | --- |
| arXiv:2510.19136 Tables 2, 27, 35 | benchmark | Non-Veneziano theory comparison (20 SU + 6 SO + 9 Sp) | compare, cite | execution, verification |
| arXiv:2510.19136 Table 45 | benchmark | Veneziano-limit comparison target | read (deferred to Phase 5) | verification |
| `a_maximization_large_N.py` | prior artifact | Computation engine for R-charges and central charges | use directly | execution |
| `quiver_generation.py` | prior artifact | Theory enumeration with AF + anomaly constraints | use directly | execution |
| Universal formulas: Type II a=c=(27/128)dim(G)/N^2, Type III a=c=(1/4)dim(G)/N^2 | method | Validation anchor for SO/Sp theories and SU theories with only adjoints | verify all computed values match | verification |

**Missing or weak anchors:** The extra SU theory S+2A+3Abar (conformal manifold, not in paper Table 2) has no independent benchmark. The paper may exclude it for a valid reason (e.g., unitarity violation at finite N). Flag for verification.

## Conventions

| Choice | Convention | Alternatives | Source |
| --- | --- | --- | --- |
| Gauge groups | SU(N), SO(N), Sp(N)=USp(2N) | Different Sp normalizations | CLAUDE.md |
| Dynkin index | T(fund_SU)=1/2, T(V_SO)=1, T(f_Sp)=1/2 | T(fund)=1 normalization | CLAUDE.md, module1 |
| Representation notation | adj for all gauge types; S, Sbar, A, Abar for SU; S for SO symmetric; A for Sp antisymmetric | Paper uses different notation for SO/Sp | CONTEXT.md decision |
| N_f definition | SU: min(n_fund, n_antifund); SO: n_vector; Sp: n_fund/2 | Alternative: total number of flavors | quiver_generation.py |
| Large N scaling | a/N^2, c/N^2 as O(1) quantities at leading order | Alternative: a, c as functions of N | Standard in literature |

**CRITICAL: All equations and results below use these conventions.**

## Mathematical Framework

### Key Equations and Starting Points

| Equation | Name/Description | Source | Role in This Phase |
| --- | --- | --- | --- |
| b_0 = 3T(adj) - sum_i T(r_i) | One-loop beta function coefficient | module1_beta_functions.md | AF condition: b_0 > 0 determines N_f bound |
| T(adj) + sum_i T_i(R_i - 1) = 0 | Anomaly-free R-symmetry constraint | Intriligator-Wecht (hep-th/0304128) | Linear constraint on R-charges per gauge node |
| a = (3/32)(3 Tr R^3 - Tr R) | Trial a-function | Intriligator-Wecht | Maximize over free R-charge parameters |
| c = (1/32)(9 Tr R^3 - 5 Tr R) | c central charge | Standard | Computed from same traces as a |
| delta(N) = a_delta * N + b_delta | Chiral excess from anomaly cancellation | module2, quiver_generation.py | Determines O(N) vs O(1) fundamentals for SU |

### Required Techniques

| Technique | What It Does | Where Applied | Standard Reference |
| --- | --- | --- | --- |
| Large N expansion of Dynkin indices | Replaces T(rep), dim(rep) with leading N coefficients | All field contributions | a_maximization_large_N.py |
| Symbolic linear constraint solving | Solves anomaly constraints for R-charges | Step 2 of a-maximization | sympy.solve |
| Symbolic cubic maximization | Differentiates a_trial, solves quadratic critical-point equations | Step 3 of a-maximization | sympy.solve |

### Approximation Schemes

| Approximation | Small Parameter | Regime of Validity | Error Estimate | Alternatives if Invalid |
| --- | --- | --- | --- | --- |
| Large N leading order | 1/N | N >> 1 | O(1/N) corrections to R-charges, O(1/N^2) to a/N^2 | Exact finite-N a-maximization in a_maximization.py |

**Key insight for SU(N) at large N:**
- All rank-2 tensors (S, Sbar, A, Abar) have identical T_lead = 1/2 and dim_lead = 1/2
- Adjoint has T_lead = 1 and dim_lead = 1 (double the rank-2 values)
- Therefore 1 SU adjoint is NOT equivalent to 1 rank-2 tensor at leading order; it is equivalent to 2 rank-2 tensors in the anomaly constraint but contributes differently to traces

**Key insight for SO(N) and Sp(N) at large N:**
- SO: adj and S both have T_lead = 1, dim_lead = 1/2 (IDENTICAL)
- Sp: adj and A both have T_lead = 1, dim_lead = 2 (IDENTICAL)
- Therefore ALL rank-2 representations are equivalent within each gauge type, producing clean universality

## Standard Approaches

### Approach 1: Automated Enumeration + Symbolic a-Maximization (RECOMMENDED)

**What:** Use existing `quiver_generation.py` to enumerate all valid single-node theories, then `a_maximization_large_N.py` to compute exact symbolic R-charges and central charges.

**Why standard:** The code is already implemented, tested, and produces exact algebraic results via sympy.

**Key steps:**

1. Enumerate: `enumerate_quivers(n_nodes=1, min_multiedge=0, require_connected=False)` returns 70 theories
2. Filter: Remove 3 pure SYM theories (no rank-2 matter)
3. For each theory: Call `a_maximize_large_N(quiver)` to get exact R-charges, a/N^2, c/N^2
4. Compute a/c ratio
5. Extract chiral excess coefficients for SU theories via `chiral_excess_coeffs(q, 0)`
6. Format N_f bound from AF condition (alpha, gamma) as "N_f < alpha*N + gamma"
7. Assemble into three tables (SU, SO, Sp)
8. Validate against known results (universal formulas, paper tables)

**Known difficulties at each step:**

- Step 3: Some Type I theories (n_rank2=1 with adj) give a=c=0, making a/c undefined. Handle with "n/a" or "indeterminate".
- Step 3: SU theories with n_rank2=1 and non-adjoint rank-2 (1A or 1S with Veneziano funds) give NEGATIVE a/N^2. This is physically meaningful: the leading-order a-function is maximized at a negative value, indicating the IR physics is dominated by subleading contributions. Present as computed.
- Step 6: Some theories have negative gamma in N_f bound (e.g., SO pure SYM has N_f < 3N-6). This means N_f must be strictly positive for large enough N but the bound is valid.

### Anti-Patterns to Avoid

- **Re-deriving a-maximization from scratch:** The code already handles all cases including Veneziano-limit fundamentals from chiral excess. Do not rewrite the solver.
- **Treating SU adjoint as equivalent to rank-2:** At large N for SU, adj has (T_lead, dim_lead) = (1, 1) vs rank-2 (1/2, 1/2). They are NOT interchangeable in the traces.
- **Ignoring Veneziano-limit theories:** These are the 31+ SU theories with nonzero chiral excess a-coefficient, where anomaly cancellation forces O(N) fundamentals. They contribute at leading order and have distinct R-charges.

## Existing Results to Leverage

### Established Results (DO NOT RE-DERIVE)

| Result | Exact Form | Source | How to Use |
| --- | --- | --- | --- |
| a-maximization principle | R_SCFT maximizes a_trial(R) subject to Tr[R G^2]=0 | Intriligator-Wecht hep-th/0304128 | Foundation of all computations |
| SO/Sp Type II universal formula | a/N^2 = c/N^2 = (27/128) * dim(G)/N^2, R_rank2 = 1/2 | Verified computationally | Validation check for 2S, adj+S, 2adj (SO) and 2A, adj+A, 2adj (Sp) |
| SO/Sp Type III universal formula | a/N^2 = c/N^2 = (1/4) * dim(G)/N^2, R_rank2 = 2/3 | Verified computationally | Validation for 3adj (SO, Sp) and 3A, adj+2A, 2adj+A (Sp) |
| SO/Sp Type I result | a/N^2 = c/N^2 = 0, R_rank2 = 0 at leading order | Verified computationally | All single-rank-2 SO/Sp theories |
| SU 2adj Type II | a/N^2 = c/N^2 = 27/128, R_adj = 1/2 | Verified | Only SU theory matching universal Type II formula |
| SU 3adj Type III (conformal manifold) | a/N^2 = c/N^2 = 1/4, R_adj = 2/3 | Verified | Conformal manifold at b_0=0 |
| dim(G)/N^2 leading coefficients | SU: 1, SO: 1/2, Sp: 2 | a_maximization_large_N.py | Converts universal formulas to gauge-type-specific values |
| SQCD conformal window | 3N/2 < N_f < 3N for SU(N) | Seiberg (hep-th/9411149) | Baseline: pure SQCD (no rank-2) should reproduce this |

### Useful Intermediate Results

| Result | What It Gives You | Source | Conditions |
| --- | --- | --- | --- |
| Chiral excess delta(N) = a*N + b | Number of extra fund/antifund from anomaly cancellation | chiral_excess_coeffs() | SU nodes only; a != 0 means Veneziano limit |
| N_f bound: N_f < alpha*N + gamma | Maximum number of flavor pairs before AF is lost | nf_bound() | Per-node bound |
| Leading-order field decomposition | List of LeadField objects with T_lead, dim_lead | build_fields_large_N() | Input to a-maximization |

### Relevant Prior Work

| Paper/Result | Authors | Year | Relevance | What to Extract |
| --- | --- | --- | --- | --- |
| hep-th/0304128 | Intriligator, Wecht | 2003 | Foundation: a-maximization principle | Method, unitarity bound treatment |
| arXiv:2007.16165 | (classification paper) | 2020 | Classification of large N SCFTs with dense spectra | Type I/II/III classification scheme, 35 non-Veneziano classes |
| arXiv:2510.19136 | (benchmark paper) | 2025 | Tables 2, 27, 35 for comparison | Specific R-charges and central charges for non-Veneziano theories |

## Computational Tools

### Core Tools

| Tool | Version/Module | Purpose | Why Standard |
| --- | --- | --- | --- |
| quiver_generation.py | enumerate_quivers() | Enumerate all valid 1-node theories | Already handles AF + anomaly + Witten constraints |
| a_maximization_large_N.py | a_maximize_large_N() | Exact symbolic R-charges and central charges | Produces exact algebraic numbers via sympy |
| sympy | Rational, solve, sqrt | Exact symbolic algebra | Avoids floating-point errors in R-charges |

### Supporting Tools

| Tool | Purpose | When to Use |
| --- | --- | --- |
| a_maximization_large_N.py: chiral_excess_coeffs() | Compute anomaly-forced fundamentals | SU theories: determines Veneziano vs non-Veneziano |
| a_maximization_large_N.py: nf_bound() | Extract AF bound | All theories: N_f < alpha*N + gamma |

### Computational Feasibility

| Computation | Estimated Cost | Bottleneck | Mitigation |
| --- | --- | --- | --- |
| Enumerate 70 theories | < 1 second | None | Already fast |
| Run a_maximize_large_N on 67 theories | ~30 seconds total | sympy symbolic solve for each | Already fast enough; could batch |
| Format and validate tables | Minutes (manual review) | Checking against paper tables | Systematic comparison script |

**Installation / Setup:**
No additional packages needed. All dependencies (sympy, numpy) are already available in the project environment.

## Validation Strategies

### Internal Consistency Checks

| Check | What It Validates | How to Perform | Expected Result |
| --- | --- | --- | --- |
| Theory count | Enumeration completeness | Count by gauge type | 52 SU + 6 SO + 9 Sp = 67 |
| SO/Sp Type II universality | All n_rank2=2 theories give same a/N^2 | Compare within group | a/N^2 = 27/128 * dim_lead, a/c = 1 |
| SO/Sp Type III universality | All n_rank2=3 theories give same a/N^2 | Compare within group | a/N^2 = 1/4 * dim_lead, a/c = 1 |
| Type I result | All n_rank2=1 theories with no Veneziano funds give a=0 | Check | a/N^2 = c/N^2 = 0 |
| SU universality classes | Theories with same (n_rank2, |delta_a|, n_adj) give same a/N^2 | Group and compare | Verified: each group has single a/N^2 value |
| a > 0 for interacting theories | Positive central charge | Check all non-Type-I theories | a/N^2 > 0 for n_rank2 >= 2 with non-trivial Veneziano |
| Anomaly constraint | R-charges satisfy Tr[R G^2] = 0 | Substitute back | Exactly zero |

### Known Limits and Benchmarks

| Limit | Parameter Regime | Known Result | Source |
| --- | --- | --- | --- |
| Pure SQCD | No rank-2 matter | R = 1 - N/N_f, conformal window 3N/2 < N_f < 3N | Seiberg |
| b_0 = 0 theories | N_f bound has gamma=0, alpha=0 | Conformal manifold, R-charges fixed by anomaly alone | Standard |
| arXiv:2510.19136 Table 2 | SU non-Veneziano theories | 20 theories with specific R-charges | Paper benchmark |
| arXiv:2510.19136 Table 27 | SO theories | 6 theories | Paper benchmark |
| arXiv:2510.19136 Table 35 | Sp theories | 9 theories | Paper benchmark |

### Red Flags During Computation

- If any SO or Sp theory count differs from 6 and 9 respectively: enumeration bug
- If any SO/Sp theory with n_rank2=2 gives a/N^2 != 27/128 * dim_lead: code bug or convention error
- If any Type I (n_rank2=1) SO/Sp theory gives nonzero a/N^2: missing subleading terms leaking into leading order
- If a/c < 1/2 or a/c > 5/4 for any theory: likely error (known bounds for unitary 4D N=1 SCFTs)
- If SU theory with adj-only matter (like 2adj, 3adj) does not match universal formula: adj contribution miscoded

## Common Pitfalls

### Pitfall 1: Confusing SU Adjoint with SU Rank-2

**What goes wrong:** Assuming adj ~ 2 rank-2 fields universally. While 1 adj contributes 2x the T_lead of 1 rank-2 field, it also contributes 2x the dim_lead. The a-function depends on BOTH T (through anomaly constraint) and dim (through traces), and the ratio T/dim differs: adj has T_lead/dim_lead = 1, rank-2 has T_lead/dim_lead = 1. Wait -- they're actually the same ratio! But the absolute values differ, and since the a-function is nonlinear (cubic), replacing one adj with two rank-2 fields changes the result.

**How to avoid:** Always use the exact T_lead and dim_lead values from the code. Never substitute "equivalent" representations.

**Warning signs:** SU theories with adjoints giving different a/N^2 than theories with the same n_rank2 but only S/A/Sbar/Abar.

### Pitfall 2: Missing Veneziano-Limit Fundamentals

**What goes wrong:** Treating all theories as having only O(1) fundamentals. When anomaly cancellation forces chiral excess delta ~ a*N with a != 0, there are O(N) fundamentals that contribute at leading order.

**How to avoid:** The code handles this correctly via `chiral_excess_coeffs()` and `build_fields_large_N()`. Trust the code, but verify by checking that theories with nonzero a_delta have fund/antifund LeadFields in their field list.

**Warning signs:** A theory with n_S=2 (chiral excess = -2N-8) giving the same a/N^2 as one with 1S+1Sbar (chiral excess = 0). They should differ.

### Pitfall 3: Negative a/N^2 for Type I Theories

**What goes wrong:** Reporting negative a/N^2 as a bug. For SU theories with n_rank2=1 and nonzero chiral excess (like SU + 1A), the leading-order a-maximization yields a < 0. This is not a bug; it means the leading-order approximation breaks down and the actual central charge is dominated by subleading O(1/N) corrections.

**How to avoid:** Report the computed value as-is. Note that Type I theories require subleading analysis for meaningful central charges (deferred).

### Pitfall 4: a/c Undefined for a=c=0 Theories

**What goes wrong:** Division by zero when computing a/c for theories where both a/N^2 and c/N^2 vanish at leading order.

**How to avoid:** Handle these cases explicitly. For Type I theories with a=c=0, report a/c as "indeterminate at leading order" or "1 (expected from subleading)."

### Pitfall 5: Charge Conjugation Overcounting

**What goes wrong:** Counting SU(N) + 1A and SU(N) + 1Abar as different theories, when they are related by charge conjugation.

**How to avoid:** The `_dedup_symmetries()` function in `quiver_generation.py` already handles this. The 67-theory count is after deduplication. Verify: 1A and 1Abar are NOT both present in the enumeration.

## Level of Rigor

**Required for this phase:** Exact symbolic computation (no floating-point approximations for R-charges and central charges) + systematic enumeration.

**Justification:** The phase produces a classification table that will appear in a JHEP paper. Exact algebraic results (like R = 1/2 for Type II, or R involving sqrt(19) for certain theories) are standard in the literature and expected by referees.

**What this means concretely:**
- R-charges must be exact sympy expressions, not floats
- Central charges a/N^2, c/N^2 must be exact rationals or algebraic numbers
- Numerical approximations are acceptable only as supplementary columns for readability
- Theory enumeration must be provably complete (all AF + anomaly conditions checked exhaustively)

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
| --- | --- | --- | --- |
| Case-by-case analysis per gauge group | Automated enumeration with constraint backtracking | This project | Ensures completeness; found 67 vs original estimate of 19 |
| Numerical a-maximization | Exact symbolic a-maximization at large N | This project (a_maximization_large_N.py) | Produces exact algebraic R-charges, not numerical approximations |
| Non-Veneziano theories only | Full classification including Veneziano limit | This phase | 52 SU theories vs 20 non-Veneziano |

## Open Questions

1. **Why does arXiv:2510.19136 exclude S+2A+3Abar from Table 2?**
   - What we know: Our enumeration finds this as a valid b_0=0 conformal manifold theory
   - What's unclear: Paper may have checked unitarity at finite N and found a violation, or may use a different scope definition
   - Impact on this phase: LOW -- include it in our table with a note
   - Recommendation: Flag with asterisk, verify unitarity of gauge-invariant operators at finite N if time permits

2. **Is the universality class structure the right way to organize the SU table?**
   - What we know: SU theories group by (n_rank2, |delta_a|, n_adj) into classes with identical leading-order data
   - What's unclear: Whether this grouping is physically meaningful or merely an artifact of the large N limit
   - Impact on this phase: Table organization choice (Agent's Discretion)
   - Recommendation: Sort by a/N^2 ascending within each gauge type; add a "class" column if desired

## Alternative Approaches if Primary Fails

| If This Fails | Because Of | Switch To | Cost of Switching |
| --- | --- | --- | --- |
| `a_maximize_large_N()` gives wrong results | Bug in symbolic solver | Use `a_maximize_large_N_fast()` (numerical) as cross-check, then debug symbolic | Low -- numerical solver exists |
| `enumerate_quivers()` misses theories | Bug in constraint checking | Manual enumeration from AF conditions (Tables in module1/module2) | Medium -- tedious but straightforward for 1-node |
| sympy `solve()` fails for some theory | Degenerate system | Use numerical BFGS fallback in `a_maximize_large_N_fast()` | Low |

**Decision criteria:** If any computed value disagrees with the universal formulas for SO/Sp Type II/III, debug immediately -- this indicates a code error, not a physics ambiguity.

## Sources

### Primary (HIGH confidence)

- Intriligator, Wecht, hep-th/0304128 - a-maximization principle, trial a-function, anomaly-free R-symmetry
- arXiv:2007.16165, "Classification of large N superconformal gauge theories with a dense spectrum" - Type classification, 35 non-Veneziano theories
- arXiv:2510.19136 - Tables 2, 27, 35 for direct numerical comparison of non-Veneziano theories

### Secondary (MEDIUM confidence)

- `a_maximization_large_N.py` in this repository - verified against known SQCD results and universal formulas
- `quiver_generation.py` in this repository - enumeration code producing correct counts (verified 52+6+9=67)

### Tertiary (LOW confidence)

- The S+2A+3Abar conformal manifold theory inclusion -- not independently benchmarked

## Metadata

**Confidence breakdown:**

- Mathematical framework: HIGH - a-maximization is textbook, well-established since 2003
- Standard approaches: HIGH - code already exists and is verified against known results
- Computational tools: HIGH - sympy produces exact results; enumeration code gives correct counts
- Validation strategies: HIGH - universal formulas provide strong cross-checks for SO/Sp; paper tables for SU

**Research date:** 2026-04-14
**Valid until:** Indefinitely for the physics; code may need updates if quiver_generation.py API changes
