# Research Gaps and Open Issues

**Analysis Date:** 2026-04-14

## Incomplete Derivations

**Module 4: IR Fixed Point Boundary Analysis (not implemented):**
- What exists: Complete theoretical writeup in `module4_ir_fixed_points.md` describing the boundary analysis method of hep-th/0502049. Pseudocode algorithm in Sec. 8 (lines 142-167). All required ingredients (`a_maximize_large_N`, `build_fields_large_N`, `_anomaly_matrix_exact`) exist in `a_maximization_large_N.py`.
- What's missing: No code implements the decoupled sub-quiver construction or the B_a evaluator. The full classification pipeline (Module 5) has not applied the B_a < 0 filter. The database `quivers.db` contains 7537 theories in 326 classes that have passed Modules 1-3 but NOT Module 4.
- Files: `module4_ir_fixed_points.md` (Sec. 11, line 216), `module5_classification.md` (Sec. 6, line 200)
- Impact: **The current classification is incomplete.** Some fraction of the 326 universality classes may not flow to non-trivial IR SCFTs. Without the B_a < 0 filter, the database includes theories that could be free or confined in the IR.
- Difficulty estimate: Moderate. The required functions exist; the main work is (1) building the decoupled sub-quiver Q_a for each node a, (2) running a-maximization on Q_a with bifundamentals treated as free flavors (R = 2/3), (3) computing B_a = Tr[R^(a) G_a^2].

**Single-node boundary analysis subtlety (unresolved in notes):**
- What exists: `module4_ir_fixed_points.md` Sec. 9 (lines 175-194) works through the SQCD sanity check and discovers that the boundary analysis gives B > 0 for the entire AF window, contradicting Seiberg's conformal window result.
- What's missing: The notes acknowledge the resolution ("the boundary analysis is most useful for multi-node quivers") but do not provide a corrected single-node criterion. The text contains "Wait --" and "Actually for SQCD:" markers (lines 96, 190) indicating unfinished reasoning.
- Files: `module4_ir_fixed_points.md` (lines 96-98, 186-194)
- Impact: Medium. For multi-node quivers (the project's focus), the boundary analysis is valid. But the single-node case reveals that the method has limitations that should be documented precisely.
- Difficulty estimate: Straightforward -- the resolution is well-known (Banks-Zaks window vs boundary analysis applicability).

**Incomplete beta function counting for N=2 necklace quivers:**
- What exists: `module5_classification.md` line 109: "Each node: b_0 = 3N - N - 2*N/2 = N > 0... wait, this needs careful counting."
- What's missing: The computation trails off. The N=2 necklace quiver b_0 was not completed.
- Files: `module5_classification.md` (line 109)
- Impact: Low. This is a specific example in the notes, not part of the computational pipeline.
- Difficulty estimate: Straightforward.

## Unchecked Limits

**Large N convergence (commented out, not systematically run):**
- Limit: N -> infinity for fixed theory content; numerical a/N^2 at finite N should converge to the exact large-N value
- Expected behavior: a(N)/N^2 -> a_exact/N^2 as N -> infinity, with O(1/N) corrections
- Current status: The function `_convergence_check()` exists in `a_maximization_large_N.py` (line 1224) but is **entirely commented out** (lines 1230-1287). No convergence checks have been executed.
- Files: `a_maximization_large_N.py` (lines 1224-1287)
- Why it matters: The entire classification rests on large-N leading-order results. Without convergence checks, there is no empirical verification that subleading corrections are small.

**a-theorem consistency (a_UV >= a_IR not checked):**
- Limit: Compare UV and IR central charges
- Expected behavior: a_UV > a_IR for all theories in the classification (Komargodski-Schwimmer a-theorem)
- Current status: Not checked. `module5_classification.md` Sec. 10 (line 282) lists this as a cross-check but it is not implemented.
- Files: `module5_classification.md` (line 282)
- Why it matters: A violation would indicate an error in the a-maximization or an incorrect theory identification.

**N=2 supersymmetry limit (not checked):**
- Limit: Setting adjoint R-charge to R = 2/3 should recover N=2 results
- Expected behavior: Known N=2 central charges
- Current status: Listed as a cross-check in `module5_classification.md` (line 279) but not implemented.
- Files: `module5_classification.md` (line 279)
- Why it matters: Provides a non-trivial consistency check on the a-maximization machinery.

**Seiberg duality cross-check (not performed):**
- Limit: Seiberg duals of classified theories should also appear in the classification
- Expected behavior: Dual pairs share the same a*, c* values
- Current status: Not checked. Listed in `module5_classification.md` (line 280).
- Files: `module5_classification.md` (line 280)
- Why it matters: Would validate both the enumeration completeness and the a-maximization results.

## Unjustified Approximations

**Leading-order large N truncation:**
- Where used: `a_maximization_large_N.py` (entire file), `two_node_db.py` (build pipeline)
- Justification status: Numerical evidence absent. The `_convergence_check` function that would provide evidence is commented out (lines 1230-1287 of `a_maximization_large_N.py`).
- What could go wrong: For theories near the boundary of the conformal window, subleading O(1/N) corrections to R-charges and central charges could change the qualitative physics (e.g., flip the sign of B_a, move R[O] across the unitarity bound 2/3).
- How to justify: Uncomment and run `_convergence_check` for representative theories. Compare a/N^2 at N=50, 100, 200 against the exact large-N value.
- Priority: High

**No superpotential assumed (W = 0):**
- Where used: All modules. `plan.md` (line 151) explicitly flags this as an open question: "the current framework has no superpotential."
- Justification status: Hand-waving. The classification assumes no superpotential, which is the generic case for theories without cubic gauge-invariant operators. However, for theories with adjoints + fundamentals, the superpotential W = Tr(Phi Q Qtilde) is marginal and generically present.
- What could go wrong: Theories that should have W != 0 (e.g., adjoint + fundamentals) would have incorrect R-charges because the constraint R[Phi] + R[Q] + R[Qtilde] = 2 is not imposed. This changes the number of free parameters in a-maximization and can shift the IR fixed point.
- How to justify: Identify which theories in the database have marginal superpotential terms. For those, redo a-maximization with the superpotential constraint. Compare results.
- Priority: High -- affects correctness of R-charges for a subset of theories.

**Equal-rank assumption (all nodes share same N, or fixed rank multipliers [1,2]):**
- Where used: `quiver_generation.py` (Quiver dataclass), `two_node_db.py` (build, lines 359-373)
- Justification status: Stated in `CLAUDE.md` ("All nodes share the same parameter N"). The mixed-rank extension to [2,1] and [1,2] multipliers was added (commit `6aeb4f5`).
- What could go wrong: Many interesting quiver gauge theories have different ranks at different nodes (e.g., "quiver tails" with decreasing ranks). These are excluded from the classification.
- How to justify: This is a scope limitation, not an approximation. Document it clearly. `module5_classification.md` Sec. 8 (line 223) lists "Different ranks" as a generalization direction.
- Priority: Low (scope, not correctness)

## Missing Cross-Checks

**Unitarity bound on gauge-invariant operators (incomplete in large-N solver):**
- What to verify: The large-N a-maximization checks gauge-invariant operators via `_gauge_invariant_ops_symbolic()` (`a_maximization_large_N.py`, line 342), but the database schema has no column for unitarity_ok. The build pipeline does not filter on unitarity -- it was previously added then removed (commit `c019005`: "Remove wrong R>=2/3 filter").
- Method: After computing R-charges, check R[O] >= 2/3 for all gauge-invariant operators. If violated, apply iterative decoupling (`a_maximize_large_N_with_decoupling` exists at line 461).
- Expected outcome: Some theories may require decoupled R-charges; their a/N^2 values would change.
- Files to modify: `two_node_db.py` (build pipeline), `a_maximization_large_N.py`
- Priority: High -- the commit history shows this was attempted and then reverted, suggesting the implementation had issues.

**Exact vs numerical a/N^2 discrepancy for 145 classes:**
- What to verify: 145 of 326 universality classes have only numerical (Mathematica NSolve) results; the exact sympy solver timed out at 30s. The numerical values have not been independently verified.
- Method: Increase sympy timeout or use a different symbolic solver. Alternatively, cross-check Mathematica NSolve results against scipy at high N.
- Expected outcome: Confirm numerical a/N^2 values are correct to within tolerance.
- Files to modify: `two_node_db.py` (Phase 3, lines 478-519)
- Priority: Medium

**220 diverged theories (class_id NULL):**
- What to verify: 220 theories returned `None` from Mathematica NSolve (optimizer diverged). These may be theories with no bounded a-maximum, or they may have numerical issues.
- Method: Investigate representative diverged theories: check if the a-function is unbounded (no maximum exists) or if the Hessian check is too strict.
- Expected outcome: Distinguish genuine "no SCFT" cases from numerical failures.
- Files to modify: `two_node_db.py` (lines 388-424, 789-798)
- Priority: Medium

## Numerical Concerns

**Clustering tolerance for universality classes:**
- Problem: Universality classes are identified by clustering a/N^2 values within tolerance TOL = 1e-5 (`two_node_db.py`, line 38). If two genuinely distinct classes have a/N^2 values closer than 1e-5, they will be erroneously merged. Conversely, numerical noise could split a single class.
- Files: `two_node_db.py` (line 38, `_cluster()` at lines 214-235)
- Symptoms: Classes with suspiciously many or few theories. Exact symbolic results for a class disagree with the numerical centroid.
- Resolution: Cross-check by comparing exact symbolic a/N^2 values across all theories in each class (where available). Consider tighter tolerance or a second-pass verification.

**Hessian eigenvalue tolerance:**
- Problem: The Mathematica batch solver (`a_maximization_large_N.py`, line 857-858) uses `NegativeDefiniteMatrixQ` to check the Hessian. The commit history shows this was tightened (commit `c019005`: "tighten Hessian tolerance to eigenvalue <= 0"), and a previous R >= 2/3 filter was removed as "wrong". The Hessian check is purely numerical at WorkingPrecision -> 30, so near-degenerate eigenvalues could be misclassified.
- Files: `a_maximization_large_N.py` (lines 856-858)
- Symptoms: Theories classified as having no maximum that actually do (false negatives contributing to the 220 diverged theories), or saddle points accepted as maxima (false positives).
- Resolution: Log the smallest Hessian eigenvalue for each theory. Flag theories where |lambda_min| < 1e-10 as borderline.

**sympy `simplify` performance:**
- Problem: The exact solver (`a_maximize_large_N`, line 430) calls `simplify()` on R-charges and central charges. For theories with 4+ free parameters, this can take > 30s and trigger the timeout in `_exact_with_timeout` (`two_node_db.py`, line 181-195). This accounts for the 144 exact-solve failures.
- Files: `a_maximization_large_N.py` (lines 642-653), `two_node_db.py` (lines 181-195)
- Symptoms: 144/326 classes lack exact symbolic results.
- Resolution: Replace `simplify` with lighter canonicalization (e.g., `nsimplify` or `cancel`). Increase timeout. Use Mathematica's symbolic solver for the hard cases.

**Signal-based timeout (SIGALRM) is process-wide:**
- Problem: `_exact_with_timeout` in `two_node_db.py` (line 181-195) uses `signal.SIGALRM` for timeouts. This is not thread-safe and will fail in multi-threaded environments. The broad `except (_Timeout, Exception)` at line 191 also swallows all errors silently.
- Files: `two_node_db.py` (lines 181-195)
- Symptoms: Silent failures during exact solve phase. Any exception (including bugs) is caught and treated as a timeout.
- Resolution: Use `multiprocessing` with a timeout instead of SIGALRM. Narrow the exception handler.

## Physical Consistency Issues

**Boundary analysis applicability for single-node sub-quivers:**
- Concern: When implementing Module 4, the decoupled sub-quiver Q_a (with node a removed) may itself be a single-node theory. As documented in `module4_ir_fixed_points.md` Sec. 9 (lines 175-194), the boundary analysis gives misleading results for single-node SQCD. The notes acknowledge this but do not provide a corrected prescription.
- Files: `module4_ir_fixed_points.md` (lines 170-194)
- Impact: The Module 4 implementation for two-node quivers requires a-maximizing a single-node sub-quiver. If that sub-quiver is outside its conformal window, the R-charges at the boundary are simply R = 2/3 (free values). The notes state this (line 171) but the edge cases need careful handling.
- Resolution path: For the two-node case, the decoupled sub-quiver is always a single-node theory. Use the known conformal window boundaries from the literature (arXiv:2007.16165, arXiv:2510.19136 referenced in `module5_classification.md` Sec. 2a) to determine whether the sub-quiver has a non-trivial fixed point. If not, set all fields to R = 2/3.

**SU(N) anomaly cancellation requires N_f fundamentals that are not enumerated:**
- Concern: For SU(N) nodes with nonzero chiral excess delta(N) = aN + b, anomaly cancellation requires O(N) extra fundamentals/antifundamentals. These are included in the large-N analysis (`build_fields_large_N`, line 209-235 of `a_maximization_large_N.py`) via the `chiral_excess_coeffs` function. However, only the O(N) piece contributes at leading order; the O(1) piece (b coefficient) is dropped.
- Files: `a_maximization_large_N.py` (lines 209-235), `quiver_generation.py` (lines 157-194)
- Impact: For theories where b != 0, the actual number of fundamentals at finite N differs from the large-N approximation. This is a controlled approximation (subleading in 1/N) but has never been quantified.
- Resolution path: Run the convergence check `_convergence_check` to verify.

**Superpotential-sensitive theories misclassified:**
- Concern: Theories with an adjoint chiral multiplet and fundamental matter (e.g., SU(N) + 1 adj + N_f flavors) generically have a marginal superpotential W = Tr(Phi Q Qtilde). This superpotential is not included in the a-maximization, which means R-charges for such theories may be incorrect. The plan.md (line 151) flags this explicitly.
- Files: `plan.md` (line 151), `a_maximization.py`, `a_maximization_large_N.py`
- Impact: Affects all theories in the database that have both adjoint matter and fundamental/bifundamental matter at the same node. The number of such theories is significant (any class containing adj + bifundamental edges).
- Resolution path: Identify affected theories. For each, add the superpotential constraint R[adj] + R[Q] + R[Qtilde] = 2 as an additional linear constraint in a-maximization.

## Missing Generalizations

**Three-node and higher quivers (not enumerated):**
- Current scope: Complete classification for 1-node and 2-node quivers only. `module5_classification.md` Sec. 6 (line 202): "k >= 3-node theories (not yet enumerated)."
- Natural extension: Run `enumerate_quivers(n_nodes=3)` and the full pipeline. The infrastructure supports it.
- Difficulty: Moderate (combinatorial explosion; `quiver_generation.py` `enumerate_quivers` works for arbitrary n_nodes but the number of candidates grows rapidly).
- Blocks: Complete classification of quiver gauge theories.

**Superpotential constraints:**
- Current scope: W = 0 assumed throughout.
- Natural extension: Add superpotential terms as linear constraints on R-charges. `module3_a_maximization.md` Sec. 4 (lines 53-56) describes the formalism.
- Difficulty: Moderate. Requires identifying which superpotential terms are marginal for each theory.
- Blocks: Correct R-charges for theories with marginal superpotentials.

## Stale or Dead Content

**Commented-out validation examples in a_maximization_large_N.py:**
- What: Five complete validation examples (lines 1230-1287) are entirely commented out. These include convergence checks for SU(N)+adj, SU(N)^2 circular, SU(N)+S+Sbar, SU(N)+S (chiral), and SU(N)+adj+S+Sbar. Also `_classify_single_node()` call at line 1291 is commented out.
- Files: `a_maximization_large_N.py` (lines 1230-1292)
- Action: Uncomment and run at least once to validate. Then move to a proper test suite.
- Risk: Without these checks, there is no regression testing for the exact large-N solver.

**"Wait" markers in module notes:**
- What: `module4_ir_fixed_points.md` contains "Wait --" (line 96) and `module5_classification.md` contains "wait, this needs careful counting" (line 109). These indicate unfinished reasoning left in the notes.
- Files: `module4_ir_fixed_points.md` (line 96), `module5_classification.md` (line 109)
- Action: Resolve the reasoning and clean up the notes.
- Risk: Low. These are documentation issues, not code bugs.

**Stale build_info in quivers.db:**
- What: The database build_info says n_exact=182 and n_exact_fail=144, but the actual universality_class table has 181 exact and 145 numeric-only. Off-by-one between build_info metadata and actual table contents.
- Files: `quivers.db` (build_info table vs universality_class table)
- Action: Investigate the discrepancy. Likely a counting bug in the build pipeline.
- Risk: Low (metadata only, does not affect results).

## Placeholder and Stub Content

**Module 4 implementation (code not started):**
- What: `module4_ir_fixed_points.md` Sec. 11 (line 216) states "Module 4 is not yet implemented in code."
- Files: No ir_fixed_points.py or equivalent exists.
- Needed for: Completing the classification pipeline (filtering theories that do not flow to non-trivial IR SCFTs).

**Three-node enumeration (infrastructure exists but not run):**
- What: `enumerate_quivers(n_nodes=3)` can be called but has not been run through the full pipeline.
- Files: `quiver_generation.py` (line 306), `module5_classification.md` (line 202)
- Needed for: Extending the classification beyond two-node quivers.

## Priority Ranking

**Critical (blocks correctness):**
1. Module 4 boundary analysis not implemented: 326 universality classes have not been filtered for non-trivial IR fixed points. The current database may contain theories that are free or confined in the IR.
2. Superpotential constraints missing: Theories with adjoint + fundamental matter have incorrect R-charges because W = Tr(Phi Q Qtilde) constraint is not imposed.

**High (blocks completeness):**
1. Unitarity decoupling not applied in database build: The R >= 2/3 filter was added then removed (commit history c019005, 7065a5c). Some theories may need decoupled R-charges.
2. Large-N convergence checks commented out: No empirical validation that subleading corrections are small.
3. 220 diverged theories uninvestigated: Unknown whether these are genuine "no SCFT" cases or numerical failures.

**Medium (improves quality):**
1. 145/326 classes lack exact symbolic results (sympy timeout at 30s).
2. Hessian eigenvalue tolerance not logged: near-degenerate cases may be misclassified.
3. Seiberg duality and a-theorem cross-checks not performed.
4. build_info metadata discrepancy (off-by-one in exact count).

**Low (nice to have):**
1. Three-node and higher quiver enumeration.
2. Clean up "Wait" markers in module notes.
3. Move commented-out validation examples to a test suite.
4. Thread-safe timeout mechanism in two_node_db.py.

---

*Gap analysis: 2026-04-14*
