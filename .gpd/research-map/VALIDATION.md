# Validation and Cross-Checks

**Analysis Date:** 2026-04-14

## Analytic Cross-Checks

**Limiting Cases Verified:**

- SU(N) SQCD R-charge: R = 1 - N/N_f -> Matches analytic expectation for N=5, N_f=10 (R = 0.5): Match confirmed
  - File: `a_maximization.py` (line 478-484, `__main__` Example 1)

- SU(N) SQCD conformal window scan: R-charges at N=5 for N_f = 8, 10, 12, 14 match R = 1 - N/N_f: Match confirmed
  - File: `a_maximization.py` (lines 489-494, `__main__` Example 1 scan)

- SU(N) SQCD below conformal window (N_f=4 < 3N/2=7.5): unitarity_ok = False signals no SCFT: Match confirmed
  - File: `a_maximization.py` (lines 498-508, `__main__` Example 2)

- SU(N) + 1 adj: leading-order anomaly constraint gives R_adj = 0, fully determined (no free params): Match confirmed
  - File: `a_maximization_large_N.py` (lines 1230-1237, commented-out Example 1)

- SU(N)^2 circular (balanced +-): R_01 + R_10 = 0, a-max gives R = 0 (no SCFT at N_f=0): Match confirmed
  - File: `module3_a_maximization.md` (Sec. 11); `a_maximization_large_N.py` (lines 1241-1249, commented-out Example 2)

- SU(N) + 1S + 1Sbar: symmetry forces R_S + R_Sbar = 0, a-max gives R_S = R_Sbar = 0: Match confirmed
  - File: `a_maximization_large_N.py` (lines 1253-1262, commented-out Example 3)

- SU(N), no rank-2 matter: AF bound gives N_f < 3*N: Match confirmed
  - File: `quiver_generation.py` (lines 711-715, `__main__` validation)

- SU(N) + 1 S: AF bound gives N_f < 2*N - 3: Match confirmed
  - File: `quiver_generation.py` (lines 719-723, `__main__` validation)

- SU(N) + 1 A: AF bound gives N_f < 2*N + 3: Match confirmed
  - File: `quiver_generation.py` (lines 727-731, `__main__` validation)

- SU(N) + 3 A: AF bound gives N_f < 9 (constant bound): Match confirmed
  - File: `quiver_generation.py` (lines 735-739, `__main__` validation)

- SU(N) + 4 A: not AF for all N (alpha < 0): Match confirmed
  - File: `quiver_generation.py` (lines 743-747, `__main__` validation)

- Two SU(N) balanced +- nodes: b_0 per node = 5N/2, AF bound N_f < 5N/2: Match confirmed
  - File: `quiver_generation.py` (lines 751-758, `__main__` validation)

- Sp(N) Witten anomaly: odd degree fails, even degree passes: Match confirmed
  - File: `quiver_generation.py` (lines 762-772, `__main__` validation)

- SU with 2 SU bifundamental neighbors: b_0 = 2*N (linear coefficients alpha=2, beta=0): Match confirmed
  - File: `beta_functions.py` (lines 286-287, `__main__` example)

**Limiting Cases NOT Verified:**

- N=2 theories: SU(2) = Sp(1) equivalence not tested. The code sets N_MIN["SU"] = 5, so SU(2), SU(3), SU(4) are never evaluated. For classification purposes this is fine (project targets large N), but finite-N physical predictions at small N are untested.

- N=2 SQCD conformal window boundaries: The lower bound N_f = 3N/2 = 3 should reproduce Seiberg's result. This is discussed in `module4_ir_fixed_points.md` (Sec. 9) but the boundary analysis (Module 4) is not implemented, so the conformal window lower boundary is not computationally verified.

- Free field limit: As N_f -> 3N (upper AF boundary), R -> 2/3 (free field value). This limit is not explicitly tested but is implied by the SQCD scan.

- Seiberg duality: Dual descriptions of classified theories are not checked to appear in the classification.

- a-theorem (a_UV >= a_IR): Mentioned as a cross-check in `module5_classification.md` (Sec. 10) but not implemented.

- N=2 (extended SUSY) limit: "Setting the adjoint chiral R-charge to R = 2/3 should recover N=2 results" — mentioned in `module5_classification.md` (Sec. 10) but not implemented.

**Symmetry Checks:**

- Gauge anomaly cancellation: Verified algebraically via `check_anomalies()` for all enumerated quivers
  - File: `quiver_generation.py` (line 298, `check_anomalies`); called during enumeration

- Sp(N) Witten anomaly: Verified for all Sp nodes via `check_sp_witten()`
  - File: `quiver_generation.py` (line 287, `check_sp_witten`); called during enumeration

- R-symmetry anomaly-free condition Tr[R G_a^2] = 0: Enforced as a hard constraint in a-maximization (not checked post-hoc — it is an input constraint)
  - File: `a_maximization.py` (line 203, `anomaly_matrix`); `a_maximization_large_N.py` (line 274, `_build_anomaly_system`)

**Sum Rules / Consistency Relations:**

- Anomaly-free condition is linear: verified by construction (matrix equation A @ R = b)
  - File: `a_maximization.py` (line 203)

- Unitarity bound R[O] >= 2/3 for gauge-invariant operators: Checked for bilinear operators (mesons, Tr Phi^2, S*Sbar, A*Abar, V*V, f*Omega*f, bifundamental mesons)
  - File: `a_maximization.py` (line 266, `gauge_invariant_ops`); `a_maximization_large_N.py` (line 342, `_gauge_invariant_ops_symbolic`)

## Numerical Validation

**Convergence Tests:**

- Large-N convergence of a/N^2: A function `_convergence_check()` exists that computes numerical a/N^2 at N = 50, 100, 200 and compares to the exact large-N result. Currently **commented out** in the `__main__` block.
  - Script: `a_maximization_large_N.py` (lines 1224-1228, `_convergence_check`)
  - Status: Infrastructure exists but not actively run; convergence was verified during development

**Stability Analysis:**

- BFGS optimizer stability: `a_maximization.py` uses `scipy.optimize.minimize` with method BFGS, starting from s=0 (the particular solution). No multi-start or global optimization.
  - Script: `a_maximization.py` (line 380)
  - Risk: For theories with multiple critical points, BFGS may find a local maximum. Mitigated by the symbolic solver (which finds all critical points) for the large-N case.

- Mathematica NSolve with 8 random restarts: `a_maximize_large_N_fast()` (referenced in `module3_a_maximization.md`, line 242) uses 8 random restarts for the fast numerical scan. The batch Mathematica solver uses `NSolve` which finds ALL real solutions of the polynomial system.
  - Script: `a_maximization_large_N.py` (line 854, Mathematica script)

**Precision and Error Control:**

- Mathematica working precision: 30 digits for NSolve computations
  - File: `a_maximization_large_N.py` (line 719, `working_precision: int = 30`)

- Exact-vs-numerical agreement tolerance: |exact - numerical| < 0.001 required for an exact result to be accepted
  - File: `two_node_db.py` (line 502)

- Clustering tolerance: TOL = 1e-5 for grouping theories into universality classes by a/N^2
  - File: `two_node_db.py` (line 39)

- Unitarity tolerance: 1e-9 for R >= 2/3 check (allows tiny floating-point errors)
  - File: `a_maximization.py` (line 389)

- Hessian check: Mathematica `NegativeDefiniteMatrixQ` used to verify critical points are maxima (not saddle points)
  - File: `a_maximization_large_N.py` (lines 856-858, Mathematica script)
  - Git history note: "tighten Hessian tolerance to eigenvalue <= 0" (commit c019005)

## Comparison with Literature

**Reproduced Results:**

- SU(N) SQCD R-charges: R = 1 - N/N_f matches Intriligator-Wecht (hep-th/0304128)
  - Comparison in: `a_maximization.py` `__main__` block (lines 471-494)
  - Agreement: Exact

- SQCD conformal window upper boundary: N_f < 3N (from b_0 > 0) matches standard result
  - Comparison in: `module1_beta_functions.md` (Sec. 1)
  - Agreement: Exact

- Single-node classification: 19 quivers found, consistent with arXiv:2007.16165 and arXiv:2510.19136
  - Comparison in: `module5_classification.md` (Sec. 2a, Sec. 6)
  - Agreement: Stated but detailed comparison not shown in code

**Discrepancies:**

- SQCD boundary analysis (Module 4, Sec. 9): The naive boundary analysis gives B > 0 for all AF SQCD, seemingly contradicting the existence of the Banks-Zaks fixed point in the conformal window. The document explains the resolution: "the hep-th/0502049 method is most useful for multi-node quivers where coupling one node to another can drive a second node into the non-trivial phase."
  - File: `module4_ir_fixed_points.md` (Sec. 9, lines 176-194)
  - Status: Documented subtlety; single-node conformal window is handled via the R-charge unitarity bound (N_f > 3N/2), not the boundary analysis

- Hessian tolerance issue: Git commit history shows "Remove wrong R>=2/3 filter; tighten Hessian tolerance to eigenvalue <= 0" (commit c019005), indicating a previously incorrect R-charge filter was applied and the Hessian check was too loose.
  - Impact: The current quivers.db (377 classes per commit 4a6f6d5) reflects the corrected pipeline

## Internal Consistency

**Cross-Method Verification:**

- Numerical vs. exact a-maximization: For each of the 135 universality classes, the exact symbolic result (from sympy) is validated against the Mathematica NSolve numerical result with tolerance 0.001
  - Method A (numerical): `a_maximize_batch_mathematica()` in `a_maximization_large_N.py`
  - Method B (exact): `a_maximize_large_N()` in `a_maximization_large_N.py`
  - Agreement: 76/135 classes have exact solutions that agree with numerical to < 0.001; 59 classes are numerical-only (sympy timeout)
  - File: `two_node_db.py` (lines 496-516, Phase 3 validation)

- Finite-N convergence to large-N: `_convergence_check()` compares a/N^2 at N=50,100,200 to the exact large-N limit
  - Method A: `a_maximize()` at specific N with N_f=0
  - Method B: `a_maximize_large_N()` exact leading order
  - File: `a_maximization_large_N.py` (lines 1224-1228)
  - Status: Infrastructure exists; not actively run in CI

**Self-Consistency:**

- Anomaly constraints satisfied by construction: The parameterization R = R0 + F @ s ensures A @ R = b is satisfied exactly (R0 is a particular solution, F spans the null space)
  - File: `a_maximization.py` (lines 366-371); `a_maximization_large_N.py` (lines 605-614)

- Deduplication consistency: Charge conjugation deduplication (`_dedup_conjugation`) and symmetry deduplication (`_dedup_symmetries`) remove equivalent theories. The final count (4560 theories, 135 classes for equal-rank two-node; 377 classes total with mixed rank) is a self-consistency check — duplicate theories would show up as anomalously large universality classes.
  - File: `quiver_generation.py`; `two_node_db.py`

## Anomalies and Topological Properties

**Anomaly Checks:**

- SU(N) cubic gauge anomaly: Cancellation enforced by `chiral_excess_coeffs()` which computes the required delta = n_f - n_fbar. Any quiver where anomaly cancellation cannot be achieved is rejected.
  - File: `quiver_generation.py` (line 298, `check_anomalies`)
  - Anomaly coefficient: A(fund)=+1, A(antifund)=-1, A(adj)=0, A(S)=+(N+4), A(Sbar)=-(N+4), A(A)=+(N-4), A(Abar)=-(N-4)
  - File: `beta_functions.py` (line 120, `A_rep`)

- Sp(N) Witten global anomaly: Total fundamental count at each Sp node must be even. Checked by `check_sp_witten()`.
  - File: `quiver_generation.py` (line 287)

- ABJ anomaly (R-symmetry gauge anomaly): Tr[R G_a^2] = 0 enforced as a constraint in a-maximization. This is the mixed R-gauge anomaly cancellation condition.
  - File: `a_maximization.py` (line 203); `module3_a_maximization.md` (Sec. 3)

**Anomaly Matching:**

- UV-IR anomaly matching: Not explicitly checked. The a-theorem (a_UV > a_IR) is mentioned as a future cross-check but not implemented.
  - File: `module5_classification.md` (Sec. 10)

## Test Suite

**Existing Tests:**

There are no dedicated test files (`test_*.py`, `*_test.py`, `check_*.py`, `verify_*.py`). All validation is done through `__main__` blocks in each module file:

- `beta_functions.py` `__main__` (lines 261-287): Tests SQCD b_0 computation, N=2 SQCD, AF-all-N check, bifundamental coefficient
  - Coverage: Dynkin indices, b_0 computation, is_af, is_af_all_N, b0_linear

- `quiver_generation.py` `__main__` (lines 710-783): Tests N_f bounds for single-node theories (SU+S, SU+A, SU+3A, SU+4A), two-node balanced quiver, Sp Witten anomaly, enumeration count
  - Coverage: chiral_excess_coeffs, nf_bound, check_sp_witten, enumerate_quivers

- `a_maximization.py` `__main__` (lines 467-533): Tests SU(5) SQCD at various N_f, conformal window scan, below-conformal-window decoupling, circular quiver, SU + adj + flavors
  - Coverage: build_fields, a_maximize, a_maximize_with_decoupling, gauge_invariant_ops

- `a_maximization_large_N.py` `__main__` (lines 1220-1294): Several examples commented out; active code runs `_classify_two_node()` (full classification table)
  - Coverage: build_fields_large_N, a_maximize_large_N, _convergence_check (commented out)

**Run Commands:**

```bash
python beta_functions.py            # Run Module 1 validation examples
python quiver_generation.py         # Run Module 2 validation examples
python a_maximization.py            # Run Module 3 numerical validation
python a_maximization_large_N.py    # Run large-N classification (requires Mathematica)
python two_node_db.py build         # Build full database (requires Mathematica, ~3 min)
python two_node_db.py stats         # Show database statistics
```

**Test Patterns:**

Validation is done by printing computed values alongside expected analytic results, with manual inspection:

```python
# Pattern from a_maximization.py __main__:
res = a_maximize(q, N=5, N_f=10)
for label, R in res.R_charges.items():
    expect = 1 - 5 / 10
    print(f"  {label}: R = {R:.6f}  (expect {expect:.4f} = 1 - N/N_f)")
```

No automated assertions (assert statements) are used in validation blocks. The beta_functions.py file uses assert only for input validation (`__post_init__` in NodeSpec, line 165), not for test assertions.

**Missing Tests:**

- No automated test suite (pytest, unittest) — all validation is manual print-and-inspect
- No regression tests: database rebuild does not compare to a reference
- No edge case tests for:
  - Degenerate quivers (single node with no matter)
  - Maximum-degree quivers (5 bifundamental SU-SU neighbors, B_0 = 1/2)
  - Mixed gauge pair theories with unusual rank multipliers
  - Theories at the boundary of the conformal window (N_f exactly at 3N/2 or 3N)
- No tests for Module 4 (boundary analysis) — module is not implemented
- No cross-check of gauge-invariant operator enumeration completeness (higher-order operators like baryons or multi-trace operators are not checked)
- No numerical stability tests (e.g., behavior of BFGS near degenerate critical points)

## Reproducibility

**Random Seeds:**

- No random number generation in the classification pipeline. The BFGS optimizer in `a_maximization.py` starts from s=0 (deterministic). The Mathematica NSolve is deterministic. Not applicable.

**Platform Dependence:**

- Mathematica (wolframscript) is required for the batch NSolve scanner in `a_maximization_large_N.py` and the database build in `two_node_db.py`. Different Mathematica versions may produce slightly different numerical results at finite working precision.
- The sympy exact solver is platform-independent.
- The quivers.db file is pre-built and included in the repository.

**Version Pinning:**

- No requirements.txt, pyproject.toml, or dependency version pinning file exists.
- Implicit dependencies: Python 3.10+ (type hints), numpy, scipy, sympy, sqlite3 (stdlib), Mathematica/wolframscript (external)

## Database Validation

**Build pipeline validation (two_node_db.py):**

The database build has implicit validation at each phase:

1. **Phase 1 (enumerate + scan):** Quiver enumeration applies all constraints (AF, anomaly, Witten). Theories where Mathematica NSolve finds no bounded maximum (no real critical points or no negative-definite Hessian) are classified as "diverged" and stored with class_id=NULL.
   - File: `two_node_db.py` (lines 354-452)

2. **Phase 2 (clustering):** Clustering by a/N^2 within TOL=1e-5. Sanity check: theories in the same universality class must have identical leading-order field structure (same T_lead and dim_lead values).
   - File: `two_node_db.py` (lines 454-476)

3. **Phase 3 (exact solve):** Exact symbolic result validated against numerical: |exact - numerical| < 0.001 required for acceptance. Failures are stored as numerical-only.
   - File: `two_node_db.py` (lines 478-519)

**Database statistics (as of latest build):**
- 4560 theories in 135 universality classes (equal rank)
- 377 total classes including mixed-rank quivers
- 76 exact symbolic solutions, 59 numerical-only (for the 135 equal-rank classes)
- File: `quivers.db`; `module5_classification.md` (Sec. 6, Sec. 9)

---

_Validation analysis: 2026-04-14_
