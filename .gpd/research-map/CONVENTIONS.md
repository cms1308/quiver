# Derivation Quality Analysis

**Analysis Date:** 2026-04-14

## Unit System and Natural Units

This project works entirely in **natural units** (hbar = c = 1). No explicit factors of hbar or c appear anywhere. All quantities are expressed in terms of the gauge group rank parameter N and flavor number N_f:

- Central charges a, c are dimensionless
- R-charges are dimensionless
- Beta function coefficient b_0 is dimensionless (coupling g is dimensionless in 4D)
- Dynkin indices T(r) are dimensionless

No SI or CGS quantities appear. The project is purely algebraic/group-theoretic — no physical dimensionful quantities (energies, lengths, etc.) are computed.

## Derivation Inventory

**One-loop beta function coefficient b_0:**
- Result: b_0 = 3 T(adj) - sum_i T(r_i)
- File: `module1_beta_functions.md` (Sec. 1, Eq. for b_0)
- Method: Standard one-loop perturbative result (textbook)
- Starting point: NSVZ exact beta function at one loop
- Key steps: Sum Dynkin indices over all chiral multiplets charged under a gauge node
- Status: Complete — implemented in `beta_functions.py`, function `compute_b0()` (line 177)

**Large-N AF condition:**
- Result: b_0(N) = alpha*N + beta > 0 for all N requires alpha > 0 or (alpha = 0 and beta >= 0)
- File: `module1_beta_functions.md` (Sec. 4-5); `beta_functions.py` (line 198, `b0_linear`)
- Method: Exact — b_0 is linear in N, so two evaluations determine the coefficients
- Starting point: b_0 at specific N values
- Status: Complete

**Gauge anomaly cancellation (SU(N)):**
- Result: N*(deg+ - deg-) + (n_S - n_Sbar)*(N+4) + (n_A - n_Abar)*(N-4) + (n_f - n_fbar) = 0
- File: `module2_quiver_generation.md` (Sec. 4a); `quiver_generation.py` (line 157, `chiral_excess_coeffs`)
- Method: Exact (cubic anomaly coefficient cancellation)
- Starting point: Anomaly coefficients A(r) per representation
- Status: Complete

**a-maximization (trial a-function):**
- Result: a = (3/32)(3 Tr R^3 - Tr R), maximized subject to Tr[R G_a^2] = 0 for all gauge nodes
- File: `module3_a_maximization.md` (Sec. 5-6); `a_maximization.py` (line 252, `a_trial`)
- Method: Constrained optimization — anomaly-free conditions are linear constraints, a is cubic in R-charges
- Starting point: Intriligator-Wecht a-theorem (hep-th/0304128)
- Key steps: (1) Build anomaly constraint matrix A*R = b, (2) Find particular solution via lstsq, (3) Null space gives free parameters, (4) Maximize cubic a over free parameters via BFGS
- Status: Complete — two implementations: numerical (`a_maximization.py`) and exact symbolic (`a_maximization_large_N.py`)

**Large-N a-maximization (leading order):**
- Result: Exact algebraic R-charges and a/N^2, c/N^2 as sympy expressions
- File: `a_maximization_large_N.py` (line 430, `a_maximize_large_N`)
- Method: Exact symbolic — at large N, anomaly constraints are linear in R, critical-point equations are quadratic, yielding closed-form algebraic solutions
- Starting point: Leading-order coefficients T_lead, dim_lead (all scale as powers of N)
- Key steps: (1) Build leading-order fields, (2) Symbolic anomaly constraints via sympy, (3) Null space, (4) Symbolic differentiation of a, (5) Solve quadratic critical-point equations exactly
- Status: Complete — produces exact results for 76/135 two-node classes; 59 classes require Mathematica NSolve (timeout at 30s due to 4+ free parameters)

**Unitarity decoupling:**
- Result: When a gauge-invariant operator O has R[O] < 2/3, pin R[O] = 2/3 and re-maximize
- File: `module3_a_maximization.md` (Sec. 9); `a_maximization.py` (line 394, `a_maximize_with_decoupling`)
- Method: Iterative — add linear constraint per violating operator, redo a-maximization until convergence
- Starting point: Kutasov-Schwimmer unitarity bound on gauge-invariant operators
- Status: Complete for both numerical and symbolic implementations

**Boundary analysis (Module 4):**
- Result: B_a = Tr[R^(a) G_a^2] < 0 for all nodes a is necessary for non-trivial IR fixed point
- File: `module4_ir_fixed_points.md` (Sec. 2-4)
- Method: Evaluate NSVZ beta function numerator at the boundary where g_a -> 0
- Status: **Not yet implemented in code** (documented only)

## Approximations Made

**Large-N leading-order approximation (a_maximization_large_N.py):**
- What is neglected: All O(1) and subleading terms in N — fund-like matter with fixed N_f, O(1) corrections to Dynkin indices and dimensions
- Justification given: "At large N, all quantities scale polynomially: T ~ N, dim ~ N^2. Dividing through gives N-independent leading coefficients." (`a_maximization_large_N.py`, line 4-7)
- Justification quality: Strong — controlled expansion in 1/N with well-defined leading order
- Parameter controlling approximation: 1/N (implicit; all nodes share the same N)
- Estimated error: O(1/N) corrections to R-charges; O(1/N^2) corrections to a/N^2
- File: `a_maximization_large_N.py` (lines 1-13)

**One-loop beta function (perturbative truncation):**
- What is neglected: Two-loop and higher corrections; anomalous dimensions gamma_i at higher loops
- Justification given: "At the IR fixed point, beta(g) = 0 and gamma_i = 3R_i - 2, recovering the anomaly-free R-symmetry condition" — the one-loop result gives the correct fixed-point condition when combined with a-maximization
- Justification quality: Strong — the NSVZ exact beta function shows that the one-loop b_0 correctly determines asymptotic freedom, and a-maximization captures the exact IR physics
- Parameter controlling approximation: g^2/(16 pi^2) (coupling constant)
- File: `module1_beta_functions.md` (Sec. 6)

**Universality class clustering tolerance:**
- What is neglected: Numerical differences in a/N^2 below 10^{-5}
- Justification given: Theories with identical leading-order field structure must have identical a/N^2 at infinite precision
- Justification quality: Adequate — tolerance chosen empirically; exact symbolic computation confirms agreement for 76/135 classes
- Parameter controlling approximation: TOL = 1e-5 (clustering threshold)
- File: `two_node_db.py` (line 39); `a_maximization_large_N.py` (line 1150)

**Exclusion of superpotential deformations:**
- What is neglected: All superpotential terms W — no cubic or higher couplings are included
- Justification given: Not explicitly justified; the classification enumerates theories without specifying W
- Justification quality: Weak — the absence of superpotential means R-charges are determined purely by anomaly-free conditions and a-maximization, but physical theories may require specific W terms for consistency
- Parameter controlling approximation: N/A
- File: `module3_a_maximization.md` (Sec. 3, mentions W constraints but does not implement them)

## Assumptions Catalog

**Explicit Assumptions:**

- All gauge nodes share the same rank parameter N (equal-rank assumption)
  - File: `CLAUDE.md` (Physics Conventions); `module2_quiver_generation.md` (Sec. 2)
  - Exception: mixed-rank quivers with rank_multipliers [2,1] and [1,2] are also enumerated in `two_node_db.py` (line 363-368)

- Inter-node matter is restricted to bifundamentals only
  - File: `module1_beta_functions.md` (Sec. 3); `CLAUDE.md` (Physics Conventions)
  - Justification: Higher representations (e.g., (fund, adj)) have dim ~ N^2 multiplicity, overwhelming the 3N gauge contribution to b_0

- SO(N) spinor representations are excluded
  - File: `module1_beta_functions.md` (Sec. 2, SO table note); `module2_quiver_generation.md` (Sec. 3)
  - Justification: Spinor dimensions grow exponentially with N, violating large-N factorization

- No superpotential (W = 0)
  - File: `plan.md` (Open question 1); implicit throughout all module docs
  - Risk: Some theories may require specific W for anomaly cancellation or to avoid runaway directions

**Implicit Assumptions:**

- Unitarity bound applies only to gauge-invariant operators, not elementary fields
  - File: `module3_a_maximization.md` (Sec. 9); `a_maximization.py` (line 389)
  - Stated explicitly in documentation, but the implementation requires careful enumeration of gauge-invariant operators; the list in `gauge_invariant_ops()` is limited to simple bilinear and meson operators
  - Risk: Higher-order gauge-invariant operators (e.g., baryons, multi-trace operators) are not checked. An unchecked operator below 2/3 would invalidate the R-charge solution.

- The a-function has a unique global maximum (not just a critical point)
  - File: `a_maximization.py` uses BFGS (local optimizer); `a_maximization_large_N.py` enumerates all critical points and picks the maximum (line 627-640)
  - Risk: BFGS in `a_maximize()` may find a local maximum; the symbolic solver is exhaustive

- All IR fixed points are interacting SCFTs (not free or confined)
  - File: `module4_ir_fixed_points.md` (Sec. 1)
  - Risk: Module 4 boundary analysis (not yet implemented) is needed to confirm this

## Mathematical Rigor Assessment

**a-maximization (a_maximization.py):**
- Rigor level: Physicist-standard
- Issues:
  - Numerical optimization via BFGS — no proof that the found critical point is the global maximum
  - Anomaly constraint satisfaction is checked to floating-point precision only (lstsq residual)
  - Unitarity tolerance of 1e-9 used for R >= 2/3 check (`a_maximization.py`, line 389)
- File: `a_maximization.py`

**a-maximization, large N (a_maximization_large_N.py):**
- Rigor level: Rigorous (for the exact symbolic solver)
- Issues:
  - sympy `solve()` may miss solutions for high-degree systems (unlikely for quadratic)
  - Hessian check in Mathematica NSolve path uses `NegativeDefiniteMatrixQ` — rigorous at working precision 30
- File: `a_maximization_large_N.py`

**Quiver enumeration (quiver_generation.py):**
- Rigor level: Rigorous — systematic enumeration with backtracking pruning
- Issues:
  - Deduplication by charge conjugation (`_dedup_conjugation`) must be verified to not miss any inequivalent theories
  - Mixed-rank deduplication (`_dedup_symmetries`) across [2,1] and [1,2] multipliers
- File: `quiver_generation.py`

## Dimensional Analysis

**Consistency Checks:**

All equations in this project are dimensionless in natural units. Specific checks:

- b_0 = 3 T(adj) - sum T(r_i): Both terms are dimensionless (Dynkin indices are pure numbers). Dimensionless — verified.
  - File: `module1_beta_functions.md` (Sec. 1)

- a = (3/32)(3 Tr R^3 - Tr R): R-charges are dimensionless, traces sum over fermion multiplicities (integers). The prefactor 3/32 is a convention. Dimensionless — verified.
  - File: `module3_a_maximization.md` (Sec. 5)

- c = (1/32)(9 Tr R^3 - 5 Tr R): Same structure as a. Dimensionless — verified.
  - File: `module3_a_maximization.md` (Sec. 7)

- Anomaly-free condition: T(adj_a) + sum_i T_{G_a}(r_i)(R_i - 1) = 0: All terms are dimensionless (Dynkin index times dimensionless R-charge shift). Dimensionless — verified.
  - File: `module3_a_maximization.md` (Sec. 3)

- Gauge anomaly: sum_i A(r_i) * dim(other groups) = 0: Anomaly coefficients A(r) are integers, dim is integer. Dimensionless — verified.
  - File: `module2_quiver_generation.md` (Sec. 4a)

**Dimensional Anomalies:**
- None detected. All equations are purely algebraic/combinatorial.

## Sign and Factor Conventions

**Sign Choices:**

- Beta function sign convention: beta(g) = -(g^3/16pi^2) * b_0, so b_0 > 0 means asymptotic freedom (coupling decreases at high energy)
  - File: `module1_beta_functions.md` (Sec. 1)
  - Consistent throughout: Yes

- Anomaly coefficient signs: A(fund) = +1, A(antifund) = -1 (positive for fundamental, negative for conjugate)
  - File: `module1_beta_functions.md` (Sec. 2, SU table); `beta_functions.py` (line 120, `A_rep`)
  - Consistent throughout: Yes

- Chiral excess delta = n_f - n_fbar (positive means more fundamentals than anti-fundamentals)
  - File: `quiver_generation.py` (line 157, `chiral_excess_coeffs`)
  - Consistent throughout: Yes

- Bifundamental edge direction: `Edge(src, dst, "+-")` means (fund_src, antifund_dst)
  - File: `quiver_generation.py` (line 56-62, `Edge` docstring)
  - "+" = fund (square), "-" = antifund (square-bar)
  - Consistent throughout: Yes

**Factor Tracking:**

- Factor of 3/32 in a-function: standard 4D N=1 convention
  - File: `module3_a_maximization.md` (Sec. 5); `a_maximization.py` (line 253)
  - Consistent between numerical and symbolic implementations: Yes

- Dynkin index normalization: T(fund_SU) = 1/2 (not 1)
  - File: `CLAUDE.md`; `module1_beta_functions.md` (Sec. 2)
  - Used consistently in all code: `beta_functions.py` (line 50), `a_maximization.py`, `a_maximization_large_N.py`

## Notation Consistency

**Gauge Group Naming:**

| Symbol in docs | Symbol in code | Meaning |
|---|---|---|
| SU(N) | `"SU"` | Special unitary group |
| SO(N) | `"SO"` | Special orthogonal group |
| Sp(N) = USp(2N) | `"Sp"` | Symplectic group; N is the rank, group order is 2N |

- File: `CLAUDE.md` (Physics Conventions); `beta_functions.py` (line 5, `GaugeType`)
- Consistent throughout: Yes

**Representation Naming:**

| Physics symbol | Code string | Gauge group |
|---|---|---|
| square (fund) | `"fund"` | SU |
| square-bar (antifund) | `"antifund"` | SU |
| adj | `"adj"` | SU, SO, Sp |
| S (rank-2 symmetric) | `"S"` | SU, SO |
| S-bar (conj. symmetric) | `"Sbar"` | SU |
| A (rank-2 antisymmetric) | `"A"` | SU, Sp |
| A-bar (conj. antisymmetric) | `"Abar"` | SU |
| V (vector) | `"V"` | SO |
| f (fundamental) | `"fund"` | Sp |

- File: `CLAUDE.md`; `beta_functions.py` (lines 19-29, `VALID_REPS`)
- CLAUDE.md prescribes: "use A for rank-2 antisymmetric, never Lambda or Lambda^2"
- Consistent throughout: Yes

**Edge representation codes:**

| Code | Meaning | Gauge pair |
|---|---|---|
| `"+-"` | (fund, antifund) | SU-SU |
| `"++"` | (fund, fund) | SU-SU |
| `"--"` | (antifund, antifund) | SU-SU |
| `"+"` | (fund, V/f) | SU-SO, SU-Sp |
| `"-"` | (antifund, V/f) | SU-SO, SU-Sp |
| `"std"` | self-conjugate | SO-SO, SO-Sp, Sp-Sp |

- File: `quiver_generation.py` (lines 40-47, `EDGE_REPS`)

**Variable naming code-to-physics mapping:**

| Code variable | Physics quantity |
|---|---|
| `b0` | One-loop beta function coefficient b_0 |
| `T_adj`, `T_rep`, `T_bifund` | Dynkin indices T(r) |
| `A_rep` | Cubic anomaly coefficient A(r) |
| `R_values`, `R_opt`, `R_charges` | R-charge assignments R[Phi_i] |
| `a_val`, `a_opt`, `a_over_N2` | Central charge a (or a/N^2 at large N) |
| `c_val`, `c_opt`, `c_over_N2` | Central charge c (or c/N^2 at large N) |
| `N_f` | Number of fundamental flavor pairs |
| `delta` | Chiral excess n_f - n_fbar |
| `alpha`, `beta` in `b0_linear` | Coefficients of b_0(N) = alpha*N + beta |
| `n_free` | Number of free parameters after imposing anomaly constraints |
| `s_syms`, `s_opt` | Free parameters in null space of anomaly matrix |

**Notation Conflicts:**
- `N` is used for both the gauge group rank parameter and (in sympy) `from sympy import ... N ...` (but the project avoids this import conflict by using specific imports)
- `alpha`, `beta` in `b0_linear` are not the beta function itself but coefficients of the linear b_0(N) expansion. Context-clear, no actual conflict.

## Sp(N) Convention

Use Sp(N) = USp(2N) throughout (rank N, fundamental dimension 2N). This is the physics convention, not the mathematics convention (Sp(2N)).

- File: `CLAUDE.md` (Physics Conventions); `module1_beta_functions.md` (Sec. 2)
- Code: `_fund_dim("Sp", N)` returns `2*N` (`beta_functions.py`, line 82)
- Consistent throughout: Yes

## Exact Arithmetic

Use `fractions.Fraction` for all Dynkin index and anomaly coefficient computations to avoid floating-point errors in the combinatorial/algebraic part of the pipeline:

- File: `beta_functions.py` (line 9, `from fractions import Fraction`; all T functions return `Fraction`)
- File: `a_maximization_large_N.py` uses `sympy.Rational` for exact symbolic computation
- File: `a_maximization.py` converts to float only at the numerical optimization stage

## N_f Parametrization Convention

Rather than fixing explicit fundamental counts, each theory is parametrized by a universal N_f:
- SU(N): N_f = min(n_f, n_fbar) — number of complete (fund, antifund) pairs
- SO(N): N_f = number of vectors
- Sp(N): N_f = half the number of fundamentals (so 2*N_f fundamentals total)

The chiral excess delta = n_f - n_fbar is fixed by anomaly cancellation and is a function of N (linear: delta(N) = a*N + b).

- File: `quiver_generation.py` (lines 8-16, docstring); `module2_quiver_generation.md` (Sec. 9)

## Minimum N Values

| Gauge group | N_MIN |
|---|---|
| SU(N) | 5 |
| SO(N) | 5 |
| Sp(N) | 3 |

- File: `beta_functions.py` (line 16, `N_MIN`)
- Note: These are chosen conservatively — SU(2) = Sp(1) and SO(3) = SU(2)/Z_2 have special properties. The minimum N values ensure the representations listed are all distinct and non-trivial.

---

_Derivation quality analysis: 2026-04-14_
