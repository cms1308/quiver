# Theoretical Frameworks

**Analysis Date:** 2026-04-14

## Physical System

**Subject:** 4D N=1 supersymmetric quiver gauge theories -- classification of asymptotically free theories with large N limit that flow to non-trivial IR superconformal fixed points.

**Scales:**

- Energy: UV asymptotically free regime (weak coupling at high energy) flowing to IR superconformal fixed point (strong coupling). No explicit energy scale; the theory is conformal at the fixed point.
- Length: Not applicable (conformal field theory at the IR fixed point has no intrinsic length scale).
- Time: Not applicable (Lorentzian 4D, but all analysis is Euclidean/algebraic).
- Dimensionless parameters:
  - N (gauge group rank parameter, shared by all nodes) -- the large N parameter
  - N_f (number of fundamental flavors per node) -- parametrizes the conformal window
  - b_0^{(a)} / N (leading beta function coefficient per node) -- controls UV behavior
  - a/N^2, c/N^2 (central charges normalized by N^2) -- IR conformal data

**Degrees of Freedom:**

- Gauge vector superfields V_a (one per quiver node a): contain gauge boson A_mu and gaugino lambda_a. Defined in `module1_beta_functions.md` (Sec. 1).
- Chiral superfields Phi_i (matter): contain complex scalar phi_i and Weyl fermion psi_i. Carry representations of one or two gauge groups. Defined in `module2_quiver_generation.md` (Sec. 3).
  - Node matter: adjoint, rank-2 symmetric (S, Sbar), rank-2 antisymmetric (A, Abar), fundamental, anti-fundamental
  - Bifundamental matter: fundamental x fundamental/anti-fundamental of two gauge groups at quiver edges

## Theoretical Framework

**Primary Framework:**

- 4D N=1 supersymmetric quantum field theory
- Formulation: Superspace / component fields. The analysis uses the NSVZ exact beta function and Intriligator-Wecht a-maximization. No Lagrangian is written explicitly; the theory is defined by gauge group content and matter representations.
- File: `README.md` (Physical Setup, lines 28-46); `CLAUDE.md` (Physics Conventions)

**Secondary/Supporting Frameworks:**

- Representation theory of SU(N), SO(N), USp(2N): Used to compute Dynkin indices, anomaly coefficients, and dimensions.
  - File: `module1_beta_functions.md` (Sec. 2); `beta_functions.py`
- Graph theory (quiver diagrams): Quivers encode the gauge theory data as directed multigraphs with labeled nodes and edges.
  - File: `module2_quiver_generation.md` (Sec. 1)
- Superconformal algebra: Fixes gaugino R-charge to R[lambda] = 1 and relates anomalous dimensions to R-charges via gamma_i = 3R_i - 2 at the fixed point.
  - File: `module3_a_maximization.md` (Sec. 2); `module1_beta_functions.md` (Sec. 6)

## Fundamental Equations

### Equation Catalog

| ID | Equation | Type | Location | Dimensions | Status | Depends On | Used By |
|----|----------|------|----------|------------|--------|------------|---------|
| EQ-001 | b_0 = 3 T(adj) - sum_i T(r_i) | Defining | `module1_beta_functions.md` (Sec. 1, boxed eq.) | dimensionless -- verified (T(r) is dimensionless in the convention used) | Postulated (standard result) | -- | EQ-002, EQ-003, EQ-008 |
| EQ-002 | mu dg/dmu = -(g^3 / 16pi^2) b_0 + O(g^5) | Derived | `module1_beta_functions.md` (Sec. 1, first eq.) | [mass]^0 (both sides dimensionless: g is dimensionless in 4D) -- verified | Standard perturbative result | EQ-001 | EQ-003 |
| EQ-003 | beta(g) = -(g^3/16pi^2) [3T(adj) - sum_i T(r_i)(1 - gamma_i)] / [1 - T(adj) g^2/(8pi^2)] | Exact (NSVZ) | `module1_beta_functions.md` (Sec. 6) | dimensionless -- verified | EQ-001 | EQ-004, EQ-008 |
| EQ-004 | T(adj_a) + sum_i T_{G_a}(r_i)(R[Phi_i] - 1) = 0 | Constraint (anomaly-free R-symmetry) | `module3_a_maximization.md` (Sec. 3, boxed eq.) | dimensionless -- verified (T and R are both dimensionless) | Derived from Tr[R G_a^2] = 0 | EQ-003 | EQ-005, EQ-008 |
| EQ-005 | a = (3/32)(3 Tr R^3 - Tr R) | Defining | `module3_a_maximization.md` (Sec. 1) | dimensionless -- verified (trace over Weyl fermions with dimensionless R-charges) | Standard 't Hooft anomaly relation | -- | EQ-006, EQ-007 |
| EQ-006 | c = (1/32)(9 Tr R^3 - 5 Tr R) | Defining | `module3_a_maximization.md` (Sec. 7) | dimensionless -- verified | Standard 't Hooft anomaly relation | -- | EQ-007 |
| EQ-007 | a_trial(R) = (3/32)[2 sum_a dim(G_a) + sum_i d_i(3(R_i-1)^3 - (R_i-1))] | Derived | `module3_a_maximization.md` (Sec. 5) | dimensionless -- verified | EQ-005 | EQ-009 |
| EQ-008 | B_a = T(adj_a) + sum_i T_{G_a}(r_i)(R_i^{(a)} - 1) | Derived | `module4_ir_fixed_points.md` (Sec. 4, Step 3) | dimensionless -- verified | EQ-004 | EQ-009 |
| EQ-009 | IR fixed point criterion: B_a < 0 for all a | Criterion | `module4_ir_fixed_points.md` (Sec. 4, Step 4) | N/A (inequality) | Cited from hep-th/0502049 | EQ-007, EQ-008 | -- |
| EQ-010 | SU(N) anomaly cancellation: N sum_b(n_{ab}^{+-} - n_{ba}^{+-} + n_{ab}^{++} - n_{ab}^{--}) + (n_f - n_fbar) + (n_S - n_Sbar)(N+4) + (n_A - n_Abar)(N-4) = 0 | Constraint | `module2_quiver_generation.md` (Sec. 4a) | dimensionless -- verified (all terms are pure numbers) | Standard gauge anomaly cancellation | -- | EQ-001 |
| EQ-011 | Sp(N) Witten anomaly: n_f^{Sp} = 0 mod 2 | Constraint | `module2_quiver_generation.md` (Sec. 4b) | N/A (discrete condition) | Standard global anomaly | -- | -- |
| EQ-012 | Tr R = sum_a dim(G_a) + sum_i d_i (R_i - 1) | Derived | `module3_a_maximization.md` (Sec. 4) | dimensionless -- verified | -- | EQ-005, EQ-006 |
| EQ-013 | Tr R^3 = sum_a dim(G_a) + sum_i d_i (R_i - 1)^3 | Derived | `module3_a_maximization.md` (Sec. 4) | dimensionless -- verified | -- | EQ-005, EQ-006 |

**Equation of Motion / Field Equations:**

Not applicable. The project does not solve classical equations of motion. The analysis is entirely algebraic: it computes one-loop beta function coefficients, anomaly cancellation conditions, and a-maximization at the conformal fixed point. The NSVZ beta function (EQ-003) provides the exact beta function, but the project uses only the one-loop coefficient (EQ-001) for UV analysis and the vanishing condition (EQ-004) for IR analysis.

**Constraints:**

- Asymptotic freedom: b_0^{(a)} > 0 for all nodes a and all valid N (EQ-001)
  - File: `module1_beta_functions.md` (Sec. 1, 5); `beta_functions.py` (`is_af_all_N`, line 241)
- Large N admissibility: b_0^{(a)} > 0 as N -> infinity. Equivalent to alpha > 0 (or alpha = 0, beta > 0) where b_0 = alpha*N + beta.
  - File: `beta_functions.py` (`b0_linear`, line 198; `is_af_all_N`, line 241)
- SU(N) cubic gauge anomaly cancellation (EQ-010)
  - File: `module2_quiver_generation.md` (Sec. 4a); `quiver_generation.py` (`chiral_excess_coeffs`)
- Sp(N) Witten global anomaly (EQ-011)
  - File: `module2_quiver_generation.md` (Sec. 4b); `quiver_generation.py` (`check_sp_witten`)
- Anomaly-free R-symmetry at each gauge node (EQ-004)
  - File: `module3_a_maximization.md` (Sec. 3)
- Unitarity bound: R[O] >= 2/3 for all gauge-invariant chiral operators O (NOT for elementary gauge-charged fields)
  - File: `module3_a_maximization.md` (Sec. 9)
- Superpotential constraint: R[W] = 2 (not currently used; no superpotential in the project)
  - File: `module3_a_maximization.md` (Sec. 3)
- Bifundamental-only inter-node matter (from large N AF)
  - File: `module1_beta_functions.md` (Sec. 3)

## Symmetries and Conservation Laws

**Exact Symmetries:**

- N=1 superconformal symmetry (at the IR fixed point): the superconformal algebra SU(2,2|1) contains the conformal group SO(4,2), Poincare supersymmetry, and the exact R-symmetry.
  - File: `module3_a_maximization.md` (Sec. 1, 2)
- Gauge symmetry: product gauge group G_1 x G_2 x ... x G_k where each G_a in {SU(N), SO(N), USp(2N)}
  - File: `README.md` (Sec. "Gauge Groups"); `module2_quiver_generation.md` (Sec. 2)
- Flavor symmetry: at each node, permutation of identical representations. For SU(N) SQCD with N_f flavors: SU(N_f)_L x SU(N_f)_R x U(1)_B x U(1)_R.
  - File: `module3_a_maximization.md` (Sec. 2, 10)

**Approximate Symmetries:**

- Large N universality: at leading order in 1/N, all rank-2 representations (S, Sbar, A, Abar) of SU(N) become indistinguishable (same T_lead = 1/2, dim_lead = N^2/2). This produces universality classes where distinct theories share the same IR physics.
  - File: `a_maximization_large_N.py` (lines 60-80); `module5_classification.md` (Sec. 9)

**Gauge Symmetries:**

- SU(N): requires cubic anomaly cancellation (EQ-010). Gauge-fixing not discussed (all analysis is gauge-invariant).
  - File: `module2_quiver_generation.md` (Sec. 4a)
- SO(N): automatically anomaly-free (real/pseudo-real representations). No additional constraint.
  - File: `module2_quiver_generation.md` (Sec. 4b)
- USp(2N) = Sp(N): pseudo-real, no cubic anomaly. Witten global anomaly requires even number of fundamentals (EQ-011).
  - File: `module2_quiver_generation.md` (Sec. 4b)

**Ward Identities / Selection Rules:**

- Tr[R G_a^2] = 0 for all gauge nodes a (ABJ anomaly cancellation for the R-symmetry): this is the key constraint (EQ-004) that determines IR R-charges.
  - File: `module3_a_maximization.md` (Sec. 3)
- At the NSVZ fixed point: gamma_i = 3R_i - 2 (superconformal Ward identity relating anomalous dimensions to R-charges).
  - File: `module1_beta_functions.md` (Sec. 6)

**Anomalies:**

- ABJ anomaly (Tr[R G_a^2]): must vanish for the R-symmetry to be non-anomalous. This is a linear constraint on R-charges, one per gauge node.
  - File: `module3_a_maximization.md` (Sec. 3)
- SU(N) cubic gauge anomaly (Tr[G_a^3]): must vanish for gauge consistency. Anomaly coefficient A(r) tabulated in `module1_beta_functions.md` (Sec. 2) and implemented in `beta_functions.py` (`A_rep`, line 119).
- Sp(N) Witten anomaly (Z_2 global anomaly): topological obstruction requiring even number of fundamental Weyl fermions.
  - File: `module2_quiver_generation.md` (Sec. 4b)
- 't Hooft anomalies kappa_RRR = Tr R^3, kappa_R = Tr R: these determine central charges a, c (EQ-005, EQ-006).
  - File: `module3_a_maximization.md` (Sec. 8)

**Topological Properties:**

- Not directly studied. The Witten anomaly is a topological (Z_2) constraint on Sp(N) matter content.

**Dualities and Correspondences:**

- Seiberg duality: referenced for SQCD conformal window. Seiberg duals of classified theories should also appear in the classification (mentioned as a cross-check).
  - File: `module5_classification.md` (Sec. 10)
- Large N universality: theories with different rank-2 matter but identical leading-order coefficients form universality classes sharing the same a/N^2 and R-charges.
  - File: `module5_classification.md` (Sec. 9); `a_maximization_large_N.py`

## Parameters and Couplings

**Fundamental Parameters:**

- N (integer >= 2 for SU, >= 3 for SO, >= 1 for Sp): shared rank parameter for all gauge nodes. The large N limit N -> infinity is the primary regime.
  - Defined in: `README.md` (Sec. "Gauge Groups"); `beta_functions.py` (`N_MIN`, line 16)
- N_f (non-negative integer): number of fundamental flavors. Parametrizes the conformal window.
  - Defined in: `module2_quiver_generation.md` (Sec. 9, "N_f Parametrization")
- Gauge group type l_V(a) in {SU, SO, Sp}: discrete label at each quiver node.
  - Defined in: `module2_quiver_generation.md` (Sec. 1)
- Edge representation type: "++" (fund, fund), "+-" (fund, antifund), "--" (antifund, antifund) for SU-SU; "+" or "-" for SU-SO/Sp; "std" for SO/Sp-SO/Sp.
  - Defined in: `quiver_generation.py` (Edge dataclass)

**Derived Quantities:**

- b_0^{(a)} = alpha*N + beta (one-loop beta function coefficient per node): linear in N.
  - Computed in: `beta_functions.py` (`compute_b0`, line 177; `b0_linear`, line 198)
- B_0^{(a)} = lim_{N->inf} b_0^{(a)}/N (leading large-N beta coefficient)
  - Computed in: `module5_classification.md` (Sec. 3)
- R_i^* (IR superconformal R-charges): determined by a-maximization subject to anomaly-free constraints.
  - Computed in: `a_maximization.py` (`a_maximize`); `a_maximization_large_N.py` (`a_maximize_large_N`)
- a^* = (3/32)(3 Tr R*^3 - Tr R*) (IR a central charge)
  - Computed in: `a_maximization.py` (`a_trial`); `a_maximization_large_N.py`
- c^* = (1/32)(9 Tr R*^3 - 5 Tr R*) (IR c central charge)
  - Computed in: `a_maximization.py` (`c_trial`); `a_maximization_large_N.py`
- delta(N) = a*N + b (chiral excess at SU(N) nodes): determined by anomaly cancellation from rank-2 and bifundamental content.
  - Computed in: `quiver_generation.py` (`chiral_excess_coeffs`)
- B_a (boundary analysis quantity Tr[R^{(a)} G_a^2]): determines IR fixed point existence.
  - Defined in: `module4_ir_fixed_points.md` (Sec. 4); NOT YET IMPLEMENTED in code.

**Dimensionless Ratios:**

- a/N^2: IR central charge normalized by N^2. Finite and N-independent at leading order. Used as the primary classifier for universality classes.
  - Computed in: `a_maximization_large_N.py`; stored in `quivers.db` (`universality_class.a_over_N2`)
- a/c: ratio of central charges. Equals 1 for N=4 SYM; generically != 1 for N=1 theories.
  - File: `module3_a_maximization.md` (Sec. 7)
- b_0/N: large N scaling of the UV beta function. Must be positive for large N AF.

## Phase Structure / Regimes

**Regimes Studied:**

- UV asymptotically free regime: g_a << 1 at high energy. One-loop beta function b_0 > 0 controls the running.
  - Applicable equations: EQ-001, EQ-002
  - File: `module1_beta_functions.md`
- IR superconformal fixed point: g_a = g_a^* with beta(g_a^*) = 0 for all a. R-charges determined by a-maximization.
  - Applicable equations: EQ-003, EQ-004, EQ-005, EQ-006, EQ-007
  - File: `module3_a_maximization.md`
- Large N limit: N -> infinity with fixed matter content labels. Leading-order analysis in 1/N.
  - Applicable equations: all, evaluated at leading order
  - File: `a_maximization_large_N.py`

**Phase Transitions / Crossovers:**

- Conformal window boundaries:
  - Upper boundary: b_0 = 0 (loss of asymptotic freedom). N_f^{upper} determined by AF condition.
  - Lower boundary: unitarity bound R[O] = 2/3 saturated, or B_a = 0 (loss of IR fixed point). N_f^{lower} determined by a-maximization + boundary analysis.
  - File: `module4_ir_fixed_points.md` (Sec. 7); `module5_classification.md` (Sec. 4a)
- Free magnetic phase (below conformal window): Seiberg duality maps to a dual description.
  - File: `module4_ir_fixed_points.md` (Sec. 9)
- Fixed point merging/annihilation: as parameters vary, fixed points can collide and disappear.
  - File: `module4_ir_fixed_points.md` (Sec. 12)

**Known Limiting Cases:**

- SU(N) SQCD: conformal window 3N/2 < N_f < 3N. R^* = 1 - N/N_f. Verified in `module3_a_maximization.md` (Sec. 10).
- Pure SU(N)^k circular quiver (no node matter, no flavors): R_Q = R_Qtilde = 0, no interacting SCFT. Verified in `module3_a_maximization.md` (Sec. 11).
- N=4 SYM (SU(N) + 1 adj, equivalent to 3 adjoints in N=1 language): a = c (not explicitly computed in project).
- Kutasov-Schwimmer models: SU(N) with one antisymmetric + flavors. Referenced in `module5_classification.md` (Sec. 4b).

## Units and Conventions

**Unit System:**

- Natural units (hbar = c = 1) throughout. No explicit factors of hbar or c appear anywhere.
- No metric signature is relevant (analysis is algebraic, not spacetime-dependent).
- All quantities in the equation catalog are dimensionless.

**Key Conventions:**

- Dynkin index normalization: T(fund_{SU(N)}) = 1/2, T(V_{SO(N)}) = 1, T(f_{Sp(N)}) = 1/2. This is stated in `CLAUDE.md` and used consistently throughout.
  - File: `CLAUDE.md` (Physics Conventions); `module1_beta_functions.md` (Sec. 2)
- USp(2N) written as Sp(N): the rank parameter N means USp(2N) has rank N, fundamental dimension 2N.
  - File: `CLAUDE.md` (Physics Conventions)
- Representation notation: A for rank-2 antisymmetric (never Lambda or Lambda^2); S for rank-2 symmetric. Sbar, Abar for conjugates.
  - File: `CLAUDE.md` (Physics Conventions)
- Code representation labels: "fund", "antifund", "adj", "S", "Sbar", "A", "Abar" for SU(N); "V", "adj", "S" for SO(N); "fund", "adj", "A" for Sp(N).
  - File: `beta_functions.py` (`VALID_REPS`, line 19)
- N_MIN: SU requires N >= 5, SO requires N >= 5, Sp requires N >= 3 (in code; differs from mathematical minimum to avoid degenerate small-N cases).
  - File: `beta_functions.py` (line 16)
- All nodes share the same N (equal-rank assumption).
  - File: `README.md`; `CLAUDE.md`
- R-charge convention: R[gaugino] = 1, R[chiral superfield Phi_i] is the variable, fermion in Phi_i has R-charge R_i - 1.
  - File: `module3_a_maximization.md` (Sec. 2)
- Unitarity bound applies to gauge-INVARIANT chiral operators only (R[O] >= 2/3), NOT to elementary gauge-charged fields.
  - File: `module3_a_maximization.md` (Sec. 9); `CLAUDE.md`

## Dynkin Index Tables

These tables are central to all computations and must be preserved exactly.

### SU(N) Representations

| Rep | Symbol | dim | T(r) | A(r) |
|-----|--------|-----|------|------|
| Fundamental | fund | N | 1/2 | +1 |
| Anti-fundamental | antifund | N | 1/2 | -1 |
| Adjoint | adj | N^2 - 1 | N | 0 |
| Rank-2 symmetric | S | N(N+1)/2 | (N+2)/2 | +(N+4) |
| Rank-2 sym conjugate | Sbar | N(N+1)/2 | (N+2)/2 | -(N+4) |
| Rank-2 antisymmetric | A | N(N-1)/2 | (N-2)/2 | +(N-4) |
| Rank-2 antisym conjugate | Abar | N(N-1)/2 | (N-2)/2 | -(N-4) |

Source: `module1_beta_functions.md` (Sec. 2); `beta_functions.py` (`T_rep`, `A_rep`)

### SO(N) Representations

| Rep | Symbol | dim | T(r) |
|-----|--------|-----|------|
| Vector | V | N | 1 |
| Adjoint | adj | N(N-1)/2 | N-2 |
| Traceless symmetric | S | N(N+1)/2 - 1 | N+2 |

Source: `module1_beta_functions.md` (Sec. 2)

### Sp(N) = USp(2N) Representations

| Rep | Symbol | dim | T(r) |
|-----|--------|-----|------|
| Fundamental | f | 2N | 1/2 |
| Adjoint (= rank-2 sym) | adj | N(2N+1) | N+1 |
| Traceless antisymmetric | A | N(2N-1) - 1 | N-1 |

Source: `module1_beta_functions.md` (Sec. 2)

### Bifundamental Contributions to b_0

| Node a | Node b | T_{G_a} per edge | T_{G_b} per edge |
|--------|--------|-------------------|-------------------|
| SU(N) | SU(N) | N/2 | N/2 |
| SU(N) | SO(N) | N/2 | N |
| SU(N) | Sp(N) | N | N/2 |
| SO(N) | SO(N) | N | N |
| SO(N) | Sp(N) | 2N | N/2 |
| Sp(N) | Sp(N) | N | N |

Source: `module1_beta_functions.md` (Sec. 4); `module2_quiver_generation.md` (Sec. 3)

### Large N Leading Coefficients

At large N, representations become degenerate:

| Gauge | Rep | T_lead = lim T/N | dim_lead = lim dim/N^2 |
|-------|-----|-------------------|------------------------|
| SU | fund/antifund | 0 | 0 |
| SU | adj | 1 | 1 |
| SU | S, Sbar, A, Abar | 1/2 | 1/2 |
| SO | V | 0 | 0 |
| SO | adj | 1 | 1/2 |
| SO | S | 1 | 1/2 |
| Sp | fund | 0 | 0 |
| Sp | adj | 1 | 2 |
| Sp | A | 1 | 2 |

Source: `a_maximization_large_N.py` (lines 42-80)

Key universality: all SU(N) rank-2 reps (S, Sbar, A, Abar) have identical T_lead = 1/2 and dim_lead = 1/2. Theories differing only in rank-2 composition are in the same universality class at leading order.

---

_Framework analysis: 2026-04-14_
