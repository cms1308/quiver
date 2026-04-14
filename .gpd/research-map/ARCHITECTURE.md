# Computational Architecture

**Analysis Date:** 2026-04-14

## Mathematical Setting

**Spaces:**

- Configuration space: discrete set of quiver gauge theories Q = (V, E, gauge_types, node_matter, rank_multipliers)
  - File: `quiver_generation.py` (lines 66-152, `Quiver` dataclass)
- R-charge space: real vector space R^{n_R} where n_R = number of distinct chiral field types
  - Constrained to an affine subspace by anomaly-free conditions (linear in R_i)
  - File: `a_maximization.py` (lines 202-234, `anomaly_matrix`)
- Central charge surface: trial a-function a(R) is a cubic polynomial on the anomaly-free subspace
  - File: `a_maximization.py` (lines 239-255, `a_trial`, `c_trial`)

**Key Mathematical Objects:**

| Object | Type | Symbol | Defined In |
|--------|------|--------|------------|
| One-loop beta coefficient | Fraction (exact rational) | b_0 | `beta_functions.py` (line 177, `compute_b0`) |
| Dynkin index | Fraction | T(r) | `beta_functions.py` (lines 35-71, `T_adj`, `T_rep`) |
| Anomaly matrix | numpy array (n_nodes x n_R) | A | `a_maximization.py` (lines 203-234, `anomaly_matrix`) |
| Trial a-function | float or sympy Expr | a(R) = (3/32)(3 Tr R^3 - Tr R) | `a_maximization.py` (line 253) |
| Trial c-function | float or sympy Expr | c(R) = (1/32)(9 Tr R^3 - 5 Tr R) | `a_maximization.py` (line 260) |
| Leading-order field descriptor | LeadField dataclass | (label, R_index, dim_lead, T_lead) | `a_maximization_large_N.py` (lines 160-167) |
| Large-N result | LargeNResult dataclass | (R_charges, a/N^2, c/N^2, unitarity_ok) | `a_maximization_large_N.py` (lines 400-425) |
| Quiver | Quiver dataclass | (gauge_types, edges, node_matter, rank_multipliers) | `quiver_generation.py` (lines 66-152) |

## Computational Pipeline

The project implements a five-module pipeline. Modules 1-3 and 5 are implemented in Python; Module 4 is not yet implemented.

**Full Pipeline:**

1. **Quiver Enumeration (Module 2):** Generate all valid 2-node quiver gauge theories
   - Input: gauge group types {SU, SO, Sp}, max multiedge count, rank multipliers
   - Output: list of `(Quiver, bounds)` tuples, deduplicated by charge conjugation and node swap
   - Script: `quiver_generation.py` (lines 412-468, `enumerate_quivers`; lines 471-534, `enumerate_quivers_mixed_rank`)
   - Deduplication: `_dedup_symmetries` (lines 614-642) removes charge-conjugate and node-swap duplicates

2. **Asymptotic Freedom Check (Module 1):** Verify b_0 > 0 for all N at each node
   - Input: `NodeSpec` (gauge_type, N, matter, bifund_neighbors)
   - Output: boolean (AF/not AF) and linear coefficients (alpha, beta) of b_0(N)
   - Script: `beta_functions.py` (lines 177-257)
   - Key: `b0_linear` evaluates b_0 at two N values and extracts exact rational slope/intercept

3. **a-Maximization (Module 3):** Determine IR R-charges by maximizing trial a-function
   - Two implementations: numerical (fixed N, N_f) and exact symbolic (large N)
   - Script: `a_maximization.py` (numerical), `a_maximization_large_N.py` (symbolic + Mathematica batch)

4. **IR Fixed Point Analysis (Module 4):** Check B_a < 0 boundary condition
   - NOT YET IMPLEMENTED
   - Specification: `module4_ir_fixed_points.md`

5. **Classification Database (Module 5):** Build SQLite DB of universality classes
   - Input: enumerated quivers from step 1
   - Output: `quivers.db` (SQLite, ~2.5 MB)
   - Script: `two_node_db.py` (line 341, `cmd_build`)
   - Three-phase build: enumerate + Mathematica NSolve scan -> cluster into classes -> exact sympy per class

**Database Build Pipeline (detailed):**

```
enumerate_quivers(n_nodes=2, max_multiedge=4)         # equal-rank
enumerate_quivers_mixed_rank([2,1])                    # mixed-rank
enumerate_quivers_mixed_rank([1,2])                    # mixed-rank
         |
         v
    _dedup_symmetries()                                # cross-set dedup
         |
         v
    a_maximize_batch_mathematica()                     # Phase 1: Mathematica NSolve
         |                                             # parallel kernels, WorkingPrecision->30
         v
    _cluster(a_vals, TOL=1e-5)                         # Phase 2: cluster by a/N^2
         |
         v
    _exact_with_timeout(rep_quiver, timeout=30s)       # Phase 3: exact sympy per class
         |
         v
    SQLite write (universality_class, theory,          # Phase 4: persist
                  morphology_class, build_info)
```

## Key Algorithms

### Exact Arithmetic Throughout Module 1-2

All beta function and anomaly calculations use `fractions.Fraction` for exact rational arithmetic. No floating-point until a-maximization optimization.
- Implementation: `beta_functions.py` (entire file uses `Fraction`)
- Significance: eliminates rounding errors in AF bounds and chiral excess coefficients

### Backtracking Enumeration with AF Pruning

Node matter enumeration uses monotonic pruning: since each additional representation decreases b_0, if count=k violates AF then all count > k do too.
- Implementation: `quiver_generation.py` (lines 308-343, `_enumerate_node_matter`)
- Complexity: effectively bounded by the small number of valid matter assignments per node (typically < 10)

### Null-Space Parameterization for Constrained Optimization

a-maximization solves `A @ R = b` (anomaly constraints) then parameterizes the solution space as `R = R0 + F @ s` where F = null_space(A).
- Numerical: `a_maximization.py` (lines 350-391) uses `np.linalg.lstsq` for R0, `scipy.linalg.null_space` for F, then `scipy.optimize.minimize(method="BFGS")` over free parameters s
- Symbolic: `a_maximization_large_N.py` (lines 548-673) uses `sympy.Matrix.nullspace()`, `sympy.solve()` on gradient equations (quadratic in s)
- Mathematica batch: `a_maximization_large_N.py` (lines 716-914) uses `LinearSolve`, `NullSpace`, `NSolve` with `WorkingPrecision->30`

### Iterative Unitarity Decoupling

When a gauge-invariant operator O has R[O] < 2/3, the constraint R[O] = 2/3 is added and a-maximization is redone. Iterates until all operators satisfy unitarity.
- Numerical: `a_maximization.py` (lines 394-462, `a_maximize_with_decoupling`)
- Symbolic: `a_maximization_large_N.py` (lines 461-503, `a_maximize_large_N_with_decoupling`)

### Universality Class Clustering

Theories are grouped into universality classes by clustering their a/N^2 values within tolerance TOL = 1e-5.
- Implementation: `two_node_db.py` (lines 214-235, `_cluster`)
- Algorithm: sort a-values, greedily merge adjacent entries within tolerance

### Morphology Classification

Theories are additionally classified by a 5-tuple morphology vector: (N_rank2_0, N_rank2_1, N_bif, N_fund_0, N_fund_1) capturing the matter content structure independent of specific representation choices.
- Implementation: `two_node_db.py` (lines 260-338)

## Solver Stack

### Python Standard Library
- `fractions.Fraction`: exact rational arithmetic for all Dynkin indices, anomaly coefficients, beta functions
- `itertools.product`: combinatorial enumeration of gauge types, edges, matter
- `dataclasses`: all data structures (NodeSpec, Quiver, Edge, ChiralField, LeadField, etc.)
- `sqlite3`: database storage and queries

### NumPy
- `np.linalg.lstsq`: particular solution of anomaly constraint system
- Array operations for R-charge vector manipulation
- Used in: `a_maximization.py`

### SciPy
- `scipy.linalg.null_space`: null space of anomaly constraint matrix (SVD-based)
- `scipy.optimize.minimize(method="BFGS")`: unconstrained optimization of trial a-function over free parameters
- Used in: `a_maximization.py`

### SymPy
- `sympy.Matrix.nullspace()`: exact rational null space computation
- `sympy.solve()`: exact solution of gradient equations (polynomial system)
- `sympy.diff()`: symbolic differentiation of cubic a-function
- `sympy.simplify()`: simplification of exact algebraic R-charges
- `sympy.Rational`: exact rational number representation
- Used in: `a_maximization_large_N.py`

### Mathematica (via wolframscript)
- `NullSpace`, `LinearSolve`: exact anomaly constraint handling
- `NSolve[..., Reals, WorkingPrecision->30]`: high-precision numerical root finding for gradient equations
- `NegativeDefiniteMatrixQ`: Hessian check to select local maxima
- `ParallelMap` + `LaunchKernels[]`: parallel processing across quivers
- Invoked as external subprocess: `subprocess.run(["wolframscript", "-file", script_path])`
- Used in: `a_maximization_large_N.py` (lines 716-914, `a_maximize_batch_mathematica`)
- Timeout: 7200s default for batch, 30s per individual exact computation

## Performance Characteristics

**Enumeration (Module 2):**
- Equal-rank 2-node quivers with max_multiedge=4: generates ~O(10^3-10^4) theories
- Mixed-rank adds [2,1] and [1,2] variants, cross-set deduplication reduces count
- Current DB: 7757 total theories (7537 converged + 220 diverged)

**Mathematica Batch NSolve (Phase 1):**
- Single wolframscript process handles entire batch (avoids per-theory startup overhead)
- Parallelized via Mathematica's `ParallelMap` + `LaunchKernels[]`
- Total batch timeout: 7200s (2 hours)
- Working precision: 30 digits

**Exact Symbolic (Phase 3):**
- Per-class timeout: 30s (configurable via `--exact-timeout`)
- Success rate: 182/326 = 55.8% of classes get exact solutions
- 144 classes fall back to numerical-only results
- Signal-based timeout via `signal.SIGALRM` (Unix only)

**Database:**
- SQLite file: ~2.5 MB for 7757 theories, 326 classes, 580 morphology classes
- Indexed on: class_id, gauge_pair, a_over_N2, morph_id

## Data Flow

```
Input parameters                    Output
─────────────                       ──────
gauge_types: [SU,SO,Sp]^n    ──>   Quiver objects
max_multiedge: int            ──>   (with AF bounds)
rank_multipliers: [int]
                                      |
                                      v
                                    Mathematica NSolve
                                    (wolframscript subprocess)
                                      |
                                      v
                                    FastNumericalResult per quiver
                                    (a/N^2, c/N^2, R-charges)
                                      |
                                      v
                                    Clustering by a/N^2 (TOL=1e-5)
                                      |
                                      v
                                    SymPy exact per class representative
                                    (LargeNResult with exact algebraic R-charges)
                                      |
                                      v
                                    quivers.db (SQLite)
                                    tables: theory, universality_class,
                                            morphology_class, build_info
```

## External Dependencies

- **Mathematica/wolframscript**: REQUIRED for batch NSolve in `a_maximize_batch_mathematica`. Without it, Phase 1 of the DB build cannot run. Must be on PATH.
- **No requirements.txt or pyproject.toml**: dependencies are implicit (numpy, scipy, sympy, plus Mathematica)
- **No virtual environment configuration**: scripts assume system Python with scientific packages installed

## Transformation Properties

**Charge Conjugation Symmetry:**
- SU representations: S <-> Sbar, A <-> Abar, fund <-> antifund, adj self-conjugate
- Edge representations: +- reverses direction, ++ <-> --
- Implementation: `quiver_generation.py` (lines 537-582, `_conjugate_matter`, `_conjugate_signature`)
- Used for deduplication: physically equivalent quivers related by charge conjugation are identified

**Node Swap Symmetry (2-node):**
- Swap node 0 <-> node 1: exchanges gauge types, matter, rank multipliers, and edge endpoints
- Implementation: `quiver_generation.py` (lines 598-611, `_node_swap_signature`)
- Combined with conjugation for full symmetry group of 2-node quivers

## Boundary and Initial Conditions

**Asymptotic Freedom Boundary:**
- b_0 > 0 at each node, for all N >= N_MIN (where N_MIN = {SU:5, SO:5, Sp:3})
- Equivalently: alpha > 0 or (alpha == 0 and beta >= 0) in b_0 = alpha*N + beta
- File: `beta_functions.py` (lines 241-256, `is_af_all_N`)

**Unitarity Bound:**
- All gauge-invariant operators must have R[O] >= 2/3
- Applied to composite operators (mesons, bilinears), NOT individual fields
- File: `a_maximization.py` (lines 266-335, `gauge_invariant_ops`)

**Witten Anomaly (Sp nodes):**
- Total number of fundamental chiral multiplets at each Sp(N) node must be even
- File: `quiver_generation.py` (lines 287-295, `check_sp_witten`)

---

*Architecture analysis: 2026-04-14*
