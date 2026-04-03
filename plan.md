# Plan: Module 3 — a-maximization (`a_maximization.py`)

## Context

Module 3 takes a quiver (from Module 2) at specific (N, N_f) and determines the exact superconformal IR R-charges by maximizing the Intriligator–Wecht a-function subject to anomaly-free constraints.

---

## Files

- **New file:** `a_maximization.py`
- **Reuse from `beta_functions.py`:** `T_adj`, `T_rep`, `T_bifund`, `N_MIN`
- **Reuse from `quiver_generation.py`:** `Quiver`, `Edge`, `chiral_excess_coeffs`

---

## 1. Helper functions: dimensions

```python
def dim_group(gauge_type: GaugeType, N: int) -> int:
    """Dimension of gauge group (number of gauginos)."""
    # SU(N): N²-1
    # SO(N): N(N-1)/2
    # Sp(N) = USp(2N): N(2N+1)

def dim_rep(gauge_type: GaugeType, rep: str, N: int) -> int:
    """Dimension of a single-node representation."""
    # SU: fund→N, antifund→N, adj→N²-1,
    #     S→N(N+1)/2, Sbar→N(N+1)/2, A→N(N-1)/2, Abar→N(N-1)/2
    # SO: V→N, adj→N(N-1)/2, S→N(N+1)/2 - 1   (traceless symmetric)
    # Sp: fund→2N, adj→N(2N+1), A→N(2N-1) - 1  (traceless antisymmetric)
```

---

## 2. ChiralField dataclass

Represents one "type" of chiral multiplet (all copies share the same R-charge):

```python
@dataclass
class ChiralField:
    label: str                           # "node0_adj", "node1_fund", "edge_0_1_pm", ...
    R_index: int                         # index into the R-charge vector
    dim: int                             # total dim (product over all gauge groups)
    T_contributions: dict[int, Fraction] # node_index → effective T at that node
```

- **node-level rep** at node `i` with count `n`:
  - `dim = n * dim_rep(g, rep, N)`
  - `T_contributions = {i: n * T_rep(g, rep, N)}`
- **bifundamental edge** (nodes `i`, `j`):
  - `dim = dim(r_i) * dim(r_j)` (from `dim_rep` at each endpoint)
  - `T_contributions = {i: T_i, j: T_j}` (from `T_bifund`)

---

## 3. `build_fields(quiver, N, N_f) -> list[ChiralField]`

Converts a Module 2 quiver into the full list of chiral fields at (N, N_f):

1. **Rank-2/adj matter** from `node_matter`: one field per rep type per node.
2. **Fund-like matter** derived from N_f + chiral excess δ:
   - SU node: `δ = chiral_excess_coeffs(quiver, node)`, then `n_f = N_f + max(δ, 0)`, `n_fbar = N_f + max(-δ, 0)`. Separate fields for fund and antifund if nonzero.
   - SO node: `n_V = N_f` vectors.
   - Sp node: `n_f = 2*N_f` fundamentals.
3. **Bifundamental edges**: one field per edge.

Each field type gets a unique `R_index`. Distinct copies of the same representation are grouped into one field (scaled dim and T).

---

## 4. Anomaly constraint matrix

```python
def build_constraint_matrix(fields, quiver, N) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns (A, b) with shape (n_nodes, n_R) and (n_nodes,) such that:
        A @ R = b
    encodes T(adj_a) + Σ_i T_a(i) * (R_i - 1) = 0 for each node a.

    A[a, i] = T_contributions[a] of field i (0 if not charged under node a)
    b[a]    = Σ_i T_a(i) - T_adj(g_a, N)
    """
```

---

## 5. Trial central charges

```python
def tr_R(R_values, fields, quiver, N) -> float:
    """Tr R = Σ_a dim(G_a) + Σ_fields dim * (R - 1)"""

def tr_R3(R_values, fields, quiver, N) -> float:
    """Tr R³ = Σ_a dim(G_a) + Σ_fields dim * (R - 1)³"""

def a_trial(R_values, fields, quiver, N) -> float:
    """a = (3/32)(3 Tr R³ - Tr R)"""

def c_trial(R_values, fields, quiver, N) -> float:
    """c = (1/32)(9 Tr R³ - 5 Tr R)"""
```

---

## 6. `a_maximize(quiver, N, N_f) -> AMaxResult`

```python
@dataclass
class AMaxResult:
    R_charges: dict[str, float]  # field_label → R-charge value
    a: float
    c: float
    unitarity_ok: bool           # all R_i >= 2/3
```

**Algorithm:**
1. `fields = build_fields(quiver, N, N_f)`
2. `A, b = build_constraint_matrix(fields, quiver, N)`
3. Particular solution `R0` via `np.linalg.lstsq(A, b)`
4. Null space `F` (columns = free flavor directions) via `scipy.linalg.null_space(A)`
5. `R(s) = R0 + F @ s`
6. Maximize `a_trial(R0 + F @ s, ...)` over `s` using `scipy.optimize.minimize`
7. Compute a*, c* at optimal s*
8. Check R_i ≥ 2/3

---

## 7. Unitarity decoupling (`a_maximize_with_decoupling`)

If any R_i < 2/3 at the solution:
- Pin that field to R = 2/3 (free field)
- Add as a fixed constraint and redo a-maximization
- Iterate until all R_i ≥ 2/3

---

## 8. Validation (`__main__`)

```python
# SU(N) SQCD at N=5, N_f=10: expect R_fund = R_antifund = 1 - N/(2*N_f) = 0.75
# Circular SU(N)² (balanced +-) at N=5: expect R_bifund = 1 (spec §11)
# SU(N) + 1 adj + N_f flavors: nontrivial maximization with one free parameter
```

---

## Open questions for the user

1. **Superpotential**: the current framework has no superpotential. If we add one later (e.g. W = Tr Φ Q Q̃), it adds extra linear constraints R[Φ] + R[Q] + R[Q̃] = 2. Should this be included now or deferred?

2. **Bifundamental R-charges**: for SU-SU "++" edges (□,□), should the two endpoint R-charges be treated as distinct or equal? (They transform under different gauge groups so are generically different.)

3. **Output format**: return a `dict[str, float]` keyed by field label, or a structured object keyed by node/edge?
