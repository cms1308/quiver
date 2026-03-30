# Module 5: Final Classification

## Goal
Combine all preceding modules to produce the complete list of 4D N=1 quiver gauge theories that are:
1. UV asymptotically free (Module 1, 2)
2. Admit a large N limit (Module 2)
3. Satisfy all consistency conditions (anomaly cancellation, Witten anomaly) (Module 2)
4. Flow to a non-trivial IR superconformal fixed point (Modules 3, 4)

---

## 1. Classification Pipeline

```
┌─────────────────────────────────────────┐
│  Enumerate quivers (Module 2)           │
│  - Graph topologies                     │
│  - Gauge group assignments              │
│  - Matter content                       │
└────────────────┬────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────┐
│  Filter: UV asymptotic freedom (Mod 1)  │
│  b_0^(a) > 0 for all nodes a           │
│  B_0^(a) > 0 at large N               │
└────────────────┬────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────┐
│  Filter: Consistency conditions (Mod 2) │
│  - SU(N) gauge anomaly = 0             │
│  - Sp(N) Witten condition              │
└────────────────┬────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────┐
│  a-maximization (Module 3)              │
│  - Compute R_i*, a*, c*                │
│  - Check unitarity bounds              │
└────────────────┬────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────┐
│  IR fixed point check (Module 4)        │
│  - Boundary analysis B_a < 0 for all a │
└────────────────┬────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────┐
│  CLASSIFIED THEORIES                    │
└─────────────────────────────────────────┘
```

---

## 2. Classification by Quiver Topology

We organize results by the structure of the quiver graph.

### 2a. Single-Node Theories

| Gauge Group | Matter Content | UV AF Condition | IR FP Condition | Notes |
|-------------|----------------|-----------------|-----------------|-------|
| SU(N) | $N_f$ fund + $\bar{N}_f$ anti-fund | $N_f < 3N$ | $N_f > \frac{3N}{2}$ | SQCD conformal window |
| SU(N) | 1 adjoint + $N_f$ fund + $\bar{N}_f$ anti-fund | $N_f < 2N$ | $N_f > N$ | $\mathcal{N}=2$ SQCD at $N_f = 2N$ marginal |
| SU(N) | 1 sym + $N_f$ fund | check AF + anomaly | from bdry analysis | Kutasov-type |
| SU(N) | 1 antisym + $N_f$ fund | check AF + anomaly | from bdry analysis | Kutasov-type |
| SO(N) | $N_f$ vectors | $N_f < 3(N-2)$ | $N_f > \frac{3(N-2)}{2}$... | Intriligator-Seiberg |
| Sp(N) | $2N_f$ fundamentals (even) | $2N_f < 3(N+1)$ | from bdry analysis | Intriligator-Seiberg |

### 2b. Two-Node Quivers

Key examples to classify:

#### SU(N) × SU(N)
- 1 bifundamental $(\square, \bar{\square})$ + 1 bifundamental $(\bar{\square}, \square)$ (circular quiver with 2 nodes)
- Extra node matter at each node

| Node 1 | Node 2 | Bifundamentals | Extra matter | Anomaly | IR FP? |
|--------|--------|----------------|--------------|---------|--------|
| SU(N) | SU(N) | $Q, \tilde{Q}$ (1 each) | — | ✓ balanced | Check $\mathcal{B}$ |
| SU(N) | SU(N) | $n$ bifund. | $n_f$ at end | balanced if $n$ symmetric | Check |

#### SU(N) × SO(N)
- Bifundamental in $(\square, V)$
- Need even number of bifundamentals if Sp is involved

#### SO(N) × Sp(N)
- Special care needed for large N scaling (SO-Sp bif. very costly)

### 2c. Linear Quivers (Chain)

$G_1 - G_2 - \cdots - G_k$ with bifundamentals on each link.

**Pure SU(N)$^k$ linear quiver:**
- End nodes: $b_0 = 5N/2$ ✓
- Interior nodes: $b_0 = 2N$ ✓
- Anomaly: each interior node has equal in/out arrows ✓
- Valid for any $k$ at large N

**With additional node matter:**
- Add adjoints, symmetric, antisymmetric at specific nodes to tune conformal window

### 2d. Circular (Necklace) Quivers

$G_1 - G_2 - \cdots - G_k - G_1$ (cyclic).

**Pure SU(N)$^k$ circular quiver (affine $A_{k-1}$ quiver):**
- All nodes: $b_0 = 2N > 0$ ✓
- Balanced at each node ✓
- Valid for any $k$ ≥ 2

At $k=1$ with 1 adjoint: $b_0 = 3N - N = 2N$ (pure $\mathcal{N}=2$ SYM, marginal not AF).

**SU(N)$^k$ with one adjoint per node (affine $\hat{A}_{k-1}$, $\mathcal{N}=2$ quiver):**
- Each node: $b_0 = 3N - N - 2 \cdot N/2 = N > 0$... wait, this needs careful counting.

### 2e. Star (Hub-and-Spoke) Quivers

One central node connected to several outer nodes.

**SU(N) hub with $k$ SU(N) spokes:**
- Hub node: $b_0 = 3N - k \cdot N/2$; requires $k \leq 5$
- Spoke nodes: $b_0 = 3N - N/2 = 5N/2 > 0$ ✓

For $k = 5$: hub has $b_0 = N/2 > 0$ ✓ (marginally AF at large N; finite N corrections matter)

---

## 3. Large N Admissible Quivers — Catalog

We organize by the **large N leading beta function coefficient** $B_0^{(a)} = \lim_{N\to\infty} b_0^{(a)}/N$.

### Requirement: $B_0^{(a)} > 0$ for all nodes

For **SU(N) nodes** with $n_b$ bifundamental neighbors (all SU-SU edges):
$$B_0 = 3 - n_{\mathrm{adj}} - \frac{n_S + n_A}{2} - \frac{n_b}{2}$$

| $n_{\mathrm{adj}}$ | $n_S + n_A$ | $n_b^{\max}$ |
|-------------------|-------------|--------------|
| 0 | 0 | 5 |
| 0 | 1 | 3 |
| 0 | 2 | 1 |
| 1 | 0 | 3 |
| 1 | 1 | 1 |
| 2 | 0 | 1 |

---

## 4. Theories Known to Have Non-Trivial IR Fixed Points

### 4a. $\mathcal{N}=1$ SQCD (SU, SO, Sp variants)
- Seiberg's conformal window: $\frac{3}{2}N < N_f < 3N$ for SU(N)

### 4b. Kutasov–Schwimmer Models
- SU(N) with one rank-2 antisymmetric + $N_f$ fundamentals
- Conformal window: $N+2 \leq N_f \leq \frac{3}{2}(N-2)$ (depends on exact conventions)

### 4c. Circular SU(N)$^k$ Quivers
- All nodes balanced, always UV-AF for any $k$
- IR analysis: $k=2$ (two nodes) is equivalent to SU(N) SQCD-like analysis

### 4d. Linear SU(N)$^k$ Quivers with Flavors
- Adding $N_f$ fundamental flavors at end nodes extends the conformal window

### 4e. Mixed Gauge Group Quivers
- SO(N)–Sp(N) pairs: special care needed due to large Dynkin indices for bifundamentals
- SU(N)–SO(N) and SU(N)–Sp(N) quivers: more natural candidates

---

## 5. Output Format for Each Classified Theory

For each valid theory, record:

```
Theory ID: [unique identifier]
Quiver type: [linear / circular / star / other]
Nodes: [(G_1, type_1), (G_2, type_2), ...]
Edges: [(a, b, rep_ab), ...]
Node matter: {a: [rep_1, rep_2, ...], ...}
Superpotential: [W = ...]

UV properties:
  b_0^(a) = [...]  for each a
  B_0^(a) = [...]  (large N coefficient)
  Anomaly cancellation: [✓/✗]

IR properties:
  R_charges: {field_i: R_i*, ...}
  a* = [...]
  c* = [...]
  a/c = [...]
  B_a = [...]  (boundary analysis values, all < 0 ✓)

Notes: [special properties, known dual, etc.]
```

---

## 6. Expected Results by Quiver Size

### 1-node theories
- Well-understood; reproduce known SQCD and Kutasov-type windows
- Several infinite families parametrized by $N_f$

### 2-node theories
- Richer structure; coupling between conformal windows of two nodes
- Expected: several families, some with enhanced symmetry in IR

### $k$-node circular quivers (affine $A_{k-1}$)
- Always UV-AF at large N
- IR analysis: need to check that all nodes are in their conformal window simultaneously
- These are related to $\mathcal{N}=2$ theories and orbifold field theories

### $k$-node linear quivers
- End nodes have larger $b_0$; interior nodes are more constrained
- Adding matter at ends shifts the conformal window boundaries

---

## 7. Special Cases and Marginal Theories

Theories with $B_0^{(a)} = 0$ for some node $a$ are **marginally AF at large N**: they are asymptotically free for any finite $N$ but the AF condition becomes marginal as $N \to \infty$.

These include:
- SU(N) with 6 bifundamental SU-SU neighbors: $B_0 = 3 - 3 = 0$
- SO(N) node with 3 SO-SU bifundamentals: $B_0 = 3 - 3 \cdot 1 = 0$ (needs careful $\mathcal{O}(1)$ analysis)

For these, the finite-N corrections to $b_0$ determine asymptotic freedom. They may or may not be included depending on whether one requires strict $b_0 > 0$ or $b_0 \geq 0$.

---

## 8. Generalization Directions

The current classification assumes all gauge groups have the same N. Natural generalizations:

1. **Different ranks**: $G_a = SU(N_a)$ with $N_a$ varying — allows for "quiver tails" with decreasing ranks
2. **Framing**: adding decoupled flavor nodes (non-gauged SU($N_f$)) at various nodes
3. **Superpotential deformations**: adding relevant/marginal operators that modify R-charges and conformal window
4. **Global symmetry structure**: classifying the flavor symmetries and anomalies of the classified theories

---

## 9. Implementation Summary

```python
def classify_quivers(max_nodes=4, N_range=(3, 20)):
    """Full classification pipeline."""

    results = []

    # Step 1: Generate candidate quivers
    candidates = generate_quivers(max_nodes)

    for quiver in candidates:
        # Step 2: Check UV asymptotic freedom
        if not check_uv_af(quiver):
            continue

        # Step 3: Check large N admissibility
        if not check_large_N_af(quiver):
            continue

        # Step 4: Check anomaly cancellation
        if not check_anomalies(quiver):
            continue

        # Step 5: a-maximization
        try:
            R_star, a_star, c_star = a_maximize(quiver)
        except NoFixedPoint:
            continue

        # Step 6: Boundary analysis for IR fixed point
        if not check_ir_fixed_point(quiver, R_star):
            continue

        results.append({
            'quiver': quiver,
            'R_charges': R_star,
            'a': a_star,
            'c': c_star,
        })

    return results
```

---

## 10. Cross-Checks and Validation

- **$\mathcal{N}=2$ theories**: Setting the adjoint chiral R-charge to $R = 2/3$ should recover $\mathcal{N}=2$ results
- **Known dualities**: Seiberg duals of classified theories should also appear in the classification (or relate to theories outside the classification via relevant deformations)
- **Unitarity**: All classified theories must satisfy $R_i^* \geq 2/3$
- **$a$-theorem consistency**: $a_{\mathrm{UV}} \geq a_{\mathrm{IR}}$ for all theories in the classification (Komargodski-Schwimmer)
