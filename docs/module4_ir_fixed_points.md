# Module 4: Non-Trivial IR Fixed Point Analysis

## Goal
For each UV-AF quiver that passed Module 2, determine whether it flows to a **non-trivial interacting IR superconformal fixed point** (as opposed to being free or trivially gapped). We follow the boundary analysis method of hep-th/0502049.

---

## 1. The Problem: Existence of an IR Fixed Point

A theory with gauge couplings $g_1, g_2, \ldots, g_k$ (one per gauge node) has a non-trivial IR fixed point if there exist **simultaneously** non-zero fixed-point couplings $g_a^* \neq 0$ for all $a$ such that all beta functions vanish:
$$\beta_{g_a}(g_1^*, \ldots, g_k^*) = 0 \quad \forall a$$

Simply having UV asymptotic freedom ($b_0^{(a)} > 0$) is necessary but **not sufficient**. One must also rule out the possibilities:
- The theory flows to a free IR fixed point (all $g_a^* = 0$)
- Some gauge couplings decouple and the theory splits into sub-sectors
- The theory is in the free magnetic phase (if a dual description exists)

---

## 2. The Boundary Analysis (hep-th/0502049)

Consider a quiver with two gauge groups $G$ and $G'$ with couplings $g$ and $g'$.

### Key Idea
The coupling space $(g, g')$ has boundaries:
- $(g = 0, g')$: gauge group $G$ is weakly coupled/free; $G'$ is a standalone gauge theory with the matter charged under $G$ treated as **free fields** with canonical R-charges
- $(g, g' = 0)$: vice versa

For a non-trivial IR fixed point $(g^*, g'^*) \neq (0,0)$ to exist, **both gauge groups must want to be interacting** near the free boundary. Concretely:

> **Criterion:** At each boundary where $g_a = 0$ (with all other couplings at their fixed-point values), the remaining beta function $\beta_{g_a}$ must be negative (i.e., $g_a$ is driven to grow). Equivalently, the one-loop Tr$[R G_a^2]$ must be negative.

---

## 3. The NSVZ Beta Function Numerator

The NSVZ beta function for coupling $g_a$ is:
$$\beta(g_a) \propto -\left[3T(\mathrm{adj}_a) - \sum_i T_{G_a}(r_i)(1-2\gamma_i)\right]$$

(negative sign so that $\beta < 0$ means the coupling runs to larger values = asymptotically free/non-trivial fixed point).

At a **boundary fixed point** where $g_a \to 0$ and $g_b \to g_b^*$:
- Matter fields charged only under $G_a$ become free: $\gamma_i = 0$, $R_i = R_i^{\mathrm{free}}$
- Matter fields charged under $G_b$ (and possibly $G_a$) sit at the $G_b$ fixed point: their anomalous dimensions are determined by $g_b^*$

The condition for $G_a$ to be driven away from $g_a = 0$ is:
$$\mathrm{Tr}[R_{\mathrm{bdy}} G_a^2] < 0$$

where $R_{\mathrm{bdy}}$ is the R-charge at the boundary fixed point:
$$\mathrm{Tr}[R_{\mathrm{bdy}} G_a^2] = T(\mathrm{adj}_a) + \sum_i T_{G_a}(r_i) (R_i^{\mathrm{bdy}} - 1)$$

---

## 4. Algorithm for the Boundary Analysis

### Setup
Let the quiver have gauge nodes $\{G_1, \ldots, G_k\}$.

### For each node $a$, evaluate $\mathrm{Tr}[R G_a^2]$ at the fixed point where $G_a$ is free:

**Step 1:** Set $g_a = 0$. The theory decomposes into:
- A free sector: all matter charged under $G_a$ and not under any $G_b$, $b\neq a$
- An interacting sector: the sub-quiver with nodes $\{G_b : b \neq a\}$ with the matter bifundamentals connected to $G_a$ now as **free flavors** of the sub-quiver

**Step 2:** Find the R-charges at the sub-quiver fixed point $(g_b^*)_{b\neq a}$ using a-maximization (Module 3). Call these R-charges $R_i^{(a)}$ (R-charges at the fixed point where $G_a$ is decoupled).

**Step 3:** Compute:
$$\mathcal{B}_a \equiv \mathrm{Tr}[R^{(a)} G_a^2] = T(\mathrm{adj}_a) + \sum_i T_{G_a}(r_i) \bigl(R_i^{(a)} - 1\bigr)$$

Note: for matter fields that are only charged under $G_a$ (not under any other node), $R_i^{(a)} = 2/3$ (free field value, since they are free at this boundary).

**Step 4:** The condition for $G_a$ to flow to a non-trivial fixed point:
$$\mathcal{B}_a < 0 \quad \Leftrightarrow \quad G_a \text{ is IR free at its decoupled fixed point} \Rightarrow \text{full coupling is needed}$$

More precisely: if $\mathcal{B}_a < 0$ at the boundary $(g_a = 0)$, the beta function for $g_a$ is **negative** at small $g_a$, driving $g_a$ to grow. This is a necessary condition for the full theory to have a non-trivial IR fixed point with $g_a^* \neq 0$.

### Non-trivial IR fixed point condition:
$$\mathcal{B}_a < 0 \quad \forall  a = 1, \ldots, k$$

If this holds for all nodes, the theory is expected to flow to a non-trivial SCFT.

---

## 5. Single-Node Special Case

For a single gauge group $G$ with $N_f$ matter fields (SQCD-like), the analysis reduces to the **conformal window**:
- UV AF boundary: $b_0 > 0$, i.e., $N_f < N_f^{\mathrm{UV}}$
- IR non-trivial FP boundary: $\mathrm{Tr}[R^{\mathrm{free}} G^2] < 0$ at $g = 0$ (all fields free)

At $g = 0$, all matter fields are free: $R_i^{\mathrm{free}} = 2/3$.

$$\mathcal{B} = T(\mathrm{adj}) + \sum_i T(r_i) (2/3 - 1) = T(\mathrm{adj}) - \frac{1}{3}\sum_i T(r_i) = \frac{b_0^{\mathrm{IR}}}{3}$$

where $b_0^{\mathrm{IR}} = 3T(\mathrm{adj}) - \sum_i T(r_i)$ is the naive one-loop beta function evaluated at $R = 2/3$.

Wait — this is just the free-field beta function, which equals $b_0$ with $T(\mathrm{adj})$ contribution replaced:
$$b_0^{\mathrm{IR}} = 3T(\mathrm{adj}) - \sum_i T(r_i) = b_0^{\mathrm{UV}}$$

That's the same as the UV beta function! The difference arises because at the IR boundary we should evaluate the **exact** (NSVZ) beta function with the IR anomalous dimensions.

**Correct statement for a single node:** The lower bound of the conformal window is where the IR R-charge hits the unitarity bound $R = 2/3$. Below this, the theory is in the free magnetic phase (for SQCD, this is Seiberg duality).

---

## 6. Two-Node Analysis: Detailed Procedure

### For a 2-node quiver $G \times G'$:

Two boundary fixed points:
- **A:** $(g = g^*, g' = 0)$ — node 0 interacting, node 1 free
- **B:** $(g = 0, g' = g'^*)$ — node 1 interacting, node 0 free

**Boundary A: $g' = 0$**

1. Build single-node sub-theory for $G$: all matter at node 0 (rank-2/adj, Veneziano fundamentals) plus bifundamentals as fundamentals of $G$.
2. A-maximize the sub-theory. Bifundamental R-charges are determined dynamically by the a-maximization, **not** fixed at $R = 2/3$. Only fields exclusively charged under $G'$ (node 1's own matter) have $R = 2/3$.
3. Compute $\mathcal{B}_{G'}$ using these R-charges:
$$\mathcal{B}_{G'} = T(\mathrm{adj}_{G'}) + \sum_{i \in \text{node 1 matter}} T_{G'}(r_i)(2/3 - 1) + \sum_{i \in \text{bif}} T_{G'}(r_i)(R_i^{(A)} - 1)$$

**Boundary B: $g = 0$**

- Symmetric analysis with nodes swapped.

**We do not impose $\mathcal{B} < 0$ at both boundaries.** Instead, we compute the condition on $N_f$ such that $\mathcal{B} \geq 0$ at each boundary separately, and store these conditions in the database.

---

## 7. Relation to the Conformal Window

The conformal window for each gauge group is bounded by UV asymptotic freedom from above:
$$N_f < N_f^{\mathrm{upper}} \quad \text{(set by } b_0 > 0\text{)}$$

The boundary analysis provides an additional constraint on $N_f$ via the condition $\mathcal{B} \geq 0$ at each boundary. However, $\mathcal{B} \leq 0$ does **not** in general give the lower bound of the conformal window (cf. the SQCD counterexample in §9). The relationship between $\mathcal{B}$ and $N_f$ can be non-monotonic in quiver theories, and the condition must be evaluated case by case.

For two-node quivers, the condition at each boundary couples the matter content of both nodes.

---

## 8. Practical Computation (Two-Node)

### Step-by-step for a 2-node quiver $Q$ with nodes 0, 1:

```
for active_node in [0, 1]:
    free_node = 1 - active_node

    # Build single-node sub-theory: active node's matter + bifundamentals
    fields = build_fields_boundary(Q, active_node)
    sub_Q = single_node_quiver(Q, active_node)

    # A-maximize the sub-theory to get R-charges
    # Bifundamental R is determined by a-maximization, NOT fixed at 2/3
    R = a_maximize(sub_Q, fields)

    # Compute B for the free node
    B = T_adj_lead(G_free) * m_free
    for field_i in Q.fields_at_node(free_node):
        if field_i is bifundamental:
            B += T_{G_free}(bif) * (R[bif] - 1)    # R from a-max
        else:
            B += T_{G_free}(r_i) * (2/3 - 1)       # free field R = 2/3

    # Store condition B >= 0 as N_f constraint
    store_B_condition(Q, boundary=active_node, B=B)
```

### Notes:
- If the sub-theory has no non-trivial fixed point (below conformal window or gaugino-only), set all R = 2/3 and compute $\mathcal{B}$ with free-field values.
- For Veneziano theories, $\mathcal{B}$ depends on $N_f$ through the a-maximized R-charges. The condition $\mathcal{B} \geq 0$ gives a constraint on $N_f$.

---

## 9. Example: SQCD Conformal Window (Sanity Check)

SU(N) with $N_f$ flavors, no superpotential.

At $g = 0$: all quarks have $R = 2/3$.

$$\mathcal{B} = N + N_f \cdot \frac{1}{2} \cdot \left(\frac{2}{3} - 1\right) + N_f \cdot \frac{1}{2} \cdot \left(\frac{2}{3} - 1\right)$$
$$= N - \frac{N_f}{3}$$

Non-trivial IR FP condition: $\mathcal{B} < 0 \Rightarrow N_f > 3N$.

But wait — UV asymptotic freedom requires $N_f < 3N$! So SQCD has $\mathcal{B} > 0$ in the entire UV-AF window, meaning no interacting fixed point from this analysis?

This contradicts Seiberg's result. The resolution: the boundary analysis $\mathcal{B} < 0$ gives the condition for the coupling to be **marginally irrelevant** at the IR free fixed point, i.e., whether the free fixed point is **stable**. For $N_f > 3N/2$, the correct analysis includes the dual (magnetic) description.

Actually for SQCD:
- $3N/2 < N_f < 3N$: interacting fixed point (Banks–Zaks for $N_f$ near $3N$)
- The correct boundary analysis for a single node uses the **exact** R-charges at $g^* \neq 0$

The hep-th/0502049 method is most useful for **multi-node quivers** where coupling one node to another can drive a second node into the non-trivial phase. In that context, the R-charges at the boundary are not simply the free values $2/3$.

---

## 10. Summary

For each two-node theory, we compute:

| Quantity | At boundary A $(g', g'=0)$ | At boundary B $(g=0, g'^*)$ |
|----------|---------------------------|------------------------------|
| Sub-theory | Node 0 + bifs as funds | Node 1 + bifs as funds |
| A-maximize | Sub-theory for node 0 | Sub-theory for node 1 |
| Compute | $\mathcal{B}_{G'}$ (free node 1) | $\mathcal{B}_G$ (free node 0) |
| Store | $N_f$ condition for $\mathcal{B}_{G'} \geq 0$ | $N_f$ condition for $\mathcal{B}_G \geq 0$ |

We do **not** impose $\mathcal{B} < 0$ at both boundaries as a filter. The conditions are computed and stored for each boundary independently.

---

## 11. Implementation Plan

Scope: **two-node theories only**. For each theory, compute $\mathcal{B}$ at both boundaries A and B, and store the $N_f$ condition for $\mathcal{B} \geq 0$ as two columns in `quivers.db`.

### Required new functions (in `a_maximization_large_N.py`):

1. `build_fields_boundary(quiver, active_node)` — build leading-order fields for single-node sub-theory at boundary. Bifundamentals included as fundamentals of the active node with T_lead only at that node.
2. `compute_B_boundary(quiver, active_node, R_charges)` — compute $\mathcal{B}$ for the free node using R-charges from step 1.
3. `boundary_analysis(quiver)` — orchestrate: build sub-theory → a-maximize → compute $\mathcal{B}$ at both boundaries.

### Database changes (in `two_node_db.py`):

- Add columns `B_cond_A TEXT` and `B_cond_B TEXT` to the `theory` table.
- Add CLI command `boundary-analysis` to compute and populate these columns.

### Existing ingredients to reuse:

- `_solve_a_max(fields, quiver)` — core symbolic a-maximization solver
- `build_fields_large_N()` — pattern for field construction
- `T_adj_lead()`, `T_rep_lead()`, `T_bifund_lead()` — scaling coefficients
- `chiral_excess_coeffs()` — Veneziano parameter extraction

---

## 12. Advanced: Fixed Point Merging and Loss

As parameters (like $N_f$ or N) vary, fixed points can:
- **Merge**: two fixed points collide and annihilate → loss of conformal window
- **Decouple**: one gauge node decouples in the IR (goes to strong coupling and confines or Higgsed)

The endpoint of the conformal window is signaled by $a_{\mathrm{SCFT}} = a_{\mathrm{free}}$ at the lower end, or $b_0 = 0$ at the upper end.
