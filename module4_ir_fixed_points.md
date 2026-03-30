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

## 6. Multi-Node Analysis: Detailed Procedure

### For a 2-node quiver $G \times G'$:

**Boundary 1: $g = 0$**

- $G$ is free; $G'$ sees bifundamental fields as $N_{\mathrm{bif}}$ free fundamentals (with $R = 2/3$)
- Compute the R-charges of all matter in the $G'$ sector using a-maximization with the bifundamentals having $R = 2/3$
- Compute $\mathcal{B}_G$ using these R-charges

**Boundary 2: $g' = 0$**

- Symmetric analysis for $G'$

**Non-trivial fixed point exists if both $\mathcal{B}_G < 0$ and $\mathcal{B}_{G'} < 0$.**

### For a $k$-node quiver:

For each node $a$, the decoupled theory is the sub-quiver with $G_a$ removed. One can also consider **simultaneous decoupling** of subsets of nodes, but typically the pairwise analysis is sufficient to determine the existence of a non-trivial fixed point.

---

## 7. Relation to the Conformal Window

The conformal window for each gauge group in the quiver is bounded by:
$$N_f^{\mathrm{lower}} < N_f < N_f^{\mathrm{upper}}$$

where:
- $N_f^{\mathrm{upper}}$: set by $b_0 > 0$ (UV AF)
- $N_f^{\mathrm{lower}}$: set by $\mathcal{B}_a < 0$ at the boundary (non-trivial IR FP)

For quiver theories, this analysis couples the conformal windows of different nodes.

---

## 8. Practical Computation

### Step-by-step for a quiver $Q$:

```
for each gauge node a in Q:
    # Build the decoupled theory Q_a = Q with G_a set free
    Q_a = Q.decouple(a)

    # All fields charged under G_a become free fields in Q_a
    # They contribute as external flavors with R = 2/3 to the remaining nodes

    # a-maximize the remaining sub-quiver Q_a to get R^{(a)}
    R_a = a_maximize(Q_a, treat_external_as_free=True)

    # Compute Tr[R^{(a)} G_a^2]
    B_a = T_adj(G_a)
    for field_i in Q.fields_charged_under(a):
        B_a += T_{G_a}(r_i) * (R_a[field_i] - 1)

    if B_a >= 0:
        # G_a does not need to become strongly coupled
        # → no non-trivial IR fixed point with G_a interacting
        Q.mark_no_nontrivial_fp()
        break

# If all B_a < 0: non-trivial IR fixed point exists
```

### Notes:
- If $Q_a$ itself has no non-trivial fixed point (i.e., the sub-quiver is in the free or confined phase), the boundary analysis is more subtle.
- For free sub-quivers, set all remaining fields to $R = 2/3$ and compute $\mathcal{B}_a$ with these values.

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

## 10. Summary of Conditions

A quiver $Q$ has a non-trivial IR superconformal fixed point if:

| Condition | Meaning |
|-----------|---------|
| $b_0^{(a)} > 0$ for all $a$ | UV asymptotically free |
| $\mathcal{B}_a < 0$ for all $a$ | Each node is driven to strong coupling |
| $R_i^* \geq 2/3$ for all chiral fields | No operator hits unitarity bound (or must decouple) |

Additionally, the existence of an IR fixed point is confirmed if:
- The a-function has a local maximum (not just a saddle)
- The theory is not equivalent to a free or confined theory in the IR

---

## 11. Advanced: Fixed Point Merging and Loss

As parameters (like $N_f$ or N) vary, fixed points can:
- **Merge**: two fixed points collide and annihilate → loss of conformal window
- **Decouple**: one gauge node decouples in the IR (goes to strong coupling and confines or Higgsed)

The endpoint of the conformal window is signaled by $a_{\mathrm{SCFT}} = a_{\mathrm{free}}$ at the lower end, or $b_0 = 0$ at the upper end.
