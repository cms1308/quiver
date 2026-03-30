# Module 1: One-Loop Beta Function Coefficients

## Goal
For each gauge node $G_a$ in a quiver, compute the one-loop beta function coefficient $b_0^{(a)}$ and determine whether the node (and the full theory) is asymptotically free at UV.

---

## 1. The One-Loop Beta Function

For a 4D N=1 supersymmetric gauge theory with gauge group $G$ and chiral superfields in representations $r_i$, the one-loop beta function is:

$$\mu \frac{dg}{d\mu} = -\frac{g^3}{16\pi^2} b_0 + \mathcal{O}(g^5)$$

$$\boxed{b_0 = 3\,T(\mathrm{adj}) - \sum_i T(r_i)}$$

where $T(r)$ is the **Dynkin index** of representation $r$, defined by:
$$\mathrm{tr}_r[T^a T^b] = T(r)\,\delta^{ab}$$

**Asymptotic freedom** requires $b_0 > 0$.

The sum runs over all chiral superfields charged under $G$. For a bifundamental field $(\square_a, \bar{\square}_b)$ shared between two nodes, it contributes to **both** $b_0^{(a)}$ and $b_0^{(b)}$.

---

## 2. Dynkin Indices

### SU(N) — normalization: $T(\square) = \frac{1}{2}$

| Representation | Symbol | Dimension | Dynkin Index $T$ | Anomaly $A$ |
|----------------|--------|-----------|-------------------|-------------|
| Fundamental | $\square$ | $N$ | $\frac{1}{2}$ | $+1$ |
| Anti-fundamental | $\bar{\square}$ | $N$ | $\frac{1}{2}$ | $-1$ |
| Adjoint | $\mathrm{adj}$ | $N^2 - 1$ | $N$ | $0$ |
| Rank-2 symmetric | $S^2$ | $\frac{N(N+1)}{2}$ | $\frac{N+2}{2}$ | $N+4$ |
| Rank-2 antisymmetric | $\Lambda^2$ | $\frac{N(N-1)}{2}$ | $\frac{N-2}{2}$ | $N-4$ |

The anomaly coefficient $A(r)$ enters gauge anomaly cancellation (see Module 2).

**Beta function for SU(N):**
$$b_0^{SU(N)} = 3N - \frac{n_f}{2} - N\,n_{\mathrm{adj}} - \frac{N+2}{2}\,n_S - \frac{N-2}{2}\,n_A - [\text{bifundamental contributions}]$$

where $n_f, n_{\mathrm{adj}}, n_S, n_A$ are the numbers of fundamentals, adjoints, symmetric, and antisymmetric chiral multiplets attached to this node.

---

### SO(N) — normalization: $T(V) = 1$

| Representation | Symbol | Dimension | Dynkin Index $T$ |
|----------------|--------|-----------|-------------------|
| Vector (fundamental) | $V$ | $N$ | $1$ |
| Adjoint (= rank-2 antisym) | $\mathrm{adj}$ | $\frac{N(N-1)}{2}$ | $N-2$ |
| Rank-2 symmetric (traceless) | $S^2$ | $\frac{N(N+1)}{2}-1$ | $N+2$ |
| Spinor | $\mathbf{s}$ | $2^{\lfloor N/2 \rfloor - 1}$ | grows exponentially |

> **Note:** Spinor representations grow exponentially with $N$ and are **excluded** in the large N classification since they violate large N factorization.

SO(N) is automatically free of gauge anomalies (no $A(r)$ constraint needed).

**Beta function for SO(N):**
$$b_0^{SO(N)} = 3(N-2) - n_V \cdot 1 - n_{\mathrm{adj}}(N-2) - n_S(N+2) - [\text{bifundamental contributions}]$$

---

### USp(2N) = Sp(N) — normalization: $T(f) = \frac{1}{2}$

The group USp(2N) has rank $N$; its fundamental representation has dimension $2N$.

| Representation | Symbol | Dimension | Dynkin Index $T$ |
|----------------|--------|-----------|-------------------|
| Fundamental | $f$ | $2N$ | $\frac{1}{2}$ |
| Adjoint (= rank-2 symmetric) | $\mathrm{adj}$ | $N(2N+1)$ | $N+1$ |
| Rank-2 antisymmetric (traceless) | $\Lambda^2$ | $N(2N-1)-1$ | $N-1$ |

> **Note:** For USp(2) $\cong$ SU(2), the antisymmetric rank-2 is a singlet and plays no role.

For Sp(N), the number of fundamental representations must be **even** to avoid the Witten anomaly (a global gauge anomaly).

**Beta function for Sp(N):**
$$b_0^{Sp(N)} = 3(N+1) - \frac{n_f}{2} - n_{\mathrm{adj}}(N+1) - n_A(N-1) - [\text{bifundamental contributions}]$$

---

## 3. Bifundamental Contributions

When two gauge nodes $G_a$ and $G_b$ are connected by a bifundamental edge, the matter field is charged under both groups. The contribution to each node's beta function:

### SU(N) × SU(N) bifundamental $(\square_a, \bar{\square}_b)$:
- To node $a$: $T_{G_a} = N \cdot \frac{1}{2} = \frac{N}{2}$  (there are $N$ fundamentals of $G_a$, one per color of $G_b$)
- To node $b$: $T_{G_b} = N \cdot \frac{1}{2} = \frac{N}{2}$

### SU(N) × SO(N) bifundamental $(\square_{SU}, V_{SO})$:
- To SU(N) node: $T_{SU} = N \cdot \frac{1}{2} = \frac{N}{2}$ (N fundamentals, one per SO(N) vector index)
- To SO(N) node: $T_{SO} = N \cdot 1 = N$ (N vectors, one per SU(N) color)

### SU(N) × Sp(N) bifundamental $(\square_{SU}, f_{Sp})$:
- To SU(N) node: $T_{SU} = 2N \cdot \frac{1}{2} = N$
- To Sp(N) node: $T_{Sp} = N \cdot \frac{1}{2} = \frac{N}{2}$

### SO(N) × Sp(N) bifundamental $(V_{SO}, f_{Sp})$:
- To SO(N) node: $T_{SO} = 2N \cdot 1 = 2N$
- To Sp(N) node: $T_{Sp} = N \cdot \frac{1}{2} = \frac{N}{2}$

---

## 4. Large N Asymptotics

At large N, the leading-N behavior of $b_0$ determines whether the theory remains AF.

| Group | $3T(\mathrm{adj})$ at large $N$ | Leading matter contribution |
|-------|-------------------------------|----------------------------|
| SU(N) | $3N$ | each $\square/S^2/\Lambda^2$: $\sim N/2$ each |
| SO(N) | $3N$ | each vector: $\sim 1$; each sym: $\sim N$ |
| Sp(N) | $3N$ | each fund: $\sim 1/2$; each $\Lambda^2$: $\sim N$ |

**Large N asymptotic freedom requires that the $\mathcal{O}(N)$ part of $b_0$ is positive:**

For **SU(N)** with $n_{\mathrm{adj}}$ adjoints, $n_S$ symmetric, $n_A$ antisymmetric, and $k$ bifundamental neighbors:
$$b_0 \approx N\!\left(3 - n_{\mathrm{adj}} - \frac{n_S + n_A}{2} - \frac{k}{2}\right) + \mathcal{O}(1)$$

Large N AF condition:
$$n_{\mathrm{adj}} + \frac{n_S + n_A}{2} + \frac{k}{2} < 3$$

where $k$ = number of bifundamental edges to other nodes (each contributing $N/2$ at large N regardless of the neighbor's gauge group).

---

## 5. Algorithm

```
for each quiver candidate Q:
    for each gauge node a with group G_a:
        b0[a] = 3 * T_adj(G_a)
        for each matter field i charged under G_a:
            b0[a] -= T_{G_a}(r_i)
        if b0[a] <= 0:
            reject Q  # not asymptotically free
    accept Q as UV-AF
```

### Input
- Quiver graph: nodes labeled by gauge group type, edges labeled by bifundamental type
- Additional node matter: adjoints, symmetric, antisymmetric at each node

### Output
- $b_0^{(a)}$ for each node
- Boolean: is the full quiver UV-asymptotically free?
- Boolean: does it remain AF at large N?

---

## 6. Simple Examples

### Example 1: SU(N) with $N_f$ fundamentals
$$b_0 = 3N - \frac{N_f}{2}$$
AF condition: $N_f < 6N$.
Large N limit: $N_f = \kappa N$ with $\kappa < 6$.

### Example 2: Linear quiver SU(N)$^k$ with bifundamentals only
Each interior node has 2 bifundamental neighbors:
$$b_0^{\mathrm{int}} = 3N - 2 \cdot \frac{N}{2} = 2N > 0 \checkmark$$
Each end node has 1 neighbor:
$$b_0^{\mathrm{end}} = 3N - \frac{N}{2} = \frac{5N}{2} > 0 \checkmark$$

### Example 3: Circular quiver SU(N)$^k$ with bifundamentals only
Each node has 2 neighbors: $b_0 = 2N > 0$ for all $k$. ✓

### Example 4: SU(N) with one adjoint chiral
$$b_0 = 3N - N = 2N > 0 \checkmark$$
(This is the $\mathcal{N}=2$ vector multiplet contribution.)

### Example 5: SU(N) with one rank-2 symmetric
$$b_0 = 3N - \frac{N+2}{2} = \frac{5N-2}{2} > 0 \checkmark$$

---

## 7. Notes on the NSVZ Exact Beta Function

At higher loops (and exactly via NSVZ):
$$\beta(g) = -\frac{g^3}{16\pi^2}\,\frac{3T(\mathrm{adj}) - \sum_i T(r_i)(1 - 2\gamma_i)}{1 - T(\mathrm{adj})\,g^2/(8\pi^2)}$$

where $\gamma_i$ is the anomalous dimension of field $i$. At the IR fixed point, $\beta(g) = 0$ (numerator = 0) and $\gamma_i = R_i - 1$ (from superconformal algebra), recovering the condition:
$$3T(\mathrm{adj}) = \sum_i T(r_i)(3 - 2R_i) \quad \Leftrightarrow \quad \mathrm{Tr}[R\,G^2] = 0$$

This links the beta function vanishing to the anomaly-free R-symmetry condition used in a-maximization (Module 3).
