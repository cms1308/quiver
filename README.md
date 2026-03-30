# Classification of 4D N=1 Asymptotically Free Quiver Gauge Theories with Large N Limit

## Project Overview

This project classifies all 4D N=1 quiver gauge theories that are:
1. **Asymptotically free** at UV
2. **Admit a large N limit** (remain AF as N → ∞)
3. **Have a non-trivial IR superconformal fixed point**

The strategy follows hep-th/0502049 and related works on a-maximization and the conformal window.

---

## Modules

| Module | Description | File |
|--------|-------------|------|
| 1 | One-loop beta function coefficients | [module1_beta_functions.md](module1_beta_functions.md) |
| 2 | Quiver diagram generation & constraints | [module2_quiver_generation.md](module2_quiver_generation.md) |
| 3 | IR R-charges via a-maximization | [module3_a_maximization.md](module3_a_maximization.md) |
| 4 | Non-trivial IR fixed point analysis | [module4_ir_fixed_points.md](module4_ir_fixed_points.md) |
| 5 | Final classification | [module5_classification.md](module5_classification.md) |

---

## Physical Setup

### Gauge Groups
Each node in the quiver carries one of:
- **SU(N)** — requires gauge anomaly cancellation
- **SO(N)** — automatically anomaly-free
- **Sp(N)** [i.e. USp(2N)] — requires even number of fundamental representations

All nodes share the **same N** (the large N parameter). Generalizations to different ranks are deferred.

### Matter Representations
Edges and node decorations can carry:

| Gauge Group | Allowed Representations |
|-------------|------------------------|
| SU(N) | fund $\square$, anti-fund $\bar{\square}$, adjoint, sym $S$, sym conjugate $\bar{S}$, antisym $A$, antisym conjugate $\bar{A}$ |
| SO(N) | vector $V$, adjoint, sym $S$ (traceless), spinor (excluded at large N) |
| Sp(N) | fund $f$, adjoint, antisym $A$ (traceless) |

Bifundamentals between nodes carry representations of both gauge groups simultaneously.

### Asymptotic Freedom Condition
$$b_0^{(a)} = 3 T(\mathrm{adj}_a) - \sum_{i} T_{G_a}(r_i) > 0 \quad \forall  a$$

### Large N Condition
$b_0^{(a)} > 0$ must hold for all N and in particular as N → ∞. This restricts the total matter content at each node to be sub-leading relative to $3T(\mathrm{adj})$.

---

## Pipeline

```
Enumerate quivers (Module 2)
        │
        ▼
Check UV asymptotic freedom (Module 1)
        │
        ▼
a-maximization → R-charges (Module 3)
        │
        ▼
Check IR fixed point existence (Module 4)
        │
        ▼
Final classification (Module 5)
```

---

## Key References

- **a-maximization**: Intriligator & Wecht, hep-th/0304128
- **IR fixed point boundary analysis**: hep-th/0502049
- **NSVZ beta function**: Novikov, Shifman, Vainshtein, Zakharov (1983)
- **Quiver gauge theories**: Douglas & Moore (1996); various reviews

---

## Status

- [ ] Module 1: Beta function computation
- [ ] Module 2: Quiver generation
- [ ] Module 3: a-maximization
- [ ] Module 4: IR fixed point analysis
- [ ] Module 5: Classification
