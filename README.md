# Classification of 4D N=1 Asymptotically Free Quiver Gauge Theories with Large N Limit

## Project Overview

This project classifies all 4D N=1 quiver gauge theories that are:
1. **Asymptotically free** at UV
2. **Admit a large N limit** (remain AF as N → ∞)
3. **Have a non-trivial IR superconformal fixed point**

The strategy follows hep-th/0502049 and related works on a-maximization and the conformal window.

---

## Repository Structure

```
quivers.db                      ← pre-built database (4560 theories, 135 classes)

# Core pipeline
beta_functions.py               ← Module 1: one-loop b_0 computation
quiver_generation.py            ← Module 2: quiver enumeration, anomaly cancellation
a_maximization.py               ← Module 3: numerical a-max at fixed (N, N_f)
a_maximization_large_N.py       ← Module 3: exact symbolic large-N a-max
two_node_db.py                  ← Module 5: SQLite DB + CLI for two-node classification

docs/                           ← module documentation (module1–5 .md)
scripts/                        ← dump scripts → paper/sections/generated/
paper/
  main.tex, references.bib
  sections/                     ← hand-written LaTeX sections
  sections/generated/           ← auto-generated TeX (run scripts/ to regenerate)
  figures/
results/                        ← analysis output (classification JSONs, tables)
archive/                        ← DB backups, one-off scripts
```

---

## Modules

| Module | Description | Documentation | Implementation |
|--------|-------------|---------------|----------------|
| 1 | One-loop beta function coefficients | [docs/module1](docs/module1_beta_functions.md) | `beta_functions.py` |
| 2 | Quiver enumeration & constraints | [docs/module2](docs/module2_quiver_generation.md) | `quiver_generation.py` |
| 3 | IR R-charges via a-maximization | [docs/module3](docs/module3_a_maximization.md) | `a_maximization.py`, `a_maximization_large_N.py` |
| 4 | Non-trivial IR fixed point analysis | [docs/module4](docs/module4_ir_fixed_points.md) | — |
| 5 | Classification & database | [docs/module5](docs/module5_classification.md) | `two_node_db.py` |

---

## Physical Setup

### Gauge Groups
Each node carries one of **SU(N)**, **SO(N)**, or **Sp(N)** [i.e. USp(2N)]. All nodes share the same N (the large-N parameter).

### Matter Representations

| Gauge Group | Representations |
|-------------|----------------|
| SU(N) | $\square$, $\bar\square$, adj, $S$, $\bar S$, $A$, $\bar A$ |
| SO(N) | $V$, adj, $S$ (traceless sym) |
| Sp(N) | $f$, adj, $A$ (traceless antisym) |

Inter-node matter is bifundamental only; higher representations (e.g. $(\square, \mathrm{adj})$) are excluded because their $\dim \sim N^2$ multiplicity destroys large-N AF.

### Pipeline

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
Classification & database (Module 5)
```

---

## Database

`quivers.db` contains the full two-node classification:
- **4560 theories** across **135 universality classes**
- Theories grouped by gauge pair (SU-SU, SU-SO, SU-Sp, SO-SO, SO-Sp, Sp-Sp)
- Non-Veneziano classes: 9 classes (SU-SU only) with exact symbolic a-maximization results
- R-charges exact (symbolic) where solvable, numerical otherwise

To query: `python3 two_node_db.py --help`

## Generating Paper Sections

The auto-generated LaTeX sections in `paper/sections/generated/` are produced by:

```bash
python3 scripts/dump_tex_tables.py          # classification_tables.tex
python3 scripts/dump_no_veneziano_tables.py # no_veneziano_tables.tex
python3 scripts/dump_nf_bcond_table.py      # nf_bcond_table.tex
```

---

## Key References

- **a-maximization**: Intriligator & Wecht, hep-th/0304128
- **IR fixed point boundary analysis**: hep-th/0502049
- **NSVZ beta function**: Novikov, Shifman, Vainshtein, Zakharov (1983)
