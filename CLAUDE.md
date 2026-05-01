# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Classification of 4D N=1 asymptotically free quiver gauge theories that admit a large N limit and flow to non-trivial IR superconformal fixed points.

## Repository Structure

```
README.md                       ← pipeline overview and key references
quivers.db                      ← pre-built database (4560 theories, 135 classes)

# Core pipeline
beta_functions.py               ← Module 1: one-loop b_0 computation
quiver_generation.py            ← Module 2: quiver enumeration, anomaly cancellation
a_maximization.py               ← Module 3: numerical a-max at fixed (N, N_f)
a_maximization_large_N.py       ← Module 3: exact symbolic large-N a-max
two_node_db.py                  ← Module 5: SQLite DB + CLI for two-node classification

docs/                           ← module documentation (module1–5 .md files)
scripts/                        ← dump scripts that generate paper/sections/generated/
paper/
  main.tex, references.bib      ← LaTeX paper
  sections/                     ← hand-written section files
  sections/generated/           ← auto-generated TeX (run scripts/ to regenerate)
  figures/
results/                        ← analysis output (classification JSONs, tables)
archive/                        ← DB backups, one-off scripts, build artifacts
```

## Physics Conventions

- **Gauge groups:** SU(N), SO(N), USp(2N) — written as Sp(N) throughout. All nodes share the same parameter N.
- **Dynkin index normalization:** $T(\square_{SU(N)}) = 1/2$, $T(V_{SO(N)}) = 1$, $T(f_{Sp(N)}) = 1/2$.
- **Representation notation:**
  - SU(N): $\square$, $\bar{\square}$, adj, $S$, $\bar{S}$, $A$, $\bar{A}$ (use $A$ for rank-2 antisymmetric, never $\Lambda$ or $\Lambda^2$)
  - SO(N): $V$, adj, $S$ (traceless symmetric)
  - Sp(N): $f$, adj, $A$ (traceless antisymmetric)
- **Inter-node matter:** bifundamentals only (fundamental/vector of each end). Higher representations like $(\square, \mathrm{adj})$ are excluded because their $\dim \sim N^2$ multiplicity destroys large N asymptotic freedom.
- **Large N asymptotics:** $T(\mathrm{adj}) \sim N$ for all three group types. Fundamentals/vectors contribute $\mathcal{O}(1)$ per field; rank-2 tensors contribute $\sim N/2$ (SU) or $\sim N$ (SO, Sp) per field; adjoints contribute $\sim N$ per field.

## Key Physics Facts to Preserve

- **SU(N) anomaly cancellation:** for pure bifundamental quivers, each SU(N) node must be balanced (equal in- and out-arrows). With rank-2 tensors: $N(\deg^+ - \deg^-) + (n_S - n_{\bar{S}})(N+4) + (n_A - n_{\bar{A}})(N-4) = 0$.
- **Sp(N) Witten anomaly:** total number of fundamental chiral multiplets at each Sp(N) node must be even.
- **Large N AF condition for SU(N) node:** $n_{\mathrm{adj}} + \frac{(n_S+n_{\bar{S}})+(n_A+n_{\bar{A}})}{2} + \frac{k}{2} < 3$, where $k$ = number of bifundamental neighbors.
- **IR fixed point criterion (hep-th/0502049):** at the boundary where $g_a \to 0$, compute $\mathcal{B}_a = \mathrm{Tr}[R^{(a)} G_a^2]$ using R-charges from a-maximizing the decoupled sub-quiver. Non-trivial IR fixed point requires $\mathcal{B}_a < 0$ for all nodes $a$.
- **a-maximization anomaly-free condition per gauge node:** $T(\mathrm{adj}_a) + \sum_i T_{G_a}(r_i)(R[\Phi_i] - 1) = 0$.
